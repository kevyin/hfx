/* -----------------------------------------------------------------------------
 *
 * Module    : Ion Series
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include "utils.h"
#include "device.h"
#include "texture.h"
#include "ion_series.h"
#include "algorithms.h"

#include <stdint.h>

#define DEBUG

#ifdef DEBUG
#define MAX_PEP_LEN 1000
#endif

/*
 * Scan a warp-sized chunk of data. Because warps execute instructions in SIMD
 * fashion, there is no need to synchronise in order to share data. The most
 * efficient algorithm is the step-efficient method of Hillis & Steele that
 * takes log(N) steps, rather than the work-efficient tree-based algorithm
 * described by Blelloch that takes 2 * log(N) steps.
 */
template <class T, bool inclusive>
static __device__ T
scan_warp(T val, volatile T* s_data)
{
    const uint32_t idx  = threadIdx.x;
    const uint32_t lane = threadIdx.x & (WARP_SIZE-1);

    /*
     * If we double the size of the s_data array and pad the bottom half with
     * zero, then we can avoid branching (although there is plenty already).
     *
     * In device emulation mode, the warp size is 1 and so sync-less operation
     * does not work.
     */
    s_data[idx] = val;                                                        __EMUSYNC;
#ifdef __DEVICE_EMULATION__
    val = (lane >=  1) ? s_data[idx -  1] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  2) ? s_data[idx -  2] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  4) ? s_data[idx -  4] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  8) ? s_data[idx -  8] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >= 16) ? s_data[idx - 16] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
#else
    if (lane >=  1) s_data[idx] = val = val + s_data[idx -  1];
    if (lane >=  2) s_data[idx] = val = val + s_data[idx -  2];
    if (lane >=  4) s_data[idx] = val = val + s_data[idx -  4];
    if (lane >=  8) s_data[idx] = val = val + s_data[idx -  8];
    if (lane >= 16) s_data[idx] = val = val + s_data[idx - 16];
#endif

    if (inclusive) return s_data[idx];
    else           return (lane > 0) ? s_data[idx - 1] : 0;
}


__inline__ __device__ static float
ionMZ(const float m, const float c)
{
    return __fdividef(m + MASS_H * c, c);
}

__inline__ __device__ static uint32_t
binMZ(const float mz)
{
    return rintf(__fdividef(mz, BIN_WIDTH_MONO));
}

__inline__ __device__ static void
addIon(uint32_t *d_spec, const uint32_t N, const int32_t x, const uint32_t y)
{
    if (0 <= x && x < N) atomicMax(&d_spec[x], y);
}


template <uint32_t charge>
__device__ void
addIonsAB(uint32_t *d_spec, const uint32_t N, const float mass)
{
    float   m;
    int32_t x;

    // A-ions
    addIon(d_spec, N, binMZ(ionMZ(mass - MASS_CO, charge)), 10);

    // B-ions
    m = ionMZ(mass, charge);
    x = binMZ(m);

    addIon(d_spec, N, x,   50);
    addIon(d_spec, N, x+1, 25); // technically, should be binMZ(m+1)
    addIon(d_spec, N, x-1, 25);

    addIon(d_spec, N, binMZ(m - __fdividef(MASS_H2O, charge)), 10);
    addIon(d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10);
}


template <uint32_t charge>
__device__ void
addIonsY(uint32_t *d_spec, const uint32_t N, const float mass)
{
    float   m = ionMZ(mass + MASS_H2O, charge);
    int32_t x = binMZ(m);

    // Y-ions
    addIon(d_spec, N, x,   50);
    addIon(d_spec, N, x+1, 25);
    addIon(d_spec, N, x-1, 25);

    addIon(d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10);
}


template <uint32_t charge>
__device__ void
addIons_k(uint32_t *d_spec, const uint32_t N, const float b_mass, const float y_mass)
{
    addIonsAB<charge>(d_spec, N, b_mass);
    addIonsY <charge>(d_spec, N, y_mass);
}


/*
 * Return the mass of an amino acid residue in atomic mass units, for the given
 * short abbreviation.
 */
template <bool UseCache>
__device__ float
getAAMass(const float *d_mass, const char aa)
{
    return fetch_x<UseCache>(aa - 'A', d_mass);
}

/*
 * Determine how much the modification shift the residual mass by
 */
float
getModResDelta(const float *d_mod_ma_mass, const uint8_t *d_mod_ma_count, const uint32_t mod_num_ma) 
{
    thrust::device_ptr<const float>   d_mod_ma_mass_th(d_mod_ma_mass);
    thrust::device_ptr<const uint8_t> d_mod_ma_count_th(d_mod_ma_count);
    //float delta = 0;
    //for (uint32_t i = 0; i < mod_num_ma; ++i)
    //{
        //delta += d_mod_ma_mass_th[i]*d_mod_ma_count_th[i];
    //}
    //return delta;
    
    return thrust::inner_product(d_mod_ma_mass_th, d_mod_ma_mass_th + mod_num_ma, d_mod_ma_count_th, (float)0);
}

__device__ bool
isToBeModded(const uint32_t *d_mpep_unrank, const uint8_t *d_mod_ma_count, uint32_t ma_idx, uint32_t ith_ma)
{
    bool res = false;
    uint32_t unrank_idx = 0;
    for (uint32_t i = 0; i < ma_idx; ++i)
    {
        unrank_idx += d_mod_ma_count[i];
    }

    for (uint32_t i = 0; i < d_mod_ma_count[ma_idx]; ++i)
    {
        if (d_mpep_unrank[unrank_idx + i] == ith_ma) 
            res = true;
    }
    return res;
}


/*
 * Generate theoretical spectra for a collection of peptide fragments. The
 * 'ions' array contains the individual amino-acid masses for the database
 * entries. We are interested in the sequences generated between the terminal
 * indices (tc,tn) of the locations specified in the 'idx' array.
 *
 * A warp of threads iterates between the (tc,tn) indices, generating the b- and
 * y-ion mass ladders. A (long) sequence of (slow) global atomic update requests
 * is subsequently issued. The input d_spec should be initially zero, and on
 * output will contain the theoretical spectra peaks in a dense (although
 * mostly zero) matrix.
 */
template <uint32_t BlockSize, uint32_t MaxCharge, bool UseCache>
__global__ static void
addModIons_core
(
    uint32_t            *d_mspec,
    uint8_t             *d_mions,
    const float         *d_residual,    // peptide residual mass
    const float         *d_mass,        // lookup table for ion character codes ['A'..'Z']
    const uint8_t       *d_ions,        // individual ion character codes (the database)
    const uint32_t      *d_tc,          // c-terminal indices
    const uint32_t      *d_tn,          // n-terminal indices
    const uint32_t      *d_mpep_idx,
    const uint32_t      *d_mpep_unrank,
    const uint32_t      num_mpep,
    const uint8_t       *d_mod_ma,
    const uint8_t       *d_mod_ma_count,
    const float         *d_mod_ma_mass,
    const uint32_t      mod_num_ma,
    const uint32_t      len_spec,
    const float         mass_delta,
    const uint32_t      total_ma_count
)
{
    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t vectorsPerBlock = BlockSize / WARP_SIZE;
    const uint32_t numVectors      = vectorsPerBlock * gridDim.x;
    const uint32_t thread_id       = BlockSize * blockIdx.x + threadIdx.x;
    const uint32_t vector_id       = thread_id / WARP_SIZE;
    const uint32_t thread_lane     = threadIdx.x & (WARP_SIZE-1);

    __shared__ volatile float s_data[BlockSize];
    // Keep a record of ith moddable acid as the pep is traversed
    __shared__ volatile uint32_t s_pep_ith_ma[MAX_MA][BlockSize];
    //__shared__ volatile uint32_t s_test[BlockSize];


    for (uint32_t row = vector_id; row < num_mpep; row += numVectors)
    {
        const uint32_t mpep_idx  = d_mpep_idx[row];
        const uint32_t row_start = d_tc[mpep_idx];
        const uint32_t row_end   = d_tn[mpep_idx];
        const float    residual  = d_residual[mpep_idx] + mass_delta;

        uint32_t       *spec     = &d_mspec[row * len_spec];
        float          b_mass;
        float          y_mass;


        s_data[threadIdx.x]      = 0;
        //s_test[threadIdx.x]      = 0;
        for (int mod = 0; mod < mod_num_ma; mod++)
        {
            s_pep_ith_ma[mod][threadIdx.x] = 0;
        }

        /*
         * Have all threads read in values for this segment, writing the
         * spectral peaks out to global memory (very, very slowly...)
         */
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            bool is_ma = false;
            uint32_t ma_idx = 0;
            /*
             * Load the ion mass, and propagate the partial scan results
             */
            // is this a modable acid
            for (int m = 0; m < mod_num_ma; m++) 
            {
                uint32_t count = 0;
                if (d_mod_ma[m] == d_ions[j]) 
                {
                    is_ma = true;
                    count = 1;
                    ma_idx = m;
                }

                if (thread_lane == 0) {
                    count += s_pep_ith_ma[m][threadIdx.x + (WARP_SIZE-1)];
                    //count += s_test[threadIdx.x + (WARP_SIZE-1)];
                }
                scan_warp<uint32_t, true>(count, s_pep_ith_ma[m]); 
            }

            uint32_t ith_ma = s_pep_ith_ma[ma_idx][threadIdx.x] - 1;

            //if (is_ma)
            bool modded = isToBeModded(d_mpep_unrank + row*total_ma_count, d_mod_ma_count, ma_idx, ith_ma);

            if (is_ma && modded)
            {
                b_mass = getAAMass<UseCache>(d_mass, d_ions[j]) + d_mod_ma_mass[ma_idx];
                d_mions[row*MAX_PEP_LEN + j - row_start] = d_ions[j] + 32;
            } else {
                b_mass = getAAMass<UseCache>(d_mass, d_ions[j]);
                d_mions[row*MAX_PEP_LEN + j - row_start] = d_ions[j];
            }
            d_mions[row*MAX_PEP_LEN + j - row_start] = (uint8_t)(((int)'0')+(ith_ma));
                

            if (thread_lane == 0)
            {
                b_mass += s_data[threadIdx.x + (WARP_SIZE-1)];
            }

            /*
             * Generate fragment mass ladder
             */
            b_mass = scan_warp<float,true>(b_mass, s_data);
            y_mass = residual - b_mass;

            if (1 <= MaxCharge) addIons_k<1>(spec, len_spec, b_mass, y_mass);
            if (2 <= MaxCharge) addIons_k<2>(spec, len_spec, b_mass, y_mass);
            if (3 <= MaxCharge) addIons_k<3>(spec, len_spec, b_mass, y_mass);
            if (4 <= MaxCharge) addIons_k<4>(spec, len_spec, b_mass, y_mass);
        }
    }
}


/*
 * Select a number of threads and blocks. Each block will have at least one full
 * warp, as required by the core kernel
 */
static void
addModIons_control(uint32_t N, uint32_t &blocks, uint32_t &threads)
{
    threads = (N < MAX_THREADS) ? max(WARP_SIZE, ceilPow2(N)) : MAX_THREADS;
    blocks  = (N + threads - 1) / threads;
    blocks  = min(blocks, MAX_BLOCKS);
}


template <uint32_t MaxCharge, bool UseCache>
static void
addModIons_dispatch
(
    uint32_t            *d_mspec,
    uint8_t             *d_mions,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_mpep_idx,
    const uint32_t      *d_mpep_unrank,
    const uint32_t      num_mpep,
    const uint8_t       *d_mod_ma,
    const uint8_t       *d_mod_ma_count,
    const float         *d_mod_ma_mass,
    const uint32_t      mod_num_ma,
    const uint32_t      len_spec,
    const float         mass_delta,
    const uint32_t      total_ma_count
)
{
    uint32_t blocks;
    uint32_t threads;

    if (UseCache)
        bind_x(d_mass);

    addModIons_control(num_mpep, blocks, threads);
    switch (threads)
    {
    //case 512: addModIons_core<512,MaxCharge,UseCache><<<blocks,threads>>>(d_mspec, d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx,, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, len_spec, mass_delta, total_ma_count); break;
    //case 256: addModIons_core<256,MaxCharge,UseCache><<<blocks,threads>>>(d_mspec, d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx,, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, len_spec, mass_delta, total_ma_count); break;
    case 128: addModIons_core<128,MaxCharge,UseCache><<<blocks,threads>>>(d_mspec, d_mions, d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, mass_delta, total_ma_count); break;
    case  64: addModIons_core< 64,MaxCharge,UseCache><<<blocks,threads>>>(d_mspec, d_mions, d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, mass_delta, total_ma_count); break;
    case  32: addModIons_core< 32,MaxCharge,UseCache><<<blocks,threads>>>(d_mspec, d_mions, d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, mass_delta, total_ma_count); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }

    if (UseCache)
      unbind_x(d_mass);
}

void addModIons
(
    uint32_t            *d_out_mspec,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_mpep_idx,
    const uint32_t      *d_mpep_unrank,
    const uint32_t      num_mpep,
    const uint8_t       *d_mod_ma,
    const uint8_t       *d_mod_ma_count,
    const float         *d_mod_ma_mass,
    const uint32_t      mod_num_ma,
    const uint32_t      max_charge,
    const uint32_t      len_spec
)
{
    thrust::device_vector<uint8_t> d_mions(MAX_PEP_LEN*num_mpep);

    float delta = getModResDelta(d_mod_ma_mass, d_mod_ma_count, mod_num_ma);

    thrust::device_ptr<const uint8_t> d_mod_ma_count_th(d_mod_ma_count);
    const uint32_t total_ma_count = thrust::reduce(d_mod_ma_count_th, d_mod_ma_count_th + mod_num_ma);

    switch (max_charge)
    {
    case 1: addModIons_dispatch<1,true>(d_out_mspec,  d_mions.data().get(), d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, delta, total_ma_count); break;
    case 2: addModIons_dispatch<2,true>(d_out_mspec, d_mions.data().get(), d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, delta, total_ma_count); break;
    case 3: addModIons_dispatch<3,true>(d_out_mspec, d_mions.data().get(), d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, delta, total_ma_count); break;
    case 4: addModIons_dispatch<4,true>(d_out_mspec, d_mions.data().get(), d_residual, d_mass, d_ions, d_tc, d_tn, d_mpep_idx, d_mpep_unrank, num_mpep, d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, len_spec, delta, total_ma_count); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }

#ifdef DEBUG
    std::cout << "Checking generated spectrums" << std::endl;
    std::cout << "delta was: " << delta << std::endl;
    // To check the spectrums above, create these modified peptides serially and call addIons to create spectrums. Then compare the two spectrums
    // NB currently only allows 1 alternative mass for an acid. lowercase is used for the modified mass
    thrust::device_vector<uint32_t> d_out_check_spec(len_spec*num_mpep);
    getSpecNonParallel(d_out_check_spec.data().get(),
                     d_mions.data().get(),
                     d_residual, d_mass, d_ions, d_tc, d_tn, 
                     d_mpep_idx, d_mpep_unrank, num_mpep,
                     d_mod_ma, d_mod_ma_count, d_mod_ma_mass, mod_num_ma, delta,
                     max_charge, len_spec);

    // compare
    thrust::device_ptr<uint32_t> d_out_mspec_th(d_out_mspec);
    if (!thrust::equal(d_out_check_spec.begin(), d_out_check_spec.end(), d_out_mspec_th)) {
        std::cerr << "Spectrums doesn't seem to be correct" << std::endl;

        uint32_t cnt = 0;
        for (uint32_t i = 0; i < num_mpep; ++i) {
            for (uint32_t j = 0; j < len_spec; ++j) {
                uint32_t pos = i*len_spec + j;
                if (d_out_check_spec[pos] != d_out_mspec_th[pos]) {
                    //std::cout << "check " << d_out_check_spec[pos] << " != " << d_out_mspec_th[pos] << std::endl;
                    ++cnt;
                    break;
                } else {
                }

            }
        }
        std::cout << "num specs not right: " << cnt << " out of " << num_mpep << std::endl;
        exit(1);
    } else {
        std::cout << "spectrum seems ok" << std:: endl;
    }
#endif
#undef DEBUG

}

