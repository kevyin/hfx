/* -----------------------------------------------------------------------------
 *
 * Module    : 
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
#include "algorithms.h"
#include "functors.h"
#include "functional.hpp"

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/transform_scan.h>
#include <thrust/copy.h>
#include <thrust/scan.h>

#include <time.h>
#include <stdio.h>

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

template <uint32_t BlockSize, uint32_t MaxMA>
__global__ static void
findModablePeptides_core
(
    uint32_t            *d_valid,
    uint32_t            *d_pep_ma_count,    // 2d array, count of each ma in each peptide
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    
    const uint8_t       *d_ma,
    const uint32_t      num_ma,

    const uint32_t      *d_mod_ma_count,

    uint32_t            *d_pep_idx, 
    uint32_t            *d_pep_mod_idx, 

    uint32_t            num_pep_total
)
{

    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t vectorsPerBlock = BlockSize / WARP_SIZE;
    const uint32_t numVectors      = vectorsPerBlock * gridDim.x;
    const uint32_t thread_id       = BlockSize * blockIdx.x + threadIdx.x;
    const uint32_t vector_id       = thread_id / WARP_SIZE;
    const uint32_t thread_lane     = threadIdx.x & (WARP_SIZE-1);

    assert(MaxMA >= num_ma);

    // Keep a count for each mod
    __shared__ volatile uint32_t s_data[MaxMA][BlockSize];

    for (uint32_t row = vector_id; row < num_pep_total; row += numVectors)
    {
        const uint32_t idx       = d_pep_idx[row];
        const uint32_t row_start = d_tc[idx];
        const uint32_t row_end   = d_tn[idx];

        const uint32_t mod_idx   = d_pep_mod_idx[row];

        for (int mod = 0; mod < num_ma; mod++)
        {
            s_data[mod][threadIdx.x] = 0;
        }

        // for each acid check if this is a modable acid, if so record the count
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            bool isModable = true; // assume modability
            for (int mod = 0; mod < num_ma; mod++) 
            {
                uint32_t count = 0;
                if (d_ma[mod] == d_ions[j])
                    count++;

                if (thread_lane == 0)
                    count += s_data[mod][threadIdx.x + (WARP_SIZE-1)];
                
                count = scan_warp<uint32_t, true>(count, s_data[mod]); 

                // if last aa in peptide
                if (j == row_end-1) { 
                    // check if not enough aa for mod
                    if (count >= d_mod_ma_count[mod_idx*num_ma + mod]) 
                    {
                        //isModable = true;
                        // record pep aa counts
                        d_pep_ma_count[row*num_ma + mod] = count;
                    } else {
                        isModable = false;
                    }
                }
            }

            if (j == row_end-1) // last aa
            {
                // @TODO ?? if modable, calculate and record total number of modified peptides else 0
                //       ?? requires choose function so may not be best place to do it?
                d_valid[row] = (isModable) ? 1 : 0;
            }
        }

    }
}

/**
 * Select a number of threads and blocks. Each block will have at least one full
 * warp, as required by the core kernel
 */
static void
findByMod_control(uint32_t N, uint32_t &blocks, uint32_t &threads)
{
    threads = (N < MAX_THREADS) ? max(WARP_SIZE, ceilPow2(N)) : MAX_THREADS;
    blocks  = (N + threads - 1) / threads;
    blocks  = min(blocks, MAX_BLOCKS);
}

/**
 * fillMatrixRow
 * Used with in conjunction with transform
 * Fills a matrix where values of a row are taken from rowValues
 * eg row values: [3, 4, 2]
 * matrix:
 * 3, 3, 3, 3, ...
 * 4, 4, 4, 4, ...
 * 2, 2, 2, 2, ...
 */
template <typename T>
struct fillMatrixRow: public thrust::unary_function<uint32_t,T>
{

    thrust::device_ptr<T> rowValues;
    uint32_t width;

    __host__ __device__
    fillMatrixRow(thrust::device_ptr<T> _values, uint32_t _w) : rowValues(_values), width(_w) {}

    __host__ __device__ bool operator() (uint32_t e)
    {
        uint32_t r = e/width;
        return rowValues[r];
    }
};

//#ifdef _DEBUG

/**
 * filterByModNonParallel
 */
//void
//filterByModNonParallel 
//(
    //thrust::host_vector<uint32_t>& h_check_valid,
    //thrust::host_vector<uint32_t>& h_check_pep_ma_count,

    //thrust::device_ptr<const uint8_t>&  d_ions_th,
    //thrust::device_ptr<const uint32_t>& d_tc_th,
    //thrust::device_ptr<const uint32_t>& d_tn_th,

    //thrust::device_ptr<const uint32_t>& d_sub_idx_th,
    //const uint32_t&                     sub_idx_length,

    //thrust::device_ptr<const uint8_t>&  d_ma_th,
    //thrust::device_ptr<const uint8_t>&  d_ma_count_th,
    //const uint32_t&                     ma_length
//)
//{

    //for (uint32_t p = 0; p < sub_idx_length; ++p) {
        //thrust::host_vector<uint32_t> ma_counts(ma_length);
        //uint32_t p_idx = d_sub_idx_th[p];

        //// count the ma's in the peptide
        //for (uint32_t a = d_tc_th[p_idx]; a < d_tn_th[p_idx]; ++a) {
            //uint8_t acid = d_ions_th[a];

            //for (uint32_t ma = 0; ma < ma_length; ++ma) {
                //if (d_ma_th[ma] == acid) {
                    //++ma_counts[ma];
                    //break;
                //}
            //}
        //}

        //bool modable = true;
        //for (uint32_t ma = 0; ma < ma_length && modable; ++ma) {
            //if (ma_counts[ma] < d_ma_count_th[ma]) {
                //modable = false;
            //}
        //}

        //if (modable) {
            //h_check_valid.push_back(p_idx);    
            //h_check_pep_ma_count.insert(h_check_pep_ma_count.end(), ma_counts.begin(), ma_counts.end());
        //}
    //}
//}

/**
 * checkFindModablePeptides
 * check the function by comparing with serial methods
 * @TODO
 */
//bool
//checkFindModablePeptides
//(
    //const uint8_t                   *d_ions,
    //const uint32_t                  *d_tc,
    //const uint32_t                  *d_tn,

    //const uint32_t                  *d_sub_idx,
    //const uint32_t                  sub_idx_length,

    //const uint8_t                   *d_ma,
    //const uint8_t                   *d_ma_count,
    //const uint32_t                  ma_length,

    //const uint32_t                  out_numValid,
    //thrust::device_ptr<uint32_t>    d_out_valid,
    //thrust::device_ptr<uint32_t>    d_out_pep_ma_count

//)
//{
    //std::cout << "checking filter by mod" << std::endl;
    
    //// initialize thrust device pointers
    //thrust::device_ptr<const uint8_t>   d_ions_th(d_ions);
    //thrust::device_ptr<const uint32_t>  d_tc_th(d_tc);
    //thrust::device_ptr<const uint32_t>  d_tn_th(d_tn);

    //thrust::device_ptr<const uint32_t>  d_sub_idx_th(d_sub_idx);
    //thrust::device_ptr<const uint8_t>   d_ma_th(d_ma);
    //thrust::device_ptr<const uint8_t>   d_ma_count_th(d_ma_count);

    //// output vectors
    //thrust::host_vector<uint32_t> h_check_valid;
    //thrust::host_vector<uint32_t> h_check_pep_ma_count;

    //filterByModNonParallel(h_check_valid, h_check_pep_ma_count, 
                           //d_ions_th, d_tc_th, d_tn_th,
                           //d_sub_idx_th, sub_idx_length,
                           //d_ma_th, d_ma_count_th, ma_length);

    //// comparison
    //thrust::device_vector<uint32_t> d_check_valid = h_check_valid;
    //thrust::device_vector<uint32_t> d_check_pep_ma_count = h_check_pep_ma_count;
    //if (h_check_valid.size() != out_numValid) {
        //std::cout << "h_check_valid.size " << h_check_valid.size() << " != " << "out_numValid " << out_numValid << std::endl;
        //exit(1);
    //}
     
    //if (!thrust::equal(d_out_valid, d_out_valid + out_numValid, d_check_valid.begin())) {
        //std::cout << "d_out_valid doesn't seem to be correct" << std::endl;
        //exit(1);

    //}

    //if (!thrust::equal(d_out_pep_ma_count, d_out_pep_ma_count + out_numValid*ma_length, d_check_pep_ma_count.begin())) {
        //std::cout << "d_out_pep_ma_count doesn't seem to be correct" << std::endl;
        //std::cout << "Printing check_valid and h_check_pep_ma_count" << std::endl;
        //for(int i = 0; i < h_check_valid.size(); i++) {
            //std::cout << "h_check_valid " << h_check_valid[i] << std::endl;
            //for (int j = 0; j < ma_length; j++) {
                //uint32_t check = h_check_pep_ma_count[i*ma_length + j];
                //uint32_t out   = d_out_pep_ma_count[i*ma_length + j];
                //if (check == out) 
                    //std::cout << check << " == " << out << " | ";
                //else
                    //std::cout << check << " != " << out << " OUT ERR" << " | ";               
                //uint32_t idx = h_check_valid[i];
                //printPeptide(idx, d_ions, d_tc, d_tn);
            //}
            //std::cout << std::endl;
        //}
        //exit(1);
    //}

        //// print d_out_valid after compaction
        ////std::cout << "Printing out_valid" << std::endl;
        ////for(int i = 0; i < out_numValid; i++) {
            ////std::cout << "d_out_valid " << d_out_valid[i] << std::endl;
            ////for (int j = 0; j < ma_length; j++) {
                ////std::cout << d_out_pep_ma_count[i*ma_length + j] << " ";
            ////}
            ////std::cout << std::endl;
        ////}
    //std::cout << "checking ok" << std::endl;
    //return true;
//}
//#endif

template <uint32_t MaxMA>
void
findModablePeptides_dispatch
(
    uint32_t            *d_valid,
    uint32_t            *d_pep_ma_count,    // 2d array, count of each ma in each peptide
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    
    const uint8_t       *d_ma,
    const uint32_t      num_ma,

    const uint32_t      *d_mod_ma_count,

    uint32_t            *d_out_pep_idx_raw, 
    uint32_t            *d_out_pep_mod_idx_raw, 

    uint32_t            num_pep_total
)
{
    uint32_t            threads;
    uint32_t            blocks;
    // control
    findByMod_control(num_pep_total, blocks, threads);

    // core
    switch (threads)
    {
    //case 512: findModablePeptides_core<512,MaxMA><<<blocks, threads>>>(d_valid, d_pep_ma_count, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;
    case 256: findModablePeptides_core<256,MaxMA><<<blocks, threads>>>(d_valid, d_pep_ma_count, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;
    case 128: findModablePeptides_core<128,MaxMA><<<blocks, threads>>>(d_valid, d_pep_ma_count, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;
    case  64: findModablePeptides_core< 64,MaxMA><<<blocks, threads>>>(d_valid, d_pep_ma_count, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;
    case  32: findModablePeptides_core< 32,MaxMA><<<blocks, threads>>>(d_valid, d_pep_ma_count, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }
}

/**
 * Binary search for begin position
 */
template <typename T>
struct findBegin : public thrust::unary_function<float,T>
{
    //thrust::device_ptr<const float> d_r;
    const float    *d_r;
    //thrust::device_ptr<const uint32_t> d_pep_idx_r_sorted;
    const uint32_t *d_pep_idx_r_sorted;
    const uint32_t num_pep;
    const float    mass;
    const float    eps;

    __host__ __device__
    findBegin(const float    *_r, 
              const uint32_t *_r_sorted,
              const uint32_t _n,
              const float    _m,
              const float    _eps) 
             : d_r(_r), d_pep_idx_r_sorted(_r_sorted), num_pep(_n), mass(_m), eps(_eps) {}

    __host__ __device__ T operator() (float delta)
    {
        const float target = mass - delta - eps;
        T min = 0;
        T max = num_pep - 1;

        while (min < max) {
            T mid = (min + max) / 2;
            const float cur = d_r[d_pep_idx_r_sorted[mid]];
            if (cur >= target) {
                max = mid;
            } else {
                min = mid + 1;    
            }
        }

        return min;
    }
};

/**
 * Binary search for end position
 */
template <typename T>
struct findEnd : public thrust::unary_function<float,T>
{
    //thrust::device_ptr<const float> d_r;
    const float    *d_r;
    const uint32_t *d_pep_idx_r_sorted;
    const uint32_t num_pep;
    const float    mass;
    const float    eps;

    __host__ __device__
    findEnd (const float     *_r, 
              const uint32_t *_r_sorted,
              const uint32_t _n,
              const float    _m,
              const float    _eps) 
             : d_r(_r), d_pep_idx_r_sorted(_r_sorted), num_pep(_n), mass(_m), eps(_eps) {}

    __host__ __device__ T operator() (float delta)
    {
        const float target = mass - delta + eps;
        T min = 0;
        T max = num_pep - 1;

        while (min < max) {
            T sum = min + max;
            T mid = (sum & 1) ? ((sum / 2) + 1) : (sum/2); // get the ceil((min+max)/2)
            const float cur = d_r[d_pep_idx_r_sorted[mid]];
            if (cur > target) {
                max = mid - 1;
            } else {
                min = mid;    
            }
        }
        if (max != 0)
            return max + 1;
        else 
            return max;
        
    }
};


/**
 * For each modification, finds peptides within range by finding the begin and end ranges to pep_idx_r_sorted
 */
uint32_t
findBeginEnd_f
(
    uint32_t            *d_begin_raw,
    uint32_t            *d_end_raw,
    //uint32_t            *d_num_pep_raw,
    uint32_t            *d_num_pep_scan_raw,

    const float         *d_r,
    const uint32_t      *d_pep_idx_r_sorted,
    const uint32_t      num_pep,

    const float         *d_mod_delta_raw,
    const uint32_t      num_mod,
    const float         mass,
    const float         eps
)
{
#ifdef _BENCH
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
#endif 

    if (num_mod == 0) return 0;
    // setup device ptrs
    thrust::device_ptr<uint32_t> d_begin(d_begin_raw);
    thrust::device_ptr<uint32_t> d_end(d_end_raw);
    //thrust::device_ptr<uint32_t> d_num_pep(d_num_pep_raw);
    thrust::device_ptr<uint32_t> d_num_pep_scan(d_num_pep_scan_raw);

    thrust::device_ptr<const float> d_mod_delta(d_mod_delta_raw);

    // Find begin and ends
    thrust::transform(d_mod_delta, d_mod_delta + num_mod, d_begin, findBegin<uint32_t>(d_r, d_pep_idx_r_sorted, num_pep, mass, eps));
    thrust::transform(d_mod_delta, d_mod_delta + num_mod, d_end, findEnd<uint32_t>(d_r, d_pep_idx_r_sorted, num_pep, mass, eps));

    // calc number of peptides in a range and scan values

    thrust::transform_exclusive_scan(thrust::make_zip_iterator(thrust::make_tuple(d_end, d_begin)),
                                     thrust::make_zip_iterator(thrust::make_tuple(d_end + num_mod, d_begin + num_mod)),
                                     d_num_pep_scan,
                                     minus_2tuple<uint32_t>(),
                                     0,
                                     thrust::plus<uint32_t>()); 
    

#ifdef _DEBUG
// check begin and ends are sane
    thrust::device_ptr<const float> d_r_th(d_r);
    thrust::device_ptr<const uint32_t> d_pep_idx_r_sorted_th(d_pep_idx_r_sorted);
    for (uint32_t i = 0; i < num_mod; ++i) {
        const float target_l = mass - d_mod_delta[i] - eps;
        const float target_u = mass - d_mod_delta[i] + eps;

        //std::cerr << "num_pep " << d_num_pep[i] << " scan " << d_num_pep_scan[i] << " beg " << d_begin[i] << " end " << d_end[i] << std::endl;
        if (d_begin[i] < d_end[i]) {
            const float beginMass = d_r_th[d_pep_idx_r_sorted_th[d_begin[i]]];
            const float lastMass = d_r_th[d_pep_idx_r_sorted_th[d_end[i] - 1]];

            //std::cout << target_l << " " << beginMass << " " << d_begin[i] << " " << d_num_pep[i] << " " << d_end[i] << " " << lastMass << " " << target_u << std::endl;
            if (!(target_l <= beginMass) ||
                !(target_u >= lastMass)) {
                std::cerr << "findBeginEnd doesn't seem to be correct (1)" << std::endl;
                exit(1);
            }
            if (d_begin[i] > 0) {
                const float beginMass_ = d_r_th[d_pep_idx_r_sorted_th[d_begin[i] - 1]];
                if (!(target_l > beginMass_)) {
                    std::cerr << "findBeginEnd doesn't seem to be correct (2)" << std::endl;
                    exit(2);
                }
            }
            if (d_end[i] < num_pep) {
                const float endMass_ = d_r_th[d_pep_idx_r_sorted_th[d_end[i]]];
                if (!(target_u < endMass_)) {
                    std::cerr << "findBeginEnd doesn't seem to be correct (3)" << std::endl;
                    exit(3);
                }
            }
        } else if (d_begin[i] > d_end[i]) {
                    std::cerr << "findBeginEnd doesn't seem to be correct (4)" << std::endl;
                    exit(4);

        }
    }
#endif

    //uint32_t total = thrust::reduce(d_num_pep, d_num_pep + num_mod);

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr << "Time elapsed for findBeginEnd_f: " << difftime(t_end,t_beg) << " %.2lf seconds\n" << std::endl;
#endif 
    // return total peptides
    uint32_t last_mod_idx = num_mod - 1;
    return d_num_pep_scan[last_mod_idx] + d_end[last_mod_idx] - d_begin[last_mod_idx];
}

/**
 * The peptides from begin to end are laid out to be search in parallel together
 * pep_idx is the idx to the original r,c,n array
 * pep_mod_idx is the idx to the modification which it can be applied to
 */
template <typename T>
struct fillPepAndModIdx : public thrust::unary_function<T, void>
{
    //thrust::device_ptr<T> d_out_pep_idx;
    T *d_out_pep_idx;
    T *d_out_pep_mod_idx;

    const T *d_pep_idx_r_sorted;
    const T *d_begin;
    const T *d_end;
    const T *d_num_pep_scan;

    __host__ __device__
    fillPepAndModIdx(T *_pi,
                     T *_pmi,
                     const T *_pirs,
                     const T *_b,
                     const T *_e,
                     const T *_nps)
                    : d_out_pep_idx(_pi), d_out_pep_mod_idx(_pmi), d_pep_idx_r_sorted(_pirs), d_begin(_b), d_end(_e), d_num_pep_scan(_nps) {}

    __host__ __device__ void operator() (T mod_idx)
    {
        T pep_idx;
        T out_idx = d_num_pep_scan[mod_idx]; // position in the out array
        for (T i = d_begin[mod_idx]; i != d_end[mod_idx]; ++i) {
            pep_idx = d_pep_idx_r_sorted[i];

            d_out_pep_idx[out_idx] = pep_idx;
            d_out_pep_mod_idx[out_idx] = mod_idx;
            ++out_idx;
        }
    }
};

/**
 * findModablePeptides
 * Check for peptides which have enough of the appropriate acids to 
 * apply a modification
 * record those that are valid and the number of each modable acid in 
 * the peptide
 */ 
/*uint32_t*/
/*findModablePeptides*/
/*(*/
    /*uint32_t            *d_out_valid_raw, // valid peptide array indices*/
    /*uint32_t            *d_out_pep_idx_raw, */
    /*uint32_t            *d_out_pep_mod_idx_raw, */
    /*uint32_t            *d_out_pep_ma_count_raw,   // 2d array, count of each ma in each peptide*/
    /*uint32_t            num_pep_total,*/

    /*const uint8_t       *d_ions,*/
    /*const uint32_t      *d_tc,*/
    /*const uint32_t      *d_tn,*/

    /*const uint32_t      *d_pep_idx_r_sorted_raw,*/

    /*const uint32_t      *d_begin_raw,*/
    /*const uint32_t      *d_end_raw,*/
    /*const uint32_t      *d_num_pep_scan_raw,*/
    /*const uint32_t      *d_mod_ma_count,*/
    /*const uint32_t      num_mod,*/

    /*const uint8_t       *d_ma,*/
    /*const uint32_t      num_ma*/
/*)*/
/*{*/
/*#ifdef _BENCH*/
    /*time_t t_beg, t_end;*/
    /*time(&t_beg);*/
/*#endif*/
    /*//std::cout << "findModablePeptides" << std::endl;*/
    /*//printGPUMemoryUsage();*/

    /*thrust::device_ptr<uint32_t> d_out_valid(d_out_valid_raw);*/

    /*// fill arrays for use later on*/
    /*thrust::counting_iterator<uint32_t> first(0);*/
    /*thrust::counting_iterator<uint32_t> last = first + num_mod;*/
    /*thrust::for_each(first, last, fillPepAndModIdx<uint32_t>(d_out_pep_idx_raw, d_out_pep_mod_idx_raw, d_pep_idx_r_sorted_raw, d_begin_raw, d_end_raw, d_num_pep_scan_raw));*/

    /*// non compacted arrays*/
    /*thrust::device_vector<uint32_t> d_valid_v(num_pep_total);*/

    /*switch (num_ma)*/
    /*{*/
    /*case 1: findModablePeptides_dispatch<1>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 2: findModablePeptides_dispatch<2>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 3: findModablePeptides_dispatch<3>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 4: findModablePeptides_dispatch<4>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 5: findModablePeptides_dispatch<5>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 6: findModablePeptides_dispatch<6>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 7: findModablePeptides_dispatch<7>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 8: findModablePeptides_dispatch<8>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 9: findModablePeptides_dispatch<9>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 10: findModablePeptides_dispatch<10>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*default:*/
        /*assert(!"Non-exhaustive patterns in match");*/
    /*}*/

    /*// compact : d_valid*/
    /*// copy if d_valid[i] > 0*/
    
    /*last = first + num_pep_total;*/
    /*thrust::device_ptr<uint32_t> d_out_valid_end =*/
        /*thrust::copy_if(first, last, d_valid_v.begin(), d_out_valid, greaterThan<const uint32_t>(0));*/

    /*const uint32_t numValid = d_out_valid_end - d_out_valid;*/

/*#ifdef _BENCH*/
    /*cudaThreadSynchronize();*/
    /*time(&t_end);*/
    /*printf ("Time elapsed for findModablePeptides: %.2lf seconds\n", difftime(t_end,t_beg));*/
/*#endif */
    /*return numValid;*/

    /*//// compact : d_pep_ma_count*/
    /*//// First expand d_valid to a 2d matrix ie same length as d_pep_ma_count*/
    /*//// eg for ma_length 3 and d_valid [1,1,0,1] to*/
    /*//// [1,1,1,*/
    /*////  1,1,1,*/
    /*////  0,0,0,*/
    /*////  1,1,1]*/
    /*//// use this array to determine whether elements in pep_ma_count should be copied*/
    /*//const uint32_t N = ma_length*sub_idx_length;*/
    /*//thrust::counting_iterator<uint32_t> first(0);*/
    /*//thrust::counting_iterator<uint32_t> last = first + N;*/
    /*//thrust::device_vector<uint32_t> d_valid_padded_th(N);*/

    /*//thrust::transform(first, last, d_valid_padded_th.begin(), fillMatrixRow<const uint32_t>(d_valid, ma_length));  */
    
    /*//thrust::device_ptr<uint32_t> d_out_pep_ma_count_end =*/
        /*//thrust::copy_if(d_pep_ma_count, d_pep_ma_count + N, d_valid_padded_th.begin(), d_out_pep_ma_count, greaterThan<const uint32_t>(0));*/

/*//#ifdef _DEBUG*/
    /*//printGPUMemoryUsage();*/

    /*//checkFindModablePeptides(d_ions, d_tc, d_tn,*/
                             /*//d_sub_idx, sub_idx_length,*/
                             /*//d_ma, d_ma_count, ma_length,*/
                             /*//numValid, d_out_valid, d_out_pep_ma_count);*/
/*//#endif*/

    /*//return numValid;*/
/*}*/

template <typename T>
struct fillPepMACount: public thrust::unary_function<T, T>
{
    const T *d_pep_idx;
    const T *d_tc;
    const T *d_tn;
    const uint8_t *d_ions;
    const uint8_t *d_ma;
    const T num_ma;  

    __host__ __device__
    fillPepMACount(const T *_pep_idx,
                   const T *_tc,
                   const T *_tn,
                   const uint8_t *_ions,
                   const uint8_t *_ma,
                   const T _num_ma)
                    : d_pep_idx(_pep_idx),
                      d_tc(_tc),
                      d_tn(_tn),
                      d_ions(_ions),
                      d_ma(_ma),
                      num_ma(_num_ma) {}

    __host__ __device__ T operator() (T i)
    {
        T idx     = i / num_ma;
        T pep_idx = d_pep_idx[idx];
        T ma_idx  = i % num_ma;
        uint8_t ma  = d_ma[ma_idx];
        T count   = 0;
        for (T a = d_tc[pep_idx] ; a < d_tn[pep_idx]; ++a) {
            uint8_t ion = d_ions[a];    
            if (ion == ma) {
                count++;
            }
        }
        return count;
    }
};

template <typename T>
struct checkModable: public thrust::unary_function<T, T>
{
    const T *d_pep_ma_count;
    const T *d_pep_mod_idx;
    const T *d_mod_ma_count;
    const T num_ma;  

    __host__ __device__
    checkModable(const T *_pep_ma_count,
                 const T *_pep_mod_idx,
                 const T *_mod_ma_count,
                 const T _num_ma)
                 : d_pep_ma_count(_pep_ma_count),
                   d_pep_mod_idx(_pep_mod_idx),
                   d_mod_ma_count(_mod_ma_count),
                   num_ma(_num_ma) {}

    __host__ __device__ T operator() (T idx)
    {
        T pep_start = num_ma * idx ;
        T mod_start = num_ma * d_pep_mod_idx[idx];
        for (T ma = 0; ma < num_ma; ++ma) {
            if (d_pep_ma_count[pep_start + ma] < d_mod_ma_count[mod_start + ma]) {
                return 0;
            }
        }
        return 1;
    }
};

// Using thrust
uint32_t
findModablePeptides
(
    uint32_t            *d_out_valid_raw, // valid peptide array indices
    uint32_t            *d_out_pep_idx_raw, 
    uint32_t            *d_out_pep_mod_idx_raw, 
    uint32_t            *d_out_pep_ma_count_raw,   // 2d array, count of each ma in each peptide
    uint32_t            num_pep_total,

    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,

    const uint32_t      *d_pep_idx_r_sorted_raw,

    const uint32_t      *d_begin_raw,
    const uint32_t      *d_end_raw,
    const uint32_t      *d_num_pep_scan_raw,
    const uint32_t      *d_mod_ma_count,
    const uint32_t      num_mod,

    const uint8_t       *d_ma,
    const uint32_t      num_ma
)
{
#ifdef _BENCH
    cudaThreadSynchronize();
    std::cerr<< "findModablePeptides" << std::endl;
    time_t t_beg, t_end;
    time(&t_beg);
    //printGPUMemoryUsage();
#endif 

    thrust::device_ptr<uint32_t> d_out_valid(d_out_valid_raw);
    thrust::device_ptr<uint32_t> d_out_pep_ma_count(d_out_pep_ma_count_raw);


    // fill arrays for use later on
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_mod;
    thrust::for_each(first, last, fillPepAndModIdx<uint32_t>(d_out_pep_idx_raw, d_out_pep_mod_idx_raw, d_pep_idx_r_sorted_raw, d_begin_raw, d_end_raw, d_num_pep_scan_raw));

    // non compacted arrays
    thrust::device_vector<uint32_t> d_valid_v(num_pep_total);
    /*thrust::device_vector<uint32_t> d_out_pep_ma_count2(num_ma*num_pep_total);*/

    // count and valid
    
    /*switch (num_ma)*/
    /*{*/
    /*case 1: findModablePeptides_dispatch<1>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 2: findModablePeptides_dispatch<2>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 3: findModablePeptides_dispatch<3>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 4: findModablePeptides_dispatch<4>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 5: findModablePeptides_dispatch<5>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 6: findModablePeptides_dispatch<6>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 7: findModablePeptides_dispatch<7>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 8: findModablePeptides_dispatch<8>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 9: findModablePeptides_dispatch<9>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*case 10: findModablePeptides_dispatch<10>(d_valid_v.data().get(), d_out_pep_ma_count_raw, d_ions, d_tc, d_tn, d_ma, num_ma, d_mod_ma_count, d_out_pep_idx_raw, d_out_pep_mod_idx_raw, num_pep_total); break;*/
    /*default:*/
        /*assert(!"Non-exhaustive patterns in match");*/
    /*}*/

    last = first + num_pep_total * num_ma;
    thrust::transform(first, last, d_out_pep_ma_count, fillPepMACount<uint32_t>(d_out_pep_idx_raw, d_tc, d_tn, d_ions, d_ma, num_ma));

    last = first + num_pep_total;
    thrust::transform(first, last, d_valid_v.begin(), checkModable<uint32_t>(d_out_pep_ma_count_raw, d_out_pep_mod_idx_raw, d_mod_ma_count, num_ma));

    // compact : d_valid
    // copy if d_valid[i] > 0
    
    last = first + num_pep_total;
    thrust::device_ptr<uint32_t> d_out_valid_end =
        thrust::copy_if(first, last, d_valid_v.begin(), d_out_valid, greaterThan<const uint32_t>(0));


    const uint32_t numValid = d_out_valid_end - d_out_valid;

    /*for (uint32_t iv = 0; iv < numValid; iv++) {*/
        /*uint32_t i = d_out_valid[iv];*/
        /*std::cout << "mpep i = " << i << " ";*/
        /*for (uint32_t ma = 0; ma < num_ma; ++ma) {*/
            /*uint32_t idx = i*num_ma + ma;*/
            /*if (d_out_pep_ma_count2[idx] != d_out_pep_ma_count[idx]) {*/
                /*std::cout << "orig " << d_out_pep_ma_count[idx] << " ";*/
                /*std::cout << "new " << d_out_pep_ma_count2[idx];*/
                /*std::cout << "DIFFFF";*/
            /*}*/
        /*}*/
        /*std::cout << std::endl;*/
    /*}  */

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr << "Time elapsed for findModablePeptides: " << difftime(t_end,t_beg) << "%.2lf seconds\n" << std::endl; 
#endif 

    return numValid;
}
