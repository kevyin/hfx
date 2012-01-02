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

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>

#define DEBUG

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

template <uint32_t BlockSize>   
__global__ static void
findModablePeptides_core
(
    uint32_t            *d_valid,
    uint32_t            *d_pep_ma_count,    // 2d array, count of each ma in each peptide

    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,

    const uint32_t      *d_sub_idx,
    const uint32_t      sub_idx_length,

    const uint8_t       *d_ma,
    const uint8_t       *d_ma_count,
    const uint32_t      ma_length
)
{

    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t vectorsPerBlock = BlockSize / WARP_SIZE;
    const uint32_t numVectors      = vectorsPerBlock * gridDim.x;
    const uint32_t thread_id       = BlockSize * blockIdx.x + threadIdx.x;
    const uint32_t vector_id       = thread_id / WARP_SIZE;
    const uint32_t thread_lane     = threadIdx.x & (WARP_SIZE-1);

    assert(MAX_MA >= ma_length);

    // Keep a count for each mod
    __shared__ volatile uint32_t s_data[MAX_MA][BlockSize];

    for (uint32_t row = vector_id; row < sub_idx_length; row += numVectors)
    {
        const uint32_t idx       = d_sub_idx[row];
        const uint32_t row_start = d_tc[idx];
        const uint32_t row_end   = d_tn[idx];

        for (int mod = 0; mod < ma_length; mod++)
        {
            s_data[mod][threadIdx.x] = 0;
        }

        // for each acid check if this is a modable acid, if so record the count
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            bool isModable = true; // assume modability
            for (int mod = 0; mod < ma_length; mod++) 
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
                    if (count >= d_ma_count[mod]) 
                    {
                        //isModable = true;
                        // record pep aa counts
                        d_pep_ma_count[row*ma_length + mod] = count;
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

/*
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

template <typename T>
struct greaterThan : public thrust::unary_function<T,bool>
{
    T bound;

    __host__ __device__
    greaterThan(T _m) : bound(_m) {}

    __host__ __device__ bool operator() (T x)
    {
        return (bound < x);
    }
};

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

#ifdef DEBUG

/**
 * filterByModNonParallel
 */
void
filterByModNonParallel 
(
    thrust::host_vector<uint32_t>& h_check_valid,
    thrust::host_vector<uint32_t>& h_check_pep_ma_count,

    thrust::device_ptr<const uint8_t>&  d_ions_th,
    thrust::device_ptr<const uint32_t>& d_tc_th,
    thrust::device_ptr<const uint32_t>& d_tn_th,

    thrust::device_ptr<const uint32_t>& d_sub_idx_th,
    const uint32_t&                     sub_idx_length,

    thrust::device_ptr<const uint8_t>&  d_ma_th,
    thrust::device_ptr<const uint8_t>&  d_ma_count_th,
    const uint32_t&                     ma_length
)
{

    for (uint32_t p = 0; p < sub_idx_length; ++p) {
        thrust::host_vector<uint32_t> ma_counts(ma_length);
        uint32_t p_idx = d_sub_idx_th[p];

        // count the ma's in the peptide
        for (uint32_t a = d_tc_th[p_idx]; a < d_tn_th[p_idx]; ++a) {
            uint8_t acid = d_ions_th[a];

            for (uint32_t ma = 0; ma < ma_length; ++ma) {
                if (d_ma_th[ma] == acid) {
                    ++ma_counts[ma];
                    break;
                }
            }
        }

        bool modable = true;
        for (uint32_t ma = 0; ma < ma_length && modable; ++ma) {
            if (ma_counts[ma] < d_ma_count_th[ma]) {
                modable = false;
            }
        }

        if (modable) {
            h_check_valid.push_back(p_idx);    
            h_check_pep_ma_count.insert(h_check_pep_ma_count.end(), ma_counts.begin(), ma_counts.end());
        }
    }
}

/**
 * checkFindModablePeptides
 * check the function by comparing with serial methods
 * @TODO
 */
bool
checkFindModablePeptides
(
    const uint8_t                   *d_ions,
    const uint32_t                  *d_tc,
    const uint32_t                  *d_tn,

    const uint32_t                  *d_sub_idx,
    const uint32_t                  sub_idx_length,

    const uint8_t                   *d_ma,
    const uint8_t                   *d_ma_count,
    const uint32_t                  ma_length,

    const uint32_t                  out_numValid,
    thrust::device_ptr<uint32_t>    d_out_valid,
    thrust::device_ptr<uint32_t>    d_out_pep_ma_count

)
{
    std::cout << "checking filter by mod" << std::endl;
    
    // initialize thrust device pointers
    thrust::device_ptr<const uint8_t>   d_ions_th(d_ions);
    thrust::device_ptr<const uint32_t>  d_tc_th(d_tc);
    thrust::device_ptr<const uint32_t>  d_tn_th(d_tn);

    thrust::device_ptr<const uint32_t>  d_sub_idx_th(d_sub_idx);
    thrust::device_ptr<const uint8_t>   d_ma_th(d_ma);
    thrust::device_ptr<const uint8_t>   d_ma_count_th(d_ma_count);

    // output vectors
    thrust::host_vector<uint32_t> h_check_valid;
    thrust::host_vector<uint32_t> h_check_pep_ma_count;

    filterByModNonParallel(h_check_valid, h_check_pep_ma_count, 
                           d_ions_th, d_tc_th, d_tn_th,
                           d_sub_idx_th, sub_idx_length,
                           d_ma_th, d_ma_count_th, ma_length);

    // comparison
    thrust::device_vector<uint32_t> d_check_valid = h_check_valid;
    thrust::device_vector<uint32_t> d_check_pep_ma_count = h_check_pep_ma_count;
    if (h_check_valid.size() != out_numValid) {
        std::cout << "h_check_valid.size " << h_check_valid.size() << " != " << "out_numValid " << out_numValid << std::endl;
        exit(1);
    }
     
    if (!thrust::equal(d_out_valid, d_out_valid + out_numValid, d_check_valid.begin())) {
        std::cout << "d_out_valid doesn't seem to be correct" << std::endl;
        exit(1);

    }

    if (!thrust::equal(d_out_pep_ma_count, d_out_pep_ma_count + out_numValid*ma_length, d_check_pep_ma_count.begin())) {
        std::cout << "d_out_pep_ma_count doesn't seem to be correct" << std::endl;
        std::cout << "Printing check_valid and h_check_pep_ma_count" << std::endl;
        for(int i = 0; i < h_check_valid.size(); i++) {
            std::cout << "h_check_valid " << h_check_valid[i] << std::endl;
            for (int j = 0; j < ma_length; j++) {
                uint32_t check = h_check_pep_ma_count[i*ma_length + j];
                uint32_t out   = d_out_pep_ma_count[i*ma_length + j];
                if (check == out) 
                    std::cout << check << " == " << out << " | ";
                else
                    std::cout << check << " != " << out << " OUT ERR" << " | ";               
                uint32_t idx = h_check_valid[i];
                printPeptide(idx, d_ions, d_tc, d_tn);
            }
            std::cout << std::endl;
        }
        exit(1);
    }

        // print d_out_valid after compaction
        //std::cout << "Printing out_valid" << std::endl;
        //for(int i = 0; i < out_numValid; i++) {
            //std::cout << "d_out_valid " << d_out_valid[i] << std::endl;
            //for (int j = 0; j < ma_length; j++) {
                //std::cout << d_out_pep_ma_count[i*ma_length + j] << " ";
            //}
            //std::cout << std::endl;
        //}
    std::cout << "checking ok" << std::endl;
    return true;
}
#endif

/**
 * findModablePeptides
 * Check for peptides which have enough of the appropriate acids to 
 * apply a modification
 * record those that are valid and the number of each modable acid in 
 * the peptide
 */ 
uint32_t
findModablePeptides
(
    uint32_t            *d_out_valid_raw,           // valid indices from sub_idx
    uint32_t            *d_out_pep_ma_count_raw,    // 2d array, count of each ma in each peptide

    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,

    const uint32_t      *d_sub_idx,
    const uint32_t      sub_idx_length,

    const uint8_t       *d_ma,
    const uint8_t       *d_ma_count,
    const uint32_t      ma_length
)
{
    uint32_t            threads;
    uint32_t            blocks;

    // non compacted arrays
    thrust::device_vector<uint32_t> d_valid_v(sub_idx_length);
    thrust::device_vector<uint32_t> d_pep_ma_count_v(sub_idx_length*ma_length);

    // control
    findByMod_control(sub_idx_length, blocks, threads);

    // core
    switch (threads)
    {
    case 128: findModablePeptides_core<128><<<blocks, threads>>>(d_valid_v.data().get(), d_pep_ma_count_v.data().get(), d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    case  64: findModablePeptides_core< 64><<<blocks, threads>>>(d_valid_v.data().get(), d_pep_ma_count_v.data().get(), d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    case  32: findModablePeptides_core< 32><<<blocks, threads>>>(d_valid_v.data().get(), d_pep_ma_count_v.data().get(), d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    }

    // compact
    
    // compact : prepare thrust ptrs 
    thrust::device_ptr<const uint32_t>  d_valid(d_valid_v.data().get());
    thrust::device_ptr<uint32_t>        d_out_valid(d_out_valid_raw);

    thrust::device_ptr<const uint32_t>  d_pep_ma_count(d_pep_ma_count_v.data().get());
    thrust::device_ptr<uint32_t>        d_out_pep_ma_count(d_out_pep_ma_count_raw);
   
    // compact : d_valid
    // copy into d_out_ if d_valid[i] > 0
    thrust::device_ptr<const uint32_t> d_sub_idx_th(d_sub_idx); 
    thrust::device_ptr<uint32_t> d_out_valid_end =
        thrust::copy_if(d_sub_idx_th, d_sub_idx_th + sub_idx_length, d_valid, d_out_valid, greaterThan<const uint32_t>(0));

    const uint32_t numValid = d_out_valid_end - d_out_valid;

    // compact : d_pep_ma_count
    // First expand d_valid to a 2d matrix ie same length as d_pep_ma_count
    // eg for ma_length 3 and d_valid [1,1,0,1] to
    // [1,1,1,
    //  1,1,1,
    //  0,0,0,
    //  1,1,1]
    // use this array to determine whether elements in pep_ma_count should be copied
    const uint32_t N = ma_length*sub_idx_length;
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + N;
    thrust::device_vector<uint32_t> d_valid_padded_th(N);

    thrust::transform(first, last, d_valid_padded_th.begin(), fillMatrixRow<const uint32_t>(d_valid, ma_length));  
    
    thrust::device_ptr<uint32_t> d_out_pep_ma_count_end =
        thrust::copy_if(d_pep_ma_count, d_pep_ma_count + N, d_valid_padded_th.begin(), d_out_pep_ma_count, greaterThan<const uint32_t>(0));

#ifdef DEBUG
    printGPUMemoryUsage();

    checkFindModablePeptides(d_ions, d_tc, d_tn,
                             d_sub_idx, sub_idx_length,
                             d_ma, d_ma_count, ma_length,
                             numValid, d_out_valid, d_out_pep_ma_count);
#endif

    return numValid;
}

