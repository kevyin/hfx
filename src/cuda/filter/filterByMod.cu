/* -----------------------------------------------------------------------------
 *
 * Module    : 
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
/*#include "texture.h"*/
/*#include "algorithms.h"*/

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include "algorithms.h"


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

    __shared__ volatile uint32_t s_data[BlockSize];

    for (uint32_t row = vector_id; row < sub_idx_length; row += numVectors)
    {
        const uint32_t idx       = d_sub_idx[row];
        const uint32_t row_start = d_tc[idx];
        const uint32_t row_end   = d_tn[idx];

        s_data[threadIdx.x] = 0;


        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            // check if acid is one that should be modded
            uint32_t count = 0;
            if (d_ma[0] == d_ions[j])
                count++;

            if (thread_lane == 0)
                count += s_data[threadIdx.x + (WARP_SIZE-1)];
            
            count = scan_warp<uint32_t, true>(count, s_data); 

            if (j == row_end-1) {
                d_valid[row] = count;
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

uint32_t
findModablePeptides
(
    uint32_t            *d_out_raw,

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
    printf("This is findModablePeptides\n");
    uint32_t            threads;
    uint32_t            blocks;
    uint32_t            *d_valid_raw;

    CUDA_SAFE_CALL( cudaMalloc((void**) &d_valid_raw, sub_idx_length* sizeof(uint32_t)) );

    // control
    findByMod_control(sub_idx_length, blocks, threads);

    // core
    switch (threads)
    {
    case 128: findModablePeptides_core<128><<<blocks, threads>>>(d_valid_raw, d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    case  64: findModablePeptides_core< 64><<<blocks, threads>>>(d_valid_raw, d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    case  32: findModablePeptides_core< 32><<<blocks, threads>>>(d_valid_raw, d_ions, d_tc, d_tn, d_sub_idx, sub_idx_length, d_ma, d_ma_count, ma_length); break;
    }

    // compact
    /*uint32_t N;*/
    /*N = compactIndices((uint32_t*)  d_out, d_valid, sub_idx_length);*/

    thrust::device_ptr<const uint32_t>  d_valid(d_valid_raw);
    thrust::device_ptr<uint32_t>        d_out(d_out_raw);

    // print before compaction
    thrust::host_vector<uint32_t> H(d_valid, d_valid + sub_idx_length);
    std::cout << "Printing results. no. in subset " << sub_idx_length << std::endl;
    for(int i = 0; i < sub_idx_length; i++) {
        std::cout << "d_valid " << H[i] << std::endl;
    }
   
    // copy if there are enough aa's to apply the mod
    thrust::device_ptr<const uint32_t> d_sub_idx_th(d_sub_idx); 
    thrust::device_ptr<uint32_t> d_out_end =
        thrust::copy_if(d_sub_idx_th, d_sub_idx_th + sub_idx_length, d_valid, d_out, greaterThan<const uint32_t>(0));


    cudaFree(d_valid_raw);
    /*cudaMemcpy(&h_total, (d_rtotal + num_nIdx -1), sizeof(uint32_t), cudaMemcpyDeviceToHost);*/


    // print d_out after compaction
    thrust::host_vector<uint32_t> H_compact(d_out, d_out_end);
    std::cout << "Printing results" << std::endl;
    for(int i = 0; i < d_out_end - d_out; i++) {
        std::cout << "d_out " << H_compact[i] << std::endl;
    }
    

    return 0;
}

