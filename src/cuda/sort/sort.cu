/* -----------------------------------------------------------------------------
 *
 * Module    : Sort
 * Copyright : (c) [2009..2011] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/device_vector.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>
#include <stdint.h>

#include "algorithms.h"
#include "functors.h"

void sort_val_f(float *d_keys_raw, uint32_t *d_vals_raw, uint32_t N)
{
    thrust::device_ptr<float>    d_keys(d_keys_raw);
    thrust::device_ptr<uint32_t> d_vals(d_vals_raw);

    thrust::sort_by_key(d_keys, d_keys + N, d_vals);
}

void sort_idx_f(float *d_keys_raw, uint32_t *d_idx_raw, uint32_t N)
{
    thrust::device_ptr<float>    d_keys(d_keys_raw);
    thrust::device_ptr<uint32_t> d_idx(d_idx_raw);

    thrust::sequence(d_idx, d_idx + N);
    thrust::sort_by_key(d_keys, d_keys + N, d_idx);
}

void sort_rf(float *d_keys_raw, uint32_t *d_vals_raw, uint32_t N)
{
#ifdef _BENCH
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
    std::cerr << "sort_rf" << std::endl;
#endif
    thrust::device_ptr<float>    d_keys(d_keys_raw);
    thrust::device_ptr<uint32_t> d_vals(d_vals_raw);

    thrust::sort_by_key(d_keys, d_keys+N, d_vals, thrust::greater<float>());

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr<< "Time elapsed for sort_rf: " << difftime(t_end,t_beg) << " seconds" << std::endl;
#endif
}

