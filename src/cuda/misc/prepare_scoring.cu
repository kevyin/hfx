/* -----------------------------------------------------------------------------
 *
 * Module    : Prepare scoring
 * Copyright : (c) [2012] Kevin Ying 
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/copy.h>

#include "algorithms.h"

using namespace thrust;


template <typename T>
struct fill_indices
{
    T       *d_out;
    const T *d_values;

    __host__ __device__ 
    fill_indices(T       *_o,
                 const T *_v) : 
                 d_out(_o), d_values(_v) {}

    template <typename Tuple>
    __host__ __device__ 
    void operator() (Tuple x)
    {
        uint32_t idx  = thrust::get<0>(x);
        uint32_t n    = thrust::get<1>(x);
        T        *out = d_out + thrust::get<2>(x);

        for (size_t i = 0; i < n; ++i)
        {
            out[i] = d_values[idx]; 
            idx++;
        }
    }
};

void prepare_scoring
(
    uint32_t        *d_out_spec_pep_idx,

    const uint32_t  *d_pep_idx_r_sorted,
    const uint32_t  *d_spec_begin_raw,
    const uint32_t  *d_spec_num_pep_raw,
    const uint32_t  num_spec
)
{

    device_ptr<const uint32_t> d_spec_begin(d_spec_begin_raw);
    device_ptr<const uint32_t> d_spec_num_pep(d_spec_num_pep_raw);
    device_vector<uint32_t> d_spec_num_pep_scan(num_spec);

    exclusive_scan(d_spec_num_pep, d_spec_num_pep + num_spec, d_spec_num_pep_scan.begin());

    for_each_n(make_zip_iterator(make_tuple(d_spec_begin, d_spec_num_pep, d_spec_num_pep_scan.begin())),
               num_spec,
               fill_indices<uint32_t>(d_out_spec_pep_idx, d_pep_idx_r_sorted));
}
