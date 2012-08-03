/* -----------------------------------------------------------------------------
 *
 * Module    : Filter
 * Copyright : (c) [2009..2011] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <stdint.h>

#include <thrust/transform_scan.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include "algorithms.h"
#include "functors.h"
#include "functional.hpp"

using namespace thrust;

template <typename T>
struct interval : public thrust::unary_function<T,bool>
{
    T min_val;
    T max_val;

    __host__ __device__
    interval(T _m, T _n) : min_val(_m), max_val(_n) {}

    __host__ __device__ bool operator() (T x)
    {
        return (min_val <= x && x <= max_val);
    }
};

uint32_t
findIndicesInRange_f
(
    const float         *d_in_raw,
    uint32_t            *d_indices_raw,
    const uint32_t      N,
    const float         min_val,
    const float         max_val
)
{
#ifdef _BENCH
    std::cerr << "findIndicesInRange_f" << std::endl;
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
#endif 

    thrust::device_ptr<const float>     d_in(d_in_raw);
    thrust::device_ptr<uint32_t>        d_indices(d_indices_raw);

    // define the sequence [0, N)
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + N;

    // compute indices of elements in range
    thrust::device_ptr<uint32_t> indices_end =
        thrust::copy_if(first, last, d_in, d_indices, interval<const float>(min_val, max_val));

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr<< "Time elapsed for findIndicesInRange_f : " << difftime(t_end,t_beg) << " seconds" << std::endl;
#endif 

    return indices_end - d_indices;
}

uint32_t
findSpecCandsByMass
(
    uint32_t            *d_spec_begin_raw,
    uint32_t            *d_spec_end_raw,
    uint32_t            *d_spec_num_pep_raw,

    const float         *d_r_raw,
    const uint32_t      *d_pep_idx_r_sorted_raw,
    const uint32_t      num_pep,

    const float         *d_spec_masses_raw,
    const uint32_t      num_spec, 
    const float         eps
)
{
    if (num_spec == 0) return 0;
    device_ptr<uint32_t> d_spec_begin(d_spec_begin_raw);
    device_ptr<uint32_t> d_spec_end(d_spec_end_raw);
    device_ptr<uint32_t> d_spec_num_pep(d_spec_num_pep_raw);

    device_ptr<const float> d_r(d_r_raw);
    device_ptr<const uint32_t> d_pep_idx_r_sorted(d_pep_idx_r_sorted_raw);
    device_ptr<const float> d_spec_masses(d_spec_masses_raw);

    // Finding begining and ends
    permutation_iterator<device_ptr<const float>, device_ptr<const uint32_t> > iter(d_r, d_pep_idx_r_sorted);

    lower_bound(iter,iter+num_pep,
                make_transform_iterator(d_spec_masses, 
                                        thrust::bind2nd(minus<float>(), eps)),
                make_transform_iterator(d_spec_masses+num_spec, 
                                        thrust::bind2nd(minus<float>(), eps)),
                d_spec_begin);

    upper_bound(iter,iter+num_pep,
                make_transform_iterator(d_spec_masses, 
                                        thrust::bind2nd(plus<float>(), eps)),
                make_transform_iterator(d_spec_masses+num_spec, 
                                        thrust::bind2nd(plus<float>(), eps)),
                d_spec_end);

    // Calculating number of peptides for every spec

    thrust::transform(d_spec_end, d_spec_end + num_spec, 
                      d_spec_begin,
                      d_spec_num_pep,
                      minus<uint32_t>());

    // return total peptides across all scans
    return thrust::reduce(d_spec_num_pep, d_spec_num_pep + num_spec);
    
}
