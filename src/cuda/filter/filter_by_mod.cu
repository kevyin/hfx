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
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/scan.h>

#include <time.h>
#include <stdio.h>

using namespace thrust;

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
    device_ptr<const float> d_r_th(d_r);
    device_ptr<const uint32_t> d_pep_idx_r_sorted_th(d_pep_idx_r_sorted);

    // Find begin and ends

    permutation_iterator<device_ptr<const float>, device_ptr<const uint32_t> > iter(d_r_th, d_pep_idx_r_sorted_th);
    lower_bound(iter,iter+num_pep,
                make_transform_iterator(d_mod_delta, thrust::bind1st(minus<float>(), mass-eps)),
                make_transform_iterator(d_mod_delta+num_mod, thrust::bind1st(minus<float>(), mass-eps)),
                d_begin);

    upper_bound(iter,iter+num_pep,
                make_transform_iterator(d_mod_delta, thrust::bind1st(minus<float>(), mass+eps)),
                make_transform_iterator(d_mod_delta+num_mod, thrust::bind1st(minus<float>(), mass+eps)),
                d_end);

    // calc number of peptides in a range and scan values

    thrust::transform_exclusive_scan(
        thrust::make_zip_iterator(thrust::make_tuple(d_end, d_begin)),
        thrust::make_zip_iterator(thrust::make_tuple(d_end + num_mod, d_begin + num_mod)),
        d_num_pep_scan,
        minus_2tuple<uint32_t>(),
        0,
        thrust::plus<uint32_t>()); 
    

#ifdef _DEBUG
// check begin and ends are sane
    thrust::device_ptr<const float> d_r_th(d_r);
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
                return false;
            }
        }
        return true;
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
    uint32_t            num_cand_total,

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


    // fill arrays mod_idx and pep_idx. will look like:
    // mod_idx  pep_idx
    // 0        a0
    // 0        .
    // 0        .
    // 0        b0
    // 1        a1
    // 1        .
    // 1        .
    // 1        b1
    // Where a{n} .. b{n} are indices to the peptides which are inrange for the modificiation 
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_mod;
    thrust::for_each(first, last, fillPepAndModIdx<uint32_t>(d_out_pep_idx_raw, d_out_pep_mod_idx_raw, d_pep_idx_r_sorted_raw, d_begin_raw, d_end_raw, d_num_pep_scan_raw));

    // non compacted arrays
    thrust::device_vector<bool> d_valid_v(num_cand_total);
    /*thrust::device_vector<uint32_t> d_out_pep_ma_count2(num_ma*num_cand_total);*/

    // count and valid
    
    last = first + num_cand_total * num_ma;
    thrust::transform(first, last, d_out_pep_ma_count, fillPepMACount<uint32_t>(d_out_pep_idx_raw, d_tc, d_tn, d_ions, d_ma, num_ma));

    last = first + num_cand_total;
    thrust::transform(first, last, d_valid_v.begin(), checkModable<uint32_t>(d_out_pep_ma_count_raw, d_out_pep_mod_idx_raw, d_mod_ma_count, num_ma));

    // compact
    //device_ptr<uint32_t> d_pep_mod_idx(d_out_pep_mod_idx_raw);
    //device_ptr<uint32_t> d_pep_idx(d_out_pep_idx_raw);
    //device_ptr<uint32_t> d_pep_ma_count(d_out_pep_ma_count_raw);


    //thrust::remove_if(make_zip_iterator(make_tuple(d_pep_mod_idx, d_pep_idx),
                                        //make_tuple(d_pep_mod_idx+num_cand_total, d_pep_idx+num_cand_total),
                                        //d_valid_v.begin(),
                                        //logical_not<bool>());
    //thrust::remove_if(d_pep_ma_count, d_pep_ma_count+num_cand_total*num_ma,
                      //d_valid_v.begin(),
                      //logical_not<bool>());

    // compact : d_valid
    // copy if d_valid[i] > 0
    
    last = first + num_cand_total;
    thrust::device_ptr<uint32_t> d_out_valid_end =
        thrust::copy_if(first, last, d_valid_v.begin(), d_out_valid, greaterThan<const uint32_t>(0));


    const uint32_t numValid = d_out_valid_end - d_out_valid;

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr << "Time elapsed for findModablePeptides: " << difftime(t_end,t_beg) << "%.2lf seconds\n" << std::endl; 
#endif 

    return numValid;
}
