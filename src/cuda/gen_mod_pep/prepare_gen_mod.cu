/* -----------------------------------------------------------------------------
 *
 * Module    : Prepare gen mods
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 *      
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/for_each.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "algorithms.h"


template <typename T>
struct prepare : public thrust::unary_function<T, void>
{
    uint32_t          *d_out_mpep_pep_idx;
    uint32_t          *d_out_mpep_pep_mod_idx;
    uint32_t          *d_out_mpep_rank;
    uint32_t          *d_out_mpep_ith_valid;
    uint32_t          *d_out_mpep_mod_ma_count_sum;

    const uint32_t    *d_mod_ma_count_sum;

    const uint32_t    *d_pep_idx;
    const uint32_t    *d_pep_mod_idx;

    const uint32_t    *d_pep_valid_idx;
    const uint32_t    *d_pep_num_mpep;
    const uint32_t    *d_pep_num_mpep_scan;

    __host__ __device__
    prepare(uint32_t*       _ompi,
            uint32_t*       _ompmi,
            uint32_t*       _omr,
            uint32_t*       _omiv,
            uint32_t*       _ommmcs,

            const uint32_t* _mmcs,

            const uint32_t* _pi,
            const uint32_t* _pmi,
            const uint32_t* _pvi,
            const uint32_t* _pnm,
            const uint32_t* _pnms) :
            d_out_mpep_pep_idx(_ompi), d_out_mpep_pep_mod_idx(_ompmi), d_out_mpep_rank(_omr), d_out_mpep_ith_valid(_omiv), d_out_mpep_mod_ma_count_sum(_ommmcs), d_mod_ma_count_sum(_mmcs), d_pep_idx(_pi), d_pep_mod_idx(_pmi), d_pep_valid_idx(_pvi), d_pep_num_mpep(_pnm), d_pep_num_mpep_scan(_pnms) {}

    __host__ __device__ void operator() (T ith_valid)
    {
        const uint32_t valid_pep   = d_pep_valid_idx[ith_valid];
        const uint32_t pep_idx     = d_pep_idx[valid_pep];
        const uint32_t pep_mod_idx = d_pep_mod_idx[valid_pep];

        const uint32_t mmmcs = d_mod_ma_count_sum[pep_mod_idx];

        const uint32_t start = d_pep_num_mpep_scan[ith_valid];

        const uint32_t num_mpep = d_pep_num_mpep[ith_valid];
        for (uint32_t i = 0; i < num_mpep; ++i)
        {
            const uint32_t mpep_idx = start + i;
            d_out_mpep_pep_idx[mpep_idx]            = pep_idx;
            d_out_mpep_pep_mod_idx[mpep_idx]        = pep_mod_idx;
            //d_out_mpep_pep_idx[mpep_idx]            = start;
            d_out_mpep_rank[mpep_idx]               = i;
            d_out_mpep_ith_valid[mpep_idx]          = ith_valid;
            d_out_mpep_mod_ma_count_sum[mpep_idx]   = mmmcs;
        }
    }
};

uint32_t
prepareGenMod
(
    uint32_t          *d_out_mpep_pep_idx,
    uint32_t          *d_out_mpep_pep_mod_idx,
    uint32_t          *d_out_mpep_rank,
    uint32_t          *d_out_mpep_ith_valid,
    uint32_t          *d_out_mpep_mod_ma_count_sum,
    uint32_t          *d_out_mpep_mod_ma_count_sum_scan,

    const uint32_t    *d_mod_ma_count_sum,

    const uint32_t    *d_pep_idx,
    const uint32_t    *d_pep_mod_idx,

    const uint32_t    *d_pep_valid_idx,
    const uint32_t    *d_pep_num_mpep,
    const uint32_t    num_pep,
    const uint32_t    num_mpep
)
{
    thrust::device_ptr<const uint32_t> d_pep_num_mpep_th(d_pep_num_mpep);
    thrust::device_ptr<const uint32_t> d_pep_valid_idx_th(d_pep_valid_idx);

    thrust::device_vector<uint32_t> d_pep_num_mpep_scan(num_pep);

    thrust::exclusive_scan(d_pep_num_mpep_th, d_pep_num_mpep_th + num_pep, d_pep_num_mpep_scan.begin());

    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_pep;
    thrust::for_each(first, last, prepare<const uint32_t>(d_out_mpep_pep_idx, d_out_mpep_pep_mod_idx, d_out_mpep_rank, d_out_mpep_ith_valid, d_out_mpep_mod_ma_count_sum, d_mod_ma_count_sum, d_pep_idx, d_pep_mod_idx, d_pep_valid_idx, d_pep_num_mpep, d_pep_num_mpep_scan.data().get())); 

    thrust::device_ptr<uint32_t> d_out_mpep_mod_ma_count_sum_scan_th(d_out_mpep_mod_ma_count_sum_scan);
    thrust::device_ptr<uint32_t> d_out_mpep_mod_ma_count_sum_th(d_out_mpep_mod_ma_count_sum);
    thrust::exclusive_scan(d_out_mpep_mod_ma_count_sum_th, d_out_mpep_mod_ma_count_sum_th + num_mpep, d_out_mpep_mod_ma_count_sum_scan_th);


    //thrust::device_ptr<uint32_t> d_out_mpep_pep_idx_th(d_out_mpep_pep_idx);
    //thrust::device_ptr<uint32_t> d_out_mpep_rank_th(d_out_mpep_rank);
    //for (uint32_t i = 0; i < num_mpep; ++i) {
        //std::cout << "pep_idx " << d_out_mpep_pep_idx_th[i] << " rank " << d_out_mpep_rank_th[i] 
        //<< " mpep_mod_ma_count_sum " << d_out_mpep_mod_ma_count_sum_th[i] << " " << d_out_mpep_mod_ma_count_sum_scan_th[i] << std::endl;        
    //}

    const uint32_t ma_total = thrust::reduce(d_out_mpep_mod_ma_count_sum_th, d_out_mpep_mod_ma_count_sum_th + num_mpep);
    //std::cout << d_out_mpep_mod_ma_count_sum_th[num_mpep - 1] << " ma_total" << ma_total << std::endl;
    //
    //thrust::device_ptr<const uint32_t> d_mpep_mod_ma_count_sum_scan_th(d_out_mpep_mod_ma_count_sum_scan);
    //for (uint32_t i = 0; i < 8; i++) {
        //std::cout << " ma_num_comb_scan " << d_mpep_mod_ma_count_sum_scan_th[i] << std::endl;
    //}
    //std::cout << " ma_total " << ma_total << " lat count " << d_out_mpep_mod_ma_count_sum_th[num_mpep-1] << std::endl;
    return ma_total;
}

