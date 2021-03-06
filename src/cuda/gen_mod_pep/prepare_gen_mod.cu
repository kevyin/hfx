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

#include <time.h>

template <typename T>
struct prepare : public thrust::unary_function<T, void>
{
    uint32_t          *d_out_mpep_pep_idx;
    uint32_t          *d_out_mpep_pep_mod_idx;
    uint32_t          *d_out_mpep_rank;
    uint32_t          *d_out_mpep_ith_cand;
    uint32_t          *d_out_mpep_mod_ma_count_sum;

    const uint32_t    *d_mod_ma_count_sum;

    const uint32_t    *d_pep_idx;
    const uint32_t    *d_pep_mod_idx;

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
            const uint32_t* _pnm,
            const uint32_t* _pnms) :
            d_out_mpep_pep_idx(_ompi), d_out_mpep_pep_mod_idx(_ompmi), d_out_mpep_rank(_omr), d_out_mpep_ith_cand(_omiv), d_out_mpep_mod_ma_count_sum(_ommmcs), d_mod_ma_count_sum(_mmcs), d_pep_idx(_pi), d_pep_mod_idx(_pmi), d_pep_num_mpep(_pnm), d_pep_num_mpep_scan(_pnms) {}

    __host__ __device__ void operator() (T idx)
    {
        const uint32_t cand_idx   = idx;
        const uint32_t pep_idx     = d_pep_idx[cand_idx];
        const uint32_t pep_mod_idx = d_pep_mod_idx[cand_idx];

        const uint32_t mmmcs = d_mod_ma_count_sum[pep_mod_idx];

        const uint32_t start = d_pep_num_mpep_scan[idx];

        const uint32_t num_mpep = d_pep_num_mpep[idx];
        for (uint32_t i = 0; i < num_mpep; ++i)
        {
            const uint32_t mpep_idx = start + i;
            d_out_mpep_pep_idx[mpep_idx]            = pep_idx;
            d_out_mpep_pep_mod_idx[mpep_idx]        = pep_mod_idx;
            //d_out_mpep_pep_idx[mpep_idx]            = start;
            d_out_mpep_rank[mpep_idx]               = i;
            d_out_mpep_ith_cand[mpep_idx]          = idx;
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
    uint32_t          *d_out_mpep_ith_cand,
    uint32_t          *d_out_mpep_mod_ma_count_sum_scan,

    const uint32_t    *d_mod_ma_count_sum,

    const uint32_t    *d_pep_idx,
    const uint32_t    *d_pep_mod_idx,

    const uint32_t    *d_pep_num_mpep,
    const uint32_t    num_cand,
    const uint32_t    num_mpep
)
{
#ifdef _BENCH
    std::cout << "prepareGenMod" << std::endl;
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
#endif 

    //printGPUMemoryUsage();

    thrust::device_ptr<const uint32_t> d_pep_num_mpep_th(d_pep_num_mpep);
    thrust::device_ptr<uint32_t> d_out_mpep_mod_ma_count_sum_scan_th(d_out_mpep_mod_ma_count_sum_scan);

    thrust::device_vector<uint32_t> d_pep_num_mpep_scan(num_cand);
    thrust::device_vector<uint32_t> d_mpep_mod_ma_count_sum(num_mpep);

    thrust::exclusive_scan(d_pep_num_mpep_th, d_pep_num_mpep_th + num_cand, d_pep_num_mpep_scan.begin());

    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_cand;
    thrust::for_each(first, last, prepare<const uint32_t>(d_out_mpep_pep_idx, d_out_mpep_pep_mod_idx, d_out_mpep_rank, d_out_mpep_ith_cand, d_mpep_mod_ma_count_sum.data().get(), d_mod_ma_count_sum, d_pep_idx, d_pep_mod_idx, d_pep_num_mpep, d_pep_num_mpep_scan.data().get())); 

    thrust::exclusive_scan(d_mpep_mod_ma_count_sum.begin(), d_mpep_mod_ma_count_sum.end(), d_out_mpep_mod_ma_count_sum_scan_th);

    const uint32_t ma_total = d_out_mpep_mod_ma_count_sum_scan_th[num_mpep-1] + d_mpep_mod_ma_count_sum[num_mpep-1];

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    printf ("Time elapsed for prepareGenMod: %.2lf seconds\n", difftime(t_end,t_beg));
#endif 
    return ma_total;
}

