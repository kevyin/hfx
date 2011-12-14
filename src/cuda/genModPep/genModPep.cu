/* -----------------------------------------------------------------------------
 *
 * Module    : Generate Modified Peptides
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <stdint.h>

#include "utils.h"
#include "algorithms.h"
#include "combinatorics.h"

#include <stdlib.h>
#include <stdio.h>

template <typename T>
struct calcNumModCand : public thrust::unary_function<T, uint32_t>
{

    thrust::device_ptr<const uint32_t>    d_pep_ma_count;
    thrust::device_ptr<const uint8_t>     d_mod_ma_count;
    const uint32_t                        mod_num_ma;

    __host__ __device__
    calcNumModCand(thrust::device_ptr<const uint32_t>    _pep_ma_c, 
                   thrust::device_ptr<const uint8_t>     _ma_c, 
                   const uint32_t                        len) : 
                   d_pep_ma_count(_pep_ma_c), d_mod_ma_count(_ma_c), mod_num_ma(len) {}

    __host__ __device__ uint32_t operator() (T pep)
    {
        thrust::device_ptr<const uint32_t> pep_ma_count (d_pep_ma_count + pep*mod_num_ma);
        uint32_t N = 1;
        for (int i = 0; i < mod_num_ma; i++)
        {
            N *= choose(pep_ma_count[i], d_mod_ma_count[i]);    
        }
        return N;
    }
};

template <typename T>
struct calcNumCombPerAcid: public thrust::unary_function<T, uint32_t>
{

    thrust::device_ptr<const uint32_t>    d_pep_ma_count;
    thrust::device_ptr<const uint8_t>     d_mod_ma_count;
    const uint32_t                        mod_num_ma;

    __host__ __device__
    calcNumCombPerAcid(thrust::device_ptr<const uint32_t>    _pep_ma_c, 
                   thrust::device_ptr<const uint8_t>     _mod_ma_c, 
                   const uint32_t                        _mod_num_ma) : 
                   d_pep_ma_count(_pep_ma_c), d_mod_ma_count(_mod_ma_c), mod_num_ma(_mod_num_ma) {}

    __host__ __device__ uint32_t operator() (T pep_ma_num_comb_idx)
    {
        uint32_t pep_ma_count_row_idx = pep_ma_num_comb_idx / mod_num_ma;
        uint32_t mod_ma_count_idx     = pep_ma_num_comb_idx % mod_num_ma;

        thrust::device_ptr<const uint32_t> pep_ma_count_row (d_pep_ma_count + mod_num_ma*pep_ma_count_row_idx);
        return choose(pep_ma_count_row[mod_ma_count_idx], d_mod_ma_count[mod_ma_count_idx]);    
    }
};

template <typename T>
struct calcNumMPepPerPep: public thrust::unary_function<T, uint32_t>
{

    thrust::device_ptr<const uint32_t>  d_pep_ma_num_comb;
    thrust::device_ptr<uint32_t>        d_pep_ma_num_comb_scan;
    const uint32_t                      mod_num_ma;

    __host__ __device__
    calcNumMPepPerPep(thrust::device_ptr<const uint32_t>    _pep_ma_num_comb, 
                      thrust::device_ptr<uint32_t>          _pep_ma_num_comb_scan,
                      const uint32_t                        _mod_num_ma) : 
                      d_pep_ma_num_comb(_pep_ma_num_comb), 
                      d_pep_ma_num_comb_scan(_pep_ma_num_comb_scan), 
                      mod_num_ma(_mod_num_ma) {}

    __host__ __device__ uint32_t operator() (T pep)
    {
        uint32_t pep_ma_num_comb_idx = pep * mod_num_ma;

        thrust::device_ptr<const uint32_t> pep_ma_num_comb_row(d_pep_ma_num_comb + pep_ma_num_comb_idx);
        thrust::device_ptr<uint32_t>       pep_ma_num_comb_scan_row(d_pep_ma_num_comb_scan + pep_ma_num_comb_idx);
        uint32_t num = 1;
        //for (uint32_t i = 0; i < mod_num_ma; ++i) 
        for (int32_t i = mod_num_ma - 1; i >= 0; --i) 
        {
            num *= pep_ma_num_comb_row[i];
            // raw
            uint32_t *tmp_raw = (pep_ma_num_comb_scan_row + i).get();
            *tmp_raw = num;
            //pep_ma_num_comb_scan_row[i] = num;
        }
        return num;
    }
};

uint32_t
calcTotalModCands
(
    uint32_t          *d_out_pep_num_mpep,
    uint32_t          *d_out_pep_ma_num_comb,
    uint32_t          *d_out_pep_ma_num_comb_scan,
    const uint32_t    nPep,                   // number of peptides

    const uint32_t    *d_pep_ma_count,        // 2d array of ma count in each pep
    const uint8_t     *d_mod_ma_count,             
    const uint32_t    mod_num_ma
) 
{
    thrust::device_ptr<uint32_t>        d_out_pep_num_mpep_th(d_out_pep_num_mpep);
    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_th(d_out_pep_ma_num_comb); 
    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_scan_th(d_out_pep_ma_num_comb_scan);
    thrust::device_ptr<const uint32_t>  d_pep_ma_count_th(d_pep_ma_count);
    thrust::device_ptr<const uint8_t>   d_mod_ma_count_th(d_mod_ma_count);
    
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + nPep*mod_num_ma;
    thrust::transform(first,last, d_out_pep_ma_num_comb_th, calcNumCombPerAcid<const uint32_t>(d_pep_ma_count_th, d_mod_ma_count_th, mod_num_ma));

    //thrust::device_vector<uint32_t> d_out_pep_ma_num_comb_scan_th(nPep*mod_num_ma);
    //thrust::counting_iterator<uint32_t> first(0);
    last = first + nPep;
    thrust::transform(first, last, d_out_pep_num_mpep_th, calcNumMPepPerPep<const uint32_t>(d_out_pep_ma_num_comb_th, d_out_pep_ma_num_comb_scan_th, mod_num_ma));

    ////print scans
    //std::cout << "ma_num_comb not scanned and scanned" << std::endl;
    //for (int i = 0; i < nPep; i++)               
    //{
        //std::cout << "pep ma num_comb and scan" << i << " ";
        //for (int j = 0; j < mod_num_ma; j++)               
        //{
            //int k = i * mod_num_ma + j;
            //uint32_t q = d_mod_ma_count_th[j];
            //std::cout << "pep_ma_num_comb " << d_out_pep_ma_num_comb_th[k] << " ";
            //std::cout << "scan " << d_out_pep_ma_num_comb_scan_th[k] << " | ";
        //}
        //std::cout << std::endl;
    //}


    // compare with calcNumModCand method
    // malloc d_tmp
    thrust::device_ptr<uint32_t> d_tmp = thrust::device_malloc<uint32_t>(nPep);
    thrust::transform(first, last, d_tmp, calcNumModCand<const uint32_t>(d_pep_ma_count_th, d_mod_ma_count_th, mod_num_ma));

    // check tmp is same as d_pep_ma_count
    if (not thrust::equal(d_out_pep_num_mpep_th, d_out_pep_num_mpep_th + nPep, d_tmp))
    {
        std::cerr << "pep_num_mpep doesn't seem to be correct" << std::endl;
        exit(1);
    }

    thrust::device_free(d_tmp);



    ////
    ////print
    //std::cout << "pep_mpep_count " << std::endl;
    //for (int i = 0; i < nPep; i++)               
    //{
        //std::cout << d_out_pep_num_mpep_th[i] << std::endl;
    //}

    uint32_t total = thrust::reduce(d_out_pep_num_mpep_th, d_out_pep_num_mpep_th + nPep);
    //std::cout << "total mpeps " << total << std::endl;

    return total;
    //return 0;
}

__host__ __device__ void
add_ma_ranks 
(
    uint32_t                            *d_out_ma_ranks_row,
    const uint32_t                      mpep_rank,
    thrust::device_ptr<const uint32_t>  d_pep_ma_num_comb_scan_row,
    const uint32_t                      mod_num_ma
)
{
    //thrust::device_vector<uint32_t> ma_ranks(mod_num_ma);
    
    uint32_t *ma_ranks_raw = d_out_ma_ranks_row;

    uint32_t x = mpep_rank;

    for (uint32_t i = 0; i < mod_num_ma; ++i)
    {
        if (i != mod_num_ma - 1)
        {
            uint32_t next_pep_ma_num_comb_scan = d_pep_ma_num_comb_scan_row[i + 1];
            if (x >= next_pep_ma_num_comb_scan)
            {
                ma_ranks_raw[i] = x / next_pep_ma_num_comb_scan;
                x = x % next_pep_ma_num_comb_scan;
            }

        } else {
            ma_ranks_raw[i] = x;
        }
    }
}

template <typename T>
struct add_mpep_unrank: public thrust::unary_function<T, uint32_t>
{
    thrust::device_ptr<uint32_t>          d_out_mpep_unrank;

    //thrust::device_ptr<const uint32_t>    d_pep_ma_num_comb;
    thrust::device_ptr<const uint32_t>    d_pep_ma_num_comb_scan;
    thrust::device_ptr<const uint32_t>    d_mpep_rank;
    thrust::device_ptr<const uint32_t>    d_mpep_rank_ith_pep;

    thrust::device_ptr<const uint32_t>    d_pep_ma_count;

    thrust::device_ptr<const uint8_t>     d_mod_ma_count;
    const uint32_t                        mod_num_ma;

    const uint32_t                        sum_mod_ma_count;


    __host__ __device__
    add_mpep_unrank (thrust::device_ptr<uint32_t>       _out_mpep_unrank, 

                     thrust::device_ptr<const uint32_t> _pep_ma_num_comb_scan,
                     thrust::device_ptr<const uint32_t> _mpep_rank,
                     thrust::device_ptr<const uint32_t> _mpep_rank_ith_pep,
                     thrust::device_ptr<const uint32_t> _pep_ma_count,
                     thrust::device_ptr<const uint8_t>  _mod_ma_count,
                     const uint32_t                     _mod_num_ma,
                     const uint32_t                     _sum_mod_ma_count) : 
                     d_out_mpep_unrank      (_out_mpep_unrank), 
                     d_pep_ma_num_comb_scan (_pep_ma_num_comb_scan),
                     d_mpep_rank            (_mpep_rank),
                     d_mpep_rank_ith_pep    (_mpep_rank_ith_pep),
                     d_pep_ma_count         (_pep_ma_count),
                     d_mod_ma_count         (_mod_ma_count),
                     mod_num_ma             (_mod_num_ma),
                     sum_mod_ma_count       (_sum_mod_ma_count) {}

    __host__ __device__ uint32_t operator() (T mpep)
    {
        uint32_t mpep_rank = d_mpep_rank[mpep];
        uint32_t ith_pep  = d_mpep_rank_ith_pep[mpep];
        // mpep_rank -> ma_ranks
        //thrust::device_vector<uint32_t> d_ma_ranks(mod_num_ma);
        uint32_t d_ma_ranks[MAX_MA];
        add_ma_ranks(d_ma_ranks, mpep_rank, d_pep_ma_num_comb_scan + d_mpep_rank_ith_pep[mpep]*mod_num_ma, mod_num_ma);
        //add_ma_ranks(d_out_mpep_unrank + mpep*sum_mod_ma_count, mpep_rank, d_pep_ma_num_comb_scan + ith_pep*mod_num_ma, mod_num_ma);

        uint32_t unrank_pos = 0;
        
        // foreach ma, unrank and add to mpep_unrank
        //add_mpep_unrank(d_out_mpep_unrank + mpep*sum_mod_ma_count);
        for (uint32_t i = 0; i < mod_num_ma; ++i)
        {
            unrankComb((d_out_mpep_unrank + mpep*sum_mod_ma_count + unrank_pos).get(), 
                       d_pep_ma_count[ith_pep*mod_num_ma + i], 
                       d_mod_ma_count[i], 
                       d_ma_ranks[i]);
            unrank_pos += d_mod_ma_count[i];
        }
        
        return mpep;
        //return num;
    }
};

void
genModCands
(                                                                     
    uint32_t        *d_out_mpep_idx,
    uint32_t        *d_out_mpep_rank,
    uint32_t        *d_out_mpep_unrank,
    const uint32_t  total,

    const uint32_t  *d_pep_idx,
    const uint32_t  *d_pep_num_mpep,
    const uint32_t  *d_pep_ma_num_comb,
    const uint32_t  *d_pep_ma_num_comb_scan,
    const uint32_t  num_pep,
   
    const uint32_t  *d_pep_ma_count,        // 2d array of ma count in each pep

    const uint8_t   *d_mod_ma_count,             
    const uint32_t  mod_num_ma
)
{
    std::cout << "genModCands" << std::endl;
    thrust::device_ptr<uint32_t> d_out_mpep_idx_th(d_out_mpep_idx);
    thrust::device_ptr<uint32_t> d_out_mpep_rank_th(d_out_mpep_rank);
    thrust::device_ptr<uint32_t> d_out_mpep_unrank_th(d_out_mpep_unrank);
   
    thrust::device_ptr<const uint32_t> d_pep_idx_th(d_pep_idx);
    thrust::device_ptr<const uint32_t> d_pep_num_mpep_th(d_pep_num_mpep);
    thrust::device_ptr<const uint32_t> d_pep_ma_num_comb_th(d_pep_ma_num_comb);
    thrust::device_ptr<const uint32_t> d_pep_ma_num_comb_scan_th(d_pep_ma_num_comb_scan);

    thrust::device_ptr<const uint32_t> d_pep_ma_count_th(d_pep_ma_count);
    thrust::device_ptr<const uint8_t>  d_mod_ma_count_th(d_mod_ma_count);



    std::vector<uint32_t> mpep_rank;
    std::vector<uint32_t> mpep_rank_ith_pep; // helper vector, corres pep for each mpep_rank
    std::vector<uint32_t> mpep_idx;
    for (uint32_t i = 0; i < num_pep; i++)
    {
        for (uint32_t k = 0; k < d_pep_num_mpep_th[i]; k++)
        {
            // add mpep_rank
            mpep_rank.push_back(k);
            // add mpep_rank
            mpep_rank_ith_pep.push_back(i);
            // add idx
            mpep_idx.push_back(d_pep_idx_th[i]);
        }
    }
    if (mpep_rank.size() != total) {
        std::cerr << "mpep_rankvector not the correct length" << std::endl;
        exit(1);
    }

    // copy to device idx
    thrust::copy(mpep_idx.begin(), mpep_idx.end(), d_out_mpep_idx_th);

    // copy to device mpep_rank
    thrust::copy(mpep_rank.begin(), mpep_rank.end(), d_out_mpep_rank_th);

    // copy to device mpep_rank_ith_pep
    thrust::device_vector<uint32_t> d_mpep_rank_ith_pep_th(total);
    thrust::copy(mpep_rank_ith_pep.begin(), mpep_rank_ith_pep.end(), d_mpep_rank_ith_pep_th.data());

    // unrank mpep_rank and generate encoding which represents which acids are to be modified 
    // this information will be stored in a 2D matrix
    // The dimension of this matrix for a peptide is:
    // pep_num_mpep * (sum mod_ma_count)  # note: only (sum mod_ma_count) is constant
    // 
    //    ie. foreach mpep_rank value there will be a vector of integers representing the ith acid to modify
    //    eg. a modification with 2 a's, 1 d and 1 c.
    //    a mpep_rank might have an unranked result of: [1, 0, 0, 3]
    //    representing:
    //    a | 1 0     # the 0th and 1st a's are modified
    //    d | 0       # the 0th d is to modified
    //    c | 3       # the 3rd c is modificed
    //
    // for all peptides, the length of the flattened matrices is:
    // sum (pep_num_mpep * (sum mod_ma_count)) = total * (sum mod_ma_count)
    //
    //
    uint32_t sum_mod_ma_count = thrust::reduce(d_mod_ma_count_th, d_mod_ma_count_th + mod_num_ma);
    //
    //thrust::device_vector<uint32_t> d_out_mpep_unrank_th(total*sum_mod_ma_count);
    //
    // Add an unrank to each mpep.
    // call kernel to each mpep + mpep_rank combination
    // index to mpep_rank and mpep_idx array
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + total;
    thrust::for_each(first, last, add_mpep_unrank<uint32_t>(
                        d_out_mpep_unrank_th, 
                        d_pep_ma_num_comb_scan_th,
                        d_out_mpep_rank_th,
                        d_mpep_rank_ith_pep_th.data(),
                        d_pep_ma_count_th,
                        d_mod_ma_count_th,
                        mod_num_ma,
                        sum_mod_ma_count)); 

    // print
    std::cout << "d_out_mpep_rank" << std::endl;
    for (uint32_t i = 0; i < total; i++) {
        std::cout << d_mpep_rank_ith_pep_th[i] << " " << d_out_mpep_idx_th[i] << " " << d_out_mpep_rank_th[i] << std::endl;
    }
    std::cout << std::endl;

    // print
    std::cout << "d_out_mpep_unrank" << std::endl;
    for (uint32_t i = 0; i < total; i++) {
        for (uint32_t j = 0; j < sum_mod_ma_count; ++j) {
            std::cout << d_out_mpep_unrank_th[i*sum_mod_ma_count + j] << " " ;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
}
