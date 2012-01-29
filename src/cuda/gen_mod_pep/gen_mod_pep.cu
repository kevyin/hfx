/* -----------------------------------------------------------------------------
 *
 * Module    : Generate Modified Peptides
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * This module aims to generate information about modified peptides
 * which mions_series.cu can then use to determine which acids in a 
 * modified peptide are to be modified when generating spectrums.
 *
 * Currently different values for one acid is not accepted.
 *
 * Intro:
 * A peptide is a sequence of amino acids.
 * eg. AIGIRRAOQAANKGHJE
 * The functions in this module applies a specific
 * modification to a peptide. Something like:
 * 3A  -18
 * 1R  +42
 *
 * The ways a modification can be applied to a peptide is
 * the ways modable acids can be chosen out of the peptide to be
 * modified.
 *
 * Following the above example, there are 4A's and 2R's in the peptide.
 * Therefore there are:
 * Choose(4, 3) * Choose(2, 1) = 8 total ways the modification can be applied.
 *
 * 
 * The strategy calls two main functions, both of which are in this module.
 * calcTotalModCands -- which counts the number of relevant acids
 * genModCands       -- uses the counts to expand into rows of information
 *                      each row represents a modified peptide which
 *                      is to be later decoded into a spectrum or other forms
 *
 * More details below:
 *
 * calcTotalModCands
 * outputs:
 *      d_out_pep_num_mpep -- The number of modified versions each 
 *          peptide will generate. 
 *    
 *          For one peptide this is:
 *          num_mpep = 1
 *          foreach i in ma:
 *              num_mpep *= Choose(D_i, d_i)
 *    
 *              where
 *               ma = list of acids in the modification being applied
 *               D_i = number of acid i in the peptide
 *               d_i = number of acid i in the modification
 *    
 *       This is used to allocate memory and for the mpep_rank, 
 *       which is defined to be the
 *       0th, 1st, 2nd ... ith modified peptide generated from
 *          a peptide. This "rank" is not technically the same as the one
 *          used in combinatorics as "combinations" of a modified peptide are
 *          really combinations of it's consituant modable amino acids.
 *     
 *      d_out_pep_ma_num_comb -- A 2D array of the number of combinations each 
 *          modable amino acid can make in a peptide. 
 *    
 *          ie. Choose(D_i, d_i) from above
 *          the dimension of this 2D array would be:
 *              num_mpep * num_ma
 *    
 *              where
 *                  num_mpep = sum d_out_pep_num_mpep
 *                  num_ma   = number of different modable acids 
 *                             in the modification
 *    
 *      d_out_pep_ma_num_comb_scan -- The above scanned by row.
 *          This array and the one above helps to determine ma_rank
 *    
 *          When decoding/unranking a modified peptide into an 
 *          acid/spectrum representation, mpep_rank is used to 
 *          determine which modable acid rank (ie combination) describes this
 *          modified peptide.
 *
 *
 * genModCands
 * each row of the last two outputs below describe a modified peptide
 * outputs:
 *      d_out_mpep_rank -- "rank" or ith (from 0) modified peptide generated
 *          from a peptide.
 *          eg if d_out_pep_num_mpep values are [4, 6, 1, 2]:
 *             [0..3 0..5 0..2 0 0 1]
 *
 *      d_out_mpep_idx -- index to the original peptide ie to d_tc, d_tn etc.
 *             same length as above array, index values will have repetition
 *
 *      d_out_mpep_unrank -- unranked ma_ranks. foreach modable acid, the ith 
 *          ones which are to be modified.
 *
 *          Stored in a 2D matrix
 *          The dimension of this is:
 *          num_mpep * (sum mod_ma_count)  
 *           
 *             ie. foreach mpep_rank value there will be a row of integers 
 *                  representing the ith acid to modify
 *             eg. a modification with 2 A's, 1 D and 1 C.
 *             a mpep_rank might have an unranked row of: [1, 0, 0, 3]
 *             representing:
 *                 --------------
 *             A   | 1 | 0     # the 0th and 1st A's are modified
 *             D   | 0 |       # the 0th D is to modified
 *             C   | 3 |       # the 3rd C is modified
 *                 -------------
 *              
 *      
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "algorithms.h"
#include "combinatorics.h"

/**
 * Given a mpep_rank 
 * determine the modable acid ranks
 */
__host__ __device__ void
add_ma_ranks 
(
    uint32_t                            *d_out_ma_ranks_row,
    const uint32_t                      mpep_rank,
    const uint32_t                      *d_pep_ma_num_comb_scan_row,
    const uint32_t                      num_ma
)
{
    uint32_t *ma_ranks_raw = d_out_ma_ranks_row;

    uint32_t x = mpep_rank;

    for (uint32_t i = 0; i < num_ma; ++i)
    {
        if (i != num_ma - 1)
        {
            uint32_t next_pep_ma_num_comb_scan = d_pep_ma_num_comb_scan_row[i + 1];
            if (x >= next_pep_ma_num_comb_scan)
            {
                ma_ranks_raw[i] = x / next_pep_ma_num_comb_scan;
                x = x % next_pep_ma_num_comb_scan;
            } else {
                ma_ranks_raw[i] = 0;
            }

        } else {
            ma_ranks_raw[i] = x;
        }
    }
}

template <typename T, uint32_t NumMA>
struct add_mpep_unrank: public thrust::unary_function<T, uint32_t>
{
    uint32_t        *d_out_mpep_unrank;

    const uint32_t  *d_mod_ma_count;

    const uint32_t  *d_mpep_ith_valid;
    const uint32_t  *d_mpep_rank;
    const uint32_t  *d_mpep_mod_ma_count_sum_scan;

    const uint32_t  *d_pep_mod_idx;
    const uint32_t  *d_pep_ma_count;

    const uint32_t  *d_pep_valid_idx;
    const uint32_t  *d_pep_ma_num_comb_scan;

    __host__ __device__
    add_mpep_unrank (uint32_t        *_omu,
                     const uint32_t  *_mmc,
                     const uint32_t  *_miv,
                     const uint32_t  *_mr,
                     const uint32_t  *_mmmcss,
                 
                     const uint32_t  *_pmi,
                     const uint32_t  *_pmc,

                     const uint32_t  *_pvi,
                     const uint32_t  *_pmncs) : 
                     d_out_mpep_unrank(_omu),
                     d_mod_ma_count(_mmc),
                     d_mpep_ith_valid(_miv),
                     d_mpep_rank(_mr),
                     d_mpep_mod_ma_count_sum_scan(_mmmcss),
                     d_pep_mod_idx(_pmi),
                     d_pep_ma_count(_pmc),

                     d_pep_valid_idx(_pvi),
                     d_pep_ma_num_comb_scan(_pmncs) {}

    __host__ __device__ uint32_t operator() (T mpep)
    {
        const uint32_t ith_valid   = d_mpep_ith_valid[mpep];
        const uint32_t ith_pep     = d_pep_valid_idx[ith_valid];
        const uint32_t mpep_rank   = d_mpep_rank[mpep];
        const uint32_t pep_mod_idx = d_pep_mod_idx[ith_pep];

        uint32_t d_ma_ranks[NumMA];
        add_ma_ranks(d_ma_ranks, mpep_rank, d_pep_ma_num_comb_scan + ith_valid*NumMA, NumMA);

        uint32_t start = d_mpep_mod_ma_count_sum_scan[mpep];

        uint32_t unrank_pos = 0;
        
        // foreach ma, unrank and add to mpep_unrank
        for (uint32_t i = 0; i < NumMA; ++i)
        {
            const uint32_t mod_ma_count = d_mod_ma_count[pep_mod_idx*NumMA + i]; 
            if (mod_ma_count != 0) {
                //unrankComb(&d_out_mpep_unrank[start + unrank_pos], 
                unrankComb((d_out_mpep_unrank + start + unrank_pos), 
                           d_pep_ma_count[ith_pep*NumMA + i], 
                           mod_ma_count,
                           d_ma_ranks[i]);
                unrank_pos += mod_ma_count;
             }
        }
        
        return mpep;
        //return num;
    }
};

void
genModCands
(                                                                     
    uint32_t        *d_out_mpep_unrank,

    const uint32_t   *d_mod_ma_count,             

    const uint32_t  *d_mpep_ith_valid,
    const uint32_t  *d_mpep_rank,
    const uint32_t  *d_mpep_mod_ma_count_sum_scan,

    const uint32_t  *d_pep_mod_idx,
    const uint32_t  *d_pep_ma_count,        // 2d array of ma count in each pep

    const uint32_t  *d_pep_valid_idx,
    const uint32_t  *d_pep_ma_num_comb_scan,

    const uint32_t  num_mpep,                        
    const uint32_t  num_ma
)
{
#ifdef _BENCH
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
#endif 
    //std::cout << "gen_mod_pep" << std::endl;
    //printGPUMemoryUsage();

/*
    std::cout << "num_mpep " << num_mpep << std::endl;

    thrust::device_ptr<const uint32_t> d_mod_ma_count_th(d_mod_ma_count);
    thrust::device_ptr<const uint32_t> d_mpep_ith_valid_th(d_mpep_ith_valid);
    thrust::device_ptr<const uint32_t> d_pep_valid_idx_th(d_pep_valid_idx);
    thrust::device_ptr<const uint32_t> d_mpep_rank_th(d_mpep_rank);
    thrust::device_ptr<const uint32_t> d_pep_mod_idx_th(d_pep_mod_idx);
    thrust::device_ptr<const uint32_t> d_pep_ma_count_th(d_pep_ma_count);
    for (uint32_t i = 0; i < num_mpep; ++i) {
        //uint32_t ith_valid = d_mpep_ith_valid_th[i];
        //uint32_t valid_idx = d_pep_valid_idx_th[ith_valid];
        //uint32_t mod_idx = d_pep_mod_idx_th[valid_idx];
        const uint32_t ith_valid   = d_mpep_ith_valid_th[i];
        const uint32_t ith_pep     = d_pep_valid_idx_th[ith_valid];
        const uint32_t mpep_rank   = d_mpep_rank_th[i];
        const uint32_t pep_mod_idx = d_pep_mod_idx_th[ith_pep];
        std::cout << i << " num_ma " << num_ma 
                << " ith_valid " << ith_valid
                << " ith_pep" << ith_pep
                << " mpep_rank " << mpep_rank
                << " pep_mod_idx " << pep_mod_idx
                << " d_pep_ma_count " << d_pep_ma_count_th[ith_pep*num_ma] << " " <<  d_pep_ma_count_th[ith_pep*num_ma + 1] 
                <<  " d_mod_ma_count " << d_mod_ma_count_th[pep_mod_idx*num_ma] << " " << d_mod_ma_count_th[pep_mod_idx*num_ma +1]
                <<  std::endl;
    }
    */

    // Add an unrank to each mpep.
    // call kernel to each mpep + mpep_rank combination
    // index to mpep_rank and mpep_idx array
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_mpep;
    switch (num_ma)
    {
    case 1: thrust::for_each(first, last, add_mpep_unrank<uint32_t,1>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 2: thrust::for_each(first, last, add_mpep_unrank<uint32_t,2>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 3: thrust::for_each(first, last, add_mpep_unrank<uint32_t,3>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 4: thrust::for_each(first, last, add_mpep_unrank<uint32_t,4>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 5: thrust::for_each(first, last, add_mpep_unrank<uint32_t,5>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 6: thrust::for_each(first, last, add_mpep_unrank<uint32_t,6>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 7: thrust::for_each(first, last, add_mpep_unrank<uint32_t,7>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 8: thrust::for_each(first, last, add_mpep_unrank<uint32_t,8>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 9: thrust::for_each(first, last, add_mpep_unrank<uint32_t,9>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    case 10: thrust::for_each(first, last, add_mpep_unrank<uint32_t,10>(d_out_mpep_unrank, d_mod_ma_count, d_mpep_ith_valid, d_mpep_rank, d_mpep_mod_ma_count_sum_scan, d_pep_mod_idx, d_pep_ma_count, d_pep_valid_idx, d_pep_ma_num_comb_scan)); break;
    default:
        assert(!"Maximum number of variable acids is 10");
    }

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    printf ("Time elapsed for genMod: %.2lf seconds\n", difftime(t_end,t_beg));
#endif 

}

