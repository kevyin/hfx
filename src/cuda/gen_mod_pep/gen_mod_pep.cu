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
 *             C   | 3 |       # the 3rd C is modificed
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
 * returns num_mpep (number of modified peptides) for a peptide
 * used for checking and comparison
 */
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

/**
 * return pep_ma_num_comb
 * the number of combinations each modable acid can make in a peptide
 */
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

/**
 * find num_mpep per peptide
 * at the same time fill pep_ma_num_comb_scan
 */
template <typename T>
struct calcNumMPepPerPep: public thrust::unary_function<T, uint32_t>
{

    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_scan;
    thrust::device_ptr<const uint32_t>  d_pep_ma_num_comb;
    const uint32_t                      mod_num_ma;

    __host__ __device__
    calcNumMPepPerPep(thrust::device_ptr<const uint32_t>    _pep_ma_num_comb, 
                      thrust::device_ptr<uint32_t>          _pep_ma_num_comb_scan,
                      const uint32_t                        _mod_num_ma) : 
                      d_pep_ma_num_comb(_pep_ma_num_comb), 
                      d_out_pep_ma_num_comb_scan(_pep_ma_num_comb_scan), 
                      mod_num_ma(_mod_num_ma) {}

    __host__ __device__ uint32_t operator() (T pep)
    {
        uint32_t pep_ma_num_comb_idx = pep * mod_num_ma;

        thrust::device_ptr<const uint32_t> pep_ma_num_comb_row(d_pep_ma_num_comb + pep_ma_num_comb_idx);
        thrust::device_ptr<uint32_t>       pep_ma_num_comb_scan_row(d_out_pep_ma_num_comb_scan + pep_ma_num_comb_idx);
        uint32_t num = 1;
        //for (uint32_t i = 0; i < mod_num_ma; ++i) 
        for (int32_t i = mod_num_ma - 1; i >= 0; --i) 
        {
            num *= pep_ma_num_comb_row[i];
            // record scan results
            uint32_t *tmp_raw = (pep_ma_num_comb_scan_row + i).get();
            *tmp_raw = num;
        }
        return num;
    }
};

/**
 * Determine and record 
 *   - the number of modified peptides each peptide will generate
 *   - the number of combinations each ma within a peptide makes
 *   - the number of combinations each ma within a peptide makes scanned
 *     return total modified candidates (excludes orig)
 */
uint32_t
calcTotalModCands
(
    uint32_t          *d_out_pep_num_mpep,              // number of modified peptides each peptide will generate
    uint32_t          *d_out_pep_ma_num_comb,           // 2D array of number of combinations each ma can make
    uint32_t          *d_out_pep_ma_num_comb_scan,      // 2D array of above scanned by peptide
    const uint32_t    nPep,                             // number of peptides (unmodified)

    const uint32_t    *d_pep_ma_count,                  // 2d array of ma count in each pep

    const uint8_t     *d_mod_ma_count,                  // the modification
    const uint32_t    mod_num_ma
) 
{
    // initialise thrust ptrs
    thrust::device_ptr<uint32_t>        d_out_pep_num_mpep_th(d_out_pep_num_mpep);
    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_th(d_out_pep_ma_num_comb); 
    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_scan_th(d_out_pep_ma_num_comb_scan);
    thrust::device_ptr<const uint32_t>  d_pep_ma_count_th(d_pep_ma_count);
    thrust::device_ptr<const uint8_t>   d_mod_ma_count_th(d_mod_ma_count);
    
    // for each ma find number of combinations possible
    // ie do a choose(number in peptide, number req for mod)
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + nPep*mod_num_ma;
    thrust::transform(first,last, d_out_pep_ma_num_comb_th, calcNumCombPerAcid<const uint32_t>(d_pep_ma_count_th, d_mod_ma_count_th, mod_num_ma));

    // now determine the number of modified peptides per peptide 
    // by multiplying the number of combinations for each ma
    // at the same time fill out the _scan array
    last = first + nPep;
    thrust::transform(first, last, d_out_pep_num_mpep_th, calcNumMPepPerPep<const uint32_t>(d_out_pep_ma_num_comb_th, d_out_pep_ma_num_comb_scan_th, mod_num_ma));


#ifdef _DEBUG
    // compare with calcNumModCand method
    // malloc d_tmp
    thrust::device_ptr<uint32_t> d_tmp = thrust::device_malloc<uint32_t>(nPep);
    thrust::transform(first, last, d_tmp, calcNumModCand<const uint32_t>(d_pep_ma_count_th, d_mod_ma_count_th, mod_num_ma));

    // check tmp is same as d_pep_ma_count
    if (not thrust::equal(d_out_pep_num_mpep_th, d_out_pep_num_mpep_th + nPep, d_tmp))
    {
        std::cerr << "pep_num_mpep doesn't seem to be correct" << std::endl;
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
        ////
        ////print
        //std::cout << "pep_mpep_count " << std::endl;
        //for (int i = 0; i < nPep; i++)               
        //{
            //std::cout << d_out_pep_num_mpep_th[i] << std::endl;
        //}
        exit(1);
    }
    thrust::device_free(d_tmp);
#endif

    uint32_t total = thrust::reduce(d_out_pep_num_mpep_th, d_out_pep_num_mpep_th + nPep);

    return total;
}

/**
 * Given a mpep_rank 
 * determine the modable acid ranks
 */
__host__ __device__ void
add_ma_ranks 
(
    uint32_t                            *d_out_ma_ranks_row,
    const uint32_t                      mpep_rank,
    thrust::device_ptr<const uint32_t>  d_pep_ma_num_comb_scan_row,
    const uint32_t                      mod_num_ma
)
{
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
            } else {
                ma_ranks_raw[i] = 0;
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
        uint32_t ith_pep   = d_mpep_rank_ith_pep[mpep];

        uint32_t d_ma_ranks[MAX_MA];
        add_ma_ranks(d_ma_ranks, mpep_rank, d_pep_ma_num_comb_scan + d_mpep_rank_ith_pep[mpep]*mod_num_ma, mod_num_ma);

        uint32_t unrank_pos = 0;
        
        // foreach ma, unrank and add to mpep_unrank
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
    thrust::device_ptr<uint32_t> d_out_mpep_idx_th(d_out_mpep_idx);
    thrust::device_ptr<uint32_t> d_out_mpep_rank_th(d_out_mpep_rank);
    thrust::device_ptr<uint32_t> d_out_mpep_unrank_th(d_out_mpep_unrank);
   
    thrust::device_ptr<const uint32_t> d_pep_idx_th(d_pep_idx);
    thrust::device_ptr<const uint32_t> d_pep_num_mpep_th(d_pep_num_mpep);
    thrust::device_ptr<const uint32_t> d_pep_ma_num_comb_th(d_pep_ma_num_comb);
    thrust::device_ptr<const uint32_t> d_pep_ma_num_comb_scan_th(d_pep_ma_num_comb_scan);

    thrust::device_ptr<const uint32_t> d_pep_ma_count_th(d_pep_ma_count);
    thrust::device_ptr<const uint8_t>  d_mod_ma_count_th(d_mod_ma_count);



    thrust::host_vector<uint32_t> h_mpep_rank;
    thrust::host_vector<uint32_t> h_mpep_rank_ith_pep; // helper vector, corres pep for each mpep_rank
    thrust::host_vector<uint32_t> h_mpep_idx;

    for (uint32_t i = 0; i < num_pep; i++)
    {
        for (uint32_t k = 0; k < d_pep_num_mpep_th[i]; k++)
        {
            // add mpep_rank
            h_mpep_rank.push_back(k);
            // add mpep_rank_ith_pep
            h_mpep_rank_ith_pep.push_back(i);
            // add idx
            h_mpep_idx.push_back(d_pep_idx_th[i]);
        }                      
    }
    if (h_mpep_rank.size() != total) {
        std::cerr << "mpep_rankvector not the correct length" << std::endl;
        exit(1);
    }

    // copy to device idx
    thrust::copy(h_mpep_idx.begin(), h_mpep_idx.end(), d_out_mpep_idx_th);

    // copy to device mpep_rank
    thrust::copy(h_mpep_rank.begin(), h_mpep_rank.end(), d_out_mpep_rank_th);

    // copy to device mpep_rank_ith_pep
    //thrust::device_vector<uint32_t> d_mpep_rank_ith_pep_th(total);
    //thrust::copy(mpep_rank_ith_pep.begin(), mpep_rank_ith_pep.end(), d_mpep_rank_ith_pep_th.data());
    thrust::device_vector<uint32_t> d_mpep_rank_ith_pep_th = h_mpep_rank_ith_pep;

    // unrank mpep_rank and generate encoding which represents which acids are to be modified 
    //
    //
    uint32_t sum_mod_ma_count = thrust::reduce(d_mod_ma_count_th, d_mod_ma_count_th + mod_num_ma);
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

#ifdef _DEBUG
    std::cout << "checking genModCands" << std::endl;

    thrust::host_vector<uint32_t> h_check_mpep_unrank_th;
    thrust::host_vector<uint32_t> h_check_mpep_rank_th;

    for (uint32_t i = 0; i < num_pep; i++) { // for each pep
        for (uint32_t r = 0; r < d_pep_num_mpep_th[i]; ++r) { // foreach mpep
            // get unranks 
            thrust::host_vector<uint32_t> unrank_tmp(sum_mod_ma_count);
            uint32_t x = r; //mpep_rank
            uint32_t idx_start = 0;
            for (uint32_t j = 0; j < mod_num_ma; ++j) { // for each ma
                uint32_t ma_rank;
                if (j != mod_num_ma - 1) {
                    uint32_t next = d_pep_ma_num_comb_scan_th[i*mod_num_ma + j + 1];
                    if(x >= next) {
                        ma_rank = x / next;
                        x = x % next;
                    } else {
                        ma_rank = 0;
                    }

                } else {
                    ma_rank = x;
                }


                unrankComb(unrank_tmp.data() + idx_start,
                           d_pep_ma_count_th[i*mod_num_ma +j],
                           d_mod_ma_count_th[j],
                           ma_rank);
                idx_start += d_mod_ma_count_th[j];
            }

            h_check_mpep_rank_th.push_back(r);
            // add unranks
            h_check_mpep_unrank_th.insert(h_check_mpep_unrank_th.end(), unrank_tmp.begin(), unrank_tmp.end());
        }
        
    }

    thrust::device_vector<uint32_t> d_check_mpep_unrank_th = h_check_mpep_unrank_th;
    if (!thrust::equal(d_check_mpep_unrank_th.begin(), d_check_mpep_unrank_th.end(), d_out_mpep_unrank_th)) {
        std::cerr << "mpep_unrank doesn't seem to be correct" << std::endl;

        // print
        for (uint32_t i = 0; i < total; i++) {
            std::cout << "ith pep, mpep_idx, out_mpep_rank, check_mpep_rank" << std::endl;
            std::cout << d_mpep_rank_ith_pep_th[i] << " " << d_out_mpep_idx_th[i] << " " << d_out_mpep_rank_th[i] << " check " << h_check_mpep_rank_th[i] << std::endl;
            std::cout << "unrank out   ";
            for (uint32_t j = 0; j < sum_mod_ma_count; ++j) {
                std::cout << d_out_mpep_unrank_th[i*sum_mod_ma_count + j] << " " ;
            }
            std::cout << std::endl;

            std::cout << "unrank check ";
            for (uint32_t j = 0; j < sum_mod_ma_count; ++j) {
                std::cout << d_check_mpep_unrank_th[i*sum_mod_ma_count + j] << " " ;
            }
            std::cout << std::endl << std::endl;
        }
        std::cout << std::endl;
        exit(1);
    }


#endif
    
}

