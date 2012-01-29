/* -----------------------------------------------------------------------------
 *
 * Module    : Calculate Modified Peptides
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 *      
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "algorithms.h"
#include "combinatorics.h"

#include <time.h>

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

    const T    *d_mod_ma_count;
    const T    *d_pep_mod_idx;
    const T    *d_pep_ma_count;
    const T    *d_pep_valid_idx;
    const T    num_mod;
    const T    num_ma;

    __host__ __device__
    calcNumCombPerAcid(const T    *_mmc, 
                       const T    *_pmi, 
                       const T    *_pmc, 
                       const T    *_pvi, 
                       const T    _num_mod,
                       const T    _num_ma) : 
                   d_mod_ma_count(_mmc), d_pep_mod_idx(_pmi), d_pep_ma_count(_pmc), d_pep_valid_idx(_pvi), num_mod(_num_mod), num_ma(_num_ma) {}

    __host__ __device__ uint32_t operator() (T idx)
    {
        uint32_t valid_idx = d_pep_valid_idx[idx / num_ma];
        uint32_t mod_idx = d_pep_mod_idx[valid_idx];
        uint32_t ma_idx  = idx % num_ma;

        return choose(d_pep_ma_count[num_ma*valid_idx + ma_idx], 
                      d_mod_ma_count[num_ma*mod_idx + ma_idx]);
    }
};

/**
 * find num_mpep per peptide
 * at the same time fill pep_ma_num_comb_scan
 */
template <typename T>
struct calcNumMPepPerPep: public thrust::unary_function<T, uint32_t>
{

    uint32_t        *d_out_pep_ma_num_comb_scan;
    const uint32_t  *d_pep_ma_num_comb;
    const uint32_t  mod_num_ma;

    __host__ __device__
    calcNumMPepPerPep(const uint32_t    *_pep_ma_num_comb, 
                      uint32_t          *_pep_ma_num_comb_scan,
                      const uint32_t                        _mod_num_ma) : 
                      d_pep_ma_num_comb(_pep_ma_num_comb), 
                      d_out_pep_ma_num_comb_scan(_pep_ma_num_comb_scan), 
                      mod_num_ma(_mod_num_ma) {}

    __host__ __device__ uint32_t operator() (T pep)
    {
        uint32_t pep_ma_num_comb_idx = pep * mod_num_ma;

        const uint32_t *pep_ma_num_comb_row = d_pep_ma_num_comb + pep_ma_num_comb_idx;
        uint32_t       *pep_ma_num_comb_scan_row = d_out_pep_ma_num_comb_scan + pep_ma_num_comb_idx;
        uint32_t num = 1;
        //for (uint32_t i = 0; i < mod_num_ma; ++i) 
        for (int32_t i = mod_num_ma - 1; i >= 0; --i) 
        {
            num *= pep_ma_num_comb_row[i];
            // record scan results
            pep_ma_num_comb_scan_row [i] = num;
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
    uint32_t          *d_out_pep_num_mpep_raw,              // number of modified peptides each peptide will generate
    uint32_t          *d_out_pep_ma_num_comb_raw,           // 2D array of number of combinations each ma can make,
    uint32_t          *d_out_pep_ma_num_comb_scan_raw,      // 2D array of above scanned by peptide
    const uint32_t    *d_mod_ma_count_raw,
    const uint32_t    *d_pep_idx_raw,
    const uint32_t    *d_pep_mod_idx_raw,
    const uint32_t    *d_pep_ma_count_raw,
    const uint32_t    *d_pep_valid_idx_raw,
    const uint32_t    nPep,                             // number of peptides (unmodified)
    const uint32_t    num_mod,
    const uint32_t    num_ma
) 
{

#ifdef _BENCH
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
#endif 

    //std::cout << "calcTotalModCands" << std::endl;
    //printGPUMemoryUsage();
    // initialise thrust ptrs
    thrust::device_ptr<uint32_t>        d_out_pep_num_mpep(d_out_pep_num_mpep_raw);
    thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb(d_out_pep_ma_num_comb_raw); 
    //thrust::device_ptr<uint32_t>        d_out_pep_ma_num_comb_scan(d_out_pep_ma_num_comb_scan_raw);
    //thrust::device_ptr<const uint32_t>   d_mod_ma_count(d_mod_ma_count_raw);
    //thrust::device_ptr<const uint32_t>  d_pep_idx(d_pep_idx_raw);
    //thrust::device_ptr<const uint32_t>  d_pep_mod_idx(d_pep_mod_idx_raw);
    //thrust::device_ptr<const uint32_t>  d_pep_ma_count(d_pep_ma_count_raw);
    //thrust::device_ptr<const uint32_t>  d_pep_valid_idx(d_pep_valid_idx_raw);
    
    // for each ma, find number of combinations possible
    // ie do a choose(number in peptide, number req for mod)
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + nPep*num_ma;
    thrust::transform(first, last, d_out_pep_ma_num_comb, calcNumCombPerAcid<uint32_t>(d_mod_ma_count_raw, d_pep_mod_idx_raw, d_pep_ma_count_raw, d_pep_valid_idx_raw, num_mod, num_ma));

    // now determine the number of modified peptides per peptide 
    // by multiplying the number of combinations for each ma
    // at the same time fill out the _scan array
    last = first + nPep;
    thrust::transform(first, last, d_out_pep_num_mpep, calcNumMPepPerPep<const uint32_t>(d_out_pep_ma_num_comb_raw, d_out_pep_ma_num_comb_scan_raw, num_ma));


//#ifdef _DEBUG
    //// compare with calcNumModCand method
    //// malloc d_tmp
    //thrust::device_ptr<uint32_t> d_tmp = thrust::device_malloc<uint32_t>(nPep);
    //thrust::transform(first, last, d_tmp, calcNumModCand<const uint32_t>(d_pep_ma_count_th, d_mod_ma_count_th, mod_num_ma));

    //// check tmp is same as d_pep_ma_count
    //if (not thrust::equal(d_out_pep_num_mpep_th, d_out_pep_num_mpep_th + nPep, d_tmp))
    //{
        //std::cerr << "pep_num_mpep doesn't seem to be correct" << std::endl;
        //////print scans
        ////std::cout << "ma_num_comb not scanned and scanned" << std::endl;
        ////for (int i = 0; i < nPep; i++)               
        ////{
            ////std::cout << "pep ma num_comb and scan" << i << " ";
            ////for (int j = 0; j < mod_num_ma; j++)               
            ////{
                ////int k = i * mod_num_ma + j;
                ////uint32_t q = d_mod_ma_count_th[j];
                ////std::cout << "pep_ma_num_comb " << d_out_pep_ma_num_comb_th[k] << " ";
                ////std::cout << "scan " << d_out_pep_ma_num_comb_scan_th[k] << " | ";
            ////}
            ////std::cout << std::endl;
        ////}
        //////
        //////print
        ////std::cout << "pep_mpep_count " << std::endl;
        ////for (int i = 0; i < nPep; i++)               
        ////{
            ////std::cout << d_out_pep_num_mpep_th[i] << std::endl;
        ////}
        //exit(1);
    //}
    //thrust::device_free(d_tmp);
//#endif

    uint32_t num_mpep = thrust::reduce(d_out_pep_num_mpep, d_out_pep_num_mpep + nPep);

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    printf ("Time elapsed for calcTotalModCands: %.2lf seconds\n", difftime(t_end,t_beg));
#endif 
    return num_mpep;
}

