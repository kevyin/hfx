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
 * find num_mpep per peptide
 * at the same time fill pep_ma_num_comb_scan
 */
template <typename T>
struct calcNumMPepPerPep: public thrust::unary_function<T, uint32_t>
{
    const T    *d_mod_ma_count;
    const T    *d_pep_mod_idx;
    const T    *d_pep_ma_count;
    uint32_t   *d_out_pep_ma_num_comb_scan;
    const T    num_ma;

    __host__ __device__
    calcNumMPepPerPep(const T    *_mmc, 
                       const T    *_pmi, 
                       const T    *_pmc, 
                       uint32_t   *_pep_ma_num_comb_scan,
                       const T    _num_ma) : 

                      d_mod_ma_count(_mmc), 
                      d_pep_mod_idx(_pmi),
                      d_pep_ma_count(_pmc),
                      d_out_pep_ma_num_comb_scan(_pep_ma_num_comb_scan), 
                      num_ma(_num_ma) {}

    __host__ __device__ uint32_t operator() (T pep_idx)
    {
        uint32_t pep_ma_num_comb_idx = pep_idx * num_ma;

        uint32_t       *pep_ma_num_comb_scan_row = d_out_pep_ma_num_comb_scan + pep_ma_num_comb_idx;
        uint32_t num = 1;
        for (int32_t i = num_ma - 1; i >= 0; --i) 
        {
            // calc number of combinations each modable acid can make in a peptide
            uint32_t mod_idx = d_pep_mod_idx[pep_idx];
            num *= choose(d_pep_ma_count[num_ma*pep_idx + i],
                          d_mod_ma_count[num_ma*mod_idx + i]);

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
    uint32_t          *d_out_pep_ma_num_comb_scan_raw,      // 2D array of above scanned by peptide
    const uint32_t    *d_mod_ma_count_raw,
    const uint32_t    *d_pep_idx_raw,
    const uint32_t    *d_pep_mod_idx_raw,
    const uint32_t    *d_pep_ma_count_raw,
    const uint32_t    num_cand_massmod,                             // number of peptides (unmodified)
    const uint32_t    num_mod,
    const uint32_t    num_ma
) 
{

#ifdef _BENCH
    cudaThreadSynchronize();
    std::cout << "calcTotalModCands" << std::endl;
    time_t t_beg, t_end;
    time(&t_beg);
#endif 

    // initialise thrust ptrs
    thrust::device_ptr<uint32_t>        d_out_pep_num_mpep(d_out_pep_num_mpep_raw);
    
    // Determine the number of modified peptides per peptide 
    // by multiplying the number of combinations for each ma
    // at the same time fill out the _scan array
    //
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_cand_massmod;

    thrust::transform(first,last, d_out_pep_num_mpep, calcNumMPepPerPep<const uint32_t>(d_mod_ma_count_raw, d_pep_mod_idx_raw, d_pep_ma_count_raw, d_out_pep_ma_num_comb_scan_raw, num_ma));


    uint32_t num_mpep = thrust::reduce(d_out_pep_num_mpep, d_out_pep_num_mpep + num_cand_massmod);

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    printf ("Time elapsed for calcTotalModCands: %.2lf seconds\n", difftime(t_end,t_beg));
#endif 
    return num_mpep;
}

