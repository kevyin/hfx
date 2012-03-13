/* -----------------------------------------------------------------------------
 *
 * Module    : Algorithms
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __ALGORITHMS_H__
#define __ALGORITHMS_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Maximum modable aa's
#define MAX_MA  5

#if 0
/*
 * Scan and segmented-scan
 * Unlike Haskell-style scans, the size of the output array will not change.
 *   Pre-scans are exclusive:   prescanl  == init . scanl
 *   Post-scans are inclusive:  postscanl == scanl1
 */
void prescanl_plusui(const uint32_t *d_in, uint32_t *d_out, const uint32_t N);
void prescanr_plusui(const uint32_t *d_in, uint32_t *d_out, const uint32_t N);
#endif

/*
 * key-value sort (in-place)
 */
void sort_rf(float *d_keys, uint32_t *d_vals, uint32_t N);

/*
 * key-value sort (in-place)
 */
void sort_val_f(float *d_keys, uint32_t *d_vals, uint32_t N);

#if 0
/*
 * Permute (32-bit payload)
 */
void permute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length);
void bpermute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length);
uint32_t compact(const void *d_in, void *d_out, const uint32_t *d_flags, const uint32_t length);
uint32_t compactIndices(uint32_t *d_out, const uint32_t *d_flags, const uint32_t length);
#endif
#if 0
/*
 * Filter
 */
uint32_t
filterInRange_f
(
    const float         *d_in,
    float               *d_out,
    const uint32_t      length,
    const float         min,
    const float         max
);
#endif

uint32_t
findIndicesInRange_f
(
    const float         *d_in,
    uint32_t            *d_out,
    const uint32_t      length,
    const float         min,
    const float         max
);

uint32_t
findBeginEnd_f
(
    uint32_t            *d_begin_raw,
    uint32_t            *d_end_raw,
    uint32_t            *d_num_pep_raw,
    uint32_t            *d_num_pep_scan_raw,

    const float         *d_r_raw,
    const uint32_t      *d_pep_idx_r_sorted_raw,
    const uint32_t      num_pep,

    const float         *d_mod_delta_raw,
    const uint32_t      num_mod,
    const float         mass_,
    const float         eps
);

uint32_t
findModablePeptides
(
    uint32_t            *d_out_pep_idx_valid_raw, // valid indices
    uint32_t            *d_out_pep_idx_raw, 
    uint32_t            *d_out_pep_mod_idx_raw, 
    uint32_t            *d_out_pep_ma_count_raw,   // 2d array, count of each ma in each peptide
    uint32_t            num_pep_total,

    const uint8_t       *d_ions_raw,
    const uint32_t      *d_tc_raw,
    const uint32_t      *d_tn_raw,

    const uint32_t      *d_pep_idx_r_sorted_raw,

    const uint32_t      *d_begin_raw,
    const uint32_t      *d_end_raw,
    const uint32_t      *d_num_pep_scan_raw,
    const uint32_t      *d_mod_ma_count_raw,
    const uint32_t      num_mod_raw,

    const uint8_t       *d_ma_raw,
    const uint32_t      num_ma
);

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
); 

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
);

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

    const uint32_t  total,
    const uint32_t  mod_num_ma
);

void addModIons
(
    uint32_t            *d_out_mspec,
    const float         *d_residual,    // peptide residual mass
    const float         *d_mass,        // lookup table for ion character codes ['A'..'Z']
    const uint8_t       *d_ions,        // individual ion character codes (the database)
    const uint32_t      *d_tc,          // c-terminal indices
    const uint32_t      *d_tn,          // n-terminal indices

    const uint32_t      *_d_mpep_pep_idx,
    const uint32_t      *_d_mpep_pep_mod_idx,
    const uint32_t      *d_mpep_unrank,
    const uint32_t      *_d_mpep_mod_ma_count_sum_scan,
    const uint32_t      len_unrank,
    const uint32_t      *_num_mpep,

    const uint32_t      *d_mod_ma_count,
    const float         *d_mod_delta,

    const uint8_t       *d_ma,
    const float         *d_ma_mass,
    const uint32_t      num_ma,

    const uint32_t      len_spec,
    const uint32_t      max_charge,

    const uint32_t      start,
    const uint32_t      nsearch
);

/*
 * Generate theoretical spectra
 */
void
addIons
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
);

/*
 * count matching acids
 */
uint32_t
countFragByRule
(
    const uint8_t       *d_in_raw,
    const uint32_t      N,
    const uint32_t      ruleNo
);

/*
 * digest proteins by ruleNo
 */
uint32_t
digestByRule
(
    const uint8_t       *d_ions_raw,
    const uint32_t      Ni,             // Number of ions
    const uint32_t      *d_iseg_raw,    // protein segments in ions
    const uint32_t      Ns,             // number of segments (proteins)
    const uint32_t      Nf,             // upper bound on protein frags
    const float         *d_mass_raw,    // amino acid masses
    const uint32_t      ruleNo,

    // result vectors
    float         *d_r_raw,
    uint32_t      *d_c_raw,
    uint32_t      *d_n_raw,
    uint32_t      *d_fs_raw
);

#if 0
void
addIons_inplace
(
    float               *d_score,
    const float         *d_spec,
    const float         *d_residual,
    const float         *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
);
#endif

/*
 * Dense matrix-vector multiplication
 */
void mvm_if(float *d_y, const uint32_t *d_A, const float *d_x, const uint32_t m, const uint32_t n);

/*
 * reduce a word32 array
 */
uint32_t sum_Word32(const uint32_t *d_array, const uint32_t len);

// Non Parallel
void
getSpecNonParallel(
    uint32_t            *d_out_check_spec, 
    const uint8_t       *d_in_mions, 
    const float         *d_residual_raw, 
    const float         *d_mass_raw, 
    const uint8_t       *d_ions_raw, 
    const uint32_t      *d_tc_raw,
    const uint32_t      *d_tn_raw,
    const uint32_t      *d_mpep_pep_idx,
    const uint32_t      *d_mpep_pep_mod_idx,
    const uint32_t      *d_mpep_unrank_raw,
    const uint32_t      *d_mpep_mod_ma_count_sum_scan,
    const uint32_t      len_urank,
    const uint32_t      num_mpep,
    const uint32_t      *d_mod_ma_count_raw,
    const float         *d_mod_delta_raw,
    const uint8_t       *d_ma_raw,
    const float         *d_ma_mass_raw,
    const uint32_t      num_ma,

    const uint32_t      max_charge,
    const uint32_t      len_spec
    
);

#ifdef __cplusplus
}
#endif
#endif
