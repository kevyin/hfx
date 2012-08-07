/* -----------------------------------------------------------------------------
 *
 * Module    : Prepare scoring
 * Copyright : (c) [2012] Kevin Ying 
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/copy.h>

#include "algorithms.h"

using namespace thrust;

void prepare_scoring
(
    uint32_t        *d_out_spec_pep_idx,

    const uint32_t  *d_pep_idx_r_sorted_raw,
    const uint32_t  *spec_begin,
    const uint32_t  *spec_num_pep,
    const uint32_t  num_spec
)
{
    device_ptr<uint32_t> d_spec_pep_idx(d_out_spec_pep_idx);
    device_ptr<const uint32_t> d_pep_idx_r_sorted(d_pep_idx_r_sorted_raw);

    uint32_t offset = 0;
    for (size_t i = 0; i < num_spec; ++i)
    {
        counting_iterator<uint32_t> idx_iter(spec_begin[i]);
        copy_n(make_permutation_iterator(d_pep_idx_r_sorted, idx_iter),
               spec_num_pep[i],
               d_spec_pep_idx + offset);
        offset += spec_num_pep[i];
    }
}
