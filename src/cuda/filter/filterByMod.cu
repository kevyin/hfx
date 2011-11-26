/* -----------------------------------------------------------------------------
 *
 * Module    : Filter
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <stdint.h>

#include "algorithms.h"

uint32_t
findModablePeptides
(
    uint32_t            *d_out,

    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,

    const uint32_t      *d_sub_idx,
    const uint32_t      sub_nIdx,

    const uint8_t       *d_ma,
    const uint8_t       *d_ma_count
)
{
    return 0;
}

