
/* -----------------------------------------------------------------------------
 *
 * Module    : Combinatorics
 * Copyright : (c) [2009..2010] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __COMB_H__
#define __COMB_H__

__host__ __device__ uint32_t 
choose(const uint32_t n, const uint32_t k) 
{
    if (k == 0 || n == k) return 1;
    if (n < k) return 0;
    uint32_t i;
    uint32_t t = 1;

    if (k < n - k) {
        for (i = n; i >= n - k + 1; --i) {
            t = t * i / (n - i + 1);
        }
    } else {
        for (i = n; i >= k + 1; --i) {
            t = t * i / (n - i + 1);
        }
    }
    return t;
}

#endif
