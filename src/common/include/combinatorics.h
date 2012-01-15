
/* -----------------------------------------------------------------------------
 *
 * Module    : Combinatorics
 * Copyright : (c) [2009..2010] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __COMB_H__
#define __COMB_H__

__host__ __device__ __inline__ uint32_t 
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

__host__ __device__ __inline__ uint32_t 
largestV(const uint32_t a, const uint32_t b, const uint32_t x) {
    uint32_t v = a - 1;
    while (choose(v, b) > x)
        --v;

    return v;
}

/*
 * ans should be of length k
 */
__host__ __device__ __inline__ void
unrankComb (uint32_t *ans, const uint32_t n, const uint32_t k, const uint32_t rank) {
    //vector<uint32_t> ans(k);

    uint32_t a = n;
    uint32_t b = k;
    //uint32_t x = choose(n,k) - 1 - rank;
    uint32_t total = choose(n,k);
    uint32_t x; 
    bool usingDual = false;

    if ((total/2) >= rank) {
        x = rank;
    } else {
        x = total - 1 - rank;
        usingDual = true;
    }

    for (uint32_t i = 0; i < k; ++i) {
        ans[i] = largestV(a,b,x);    
        x = x - choose(ans[i],b);
        a = ans[i];
        b = b - 1;
    }

    if (usingDual) {
        for (uint32_t i = 0; i < k; ++i) {
            ans[i] = (n - 1) - ans[i];
        }
    }

    //return ans;
}    

#endif
