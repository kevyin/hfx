
/* -----------------------------------------------------------------------------
 *
 * Module    : Combinatorics
 * Copyright : (c) [2009..2010] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __COMB_H__
#define __COMB_H__

  __constant__  const uint32_t COMB_42_11_ [43][12] = {
        {1,0,0,0,0,0,0,0,0,0,0,0},
        {1,1,0,0,0,0,0,0,0,0,0,0},
        {1,2,1,0,0,0,0,0,0,0,0,0},
        {1,3,3,1,0,0,0,0,0,0,0,0},
        {1,4,6,4,1,0,0,0,0,0,0,0},
        {1,5,10,10,5,1,0,0,0,0,0,0},
        {1,6,15,20,15,6,1,0,0,0,0,0},
        {1,7,21,35,35,21,7,1,0,0,0,0},
        {1,8,28,56,70,56,28,8,1,0,0,0},
        {1,9,36,84,126,126,84,36,9,1,0,0},
        {1,10,45,120,210,252,210,120,45,10,1,0},
        {1,11,55,165,330,462,462,330,165,55,11,1},
        {1,12,66,220,495,792,924,792,495,220,66,12},
        {1,13,78,286,715,1287,1716,1716,1287,715,286,78},
        {1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364},
        {1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365},
        {1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368},
        {1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376},
        {1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824},
        {1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582},
        {1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960},
        {1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716},
        {1,22,231,1540,7315,26334,74613,170544,319770,497420,646646,705432},
        {1,23,253,1771,8855,33649,100947,245157,490314,817190,1144066,1352078},
        {1,24,276,2024,10626,42504,134596,346104,735471,1307504,1961256,2496144},
        {1,25,300,2300,12650,53130,177100,480700,1081575,2042975,3268760,4457400},
        {1,26,325,2600,14950,65780,230230,657800,1562275,3124550,5311735,7726160},
        {1,27,351,2925,17550,80730,296010,888030,2220075,4686825,8436285,13037895},
        {1,28,378,3276,20475,98280,376740,1184040,3108105,6906900,13123110,21474180},
        {1,29,406,3654,23751,118755,475020,1560780,4292145,10015005,20030010,34597290},
        {1,30,435,4060,27405,142506,593775,2035800,5852925,14307150,30045015,54627300},
        {1,31,465,4495,31465,169911,736281,2629575,7888725,20160075,44352165,84672315},
        {1,32,496,4960,35960,201376,906192,3365856,10518300,28048800,64512240,129024480},
        {1,33,528,5456,40920,237336,1107568,4272048,13884156,38567100,92561040,193536720},
        {1,34,561,5984,46376,278256,1344904,5379616,18156204,52451256,131128140,286097760},
        {1,35,595,6545,52360,324632,1623160,6724520,23535820,70607460,183579396,417225900},
        {1,36,630,7140,58905,376992,1947792,8347680,30260340,94143280,254186856,600805296},
        {1,37,666,7770,66045,435897,2324784,10295472,38608020,124403620,348330136,854992152},
        {1,38,703,8436,73815,501942,2760681,12620256,48903492,163011640,472733756,1203322288},
        {1,39,741,9139,82251,575757,3262623,15380937,61523748,211915132,635745396,1676056044},
        {1,40,780,9880,91390,658008,3838380,18643560,76904685,273438880,847660528,2311801440},
        {1,41,820,10660,101270,749398,4496388,22481940,95548245,350343565,1121099408,3159461968},
        {1,42,861,11480,111930,850668,5245786,26978328,118030185,445891810,1471442973,4280561376}
    };

__host__ __device__ __inline__ uint32_t 
_choose(const uint32_t n, const uint32_t k) 
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
choose(const uint32_t n, const uint32_t k) 
{
    assert (n <= 42);
    assert (k <= 11);
    if (n > 42 || k > 11) return _choose(n,k);

    return COMB_42_11_ [n][k];
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
