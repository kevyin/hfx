/* -----------------------------------------------------------------------------
 *
 * Module    : Utils
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <thrust/device_vector.h>

/* Taken from NVIDIA CUDA C Programming Guide 4.0 : B.14.4
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    #define printf(f, ...) ((void) (f, __VA_ARGS__),0)
#endif
*/

/*
 * Core assert function. Don't let this escape...
 */
#if defined(__CUDACC__) || !defined(__DEVICE_EMULATION__)
#define __assert(e, file, line) ((void)0)
#else
#define __assert(e, file, line) \
    ((void) fprintf (stderr, "%s:%u: failed assertion `%s'\n", file, line, e), abort())
#endif

/*
 * Test the given expression, and abort the program if it evaluates to false.
 * Only available in debug mode.
 */
#ifndef _DEBUG
#define assert(e)               ((void)0)
#else
#define assert(e)  \
    ((void) ((e) ? (void(0)) : __assert (#e, __FILE__, __LINE__)))
#endif


/*
 * Macro to insert __syncthreads() in device emulation mode
 *
 * In emulation mode, this is required inside warp-synchronous code, to recreate
 * the behaviour of the warp threads executing in lock-step on the GPU.
 */
#ifdef __DEVICE_EMULATION__
#define __EMUSYNC               __syncthreads()
#else
#define __EMUSYNC
#endif


/*
 * Check the return status of CUDA API calls, and abort with an appropriate
 * error string on failure.
 */
#define CUDA_SAFE_CALL_NO_SYNC(call)                                           \
    do {                                                                       \
        cudaError err = call;                                                  \
        if(cudaSuccess != err) {                                               \
            const char *str = cudaGetErrorString(err);                         \
            __assert(str, __FILE__, __LINE__);                                 \
        }                                                                      \
    } while (0)

#define CUDA_SAFE_CALL(call)                                                   \
    do {                                                                       \
        CUDA_SAFE_CALL_NO_SYNC(call);                                          \
        CUDA_SAFE_CALL_NO_SYNC(cudaThreadSynchronize());                       \
    } while (0)

/*
 * get and set acid chars with range('A'-'Z') and mods with range(0-7)
 * Expects data to be packed into 8 bit datatype eg uint8_t
 * (from least significant to most significant) 
 *  The first 5 bits are for the numbers 0-25 corresponding to the alphabet
 *  The last 3 bits are for the numbers 0-7 and represent a modification
 * 
 * The acid character is a letter representing the amino acid 'A' - 'Z'
 * The acid character idx is a number 0-25 representing the above letter (acid)
 *
 * The acid mod is a number 0-7
 *      where 0 is no modification
 *      1-7 is the 1st, 2nd .. 7th modification
 */
#define GET_ACID_CHAR(i)        ((i & 0x1FU) + 'A')
#define GET_ACID_CHAR_IDX(i)    ((i & 0x1FU))
#define SET_ACID_CHAR(i,c)      i = (i & 0xE0U) | (c - 'A'); 
#define SET_ACID_CHAR_IDX(i,ci) i = (i & 0xE0U) | ci; 

#define GET_ACID_MOD(i)         ((i & 0xE0U) >> 5)
#define SET_ACID_MOD(i,m)       i = (i & 0x1FU) | m << 5; 

#define SET_ACID_CHAR_MOD(i,c,m) i = (c - 'A') | m << 5;
#define SET_ACID_CHAR_IDX_MOD(i,ci,m) i = ci | m << 5;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Determine if the input is a power of two
 */
inline bool
isPow2(unsigned int x)
{
    return ((x&(x-1)) == 0);
}


/*
 * Compute the next highest power of two
 */
inline unsigned int
ceilPow2(unsigned int x)
{
#if 0
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
#endif

    return (isPow2(x)) ? x : 1u << (int) ceil(log2((double)x));
}


/*
 * Compute the next lowest power of two
 */
inline unsigned int
floorPow2(unsigned int x)
{
#if 0
    float nf = (float) n;
    return 1 << (((*(int*)&nf) >> 23) - 127);
#endif

    int exp;
    frexp((double)x, &exp);
    return 1 << (exp - 1);
}

inline
void printPeptide(const uint32_t& idx, const uint8_t* d_ions_raw, const uint32_t* d_tc_raw, const uint32_t* d_tn_raw) {
    thrust::device_ptr<const uint8_t>     d_ions(d_ions_raw);
    thrust::device_ptr<const uint32_t>    d_tc(d_tc_raw);
    thrust::device_ptr<const uint32_t>    d_tn(d_tn_raw);

    uint32_t begin = d_tc[idx];
    uint32_t end   = d_tn[idx];
    for (uint32_t i = begin; i < end; ++i) {
        uint8_t a = d_ions[i];
        printf("%c ", a);
    }
    printf("\n");


}

inline
void printGPUMemoryUsage ()
{
    // show memory usage of GPU
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
    if ( cudaSuccess != cuda_status ){
        printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
        exit(1);
    }

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

}

#undef __asert

#ifdef __cplusplus
}
#endif
#endif

