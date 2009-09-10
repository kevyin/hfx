/*
 * Module    : Kernels
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Functions which are internally implemented using CUDA, public interface.
 */

#ifndef __KERNELS_H__
#define __KERNELS_H__

#ifdef __cplusplus
extern "C" {
#endif

// -----------------------------------------------------------------------------
// Ion Series
// -----------------------------------------------------------------------------
void
k_buildXCorrSpecThry
(
    int                 charge,
    float               *b_ions,
    float               *y_ions,
    int                 *out,
    unsigned int        len_ions,
    unsigned int        len_spec
);


// -----------------------------------------------------------------------------
// Prelude
// -----------------------------------------------------------------------------
int  reducePlusf(float *xs, float N);

void scanl1Plusi(int *in, int *out, int N);
void scanr1Plusi(int *in, int *out, int N);

void zipWithPlusif(int *xs, float *ys, float *zs, int N);


#ifdef __cplusplus
}
#endif
#endif
/*
 * vim: syn=cuda
 */
