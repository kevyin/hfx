/* -----------------------------------------------------------------------------
 *
 * Module    : MVM
 * Copyright : (c) [2012] Kevin Ying
 * License   : BSD
 *
 * Matrix vector multiplication using cublas
 * 
 * ---------------------------------------------------------------------------*/

#include "algorithms.h"

#include <stdint.h>
#include <thrust/device_vector.h>
#include "cublas_v2.h"

/* -----------------------------------------------------------------------------
 * Instances
 * ---------------------------------------------------------------------------*/

void
mvm_ff(cublasHandle_t handle, float *d_y, const float *d_A, const float *d_x, const uint32_t m, const uint32_t n)
{
#ifdef _BENCH
    cudaThreadSynchronize();
    time_t t_beg, t_end;
    time(&t_beg);
    std::cerr << "mvm_ff" << std::endl;
#endif
    /*cublasHandle_t handle;*/
    /*cublasCreate(&handle);*/
    float alpha = 1.0;
    float beta = 0.0;
    // Because cublas uses col major storage (as opposed to row major) swap row and col values and use CUBLAS_OP_T 
    /*thrust::device_ptr<const float> d_A_th(d_A);*/
    /*for (int i = 0; i < m; i++) {*/
        /*for (int j = 0; j < n; j++) {*/
            /*std::cerr << d_A_th[i*n + j] << " ";*/
        /*}*/
        /*std::cerr << std::endl;*/
    /*}*/
    //std::cerr << "m " << m << "n " << n << "alpha " << alpha << "beta " << beta << std::endl;
    cublasStatus_t status = cublasSgemv(handle, CUBLAS_OP_T, n, m, &alpha, d_A, n, d_x, 1, &beta, d_y, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "CUBLAS FAILURE" << std::endl;
    }

    /*cublasDestroy(handle);*/

#ifdef _BENCH
    cudaThreadSynchronize();
    time(&t_end);
    std::cerr<< "Time elapsed for mvm_ff: " << difftime(t_end,t_beg) << " seconds" << std::endl;
#endif

}
