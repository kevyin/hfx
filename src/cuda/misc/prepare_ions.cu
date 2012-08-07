/* -----------------------------------------------------------------------------
 *
 * Module    : Prepare Ions 
 * Copyright : (c) [2012] Kevin Ying 
 * License   : BSD
 *
 * Converts the 8 bits used to represent ions (originally uint8_t characters) 
 * to pack in modfication information
 * The getter and setter macros eg SET_ACID_CHAR_MOD allow access to the 
 * acid or modifcation idx that acid can be applied to.
 * Details on the macros and how the 8 bits are separated is in utils.h
 * ---------------------------------------------------------------------------*/

#include <stdint.h>

#include "utils.h"
#include "device.h"
#include "algorithms.h"
#include "prepare_ions.h"

#ifdef _TEST
#include "tests.h"
#include <thrust/for_each.h>
#include <thrust/device_vector.h>
#endif



__device__ __inline__ uint8_t 
check_ion_for_ma (uint8_t ion, const uint8_t* d_ma, const uint32_t num_ma)
{
    for (size_t j = 0; j < num_ma; ++j)
    {
        if (d_ma[j] == ion) 
        {
            return j+1;
        }
    }
    return 0;
}

__device__ __inline__ void
set_ions (uint8_t* ions, const uint8_t* d_ma, const uint32_t num_ma, const uint32_t thread_idx)
{
    #pragma unroll
    for (uint32_t i = 0; i < NUM_PER_THREAD; ++i)
    {
        //uint8_t char_idx = ions[i] - 'A';
        uint8_t mod = check_ion_for_ma(ions[i], d_ma, num_ma);
        SET_ACID_CHAR_MOD(ions[i],ions[i],mod);
    }
}

template <bool UseCache>
__global__ static void
prepare_ions_core
(
    uint8_t         *d_ions,    // Flattened array of ions range('A'-'Z')
    const uint32_t  num_ions,
    const uint8_t   *d_ma,      // modable acids range('A'-'Z')
    const uint32_t  num_ma
)
{
    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t thread_id       = blockDim.x * blockIdx.x + threadIdx.x;
    uint8_t *ions;

    if (thread_id < num_ions / NUM_PER_THREAD)
    {
        ions = d_ions + (thread_id * NUM_PER_THREAD);
        set_ions(ions,d_ma,num_ma,thread_id);
    } 

    // set remaining elements with the first block
    // remaining elements are total number - number done
    uint32_t done = (gridDim.x * blockDim.x * NUM_PER_THREAD);
    uint32_t remaining = (num_ions - done); 
    assert(remaining < blockDim.x * NUM_PER_THREAD); // otherwise too few blocks assigned

    if (thread_id < (remaining / NUM_PER_THREAD))
    {
        ions = d_ions + done + (thread_id * NUM_PER_THREAD);
        set_ions(ions,d_ma,num_ma,thread_id);
    }
    uint32_t remainder = num_ions % NUM_PER_THREAD;
    if (remainder && thread_id == (remaining / NUM_PER_THREAD)) 
    {
        ions = d_ions + done + (thread_id * NUM_PER_THREAD);
        for (uint32_t i = 0; i < remainder; ++i) 
        {
            uint8_t mod = check_ion_for_ma(ions[i], d_ma,num_ma);
            SET_ACID_CHAR_MOD(ions[i],ions[i],mod);
        }
    }
}

static void
prepare_ions_control(uint32_t N, uint32_t &blocks, uint32_t &threads)
{
    threads = (N < MAX_THREADS) ? max(WARP_SIZE, ceilPow2(N)) : MAX_THREADS;
    blocks = (N / (threads * NUM_PER_THREAD));
}

void prepare_ions
(
    uint8_t* d_ions, 
    uint32_t N, 
    uint8_t* d_ma, 
    uint32_t num_ma
) 
{
#ifdef _TEST
    std::cerr << "Testing: prepare_ions" << std::endl;
    thrust::device_vector<uint8_t> d_ions_th = prepare_ions_test(d_ions, N, d_ma, num_ma);
    thrust::device_ptr<uint8_t> d_ions_(d_ions);
#endif
    uint32_t blocks;
    uint32_t threads;

    // 4 elements at a time
    prepare_ions_control(N,blocks,threads);
    prepare_ions_core<false><<<blocks,threads>>>(d_ions, N, d_ma, num_ma);

#ifdef _TEST
    bool eq = equal(d_ions_th.begin(), d_ions_th.end(), thrust::device_ptr<uint8_t>(d_ions));
    if (!eq) {
        std::cerr << "TEST: prepare_ions failed" << std::endl;
        exit(1);
    }
#endif
}
