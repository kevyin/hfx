/* -----------------------------------------------------------------------------
 *
 * Module    : Prepare Ions tests
 * Copyright : (c) [2012] Kevin Ying 
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/for_each.h>
#include <stdint.h>

#include "utils.h"
#include "device.h"
#include "tests.h"
#include "algorithms.h"
#include "prepare_ions.h"

using namespace thrust;

struct prepare_ion : public thrust::unary_function<uint8_t, uint8_t>
{
    uint8_t* d_ma;
    uint32_t num_ma;

    __host__ __device__
    prepare_ion(uint8_t* _m, uint32_t _nm) : d_ma(_m), num_ma(_nm) {}

    __host__ __device__
    uint8_t operator() (uint8_t orig_ion)
    {
        uint8_t res = 0;
        uint8_t mod = 0;
        for (size_t i = 0; i < num_ma; ++i)
        {
            if (orig_ion == d_ma[i])
            {
                mod = i+1;
            }
        }
        SET_ACID_CHAR(res, orig_ion);
        SET_ACID_MOD(res, mod);
        return res;
    }
};


device_vector<uint8_t> prepare_ions_test(uint8_t* _d_ions, uint32_t N, uint8_t* d_ma, uint32_t num_ma) 
{
    device_ptr<uint8_t> d_ions(_d_ions);
    device_vector<uint8_t> d_res(N);

    thrust::transform(d_ions, d_ions+N, d_res.begin(), prepare_ion(d_ma, num_ma));
                                    
    return d_res;
}

