/* -----------------------------------------------------------------------------
 *
 * Module    : Tests
 * Copyright : (c) [2012] Kevin Ying 
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __TESTS_H__
#define __TESTS_H__

#include <stdint.h>

thrust::device_vector<uint8_t> prepare_ions_test(uint8_t* _d_ions, uint32_t N, uint8_t* d_ma, uint32_t num_ma); 

#endif
