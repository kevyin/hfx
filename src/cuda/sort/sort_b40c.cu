/* -----------------------------------------------------------------------------
 *
 * Module    : Sort using back40computing
 * Modifications by Kevin Ying
 * 
 * Original code:
 * Copyright 2010-2012 Duane Merrill
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a scan of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 * 
 * For more information, see our Google Code project site: 
 * http://code.google.com/p/back40computing/
 * 
 * ---------------------------------------------------------------------------*/

#include <stdint.h> 
#include <stdio.h> 
// Sorting includes
#include <b40c/radix_sort/enactor.cuh>
#include <b40c/util/multiple_buffering.cuh>

#include "algorithms.h"

void sort_b40c_f(float *d_keys, uint32_t *d_vals, uint32_t N)
{
    b40c::radix_sort::Enactor enactor;
    b40c::util::DoubleBuffer<float,uint32_t> sort_storage(d_keys,d_vals);

    // use SMALL_SIZE for sorts with < 1M elements
    enactor.Sort(sort_storage, N);
}
