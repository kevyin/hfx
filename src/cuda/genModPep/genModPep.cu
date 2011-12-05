/* -----------------------------------------------------------------------------
 *
 * Module    : Generate Modified Peptides
 * Copyright : (c) [2009..2011] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <stdint.h>

#include "utils.h"
#include "algorithms.h"
#include "combinatorics.h"

template <typename T>
struct calcNumModCand : public thrust::unary_function<T, uint32_t>
{

    thrust::device_ptr<const uint32_t>    d_pep_ma_count;
    thrust::device_ptr<const uint8_t>     d_ma_count;
    const uint32_t                        ma_len;

    __host__ __device__
    calcNumModCand(thrust::device_ptr<const uint32_t>    _pep_ma_c, 
                   thrust::device_ptr<const uint8_t>     _ma_c, 
                   const uint32_t                        len) : 
                   d_pep_ma_count(_pep_ma_c), d_ma_count(_ma_c), ma_len(len) {}

    __host__ __device__ uint32_t operator() (T pep)
    {
        thrust::device_ptr<const uint32_t> pep_ma_count (d_pep_ma_count + pep*ma_len);
        uint32_t N = 1;
        for (int i = 0; i < ma_len; i++)
        {
            N *= choose(pep_ma_count[i], d_ma_count[i]);    
        }
        return N;
    }
};

uint32_t
calcTotalModCands
(
    uint32_t          *d_out_pep_mpep_count,
    const uint32_t    nPep,                   // number of peptides
    const uint32_t    *d_pep_ma_count,        // 2d array of ma count in each pep
    const uint8_t     *d_ma_count,             
    const uint32_t    ma_length
) 
{
    thrust::device_ptr<uint32_t>        d_out_pep_mpep_count_th(d_out_pep_mpep_count);
    thrust::device_ptr<const uint32_t>  d_pep_ma_count_th(d_pep_ma_count);
    thrust::device_ptr<const uint8_t>   d_ma_count_th(d_ma_count);
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + nPep;

    thrust::transform(first, last, d_out_pep_mpep_count_th, calcNumModCand<const uint32_t>(d_pep_ma_count_th, d_ma_count_th, ma_length));
    
    // print
    /*std::cout << "numModCand " << std::endl;*/
    /*for (int i = 0; i < nPep; i++)*/
    /*{*/
        /*std::cout << d_out_pep_mpep_count_th[i] << std::endl;*/
    /*}*/

    uint32_t total = thrust::reduce(d_out_pep_mpep_count_th, d_out_pep_mpep_count_th + nPep);
    std::cout << "total " << total << std::endl;

    return total;
}

void
genModCands
(                                                                     
    uint32_t        *d_out_mpep_idx,
    uint32_t        *d_out_mpep_mcomb,
    const uint32_t  total,
    const uint32_t  *d_pep_idx,
    const uint32_t  *d_pep_mpep_count,
    const uint32_t  num_pep
)
{
    std::cout << "genModCands" << std::endl;
    thrust::device_ptr<uint32_t> d_out_mpep_mcomb_th(d_out_mpep_mcomb);
    thrust::device_ptr<uint32_t> d_out_mpep_idx_th(d_out_mpep_idx);
    thrust::device_ptr<const uint32_t> d_pep_idx_th(d_pep_idx);
    thrust::device_ptr<const uint32_t> d_pep_mpep_count_th(d_pep_mpep_count);
    std::vector<uint32_t> mcomb;
    std::vector<uint32_t> mpep_idx;
    for (uint32_t i = 0; i < num_pep; i++)
    {
        for (uint32_t k = 0; k < d_pep_mpep_count_th[i]; k++)
        {
            mcomb.push_back(k);
            mpep_idx.push_back(d_pep_idx_th[i]);
        }
    }
    assert (mcomb.size() == total);

    thrust::copy(mcomb.begin(), mcomb.end(), d_out_mpep_mcomb_th);
    thrust::copy(mpep_idx.begin(), mpep_idx.end(), d_out_mpep_idx_th);

    // print
    std::cout << "d_out_mpep_mcomb" << std::endl;
    for (uint32_t i = 0; i < total; i++) {
        std::cout << d_out_mpep_idx_th[i] << " " << d_out_mpep_mcomb_th[i] << std::endl;
    }
    std::cout << std::endl;
    
}
