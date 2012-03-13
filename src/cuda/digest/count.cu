/* -----------------------------------------------------------------------------
 *
 * Module    : Count
 * Copyright : (c) [2009..2012] Kevin Ying
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/unique.h>
#include <thrust/count.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "algorithms.h"
#include "texture.h"
#include "utils.h"

/* Helper functions and function objects
 * ---------------------------------------------------------------------------*/
template <typename T>
struct greaterThan : public thrust::unary_function<T,bool>
{
    T bound;

    __host__ __device__
    greaterThan(T _m) : bound(_m) {}

    __host__ __device__ bool operator() (T x)
    {
        return (bound < x);
    }
};

template <typename T>
struct increment: public thrust::unary_function<T,T>
{
    __host__ __device__
    increment() {}

    __host__ __device__ T operator() (T x)
    {
        return (++x);
    }
};


template <typename T>
struct targets: public thrust::unary_function<T, bool>
{
    const uint32_t rule;

    __host__ __device__
    targets (const uint32_t _r) : rule(_r) {}

    __host__ __device__ bool operator () (T elem)
    {
      switch (rule) {
        case 1:
          switch (elem) {
            case 'K': return true;
            case 'R': return true;
            default : return false;
          }
        case 2:
          switch (elem) {
            case 'F': return true;
            case 'W': return true;
            case 'Y': return true;
            default : return false;
          }
        case 3: 
          switch (elem) {
            case 'R': return true;
            default : return false; 
          } 
        case 4:
          switch (elem) {
            case 'M': return true;
            default : return false; 
          } 
        case 5: 
          switch (elem) {
            case 'W': return true;
            default : return false; 
          } 
        case 6: 
          switch (elem) {
            case 'P': return true;
            default : return false; 
          } 
        
        case 7:
          switch (elem) {
            case 'E': return true;
            default : return false; 
          } 
        case 8:
          switch (elem) {
            case 'K': return true;
            default : return false; 
          } 
        case 9:
          switch (elem) {
            case 'R': return true;
            default : return false; 
          } 
        case 10:
          switch (elem) {
            case 'D': return true;
            default : return false; 
          } 
        case 11:
          switch (elem) {
            case 'F': return true;
            case 'W': return true;
            case 'Y': return true;
            case 'L': return true;
            default : return false;
          }
        case 12:
          switch (elem) {
            case 'A': return true;
            case 'L': return true;
            case 'I': return true;
            case 'V': return true;
            default : return false;
          }
        case 13:
          switch (elem) {
            case 'A': return true;
            case 'L': return true;
            case 'I': return true;
            case 'V': return true;
            case 'K': return true;
            case 'R': return true;
            case 'F': return true;
            case 'W': return true;
            case 'Y': return true;
            default : return false;
          }
        default:
          assert(!"Digestion rule not recognised");
      }
      return false;
    }
};

template <typename T_in, typename T_out>
struct calcResMass: public thrust::unary_function<T_in,T_out>
{
    const uint8_t   *d_ions;
    const float     *d_mass;
    const uint32_t  *d_c;
    const uint32_t  *d_n;

    __host__ __device__
    calcResMass(const uint8_t   *_di,
                const float     *_dm,
                const uint32_t  *_dc,
                const uint32_t  *_dn) : 
                d_ions(_di), d_mass(_dm), d_c(_dc), d_n(_dn) {}

    __host__ __device__ T_out operator() (T_in i)
    {
        float mass = 0;
        uint32_t beg = d_c[i];
        uint32_t end = d_n[i];
        for (uint32_t k = beg; k <= end; k++) {
            mass += d_mass[d_ions[k] - 'A'];
        }
        return mass;
    }
};

/* Functions
 * ---------------------------------------------------------------------------*/

uint32_t
countFragByRule
(
    const uint8_t       *d_in_raw,
    const uint32_t      N,
    const uint32_t      ruleNo
)
{

    thrust::device_ptr<const uint8_t>     d_in(d_in_raw);
    uint32_t count = 0;
    count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(ruleNo));
    //switch (ruleNo) {
      //case 0: count = 0; break;
      //case 1: count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(1)); break;
      //case 2: count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(2)); break;
      //case 3: count = thrust::count(d_in, d_in + N, 'R'); break;
      //case 4: count = thrust::count(d_in, d_in + N, 'M'); break;
      //case 5: count = thrust::count(d_in, d_in + N, 'W'); break;
      //case 6: count = thrust::count(d_in, d_in + N, 'P'); break;
      //case 7: count = thrust::count(d_in, d_in + N, 'E'); break;
      //case 8: count = thrust::count(d_in, d_in + N, 'K'); break;
      //case 9: count = thrust::count(d_in, d_in + N, 'R'); break;
      //case 10: count = thrust::count(d_in, d_in + N, 'D'); break;
      //case 11: count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(11)); break;
      //case 12: count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(12)); break;
      //case 13: count = thrust::count_if(d_in, d_in + N, targets<uint8_t>(13)); break;
      //default : 
          //std::cerr << "Digestion rule not recognised\n";
          //exit(1);
    //}
    return count;
}


uint32_t
digestByRule
(
    const uint8_t       *d_ions_raw,
    const uint32_t      Ni,             // Number of ions
    const uint32_t      *d_iseg_raw,    // protein segments in ions
    const uint32_t      Ns,             // number of segments (proteins)
    const uint32_t      Nf,             // upper bound on protein frags
    const float         *d_mass_raw,    // amino acid masses
    const uint32_t      ruleNo,

    // result vectors
    float         *d_r_raw,
    uint32_t      *d_c_raw,
    uint32_t      *d_n_raw,
    uint32_t      *d_fs_raw
)
{
    // init
    thrust::device_ptr<const uint8_t>  d_ions(d_ions_raw);
    thrust::device_ptr<const uint32_t> d_iseg(d_iseg_raw);
    thrust::device_ptr<uint32_t> d_fs(d_fs_raw);
    thrust::device_ptr<float> d_r(d_r_raw);
    thrust::device_ptr<uint32_t> d_c(d_c_raw);
    thrust::device_ptr<uint32_t> d_n(d_n_raw);
    
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last;

    // Allocate array for separation sites indicating the begining of fragments
    // The first part are sites from different proteins
    // the second part are cleavage sites in accordance to the ruleNo
    thrust::device_vector<uint32_t> d_sites(Nf+Ns); // leave extra room for seqment information
    thrust::copy(d_iseg, d_iseg + Ns, d_sites.begin());

    // find cleavage sites
    last = first + Ni;
    thrust::device_ptr<uint32_t> d_sites_end = 
        thrust::copy_if(first, last, d_ions, d_sites.data() + Ns, targets<uint8_t>(ruleNo));

    // sort separation points by key
    // important that sort is stable as unique takes the first common element,
    // and we want to always keep the protein segment points over the cleavage point
    //
    // To keep note which sites are from proteins, keep an array of markers
    thrust::device_vector<uint8_t> d_seg_markers(thrust::distance(d_sites.data(), d_sites_end));
    thrust::fill_n(d_seg_markers.begin(), Ns, 1);
    thrust::stable_sort_by_key(d_sites.data(), d_sites_end, d_seg_markers.data());
    thrust::pair< thrust::device_ptr<uint32_t>, thrust::device_ptr<uint8_t> > new_ends = 
        thrust::unique_by_key(d_sites.data(), d_sites_end, d_seg_markers.data());

    uint32_t seg_cnt = thrust::count(d_seg_markers.begin(), d_seg_markers.end(), 1);
    //if (seg_cnt != Ns) {
      //std::cerr << "ERROR1: lost some segments! " << seg_cnt << " != " << Ns << std::endl;
      //exit(1);
    //}


    // fill c n and then r
    uint32_t Nf_new;
    thrust::device_ptr<uint32_t> d_c_end;
    thrust::device_ptr<uint32_t> d_n_end;
    if (d_sites[0] == 0) {
      Nf_new = new_ends.first - d_sites.data() - 1;
      d_n_end = thrust::copy(d_sites.data() + 1, new_ends.first, d_n);
      d_c_end = thrust::copy(thrust::make_transform_iterator(d_sites.data() + 1, increment<uint32_t>()),
                             thrust::make_transform_iterator(new_ends.first, increment<uint32_t>()),
                             d_c + 1);
    } else {
    //Nf_new = thrust::distance(d_sites.data(), new_ends.first);
      Nf_new = new_ends.first - d_sites.data();
      d_n_end = thrust::copy(d_sites.data(), new_ends.first, d_n);
      d_c_end = thrust::copy(thrust::make_transform_iterator(d_sites.data(), increment<uint32_t>()),
                             thrust::make_transform_iterator(new_ends.first, increment<uint32_t>()),
                             d_c + 1);
    }
    d_c[0] = 0;
    d_n[Nf_new-1] = Ni-1;

    std::cout << "rule " << ruleNo << std::endl;
    std::cout << d_ions[0] << std::endl;
    std::cout << d_c[0] << std::endl;
    std::cout << d_n[0] << std::endl;
    std::cout << d_c[1] << std::endl;
    std::cout << d_n[1] << std::endl;
    std::cout << d_c[Nf_new-1] << std::endl;
    std::cout << d_n[Nf_new-1] << std::endl;
    std::cout << d_c[Nf_new-2] << std::endl;
    std::cout << d_n[Nf_new-2] << std::endl;
    last = first + Nf_new;
    thrust::transform(first, last, d_r, calcResMass<uint32_t,float>(d_ions_raw, d_mass_raw, d_c_raw, d_n_raw));

    std::cout << "here3.5" << std::endl;
    // fill protein segment starting locations in the list fragments.
    thrust::device_ptr<uint32_t> d_fs_end = 
      thrust::copy_if(first, last, d_seg_markers.data(), d_fs, greaterThan<uint8_t>(0)); 

    std::cout << "here4" << std::endl;
    //if (thrust::distance(d_fs, d_fs_end) != Ns) {
      //std::cerr << "ERROR2: lost some segments!" << std::endl;
      //exit(2);
    //}
    return Nf_new;
}
