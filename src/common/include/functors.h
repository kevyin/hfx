#ifndef __FUNCTORS_H__ 
#define __FUNCTORS_H__

#include <thrust/tuple.h>
#include <thrust/functional.h>

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
struct minus_2tuple
{

    template <typename Tuple>
    __host__ __device__ 
    T operator() (Tuple x)
    {
        return thrust::get<0>(x) - thrust::get<1>(x);
    }
};

#endif
