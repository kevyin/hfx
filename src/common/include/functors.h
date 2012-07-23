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

    __host__ __device__ 
    bool operator() (T x)
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

template <typename T>
struct mat_wider: public thrust::unary_function<size_t,T>
{
    T* col;
    size_t wide;

    __host__ __device__
    mat_wider(T* _c, size_t _w) : col(_c), wide(_w) {}

    __host__ __device__ 
    T operator() (size_t i)
    {
        return col[i/wide];
    }
};

#endif
