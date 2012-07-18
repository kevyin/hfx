/*
 *  Copyright 2009 Advanced Industrial Science and Technology (AIST).
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file
 *  \brief
 */
#pragma once

#include <thrust/functional.h>

namespace thrust
{

  // --------------------------------------------------------------------------
  // binder1st, bind1st
  // --------------------------------------------------------------------------
  template <typename Operation>
  class binder1st
    : public thrust::unary_function<typename binary_traits<Operation>::second_argument_type,
                                    typename binary_traits<Operation>::result_type>
  {
  protected:
    typename binary_traits<Operation>::function_type m_op;
    typename binary_traits<Operation>::first_argument_type m_value;

  public:
    __host__ __device__
    binder1st( typename binary_traits<Operation>::param_type op,
               typename binary_traits<Operation>::first_argument_type value)
      : m_op(op), m_value(value)
    {}

    __host__ __device__
    typename binary_traits<Operation>::result_type
    operator()( typename binary_traits<Operation>::second_argument_type y) const
    {
      return m_op( m_value, y);
    }
  };
  
  template <typename Operation>
  __host__ __device__
  inline binder1st<Operation>
  bind1st( const Operation& op,
           typename binary_traits<Operation>::first_argument_type x)
  {
    return binder1st<Operation>( static_cast<typename binary_traits<Operation>::param_type>(op), 
                                 x);
  }
  
  template <typename Operation>
  __host__ __device__
  inline binder1st<Operation>
  bind1st( Operation& op, 
           typename binary_traits<Operation>::first_argument_type x)
  {
    return binder1st<Operation>( op, x);
  }
  
  // --------------------------------------------------------------------------
  // binder2nd, bind2nd
  // --------------------------------------------------------------------------
  template <typename Operation>
  class binder2nd
    : public thrust::unary_function<typename binary_traits<Operation>::first_argument_type,
                                    typename binary_traits<Operation>::result_type>
  {
  protected:
    typename binary_traits<Operation>::function_type m_op;
    typename binary_traits<Operation>::second_argument_type m_value;
    
  public:
    __host__ __device__
    binder2nd( typename binary_traits<Operation>::param_type op,
               typename binary_traits<Operation>::second_argument_type value)
      : m_op(op), m_value(value)
    {}
    
    __host__ __device__
    typename binary_traits<Operation>::result_type
    operator()( typename binary_traits<Operation>::first_argument_type x) const
    {
      return m_op( x, m_value);
    }
  };
  
  template <typename Operation>
  __host__ __device__
  inline binder2nd<Operation>
  bind2nd( const Operation& op,
           typename binary_traits<Operation>::second_argument_type x)
  {
    return binder2nd<Operation>( static_cast<typename binary_traits<Operation>::param_type>(op),
                                 x);
  }
  
  template <typename Operation>
  __host__ __device__
  inline binder2nd<Operation>
  bind2nd( Operation& op,
           typename binary_traits<Operation>::second_argument_type x)
  {
    return binder2nd<Operation>( op, x);
  }
  
} // end namespace thrust
