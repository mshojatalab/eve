//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>
#include <eve/module/math.hpp>
// #include <eve/module/complex/regular/complex.hpp>
// #include <eve/module/complex/cmplx.hpp>

namespace eve
{
  template<floating_scalar_value Type> struct complex;

  namespace detail
  {
    template<floating_real_value V>
    EVE_FORCEINLINE auto log_(EVE_SUPPORTS(cpu_),
                                         cmplx_type const &,
                                         V v) noexcept
    {
     if constexpr(scalar_value<V>)
      {
        using c_t = complex<decltype(v)>;
        return c_t{log_abs(v), arg(v)};
      }
      else
      {
        using elt_t = element_type_t<V>;
        using c_t = eve::wide<eve::complex<elt_t>, eve::cardinal_t<V>>;
        return c_t{log_abs(v), arg(v)};
      }
    }
  }
}
