//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>
#include <eve/module/core.hpp>

namespace eve
{

  namespace detail
  {
    template<floating_real_value V>
    EVE_FORCEINLINE auto sqrt_(EVE_SUPPORTS(cpu_), cmplx_type const &, V v) noexcept
    {
     if constexpr(scalar_value<V>)
      {
        using c_t = complex<decltype(v)>;
        return if_else(is_positive(v), c_t{sqrt(v), zero(as(v))}, c_t{zero(as(v)), sqrt(-v)});
      }
      else
      {
        using elt_t = element_type_t<V>;
        using c_t = eve::wide<eve::complex<elt_t>, eve::cardinal_t<V>>;
        return if_else(is_positive(v), c_t{sqrt(v), zero(as(v))}, c_t{zero(as(v)), sqrt(-v)});
      }
    }
  }
}
