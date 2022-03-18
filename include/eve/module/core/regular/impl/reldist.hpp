//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/module/core/regular/if_else.hpp>
#include <eve/module/core/regular/all.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/module/core/regular/dist.hpp>
#include <eve/concept/value.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/skeleton_calls.hpp>

namespace eve::detail
{
  // -----------------------------------------------------------------------------------------------
  // regular case
  template<real_value T, real_value U>
  EVE_FORCEINLINE auto reldist_(EVE_SUPPORTS(cpu_)
                            , T const &a
                            , U const &b) noexcept
  requires compatible_values<T, U>
  {
    return arithmetic_call(reldist, a, b);
  }

  template<real_value T>
  EVE_FORCEINLINE T reldist_(EVE_SUPPORTS(cpu_)
                         , T const &a
                         , T const &b) noexcept
  requires has_native_abi_v<T>
  {
    if constexpr(integral_value<T>)
    {
      return dist(a, b);
    }
    else if constexpr(scalar_value<T>)
    {
      if((a == b) || (is_nan(a) && is_nan(b)))
        return 0.;
      else if(is_infinite(a) || is_infinite(b) || is_nan(a) || is_nan(b))
        return inf(as(a));
      else return 100. * (dist(a, b) / max(T(1), abs(a), abs(b)));
    }
    else if constexpr(simd_value<T>)
    {
       auto rd = 100. * (dist(a, b) / max(T(1), abs(a), abs(b)));
       rd = if_else((a == b) || (is_nan(a) && is_nan(b)), zero, rd);
       return if_else(is_infinite(a) || is_infinite(b) || is_nan(a) || is_nan(b), inf(as(a)), rd);
    }
  }
}
