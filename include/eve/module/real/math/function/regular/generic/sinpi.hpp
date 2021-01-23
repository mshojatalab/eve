//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright 2020 Joel FALCOU
  Copyright 2020 Jean-Thierry LAPRESTE

  Licensed under the MIT License <http://opensource.org/licenses/MIT>.
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/value.hpp>
#include <eve/constant/maxflint.hpp>
#include <eve/constant/one.hpp>
#include <eve/constant/pi.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/function/is_greater.hpp>
#include <eve/function/is_not_finite.hpp>
#include <eve/function/is_not_less_equal.hpp>
#include <eve/function/is_odd.hpp>
#include <eve/function/sin.hpp>
#include <eve/function/trigo_tags.hpp>
#include <eve/module/real/math/detail/generic/rem2.hpp>
#include <eve/module/real/math/detail/generic/trig_finalize.hpp>

#include <type_traits>

namespace eve::detail
{
  template<floating_real_value T>
  EVE_FORCEINLINE constexpr auto sinpi_(EVE_SUPPORTS(cpu_), restricted_type const &, T x) noexcept
  {
    if constexpr( has_native_abi_v<T> )
    {
      return eve::restricted(sin)(x * pi(eve::as<T>()));
    }
    else
      return apply_over(restricted(sinpi), x);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // medium,  big,  small
  template<decorator D, floating_real_value T>
  EVE_FORCEINLINE constexpr auto sinpi_(EVE_SUPPORTS(cpu_), D const &, T a0) noexcept
  {
    if constexpr( has_native_abi_v<T> )
    {
      if constexpr( scalar_value<T> )
        if( is_not_finite(a0) )
          return nan(eve::as<T>());
      auto x = eve::abs(a0);
      if constexpr( scalar_value<T> )
      {
        if( x > maxflint(eve::as<T>()) )
          return T(0);
      }
      else
      {
        x = if_else(is_not_finite(x), eve::allbits, x); // nan or Inf input
        x = if_else(is_greater(x, maxflint(eve::as(x))), eve::zero, x);
      }
      auto [fn, xr, dxr] = rem2(x);
      return sin_finalize(bitofsign(a0), fn, xr, dxr);
    }
    else
      return apply_over(D()(sinpi), a0);
  }

  template<floating_real_value T>
  EVE_FORCEINLINE constexpr auto sinpi_(EVE_SUPPORTS(cpu_), T const &a0) noexcept
  {
    if constexpr( has_native_abi_v<T> )
    {
      auto x = abs(a0);
      if( eve::all(eve::abs(x) <= T(0.25)) )
        return restricted(sinpi)(a0);
      else
        return big(sinpi)(a0);
    }
    else
      return apply_over(sinpi, a0);
  }

}
