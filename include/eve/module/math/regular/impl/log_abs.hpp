//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/module/core.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math/regular/log.hpp>

namespace eve::detail
{
  template<typename T, decorator D> //TODO change typename when complex are value or satisfy another concept
  EVE_FORCEINLINE constexpr auto log_abs_(EVE_SUPPORTS(cpu_), D const &, T x) noexcept
  requires(is_one_of<D>(types<regular_type, pedantic_type> {}))
  {
    return D()(eve::log)(eve::abs(x));
  }

  template<typename T> //TODO change typename when complex are value or satisfy another concept
  EVE_FORCEINLINE constexpr auto log_abs_(EVE_SUPPORTS(cpu_), T x) noexcept
  {
    return eve::log(eve::abs(x));;
  }
}
