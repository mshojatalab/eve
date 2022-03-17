//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>
#include <eve/module/complex/regular/complex.hpp>
#include <eve/module/complex/cmplx.hpp>

namespace eve
{
  template<floating_scalar_value Type> struct complex;

  namespace detail
  {
    template<floating_real_value V> EVE_FORCEINLINE eve::complex<V> log_(EVE_SUPPORTS(cpu_),
                                                                    cmplx_type const &,
                                                                    V v) noexcept
    {
      return {log_abs(v), arg(v)};
    }
  }
}
