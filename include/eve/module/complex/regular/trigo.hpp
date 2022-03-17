//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once
#include <eve/module/complex.hpp>

namespace eve
{
  template<like<complex> Z>
  EVE_FORCEINLINE auto tagged_dispatch( eve::tag::cosh_, Z const& a0 ) noexcept
  {
    using v_t = Z::value_type;
    auto ra0 = real(a0);
    auto ia0 = imag(a0);
    auto [s, c]   = sincos(ia0);
    auto [sh, ch] = sinhcosh(ra0);

    auto r = c*ch;
    auto i = s*sh;
    i = if_else(is_imag(a0) || is_real(a0), Z(0), i);
    auto res = Z(r, i);
    if (eve::none(is_invalid(a0))) return res;
    res = if_else(is_inf(ra0) && is_invalid(ia0)), Z(inf(as(ra0), nan(as(ra0)), res);
    res = if_else(is_nan(ra0) && is_inf(ia0)),     Z(nan(as(ra0), nan(as(ra0)), res);
    return res;
  }
