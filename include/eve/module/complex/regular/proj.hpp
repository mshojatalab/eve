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

namespace eve
{
  //================================================================================================
  //! @addtogroup complex
  //! @{
  //! @var proj
  //!
  //! @brief Callable object computing proj(x).
  //!
  //! **Required header:** `#include <eve/module/complex.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the computation of proj(x)                                |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  auto operator()(value auto x) const noexcept;
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`z`:   [value](@ref eve::value).
  //!
  //! **Return value**
  //! complex type associted to z
  //!
  //! For most z, std::proj(z)==z, but all complex infinities, even the numbers where one component
  //! is infinite and the other is NaN, become positive real infinity, (inf, 0.0) or (inf, -0.0).
  //! The sign of the imaginary (zero) component is the sign imag(z).
  //!
  //! #### Example
  //!
  //! @godbolt{doc/complex/proj.cpp}
  //!
  //!  @}
  //================================================================================================

  namespace tag { struct proj_; }
  template<> struct supports_conditional<tag::proj_> : std::false_type {};

  EVE_MAKE_CALLABLE(proj_, proj);

  namespace detail
  {
    template<floating_value V> EVE_FORCEINLINE auto proj_(EVE_SUPPORTS(cpu_), V v) noexcept
    {
     if constexpr(scalar_value<V>)
     {
       using c_t = eve::complex < V>;
       return if_else(is_infinite(v), c_t{inf(as(v)), zero(as(v))}, c_t(v));
     }
     else
     {
       using elt_t = element_type_t<V>;
       using c_t = eve::wide<eve::complex<elt_t>, eve::cardinal_t<V>>;
       return if_else(is_infinite(v), c_t{inf(as(v)), zero(as(v))}, c_t(v));
     }
    }
  }
}
