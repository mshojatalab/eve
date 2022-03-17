//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>

namespace eve
{
  //================================================================================================
  //! @addtogroup complex
  //! @{
  //! @var mul_i
  //!
  //! @brief Callable object computing mul_i part of values.
  //!
  //! **Required header:** `#include <eve/module/complex.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the computation of mul_i part                               |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  auto operator()(value auto x) const noexcept;
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`x`:   [value](@ref eve::value).
  //!
  //! **Return value**
  //! `complex < decltype(x)> equal to x multiplied by i (i*i = -1).
  //!
  //! #### Example
  //!
  //! @godbolt{doc/complex/mul_i.cpp}
  //!
  //!  @}
  //================================================================================================

  namespace tag { struct mul_i_; }
  template<> struct supports_conditional<tag::mul_i_> : std::false_type {};

  EVE_MAKE_CALLABLE(mul_i_, mul_i);

  namespace detail
  {
    template<floating_value V> EVE_FORCEINLINE
    auto mul_i_(EVE_SUPPORTS(cpu_), V v) noexcept
    {
      if constexpr(scalar_value<V>)
      {
        using c_t = eve::complex < V>;
        return c_t{zero(as(v)), v};
      }
      else
      {
        using elt_t = element_type_t<V>;
        using c_t = eve::wide<eve::complex<elt_t>, eve::cardinal_t<V>>;
        return c_t{zero(as(v)), v};
      }
    }
  }
}
