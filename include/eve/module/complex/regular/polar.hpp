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
  //! @addtogroup math
  //! @{
  //! @var polar
  //!
  //! @brief Callable object computing the polarugate value.
  //!
  //! **Required header:** `#include <eve/module/math.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the computation of a complex from a module argument pair   |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  auto operator()(floating_value auto rho, auto theta) const noexcept;
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`rho`:  the modulus of the complex
  //!`theta` the argument in radian of the complex
  //!
  //! **Return value**
  //!
  //! the complex value
  //!
  //! ---
  //!
  //! #### Example
  //!
  //! @godbolt{doc/complex/polar.cpp}
  //!
  //!  @}
  //================================================================================================
  EVE_MAKE_CALLABLE(polar_, polar);

  namespace detail
  {
    template<floating_value V,  floating_value U> EVE_FORCEINLINE auto polar_(EVE_SUPPORTS(cpu_), V rho, U theta) noexcept
    {
      return arithmetic_call(polar, rho, theta);
    }

    template<floating_value U> EVE_FORCEINLINE auto polar_(EVE_SUPPORTS(cpu_), U rho, U theta) noexcept
    {
      using elt_t = element_type_t<U>;
      using c_t = eve::wide<eve::complex<elt_t>, eve::cardinal_t<U>>;
      auto [s, c] = sincos(theta);
      return c_t{rho*c, rho*s};
    }
  }
}
