//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Project Contributors
  SPDX-License-Identifier: BSL-1.0
*/
//==================================================================================================
#pragma once

#include <eve/detail/overload.hpp>

namespace eve
{
//================================================================================================
//! @addtogroup math_invtrig
//! @{
//! @var asin
//!
//! @brief Callable object computing the arc sine.
//!
//!   **Defined in Header**
//!
//!   @code
//!   #include <eve/module/math.hpp>
//!   @endcode
//!
//!
//!   @groupheader{Callable Signatures}
//!
//!   @code
//!   namespace eve
//!   {
//!      template< eve::floating_value T >
//!      T asin(T x) noexcept;
//!   }
//!   @endcode
//!
//! **Parameters**
//!
//!   * `x`:   [floating real value](@ref eve::floating_ordered_value).
//!
//!  **Return value**
//!
//!   *  Returns the [elementwise](@ref glossary_elementwise) arc sine of the
//!      input in the range \f$[-\frac\pi2 , \frac\pi2]\f$.
//!
//!      In particular:
//!
//!      * If the element is \f$1\f$, \f$+0\f$ is returned.
//!      * If the element \f$|x| > 1\f$, `NaN` is returned.
//!      * If the element is a `Nan`, `NaN` is returned.
//!
//!      special cases are handled as if the operation was implemented by \f$-i \mathrm{asinh}(i z)\f$
//!
//!  @groupheader{Example}
//!
//!  @godbolt{doc/math/regular/asin.cpp}
//!
//!  @}
//================================================================================================
namespace tag
{
  struct asin_;
}

template<> struct supports_optimized_conversion<tag::asin_> : std::true_type
{};

EVE_MAKE_CALLABLE(asin_, asin);
}

#include <eve/module/math/regular/impl/asin.hpp>
