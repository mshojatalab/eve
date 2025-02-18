//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Project Contributors
  SPDX-License-Identifier: BSL-1.0
*/
//==================================================================================================
#pragma once

#include <eve/arch.hpp>
#include <eve/detail/overload.hpp>

namespace eve
{
//================================================================================================
//! @addtogroup math_exp
//! @{
//! @var expmx2
//!
//! @brief Callable object computing \f$e^{-x^2}\f$.
//!
//!   **Defined in Header**
//!
//!   @code
//!   #include <eve/module/math.hpp>
//!   @endcode
//!
//!   @groupheader{Callable Signatures}
//!
//!   @code
//!   namespace eve
//!   {
//!      template< value T>
//!      T expmx2(T x) noexcept;
//!   }
//!   @endcode
//!
//! **Parameters**
//!
//!   *  `x`:   [floating value](@ref eve::floating_value).
//!
//! **Parameters**
//!
//!    `x`:   [floating real value](@ref eve::value).
//!
//! **Return value**
//!
//!   Returns the [elementwise](@ref glossary_elementwise) exponential of minus the square of `x`
//!
//!  @groupheader{Example}
//!
//!  @godbolt{doc/math/regular/expmx2.cpp}
//!
//!  @groupheader{Semantic Modifiers}
//!
//!   * Masked Call
//!
//!     The call `eve::expmx2[mask](x, ...)` provides a masked version of `eve::expmx2` which is
//!     equivalent to `if_else (mask, expmx2(x, ...), x)`.
//!
//!      **Example**
//!
//!        @godbolt{doc/math/masked/expmx2.cpp}
//!  @}
//================================================================================================
EVE_MAKE_CALLABLE(expmx2_, expmx2);
}

#include <eve/module/math/regular/impl/expmx2.hpp>
