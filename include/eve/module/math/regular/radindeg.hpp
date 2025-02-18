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
//! @addtogroup math_trig
//! @{
//! @var radindeg
//!
//! @brief Callable object multiplying the input by \f$180/\pi\f$.
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
//!      template< eve::floating_value T >
//!      T radindeg(T x) noexcept;
//!   }
//!   @endcode
//!
//! **Parameters**
//!
//!`x`:   [floating value](@ref eve::floating_value).
//!
//! **Return value**
//!
//! Returns the [elementwise](@ref glossary_elementwise) the radian input converted in degree.
//!
//!  @groupheader{Example}
//!
//!  @godbolt{doc/math/regular/radindeg.cpp}
//!  @}
//================================================================================================
EVE_MAKE_CALLABLE(radindeg_, radindeg);
}

#include <eve/module/math/regular/impl/radindeg.hpp>
