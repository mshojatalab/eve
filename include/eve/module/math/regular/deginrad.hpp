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
//! @var deginrad
//!
//! @brief Callable object multiplying the input by \f$\pi/180\f$.
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
//!      T deginrad(T x) noexcept;
//!   }
//!   @endcode
//!
//! **Parameters**
//!
//!`x`:   [floating value](@ref eve::floating_value).
//!
//! **Return value**
//!
//! Returns the [elementwise](@ref glossary_elementwise) the degree input converted in radian.
//!
//!  @groupheader{Example}
//!
//!  @godbolt{doc/math/regular/deginrad.cpp}
//!  @}
//================================================================================================

EVE_MAKE_CALLABLE(deginrad_, deginrad);
}

#include <eve/module/math/regular/impl/deginrad.hpp>
