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
  //! @addtogroup core
  //! @{
  //! @var reldist
  //!
  //! @brief Callable object computing the reldistt operation.
  //!
  //! **Required header:** `#include <eve/module/core.hpp>`
  //!
  //! #### Members Functions
  //!
  //! | Member       | Effect                                                     |
  //! |:-------------|:-----------------------------------------------------------|
  //! | `operator()` | the reldistt operation   |
  //!
  //! ---
  //!
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  //!  template< value T, value U > auto operator()( T x, U y ) const noexcept requires compatible< T, U >;
  //!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //!
  //! **Parameters**
  //!
  //!`x`, `y`:   [values](@ref eve::value).
  //!
  //! **Return value**
  //!
  //!computes [elementwise](@ref glossary_elementwise) the 'units in the last place' distance betwween `x` and `y`.
  //!This is semantically equivalent to:`
  //!
  //!* if is_ordered(x,y), nb_values(x,y)/2.0 is returned
  //!* otherwise a `Nan` is returned
  //!
  //! ---
  //!
  //! #### Supported decorators
  //!
  //!  no decorators are supported
  //!
  //! #### Example
  //!
  //! @godbolt{doc/core/reldist.cpp}
  //!
  //!  @}
  //================================================================================================

  namespace tag { struct reldist_; }
  template<> struct supports_conditional<tag::reldist_> : std::false_type {};

  EVE_MAKE_CALLABLE(reldist_, reldist);
}

#include <eve/module/core/regular/impl/reldist.hpp>
