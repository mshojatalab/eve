//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/concept/vectorizable.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>
#include <eve/product_type.hpp>
#include <eve/detail/abi.hpp>

namespace eve
{
  //================================================================================================
  //! @addtogroup simd_types
  //! @{
  //================================================================================================
  //! @brief SIMD-compatible representation of complex numbers
  //!
  //! **Required header:** `#include <eve/module/complex.hpp>`
  //!
  //! eve::complex is structure representing complex number and mean to be used in conjunction with
  //! eve::wide.
  //!
  //! @tparam Type  Underlying floating point type
  //================================================================================================
  template<floating_scalar_value Type>
  struct complex : struct_support<complex<Type>, Type, Type>
  {
    using eve_disable_ordering = void;
    using parent = struct_support<complex<Type>, Type, Type>;

    /// Underlying type
    using value_type = Type;

    /// Default constructor
    EVE_FORCEINLINE complex(Type r = 0, Type i = 0)  noexcept : parent{r,i} {}

    //==============================================================================================
    friend std::ostream& operator<<(std::ostream& os, like<complex> auto const& z)
    {
      return os << real(z) << std::showpos << imag(z) << "i" << std::noshowpos;
    }

    //==============================================================================================
    //  Real/Imag parts as functions
    //==============================================================================================

    /// Retrieve the real part of the current complex number
    EVE_FORCEINLINE friend
    decltype(auto) tagged_dispatch( eve::tag::real_, like<complex> auto&& z )
    {
      return get<0>(EVE_FWD(z));
    }

    /// Retrieve the imaginary part of the current complex number
    EVE_FORCEINLINE friend
    decltype(auto) tagged_dispatch( eve::tag::imag_, like<complex> auto&& z )
    {
      return get<1>(EVE_FWD(z));
    }

    //==============================================================================================
    //  Operators
    //==============================================================================================
    EVE_FORCEINLINE friend auto& operator+= ( like<complex> auto& self
                                            , like<complex> auto const& other
                                            ) noexcept
    {
      real(self) += real(other);
      imag(self) += imag(other);
      return self;
    }

    EVE_FORCEINLINE friend auto& operator-= ( like<complex> auto& self
                                            , like<complex> auto const& other
                                            ) noexcept
    {
      real(self) -= real(other);
      imag(self) -= imag(other);
      return self;
    }

    template<like<complex> Z> EVE_FORCEINLINE friend auto operator-(Z z) noexcept
    {
      return Z{-real(z), -imag(z)};
    }

    //==============================================================================================
    //  Unary functions
    //==============================================================================================
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::abs_
                                                , like<complex> auto const& z
                                                ) noexcept
    {
      return hypot(real(z), imag(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::conj_, Z const& z ) noexcept
    {
      return Z{real(z), -imag(z)};
    }
  };

  //================================================================================================
  //! @}
  //================================================================================================
}
