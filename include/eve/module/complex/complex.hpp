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
#include <eve/module/complex.hpp>
#include <eve/product_type.hpp>
#include <eve/detail/abi.hpp>
#include <complex>

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

    ///  constructor from std:complex
    EVE_FORCEINLINE complex(std::complex<Type> c)  noexcept : parent{c.real(), c.imag()} {}

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

    template<like<complex> Z> EVE_FORCEINLINE friend auto operator+(Z z) noexcept
    {
      return z;
    }

    EVE_FORCEINLINE friend auto& operator*= ( like<complex> auto & self
                                            , like<complex> auto const & other
                                            ) noexcept
    {
      auto a = real(self);
      auto b = imag(self);
      auto c = real(other);
      auto d = imag(other);
      self = {diff_of_prod(a, c, b, d), sum_of_prod(a, d, b, c)};
      // still need some limit values treatments for inf nan
      return self;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto& operator *= ( Z & self
                                            , O const & other
                                            ) noexcept
    {
      auto a = real(self);
      auto b = imag(self);
      return self = Z{ a*other, b*other};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto operator/= ( Z & self
                                           , like<complex> auto const& other
                                           ) noexcept
    {
//       this can be a raw version
//       auto co = conj(other);
//       auto fact = sqr_abs(other);
//       auto rr = real(co)/fact;
//       auto ii = imag(co)/fact;
//       return self *= Z(rr, ii);

      auto rr =  eve::abs(real(self));
      auto ii =  eve::abs(imag(self));
      auto e =  -if_else((rr < ii), exponent(ii), exponent(rr));
      auto oother(eve::ldexp(other, e));
      auto denom =  sqr_abs(oother);
      self *= conj(oother);
      self /= denom;
      self = ldexp(self, e);
      return self;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto operator/= ( Z & self
                                           , O const& other
                                           ) noexcept
    {
      auto a = real(self);
      auto b = imag(self);
      return self = Z{ a/other, b/other};
    }

    //==============================================================================================
    //  Unary functions
    //==============================================================================================
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::abs_
                                                , like<complex> auto const& z
                                                ) noexcept
    {
      return eve::hypot(real(z), imag(z));
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::sqr_abs_
                                                , like<complex> auto const& z
                                                ) noexcept
    {
      return eve::sqr(real(z))+eve::sqr(imag(z));
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::arg_
                                                , like<complex> auto const& z
                                                ) noexcept
    {
      return eve::atan2(imag(z), real(z) );
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::conj_, Z const& z ) noexcept
    {
      return Z{real(z), -imag(z)};
    }

    //==============================================================================================
    // predicates
    //==============================================================================================
    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_imag_, Z const& z ) noexcept
    {
      return is_eqz(real(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_real_, Z const& z ) noexcept
    {
      return is_eqz(imag(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_not_finite_, Z const& z ) noexcept
    {
      return is_not_finite(imag(z)) || is_not_finite(real(z));
    }

    //==============================================================================================
    // trigo
    //==============================================================================================
    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cosh_, Z const& z ) noexcept {
//      using v_t = typename Z::value_type;
      auto rz = real(z);
      auto iz = imag(z);
      auto [s, c]   = sincos(iz);
      auto [sh, ch] = sinhcosh(rz);

      auto r = c*ch;
      auto i = s*sh;
      i = if_else(is_imag(z) || is_real(z), zero, i);
      auto res = Z(r, i);
      if (eve::none(is_not_finite(z))) return res;
      res = if_else(is_infinite(rz) && is_not_finite(iz)), Z(inf(as(rz), nan(as(rz))), res);
//       res = if_else(is_nan(rz) && is_inf(iz)),             Z(nan(as(rz), nan(as(rz))), res);
      return res;
    }
//     EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cos_, Z const& z ) noexcept;
//     EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sinh_, Z const& z ) noexcept;
//     EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sin_, Z const& z ) noexcept;


    //==============================================================================================
    //  Binary functions
    //==============================================================================================
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::ulpdist_
                                                , like<complex> auto const& z1
                                                , like<complex> auto const& z2
                                                ) noexcept
    {
      return eve::max( eve::ulpdist(real(z1), real(z2))
                     , eve::ulpdist(imag(z1), imag(z2)));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::ldexp_
                                                , Z const& z1
                                                , integral_value auto n
                                                ) noexcept
    {
      return Z{ ldexp(real(z1), n), ldexp(imag(z1), n)};
    }
  };

  //================================================================================================
  //! @}
  //================================================================================================
}

//#include <eve/module/complex/regular/trigo.hpp>
