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
#include <eve/module/complex/regular/complex.hpp>
#include <eve/product_type.hpp>
#include <eve/detail/abi.hpp>
#include <complex>
#include <iostream>
#include <iomanip>

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
      return os << std::setprecision(20) << real(z) << std::showpos << imag(z) << "i" << std::noshowpos;
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
    // predicates
    //==============================================================================================

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_equal_
                                               , Z const& z1
                                               , Z const& z2) noexcept
    {
      auto [z1r, z1i] = z1;
      auto [z2r, z2i] = z2;
      return is_equal(z1r, z2r) && is_equal(z1i, z2i);
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_equal_
                                               , numeric_type const &
                                               , Z const& z1
                                               , Z const& z2) noexcept
    {
      auto [z1r, z1i] = z1;
      auto [z2r, z2i] = z2;
      return numeric(is_equal)(z1r, z2r) && numeric(is_equal)(z1i, z2i);
    }

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
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_finite_, Z const& z ) noexcept
    {
      return is_finite(imag(z)) && is_finite(real(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_eqz_, Z const& z ) noexcept
    {
      return is_eqz(imag(z)) && is_eqz(real(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_nez_, Z const& z ) noexcept
    {
      return is_nez(imag(z)) || is_nez(real(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_not_finite_, Z const& z ) noexcept
    {
      return is_not_finite(imag(z)) || is_not_finite(real(z));
    }

     template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_infinite_, Z const& z ) noexcept
    {
      return is_infinite(imag(z)) || is_infinite(real(z));
    }

   template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::is_nan_, Z const& z ) noexcept
    {
      return is_nan(imag(z)) || is_nan(real(z));
    }

    //==============================================================================================
    //  Operators
    //==============================================================================================
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // add/sub
    ////////////////////////////////////////////////////////////////////////////////////////////////

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

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::sub_
                                                , like<complex> auto z1
                                                , like<complex> auto z2
                                                ) noexcept
    {
      return z1 -= z2;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::sub_
                                                , Z const  self
                                                , O const & other
                                                ) noexcept
    {
      auto [a, b] = self;
      return Z{ a-other, b};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::sub_
                                                , O const  a0
                                                , Z const & a1
                                                ) noexcept
    {
      auto [a, b] = a1;
      return Z{ a0-a, -b};
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::add_
                                                , like<complex> auto z1
                                                , like<complex> auto z2
                                                ) noexcept
    {
      return z1 += z2;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::add_
                                                , Z const  self
                                                , O const & other
                                                ) noexcept
    {
      auto [a, b] = self;
      return Z{ a+other, b};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::add_
                                                , O a0
                                                , Z const & a1
                                                ) noexcept
    {
      auto [a, b] = a1;
      return Z{ a+a0, b};
    }

    template<like<complex> Z> EVE_FORCEINLINE friend auto operator-(Z z) noexcept
    {
      return Z{-real(z), -imag(z)};
    }

    template<like<complex> Z> EVE_FORCEINLINE friend auto operator+(Z z) noexcept
    {
      return z;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // multiplies
    ////////////////////////////////////////////////////////////////////////////////////////////////

    EVE_FORCEINLINE friend auto& operator*= ( like<complex> auto & self
                                            , like<complex> auto const & other
                                            ) noexcept
    {
      auto [a, b] = self;
      auto [c, d] = other;
      self = {diff_of_prod(a, c, b, d), sum_of_prod(a, d, b, c)};
      return self;
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::mul_
                                                , like<complex> auto z1
                                                , like<complex> auto z2
                                                ) noexcept
    {
      return z1 *= z2;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto& operator *= ( Z & self
                                            , O const & other
                                            ) noexcept
    {
      auto [a, b] = self;
      return self = Z{ a*other, b*other};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::mul_
                                                , Z a0
                                                , O a1
                                                ) noexcept
    {
      auto [a, b] = a0;
      return Z{ a*a1, b*a1};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::mul_
                                                , O a0
                                                , Z a1
                                                ) noexcept
    {
      auto [a, b] = a1;
      return Z{a*a0, b*a0};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto operator*( Z a0
                                         , O a1
                                         ) noexcept
    {
      return mul(a0, a1);
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto  operator*( O a0
                                          , Z a1
                                          ) noexcept
    {
      return mul(a1, a0);
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::mul_
                                                , pedantic_type const &
                                                , like<complex> auto const & a0
                                                , like<complex> auto const & a1
                                                ) noexcept
    {
      auto [a, b] = a0;
      auto [c, d] = a1;
      using c_t = decltype(a0*a1);
      c_t r{diff_of_prod(a, c, b, d), sum_of_prod(a, d, b, c)};
      auto [rr, ri] = r;
      auto test = is_finite(ri) && is_finite(rr);//is_finite(r);
      if (eve::all(test)) return r;
      auto cur  = eve::logical_andnot(is_real(a0), test);
      if (eve::any(cur))
      {
        r = if_else(cur, eve::mul(a1, a), r);
        test = logical_or(test, cur);
        if (eve::all(test)) return r;
      }
      cur = eve::logical_andnot(is_imag(a0), test);
      if (eve::any(cur))
      {
        r = if_else(cur, eve::mul_i(eve::mul(a1, b)), r);
        test = logical_or(test, cur);
        if (eve::all(test)) return r;
      }
      cur = eve::logical_andnot(is_real(a1), test);
      if (eve::any(cur))
      {
        r = if_else(cur, eve::mul(a0, c), r);
        test = logical_or(test, cur);
        if (eve::all(test)) return r;
      }
      cur = eve::logical_andnot(is_imag(a1), test);
      if (eve::any(cur))
      {
        r = if_else(cur, eve::mul_i(eve::mul(a0, d)), r);
        test = logical_or(test, cur);
        if (eve::all(test)) return r;
      }
      return r;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // divide
    ////////////////////////////////////////////////////////////////////////////////////////////////

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto operator/= ( Z & self
                                           , like<complex> auto const& other
                                           ) noexcept
    {
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
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::div_
                                                , Z a0
                                                , O a1
                                                ) noexcept
    {
      auto [a, b] = a0;
      return Z{a/a1, b/a1};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto tagged_dispatch (  eve::tag::div_
                                                , O a0
                                                , Z a1
                                                ) noexcept
    {
      return mul(rec(a1), a0);
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::div_
                                                , like<complex> auto a0
                                                , like<complex> auto a1
                                                ) noexcept
    {
      return a0 /= a1;
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::div_
                                                , pedantic_type const &
                                                , like<complex> auto a0
                                                , like<complex> auto a1
                                                ) noexcept
    {
      auto rr = eve::abs(real(a0));
      auto ii = eve::abs(imag(a1));
      auto e = -if_else((rr < ii), exponent(ii), exponent(rr));
      auto aa1 = ldexp(a1, e);
      auto denom =  sqr_abs(aa1);
      auto num = mul(a0, conj(aa1));
      auto r = ldexp(div(num, denom), e);
      if (eve::all(is_finite(r))) return r;
      auto a1r = real(a1);
      r = if_else(is_eqz(denom),   mul(a0, copysign(inf(as(rr)), a1r)), r);
      r = if_else(is_infinite(a1), mul(a0, rec(copysign(denom, a1r))),  r);
      return r;
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto operator/= ( Z & self
                                           , O const& other
                                           ) noexcept
    {
      auto [a, b] = self;
      return self = Z{ a/other, b/other};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::rec_, Z const& z ) noexcept
    {
      auto [a, b] = z;
      auto n2 = sqr_abs(z);
      return Z{a/n2, b/n2};
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto operator/( Z a0
                                         , O a1
                                         ) noexcept
    {
      return div(a0, a1);
    }

    template<like<complex> Z, floating_real_value O>
    EVE_FORCEINLINE friend auto  operator/( O a0
                                          , Z a1
                                          ) noexcept
    {
      return div(a0, a1);
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
//      return eve::sqr(real(z))+eve::sqr(imag(z));
      auto [zr, zi] = z;
      return sum_of_prod(zr, zr, zi, zi);
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::sqr_
                                                , Z const& z
                                                ) noexcept
    {
      auto [zr, zi] = z;
      return Z{diff_of_prod(zr, zr, zi, zi), 2*zr*zi};
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

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::mul_i_, Z const& z ) noexcept
    {
      return Z{-imag(z), real(z)};
    }


    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::mul_mi_, Z const& z ) noexcept
    {
      return Z{imag(z), -real(z)};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sqrt_
                                               , pedantic_type const &
                                               , Z const& z ) noexcept
    {
      //always compute the sqrt of the complex with positive imaginary part
      //then conjugate if necessary
      auto [rz, iz] = z;
      auto negimag = is_ltz(iz);
      auto x = abs(rz);
      auto y = abs(iz);
      auto iaz = if_else(negimag, -iz, iz); // always >= 0 or -Nan
      auto gtxy = (x > y);
      auto gezrz = is_gez(rz);
      auto r = if_else(gtxy, y/x, x/y);
      auto rr= sqrt(inc(sqr(r)));
      auto sqrtx = sqrt(x);
      auto w = if_else(gtxy,
                       sqrtx*sqrt(half(as(r))*inc(rr)),
                       sqrt(y)*sqrt(half(as(r))*(r+rr)));
      auto is_real_z = is_real(z);
      auto res = Z(if_else(is_real_z, sqrtx, w)
                  , if_else(is_real_z, zero, iaz*half(as(r))/w));
      res = if_else(gezrz, res, Z(imag(res), real(res)));
      auto infty = inf(as(iaz));
      res = if_else(iaz == infty,  Z{infty, infty}, z);
      res = if_else(logical_andnot(rz == -infty, is_nan(iz)), Z{if_else(iaz == infty, iaz, zero), infty}, res);
      res = if_else(is_eqz(iz) && is_nan(rz), z, res);
      res = if_else(is_nan(z), Z(nan(as(iz)), nan(as(iz))), res);
      res = if_else(is_eqz(iz) && is_gez(rz), Z(rz), res);
      return if_else(negimag, conj(res), res);
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sqrt_
                                               , Z const& z ) noexcept
    {
      //always compute the sqrt of the complex with positive imaginary part
      //then conjugate if necessary
      auto [rz, iz] = z;
      auto negimag = is_ltz(iz);
      auto x = abs(rz);
      auto y = abs(iz);
      auto iaz = if_else(negimag, -iz, iz); // always >= 0 or -Nan
      auto gtxy = (x > y);
      auto gezrz = is_gez(rz);
      auto r = if_else(gtxy, y/x, x/y);
      auto rr= sqrt(inc(sqr(r)));
      auto sqrtx = sqrt(x);
      auto w = if_else(gtxy,
                       sqrtx*sqrt(half(as(r))*inc(rr)),
                       sqrt(y)*sqrt(half(as(r))*(r+rr)));
      auto is_real_z = is_real(z);
      auto res = Z(if_else(is_real_z, sqrtx, w)
                  , if_else(is_real_z, zero, iaz*half(as(r))/w));
      res = if_else(gezrz, res, Z(imag(res), real(res)));
      return if_else(negimag, conj(res), res);
    }


//     template<like<complex> Z>
//     EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sqrt_, Z const& z ) noexcept
//     {
//       return regular(sqrt)(z);
//     }

    //==============================================================================================
    // trigo
    //==============================================================================================
    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cosh_, Z const& z ) noexcept
    {
      auto [rz, iz] = z;
      auto [s, c]   = sincos(iz);
      auto [sh, ch] = sinhcosh(rz);
      auto r = c*ch;
      auto i = s*sh;
      i = if_else(is_imag(z) || is_real(z), zero, i);
      return Z(r, i);
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cosh_, pedantic_type const &, Z const& z ) noexcept
    {
      auto res = cosh(z);
      if (eve::none(is_not_finite(z))) return res;
      auto [rz, iz] = z;
      res = if_else(is_infinite(rz) && is_not_finite(iz), Z(inf(as(rz)), nan(as(rz))), res);
      res = if_else(is_nan(rz) && is_infinite(iz),        Z(nan(as(rz)), nan(as(rz))), res);
      return res;
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cos_, Z const& z ) noexcept
    {
     return cosh(mul_i(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::cos_, pedantic_type const &, Z const& z ) noexcept
    {
      return pedantic(cosh)(mul_i(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sinh_, Z const& z ) noexcept
    {
      auto [rz, iz] = z;
      auto [s, c]   = sincos(iz);
      auto [sh, ch] = sinhcosh(rz);
      return Z{c*sh, s*ch};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sinh_, pedantic_type const &, Z const& z ) noexcept
    {
      auto [r, i] = sinh(z);
      if (none(is_not_finite(z))) return Z{r, i};
      auto [zr, zi] = z;
      r = if_else(is_infinite(zr) && is_not_finite(zi), zr, r);
      i = if_else(is_infinite(zr) && is_nan(zi), nan(as(zr)), i);
      r = if_else(is_nan(zr), zr, r);
      i = if_else(is_nan(zr), zr, i);
      i = if_else(is_real(z), zero, i);
      r = if_else(is_imag(z), zero, r);
      return Z{ r, i};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sin_, Z const& z ) noexcept
    {
      return mul_mi(sinh(mul_i(z)));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::sin_, pedantic_type const &, Z const& z ) noexcept
    {
      return mul_mi(pedantic(sinh)(mul_i(z)));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::tan_, Z const& z ) noexcept
    //tan( x + iy ) = (sin2x + i sinh(2y)/(cos2x+cosh2y).
    {
      auto zz =  z+z;
      auto [rz, iz] = zz;
      auto [s, c]   = sincos(rz);
      auto [sh, ch] = sinhcosh(iz);
      auto tmp = c+ch;
      auto r = if_else(is_imag(z), zero, s/tmp);
      auto i = if_else(is_real(z), zero, sh/tmp);
      return Z{r, i};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::tanh_, Z const& z ) noexcept
    //tanh( x + iy ) = (sinh2x + i sin(2y)/(cosh2x+cos2y).
    {
      auto zz =  z+z;
      auto [rz, iz] = zz;
      auto [s, c]   = sincos(iz);
      auto [sh, ch] = sinhcosh(rz);
      auto tmp = c+ch;
      auto r = if_else(is_imag(z), zero, sh/tmp);
      auto i = if_else(is_real(z), zero, s/tmp);
      return Z{r, i};
    }

    //==============================================================================================
    //  exponential
    //==============================================================================================
    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::log_, Z const& z ) noexcept
    {
      auto argz = eve::arg(z);
      auto absz = if_else(is_nan(real(z)) && (inf(as(argz)) == imag(z)), inf(as(argz)), abs(z));
      return Z{log(absz), argz};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::log10_, Z const& z ) noexcept
    {
      auto a = if_else(is_real(z) && is_nan(imag(z)), zero, arg(z)) ;
      return Z{log10(abs(z)), a*invlog_10(as(real(z)))};
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::exp_, Z const& z ) noexcept
    {
      auto [rz, iz] = z;
      auto [s, c]   = sincos(iz);
      auto rho = exp(rz);
      return if_else(is_real(z) || rz == minf(as(rz)),
                     Z{rho},
                     Z{rho*c, rho*s});
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::exp_i_, Z const& z ) noexcept
    {
      return exp(mul_i(z));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::exp_ipi_, Z const& z ) noexcept
    {
      auto [rz, iz] = z;
      auto [s, c]   = sinpicospi(iz);
      auto rho = exp(pi(as(iz))*iz);
      return if_else(is_imag(z) || iz == inf(as(rz)),
                     Z{rho},
                     Z{rho*c, rho*s});
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch( eve::tag::proj_
                                               , Z const& z ) noexcept
    {
      return if_else(is_infinite(z), Z(inf(as(real(z))), copysign(zero(as(real(z))), imag(z))), z);
    }

    //==============================================================================================
    //  Binary functions
    //==============================================================================================
    EVE_FORCEINLINE friend auto operator == ( like<complex> auto const & a
                                             , like<complex> auto const & b
                                            ) noexcept
    {
      return (real(a) == real(b)) &&  (imag(a) == imag(b));
    }

    EVE_FORCEINLINE friend auto operator != ( like<complex> auto const & a
                                             , like<complex> auto const & b
                                            ) noexcept
    {
      return (real(a) != real(b)) || (imag(a) != imag(b));
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::ulpdist_
                                                , like<complex> auto const& z1
                                                , like<complex> auto const& z2
                                                ) noexcept
    {
      return eve::max( eve::ulpdist(real(z1), real(z2))
                     , eve::ulpdist(imag(z1), imag(z2)));
    }

    template < typename T >
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::ulpdist_
                                                , like<complex> auto const& z1
                                                , std::complex<T> const& z2
                                                ) noexcept
    {
      return eve::max( eve::ulpdist(real(z1), z2.real())
                     , eve::ulpdist(imag(z1), z2.imag()));
    }

    template<like<complex> Z>
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::ldexp_
                                                , Z const& z1
                                                , integral_value auto n
                                                ) noexcept
    {
      return Z{ ldexp(real(z1), n), ldexp(imag(z1), n)};
    }

    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::pow_
                                                 , like<complex> auto const & a
                                                 , like<complex> auto const & b
                                                 ) noexcept
    {
      return if_else(is_eqz(a)&&is_eqz(b), one, exp(a*log(b)));
    }

    template<like<complex> Z, floating_real_value T>
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::pow_
                                                , Z const & a0
                                                , T const & a1
                                                ) noexcept
    {
      auto t = arg(a0);
      auto a = abs(a0);
      return if_else(is_eqz(a)&&is_eqz(a1), one, pow(a, a1)*exp_i(t*a1));
    }

    template<like<complex> Z, floating_real_value T>
    EVE_FORCEINLINE friend auto tagged_dispatch ( eve::tag::pow_
                                                , T const & a0
                                                , Z const & a1
                                                ) noexcept
    {
      return if_else(is_eqz(a0)&&is_eqz(a1), one, exp(a1*cmplx(log)(a0)));
    }
  };

  //================================================================================================
  //! @}
  //================================================================================================
}
