//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include <eve/module/complex.hpp>

EVE_TEST_TYPES( "Check complex::operator-", eve::test::scalar::ieee_reals)
<typename T>(eve::as<T>)
{
  using w_t = eve::wide<T>;

  auto fill_r = [](auto i, auto) { return T(1+i); };
  auto fill_i = [](auto i, auto) { return T(1)/(1+i); };
  auto fill   = [](auto i, auto) { return eve::complex<T>(T(1+i),T(1)/(1+i)); };

  eve::complex<T>             z_s(T{1.234}, T{5.678});
  eve::wide<eve::complex<T>>  z_v(fill);

  TTS_EQUAL( eve::real(-z_s), T{-1.234});
  TTS_EQUAL( eve::imag(-z_s), T{-5.678});

  TTS_EQUAL( eve::real(-z_v), -w_t{fill_r});
  TTS_EQUAL( eve::imag(-z_v), -w_t{fill_i});
};

EVE_TEST_TYPES( "Check complex::operator-", eve::test::scalar::ieee_reals)
<typename T>(eve::as<T>)
{
  using w_t = eve::wide<T>;

  auto fill_r = [](auto i, auto) { return T(1+i); };
  auto fill_i = [](auto i, auto) { return T(1)/(1+i); };
  auto fill   = [](auto i, auto) { return eve::complex<T>(T(1+i),T(1)/(1+i)); };

  eve::complex<T>             z_s1(T{1.234}, T{5.678}), z_s2(T{2.468}, T{1.357});
  eve::wide<eve::complex<T>>  z_v1(fill), z_v2(fill_i,fill_r);

  TTS_EQUAL(eve::real(z_s1-z_s2), T{-1.234});
  TTS_EQUAL(eve::imag(z_s1-z_s2), T{ 4.321});

  TTS_EQUAL( eve::real(z_v1-z_s2), w_t{[&](auto i, auto c) { return fill_r(i,c) - eve::real(z_s2); } });
  TTS_EQUAL( eve::imag(z_v1-z_s2), w_t{[&](auto i, auto c) { return fill_i(i,c) - eve::imag(z_s2); } });

  TTS_EQUAL( eve::real(z_v2-z_s1), w_t{[&](auto i, auto c) { return fill_i(i,c) - eve::real(z_s1); } });
  TTS_EQUAL( eve::imag(z_v2-z_s1), w_t{[&](auto i, auto c) { return fill_r(i,c) - eve::imag(z_s1); } });

  TTS_EQUAL( eve::real(z_v1-z_v2), w_t{[&](auto i, auto c) { return fill_r(i,c) - fill_i(i,c); } });
  TTS_EQUAL( eve::imag(z_v1-z_v2), w_t{[&](auto i, auto c) { return fill_i(i,c) - fill_r(i,c); } });
};

EVE_TEST_TYPES( "Check complex::operator+", eve::test::scalar::ieee_reals)
<typename T>(eve::as<T>)
{
  using w_t = eve::wide<T>;

  auto fill_r = [](auto i, auto) { return T(1+i); };
  auto fill_i = [](auto i, auto) { return T(1)/(1+i); };
  auto fill   = [](auto i, auto) { return eve::complex<T>(T(1+i),T(1)/(1+i)); };

  eve::complex<T>             z_s1(T{1.234}, T{5.678}), z_s2(T{2.468}, T{1.357});
  eve::wide<eve::complex<T>>  z_v1(fill), z_v2(fill_i,fill_r);

  TTS_EQUAL(eve::real(z_s1+z_s2), T{3.702});
  TTS_EQUAL(eve::imag(z_s1+z_s2), T{7.035});

  TTS_EQUAL( eve::real(z_v1+z_s2), w_t{[&](auto i, auto c) { return fill_r(i,c) + eve::real(z_s2); } });
  TTS_EQUAL( eve::imag(z_v1+z_s2), w_t{[&](auto i, auto c) { return fill_i(i,c) + eve::imag(z_s2); } });

  TTS_EQUAL( eve::real(z_v2+z_s1), w_t{[&](auto i, auto c) { return fill_i(i,c) + eve::real(z_s1); } });
  TTS_EQUAL( eve::imag(z_v2+z_s1), w_t{[&](auto i, auto c) { return fill_r(i,c) + eve::imag(z_s1); } });

  TTS_EQUAL( eve::real(z_v1+z_v2), w_t{[&](auto i, auto c) { return fill_r(i,c) + fill_i(i,c); } });
  TTS_EQUAL( eve::imag(z_v1+z_v2), w_t{[&](auto i, auto c) { return fill_i(i,c) + fill_r(i,c); } });
};
