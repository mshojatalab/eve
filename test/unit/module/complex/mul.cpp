//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include <eve/module/complex.hpp>
#include <complex>

EVE_TEST_TYPES( "Check complex::operator+", eve::test::scalar::ieee_reals)
<typename T>(eve::as<T>)
{
  using w_t = eve::wide<T>;

//   auto fill   = [](auto e, auto) { return eve::complex<T>(T(1+e)-5/T(1-e),5*T(1+e)+T(1)/(1-e)); };

//  eve::wide<eve::complex<T>>  z_v1(fill), z_v2(fill_i,fill_r);

  std::complex<T>             c_s1(T{1}, T{5}), c_s2(T{2}, T{1}),  c_m(T{-3}, T{11});
//  eve::complex<T>             z_s1(T{1}, T{5}), z_s2(T{2}, T{1}),  z_m(T{-3}, T{11});
  eve::complex<T>             z_s1(c_s1), z_s2(c_s2), z_m(c_m);
  TTS_EQUAL(std::real(c_s1*c_s2), std::real(c_m));
  TTS_EQUAL(std::imag(c_s1*c_s2), std::imag(c_m));

  TTS_EQUAL(eve::real(z_s1*z_s2), eve::real(z_m));
  TTS_EQUAL(eve::imag(z_s1*z_s2), eve::imag(z_m));

  auto fill_r = [](auto e, auto) { return T(1+e); };
  auto fill_i = [](auto e, auto) { return T(1)/T(1+3*e); };

  using w_t = eve::wide<T>;
  using cw_t = eve::wide<eve::complex<T>>;
  w_t wr(fill_r);
  w_t wi(fill_i);
//   std::cout <<  "wr " << wr << std::endl;
//   std::cout <<  "wi " << wi << std::endl;

  eve::wide<eve::complex<T>> z(wr, wi);
  eve::wide<eve::complex<T>> zbar(wr, -wi);
  auto std_mul =  [](eve::complex<T> a,  eve::complex<T> b) -> eve::complex<T> {
    return std::complex<T>(eve::real(a), eve::imag(a))*std::complex<T>(eve::real(b), eve::imag(b));
  };
  auto my_mul = [std_mul](cw_t a,  cw_t b) -> cw_t { return map(std_mul, a, b); };
//   std::cout <<  "z       " << z << std::endl;
//   std::cout <<  "zbar    " << z << std::endl;
//   std::cout <<  "z*zbar  " << z*zbar << std::endl;
//   std::cout <<  "z*zbar  " << my_mul(z, zbar);

  TTS_EXPECT(eve::all(eve::is_equal(z*zbar, my_mul(z, zbar))));
  TTS_EXPECT(eve::all(eve::ulpdist( z*zbar, my_mul(z, zbar)) < 0.5));
  std::cout << eve::ulpdist( z*zbar, my_mul(z, zbar)) << std::endl;
  std::cout << z*zbar  << std::endl;
  std::cout << my_mul(z, zbar) << std::endl;

//   TTS_EQUAL( eve::real(z_v1*z_s2), w_t{[&](auto i, auto c) { return eve::real(fill(i,c)); } });
//   TTS_EQUAL( eve::imag(z_v1*z_s2), w_t{[&](auto i, auto c) { return eve::imag(fill(i,c)); } });

//   TTS_EQUAL( eve::real(z_v2*z_s1), w_t{[&](auto i, auto c) { return fill_i(i,c) * eve::real(z_s1); } });
//   TTS_EQUAL( eve::imag(z_v2*z_s1), w_t{[&](auto i, auto c) { return fill_r(i,c) * eve::imag(z_s1); } });

//   TTS_EQUAL( eve::real(z_v1*z_v2), w_t{[&](auto i, auto c) { return fill_r(i,c) * fill_i(i,c); } });
//   TTS_EQUAL( eve::imag(z_v1*z_v2), w_t{[&](auto i, auto c) { return fill_i(i,c) * fill_r(i,c); } });
};
// EVE_TEST_TYPES( "Check complex::operator*", eve::test::scalar::ieee_reals)
//   <typename T>(eve::as<T>)
// {
//   using w_t = eve::wide<T>;

//   auto r0  = [](auto e, auto) -> eve::complex<T>{ return eve::complex<T>(T(1+e),T(e)); };
//   auto r1  = [](auto e, auto) -> eve::complex<T>{ return eve::complex<T>(T(1),T(e)); };
//   auto r2  = [](auto e, auto) -> eve::complex<T>{ return eve::complex<T>(T(1+e-e*e), T(2*e+e*e)) ; };
//   auto m = r0*r1;

//   TTS_EQUAL( eve::real(r2), eve::real(m));
//   TTS_EQUAL( eve::imag(r2), eve::imag(m));

// };
