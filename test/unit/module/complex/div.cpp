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

  std::complex<T> c_s1(T{1}, T{5}), c_s2(T{2}, T{1}),  c_m(T{1.4}, T{1.8});
  eve::complex<T> z_s1(c_s1), z_s2(c_s2), z_m(c_m);
  TTS_EQUAL(std::real(c_s1/c_s2), std::real(c_m));
  TTS_EQUAL(std::imag(c_s1/c_s2), std::imag(c_m));

  TTS_EXPECT(eve::ulpdist(eve::real(z_s1/z_s2), eve::real(z_m)) <= 0.5);
  TTS_EXPECT(eve::ulpdist(eve::imag(z_s1/z_s2), eve::imag(z_m)) <= 0.5);

  std::cout << std::setprecision(15) << z_s1/z_s2 << std::endl;
  std::cout << std::setprecision(15) << c_s1/c_s2 << std::endl;
  auto fill_r = [](auto e, auto) { return T(1+e); };
  auto fill_i = [](auto e, auto) { return T(1)/T(1+3*e); };
  auto fill_r1 = [](auto e, auto) { return T(1-e/2); };
  auto fill_i1 = [](auto e, auto) { return T(-1+3.5*e); };

  using w_t = eve::wide<T>;
  using cw_t = eve::wide<eve::complex<T>>;
  w_t wr(fill_r);
  w_t wi(fill_i);
  w_t wr1(fill_r1);
  w_t wi1(fill_i1);
//   std::cout <<  "wr " << wr << std::endl;
//   std::cout <<  "wi " << wi << std::endl;

  eve::wide<eve::complex<T>> z(wr, wi);
  eve::wide<eve::complex<T>> zbar(wr, -wi);
  eve::wide<eve::complex<T>> z1(wr1, wi1);
  auto std_div =  [](eve::complex<T> a,  eve::complex<T> b) -> eve::complex<T> {
    return std::complex<T>(eve::real(a), eve::imag(a))/std::complex<T>(eve::real(b), eve::imag(b));
  };
  auto my_div = [std_div](cw_t a,  cw_t b) -> cw_t { return map(std_div, a, b); };
//   std::cout <<  "z       " << z << std::endl;
//   std::cout <<  "zbar    " << z << std::endl;
//   std::cout <<  "z/zbar  " << z/zbar << std::endl;
//   std::cout <<  "z/zbar  " << my_div(z, zbar);

  TTS_EXPECT(eve::all(eve::ulpdist( z/zbar, my_div(z, zbar)) <= 2));
  TTS_EXPECT(eve::all(eve::ulpdist( z/z1,   my_div(z, z1)  ) <= 2));
};
