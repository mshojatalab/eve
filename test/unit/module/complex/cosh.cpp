//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include <eve/module/complex.hpp>
#include <eve/module/math.hpp>
#include <eve/module/core.hpp>
#include <complex>

EVE_TEST( "Check behavior of cosh on scalar"
        , eve::test::scalar::ieee_reals
        , eve::test::generate(eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using c_t = std::complex<e_t>;
  using z_t = eve::complex<e_t>;
  for(auto e : a0)
  {
    for(auto f : a1)
    {
      TTS_EXPECT( eve::all(eve::ulpdist(eve::cosh(eve::complex<e_t>(e, f)),  z_t(std::cosh(c_t(e, f)))) <= 2));
      TTS_EXPECT( eve::all(eve::ulpdist(eve::pedantic(eve::cosh)(eve::complex<e_t>(e, f)),  z_t(std::cosh(c_t(e, f)))) <= 2));
//       TTS_ULP_EQUAL(eve::cosh(eve::complex<e_t>(e, f)),  z_t(std::cosh(c_t(e, f))), 2.0);
    }
  }
};

EVE_TEST( "Check behavior of cosh on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using ce_t = eve::complex<e_t>;
  using z_t = eve::wide<eve::complex<e_t>, typename T::cardinal_type>;
  using c_t = std::complex<e_t>;
  auto a = eve::cosh(z_t{a0,a1});
  z_t b;
  for(int i = 0; i !=  eve::cardinal_v<T>; ++i)
  {
    ce_t z(std::cosh(c_t(a0.get(i), a1.get(i))));
    //   auto b(map([](auto e,  auto f){ return z_t(std::cosh(c_t(e, f))); }, a0, a1));
    b.set(i, z);
  }
  TTS_EXPECT( eve::all(eve::ulpdist(a, b) <= 2.0) );


  z_t x = eve::if_else(a0 < 1, a, b);
  std::cout << "x  " << x << std::endl;
  std::cout << "a0 " << a0 << std::endl;
  std::cout << "a1 " << a1 << std::endl;
  std::cout << "a  " << a  << std::endl;
  std::cout << "b  " << b  << std::endl;
};


EVE_TEST( "Check behavior of pedantic(cosh) on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using ce_t = eve::complex<e_t>;
  using z_t = eve::wide<eve::complex<e_t>, typename T::cardinal_type>;
  using c_t = std::complex<e_t>;
  auto a = eve::pedantic(eve::cosh)(z_t{a0,a1});
  z_t b;
  for(int i = 0; i !=  eve::cardinal_v<T>; ++i)
  {
    ce_t z(std::cosh(c_t(a0.get(i), a1.get(i))));
    //   auto b(map([](auto e,  auto f){ return z_t(std::cosh(c_t(e, f))); }, a0, a1));
    b.set(i, z);
  }
  TTS_EXPECT( eve::all(eve::ulpdist(a, b) <= 2.0) );


  z_t x = eve::if_else(a0 < 1, a, b);
  std::cout << "x  " << x << std::endl;
  std::cout << "a0 " << a0 << std::endl;
  std::cout << "a1 " << a1 << std::endl;
  std::cout << "a  " << a  << std::endl;
  std::cout << "b  " << b  << std::endl;
};
