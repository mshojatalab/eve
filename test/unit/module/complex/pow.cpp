//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include "measures.hpp"
#include <eve/module/complex.hpp>
#include <eve/module/math.hpp>
#include <complex>

EVE_TEST( "Check behavior of pow on scalar"
        , eve::test::scalar::ieee_reals
        , eve::test::generate( eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1, T const& a2, T const& a3 )
{
  using e_t = typename T::value_type;
  using c_t = eve::complex<e_t>;
  for(auto e : a0)
  {
    for(auto f : a1)
    {
      for(auto g : a2)
      {
        for(auto h : a3)
        {
          auto evep = eve::pow(c_t(e, f), c_t(g, h));
          auto stdpd = std::pow(std::complex<double>(e, f), std::complex<double>(g, h));
          auto stdp  = c_t(stdpd.real(), stdpd.imag());
          TTS_ULP_EQUAL( evep, stdp, 2.0e50);
        }
      }
    }
  }
};

EVE_TEST( "Check behavior of pow on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate( eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1, T const& a2, T const& a3)
{
  using e_t = typename T::value_type;
  using c_t = eve::complex<e_t>;
  using z_t = eve::wide<c_t, typename T::cardinal_type>;

  TTS_ULP_EQUAL( eve::pow(z_t{a0,a1}, z_t{a2,a3}), eve::exp(eve::mul(eve::log(z_t{a0,a1}), z_t{a2,a3}) ), 100.0);
};
