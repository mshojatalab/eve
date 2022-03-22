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
#include <eve/module/core.hpp>
#include <complex>

EVE_TEST( "Check behavior of mul_i on scalar"
        , eve::test::scalar::ieee_reals
        , eve::test::generate( eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using z_t = eve::complex<e_t>;
  for(auto e : a0)
  {
    for(auto f : a1)
    {
      TTS_ULP_EQUAL(eve::mul_i(z_t(e, f)),  z_t(-f, e), 2);
      TTS_ULP_EQUAL(eve::mul_i(z_t(e, f)),  z_t(-f, e), 2.0);
    }
  }
};

EVE_TEST( "Check behavior of mul_i on wide"
        , eve::test::simd::ieee_reals
        , eve::test::generate(eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using z_t = eve::wide<eve::complex<e_t>, typename T::cardinal_type>;
  TTS_ULP_EQUAL(eve::mul_i(z_t{a0,a1}), z_t(-a1, a0), 0);
  TTS_ULP_EQUAL(eve::mul(z_t(a0,a1), eve::i), z_t(-a1, a0), 0);
  TTS_ULP_EQUAL(eve::mul(eve::i, z_t(a0,a1)), z_t(-a1, a0), 0);
  TTS_ULP_EQUAL(         eve::i* z_t(a0,a1) , z_t(-a1, a0), 0);
  TTS_ULP_EQUAL(z_t(a0,a1)*eve::i, z_t(-a1, a0), 0);
};
