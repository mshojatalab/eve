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

EVE_TEST( "Check behavior of tanh on scalar"
        , eve::test::scalar::ieee_reals
        , eve::test::generate( eve::test::randoms(-10, 10)
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
      TTS_ULP_EQUAL(eve::tanh(eve::complex<e_t>(e, f)),  z_t(std::tanh(c_t(e, f))), 3.0);
      TTS_ULP_EQUAL(eve::tanh(eve::complex<e_t>(e, f)),  z_t(std::tanh(c_t(e, f))), 3.0);
    }
  }
};

EVE_TEST( "Check behavior of tanh on wide"
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
  auto std_ch = [](auto x, auto y){return std::tanh(c_t(x, y)); };
  auto init_with_std = [std_ch](auto a0,  auto a1){
    z_t b;
    for(int i = 0; i !=  eve::cardinal_v<T>; ++i)
    {
      ce_t z(std_ch(a0.get(i), a1.get(i)));
      b.set(i, z);
    }
    return b;
  };
  TTS_ULP_EQUAL(eve::tanh(z_t{a0,a1}), init_with_std(a0, a1), 3);
};

EVE_TEST_TYPES( "Check return types of eve::abs", eve::test::scalar::ieee_reals)
  <typename T>(eve::as<T>)
{
  using e_t = eve::element_type_t<T>;
  using c_t = eve::complex<e_t>;
  using eve::as;

  const int N = 19;
  std::array<c_t, N> inputs =
    { c_t(eve::zero(as<T>()),eve::zero(as<T>())), //0
      c_t(eve::inf(as<T>()),eve::zero(as<T>())),  //1
      c_t(eve::minf(as<T>()),eve::zero(as<T>())), //2
      c_t(eve::nan(as<T>()),eve::zero(as<T>())),  //3
      c_t(eve::zero(as<T>()),eve::inf(as<T>())),  //4
      c_t(eve::inf(as<T>()),eve::inf(as<T>())),   //5
      c_t(eve::minf(as<T>()),eve::inf(as<T>())),  //6
      c_t(eve::nan(as<T>()),eve::inf(as<T>())),   //7
      c_t(eve::zero(as<T>()),eve::minf(as<T>())), //8
      c_t(eve::inf(as<T>()),eve::minf(as<T>())),  //9
      c_t(eve::minf(as<T>()),eve::minf(as<T>())), //10
      c_t(eve::nan(as<T>()),eve::minf(as<T>())),  //11
      c_t(eve::zero(as<T>()),eve::nan(as<T>())),  //12
      c_t(eve::inf(as<T>()),eve::nan(as<T>())),   //13
      c_t(eve::minf(as<T>()),eve::nan(as<T>())),  //14
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),   //15
      c_t(eve::inf(as<T>()),eve::pi (as<T>())),   //16
      c_t(eve::minf(as<T>()),eve::pi (as<T>())),  //17
      c_t(eve::nan(as<T>()),eve::pi (as<T>())),   //18
    };

  std::array<c_t, N>  results =
    { c_t(eve::zero(as<T>()),eve::zero(as<T>())),//0
      c_t(eve::nan(as<T>()),eve::zero(as<T>())), //1
      c_t(eve::nan(as<T>()),eve::zero(as<T>())), //2
      c_t(eve::nan(as<T>()),eve::zero(as<T>())), //3
      c_t(eve::zero(as<T>()),eve::nan(as<T>())), //4
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //5
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //6
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //7
      c_t(eve::zero(as<T>()),eve::nan(as<T>())), //8
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //9
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //10
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //11
      c_t(eve::zero(as<T>()),eve::nan(as<T>())), //12
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //13
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //14
      c_t(eve::nan(as<T>()),eve::nan(as<T>())),  //15
      c_t(eve::nan(as<T>()),eve::zero(as<T>())), //16
      c_t(eve::nan(as<T>()),eve::zero(as<T>())), //17
      c_t(eve::nan(as<T>()),eve::nan (as<T>())), //18
    };

  for(int i=0; i < N; i++)
  {
    std::cout << i << std::endl;
    TTS_ULP_EQUAL((eve::tanh)(inputs[i]), results[i], 1);
    TTS_ULP_EQUAL((eve::tanh)(-inputs[i]), -results[i], 1);
  }
};
