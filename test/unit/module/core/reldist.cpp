//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include <eve/module/core.hpp>

//==================================================================================================
//== Types tests
//==================================================================================================
EVE_TEST_TYPES( "Check return types of reldist"
              , eve::test::simd::ieee_reals
        )
<typename T>(eve::as<T>)
{
  using v_t = eve::element_type_t<T>;

  //regular
  TTS_EXPR_IS( eve::reldist(T(), T()  ) , T);
  TTS_EXPR_IS( eve::reldist(T(), v_t()) , T);
  TTS_EXPR_IS( eve::reldist(v_t(), T()) , T);
  TTS_EXPR_IS( eve::reldist(v_t(), v_t()) , v_t);
};

EVE_TEST_TYPES( "Check return types of reldist"
        , eve::test::simd::ieee_reals
        )
<typename T>(eve::as<T>)
{
  using eve::reldist;
  using eve::as;
  TTS_EXPR_IS( reldist(T(), T()), T);

  if constexpr(eve::floating_value<T>)
  {
    if constexpr(eve::platform::supports_invalids)
    {
      TTS_EQUAL ( reldist(eve::inf(eve::as<T>()), eve::inf(eve::as<T>()))   , eve::inf(eve::as<T>()) );
      TTS_EQUAL ( reldist(eve::minf(eve::as<T>()), eve::minf(eve::as<T>())) , eve::inf(eve::as<T>()) );
      TTS_EQUAL ( reldist(eve::nan(eve::as<T>()), eve::nan(eve::as<T>()))   , eve::inf(eve::as<T>()) );
    }

    TTS_ULP_EQUAL( reldist(T(1), eve::inc(eve::eps(as<T>())))      , T(100)*eve::eps(as<T>()), 2);
    TTS_ULP_EQUAL( reldist(T(1), T(-eve::dec(eve::eps(as<T>()))))  , T(100)*eve::eps(as<T>()), 2);
    TTS_ULP_EQUAL( reldist(T(1), T(-eve::dec(eve::eps(as<T>())/2))), T(100)*eve::eps(as<T>())/2, 2);
  }
};
