//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Project Contributors
  SPDX-License-Identifier: BSL-1.0
*/
//==================================================================================================
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>
#include <cmath>

int main()
{
  auto lmin = EVE_VALUE(eve::valmin(eve::as<EVE_VALUE>()));
  auto lmax = EVE_VALUE(eve::valmax(eve::as<EVE_VALUE>()));

  auto arg0 = eve::bench::random_<EVE_VALUE>(lmin,lmax);
  auto std__csc = [](auto x){return 1/std::sin(x);};

  eve::bench::experiment xp;
  run<EVE_VALUE>(EVE_NAME(std__csc) , xp, std__csc , arg0);
  run<EVE_VALUE>(EVE_NAME(csc) , xp, eve::csc , arg0);
  run<EVE_TYPE> (EVE_NAME(csc) , xp, eve::csc , arg0);

}
