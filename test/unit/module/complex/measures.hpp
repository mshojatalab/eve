//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/wide.hpp>
#include <eve/arch/top_bits.hpp>
#include <eve/module/core.hpp>
#include <eve/module/complex.hpp>


namespace eve
{
  template<typename T, typename N>
  inline bool compare_equal(wide<eve::complex<T>, N> const &l, wide<eve::complex<T>, N> const &r)
  {
    return eve::all(is_equal(l, r));
  }
  template<typename T, typename N>
  inline bool compare_equal(complex<T> const &l, eve::complex<T> const &r)
  {
    return is_equal(l, r);
  }

//   template<typename T>
//   inline bool compare_equal(logical<T> const &l, logical<T> const &r)
//   {
//     if constexpr(eve::simd_value<T>)  return l.bitmap() == r.bitmap();
//     else                              return l == r;
//   }

//   template<typename T>
//   inline std::string to_string(logical<T> const &l)
//   {
//     std::ostringstream str;
//     str << l;
//     return str.str();
//   }

//   template<typename T>
//   inline std::string to_string(top_bits<T> const &l)
//   {
//     std::ostringstream str;
//     str << l;
//     return str.str();
//   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//==  complex
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace tts
{
  template<typename T, typename N>
  inline double ulp_distance(eve::wide<eve::complex<T>, N> const &l, eve::wide<eve::complex<T>, N> const &r)
  {
    std::cout << "icitte ulp" << std::endl;
    return eve::maximum(eve::ulpdist(l, r));
  }

  template<typename T, typename N>
  inline double relative_distance(eve::wide<eve::complex<T>, N> const &l, eve::wide<eve::complex<T>, N> const &r)
  {
    std::cout << "icitte rel" << std::endl;
    return eve::maximum(eve::reldist(l, r));
  }

  template<typename T, typename N>
  inline double absolute_distance(eve::wide<eve::complex<T>, N> const &l, eve::wide<eve::complex<T>, N> const &r)
  {
    std::cout << "icitte abs" << std::endl;
    return eve::maximum(eve::dist(l, r));
  }

  /////////////
  //==  scalar
  /////////////

  template<typename T, typename N>
  inline double ulp_distance(eve::complex<T> const &l, eve::complex<T> const &r)
  {
    std::cout << "icitte s ulp" << std::endl;
    return eve::maximum(eve::ulpdist(l, r));
  }

  template<typename T, typename N>
  inline double relative_distance(eve::complex<T> const &l, eve::complex<T> const &r)
  {
    std::cout << "icitte s rel" << std::endl;
    return eve::maximum(eve::reldist(l, r));
  }


  template<typename T, typename N>
  inline double absolute_distance(eve::complex<T> const &l, eve::complex<T> const &r)
  {
    std::cout << "icitte s abs" << std::endl;
    return eve::maximum(eve::dist(l, r));
  }

}
