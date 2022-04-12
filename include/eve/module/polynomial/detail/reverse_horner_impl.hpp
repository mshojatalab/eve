//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once

#include <eve/module/core.hpp>
#include <eve/module/polynomial/regular/horner.hpp>
#include <eve/module/polynomial/pedantic/horner.hpp>
#include <eve/module/polynomial/numeric/horner.hpp>
#include <eve/concept/range.hpp>
#include <concepts>
#include <iterator>
#include <initializer_list>



namespace eve::detail
{
  //================================================================================================
  //== N+ 1  parameters (((..(a*x+b)*x+c)*x + ..)..)
  //================================================================================================
 //  //==  N = 0
//   template<decorator D,  value T0>
//   EVE_FORCEINLINE constexpr auto reverse_horner_impl(D const &, T0 const &) noexcept
//   {
//     return T0(0);
//   }

//   //==  N = 1
//   template<decorator D, value T0, value T1>
//   EVE_FORCEINLINE constexpr auto reverse_horner_impl(D const &
//                                             , T0 const &, T1 const &a) noexcept
//   requires compatible_values<T0, T1>
//   {
//     using r_t = common_compatible_t<T0, T1>;
//     return r_t(a);
//   }

//   //==  N = 2
//   template<decorator D, value T0, value T1, value T2>
//   EVE_FORCEINLINE constexpr auto reverse_horner_impl(D const &
//                                             , T0 const &x, T1 const &a, T2 const &b) noexcept
//   requires compatible_values<T0, T1> &&compatible_values<T1, T2>
//   {
//     using r_t = common_compatible_t<T0, T1, T2>;
//     return D()(fma)(r_t(x), b, a);
//   }

//   //==  N >= 3
//   template<decorator D,
//            value T0,
//            value T1,
//            value T2,
//            value ...Ts>
//   EVE_FORCEINLINE constexpr
//   auto reverse_horner_impl(D const & d
//                           , T0 xx, T1 a, T2 b, Ts... args) noexcept
//   {
//     using r_t = common_compatible_t<T0, T1, T2, Ts...>;
//     auto x =  r_t(xx);
//     std::array<r_t, sizeof...(args)+2> t{};
//     int i = t.size()-1;
//     t[i] = a;
//     t[--i] = b;
//     [[maybe_unused]]auto cpr = [&t, &i](auto c){t[--i] = r_t(c); return c; };
//     ((cpr(args)), ...);
//     return d(horner)(x, t);

//   }

  template<decorator D, value T0, value ...Cs>
  EVE_FORCEINLINE constexpr
  auto reverse_horner_impl(D const & d, T0 xx, Cs... cs) noexcept
  {
    using r_t = common_compatible_t<T0, Cs...>;
    auto x =  r_t(xx);
    std::array<r_t, sizeof...(cs)> c{r_t(cs)...};
    return d(reverse_horner)(x, c);
  }


  //================================================================================================
  //== Reverse_Horner with ranges
  //================================================================================================
  template<decorator D, value T0, std::ranges::bidirectional_range R>
  EVE_FORCEINLINE constexpr auto reverse_horner_impl(D const & d, T0 xx, R const & r) noexcept
  requires (compatible_values<T0, typename R::value_type> && (!simd_value<R>))
  {
    using r_t = common_compatible_t<T0, typename R::value_type>;
    auto x =  r_t(xx);
    auto cur =  std::end(r);
    auto first  =  std::begin(r);
    if (first == cur) return r_t(0);
    else  if (std::distance(cur, first) == 1) return r_t(*cur);
    else
    {
      first--;
      auto dfma = d(fma);
      auto that = r_t(0);
      auto step = [&](auto that, auto arg){
        return dfma(x, that, arg);
      };
      for (--cur; cur != first; --cur)
        that = step(that, *cur);
      return that;
    }
  }
}
