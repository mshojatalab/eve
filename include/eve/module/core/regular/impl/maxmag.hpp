//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Project Contributors
  SPDX-License-Identifier: BSL-1.0
*/
//==================================================================================================
#pragma once

#include <eve/traits/common_value.hpp>
#include <eve/concept/value.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/detail/skeleton_calls.hpp>
#include <eve/module/core/regular/abs.hpp>
#include <eve/module/core/regular/all.hpp>
#include <eve/module/core/regular/if_else.hpp>
#include <eve/module/core/regular/is_not_greater_equal.hpp>
#include <eve/module/core/regular/max.hpp>

#include <type_traits>

namespace eve::detail
{
template<value T, value U>
EVE_FORCEINLINE auto
maxmag_(EVE_SUPPORTS(cpu_), T const& a, U const& b) noexcept
-> common_value_t<T, U>
{
  return arithmetic_call(maxmag, a, b);
}

template<value T>
EVE_FORCEINLINE auto
maxmag_(EVE_SUPPORTS(cpu_), T const& a, T const& b) noexcept
{
  auto aa = eve::abs(a);
  auto bb = eve::abs(b);
  if constexpr( simd_value<T> )
  {
    auto tmp = if_else(is_not_greater_equal(aa, bb), b, eve::max(a, b));
    return if_else(is_not_greater_equal(bb, aa), a, tmp);
  }
  else { return aa < bb ? b : bb < aa ? a : eve::max(a, b); }
}

//================================================================================================
// Masked case
//================================================================================================
template<conditional_expr C, value U, value... V>
EVE_FORCEINLINE auto
maxmag_(EVE_SUPPORTS(cpu_),
        C const& cond,
        U  t,
        V ...f) noexcept
{
  return mask_op(cond, eve::maxmag, t, f...);
}

//================================================================================================
// N parameters
//================================================================================================
template<decorator D, value T0, value T1, value... Ts>
auto
maxmag_(EVE_SUPPORTS(cpu_), D const&, T0 a0, T1 a1, Ts... args) noexcept
-> common_value_t<T0, T1, Ts...>
{
  using r_t = common_value_t<T0, T1, Ts...>;
  r_t that(D()(maxmag)(r_t(a0), r_t(a1)));
  ((that = D()(maxmag)(that, r_t(args))), ...);
  return that;
}

template<value T0, value T1, value... Ts>
auto
maxmag_(EVE_SUPPORTS(cpu_), T0 a0, T1 a1, Ts... args) noexcept
-> common_value_t<T0, T1, Ts...>
{
  using r_t = common_value_t<T0, T1, Ts...>;
  r_t that(maxmag(r_t(a0), r_t(a1)));
  ((that = maxmag(that, r_t(args))), ...);
  return that;
}

//================================================================================================
// tuples
//================================================================================================
template<kumi::non_empty_product_type Ts>
auto
maxmag_(EVE_SUPPORTS(cpu_), Ts tup) noexcept
{
  if constexpr( kumi::size_v<Ts> == 1) return get<0>(tup);
  else return kumi::apply( [&](auto... m) { return maxmag(m...); }, tup);
}

template<decorator D, kumi::non_empty_product_type Ts>
auto
maxmag_(EVE_SUPPORTS(cpu_), D const & d, Ts tup) noexcept
{
  if constexpr( kumi::size_v<Ts> == 1) return get<0>(tup);
  else return kumi::apply( [&](auto... m) { return d(maxmag)(m...); }, tup);
}

}
