//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#pragma once
// **=======================================================
// helper file to include all functions
#include <eve/function/abs.hpp>
#include <eve/function/acosd.hpp>
#include <eve/function/acosh.hpp>
#include <eve/function/acos.hpp>
#include <eve/function/acospi.hpp>
#include <eve/function/acotd.hpp>
#include <eve/function/acoth.hpp>
#include <eve/function/acot.hpp>
#include <eve/function/acotpi.hpp>
#include <eve/function/acscd.hpp>
#include <eve/function/acsch.hpp>
#include <eve/function/acsc.hpp>
#include <eve/function/acscpi.hpp>
#include <eve/function/add.hpp>
#include <eve/function/all.hpp>
#include <eve/function/any.hpp>
#include <eve/function/arg.hpp>
#include <eve/function/asecd.hpp>
#include <eve/function/asech.hpp>
#include <eve/function/asec.hpp>
#include <eve/function/asecpi.hpp>
#include <eve/function/asind.hpp>
#include <eve/function/asinh.hpp>
#include <eve/function/asin.hpp>
#include <eve/function/asinpi.hpp>
#include <eve/function/atan2d.hpp>
#include <eve/function/atan2.hpp>
#include <eve/function/atan2pi.hpp>
#include <eve/function/atand.hpp>
#include <eve/function/atanh.hpp>
#include <eve/function/atan.hpp>
#include <eve/function/atanpi.hpp>
#include <eve/function/average.hpp>
#include <eve/function/binarize.hpp>
#include <eve/function/binarize_not.hpp>
#include <eve/function/bit_and.hpp>
#include <eve/function/bit_andnot.hpp>
#include <eve/function/bit_cast.hpp>
#include <eve/function/bit_ceil.hpp>
#include <eve/function/bit_floor.hpp>
#include <eve/function/bit_mask.hpp>
#include <eve/function/bit_notand.hpp>
#include <eve/function/bit_not.hpp>
#include <eve/function/bit_notor.hpp>
#include <eve/function/bitofsign.hpp>
#include <eve/function/bit_or.hpp>
#include <eve/function/bit_ornot.hpp>
#include <eve/function/bit_select.hpp>
#include <eve/function/bit_shl.hpp>
#include <eve/function/bit_shr.hpp>
#include <eve/function/bit_width.hpp>
#include <eve/function/bit_xor.hpp>
#include <eve/function/broadcast.hpp>
#include <eve/function/cbrt.hpp>
#include <eve/function/ceil.hpp>
#include <eve/function/clamp.hpp>
#include <eve/function/combine.hpp>
#include <eve/function/conj.hpp>
#include <eve/function/convert.hpp>
#include <eve/function/copysign.hpp>
#include <eve/function/cosd.hpp>
#include <eve/function/cosh.hpp>
#include <eve/function/cos.hpp>
#include <eve/function/cospi.hpp>
#include <eve/function/cotd.hpp>
#include <eve/function/coth.hpp>
#include <eve/function/cot.hpp>
#include <eve/function/cotpi.hpp>
#include <eve/function/countl_one.hpp>
#include <eve/function/countl_zero.hpp>
#include <eve/function/countr_one.hpp>
#include <eve/function/countr_zero.hpp>
#include <eve/function/cscd.hpp>
#include <eve/function/csch.hpp>
#include <eve/function/csc.hpp>
#include <eve/function/cscpi.hpp>
#include <eve/function/dec.hpp>
#include <eve/function/deginrad.hpp>
#include <eve/function/dist.hpp>
#include <eve/function/div_180.hpp>
#include <eve/function/div.hpp>
#include <eve/function/eps.hpp>
#include <eve/function/exp10.hpp>
#include <eve/function/exp2.hpp>
#include <eve/function/exp.hpp>
#include <eve/function/expm1.hpp>
#include <eve/function/exponent.hpp>
#include <eve/function/fdim.hpp>
#include <eve/function/firstbitset.hpp>
#include <eve/function/firstbitunset.hpp>
#include <eve/function/floor.hpp>
#include <eve/function/fma.hpp>
#include <eve/function/fms.hpp>
#include <eve/function/fnma.hpp>
#include <eve/function/fnms.hpp>
#include <eve/function/frac.hpp>
#include <eve/function/frexp.hpp>
#include <eve/function/gather.hpp>
#include <eve/function/gcd.hpp>
#include <eve/function/has_single_bit.hpp>
#include <eve/function/horizontal.hpp>
#include <eve/function/hyperbolic.hpp>
#include <eve/function/hypot.hpp>
#include <eve/function/if_else.hpp>
#include <eve/function/ifnot_else.hpp>
#include <eve/function/ifrexp.hpp>
#include <eve/function/inc.hpp>
#include <eve/function/inverse_hyperbolic.hpp>
#include <eve/function/inverse_trigo.hpp>
#include <eve/function/is_denormal.hpp>
#include <eve/function/is_equal.hpp>
#include <eve/function/is_eqz.hpp>
#include <eve/function/is_even.hpp>
#include <eve/function/is_finite.hpp>
#include <eve/function/is_flint.hpp>
#include <eve/function/is_gez.hpp>
#include <eve/function/is_greater_equal.hpp>
#include <eve/function/is_greater.hpp>
#include <eve/function/is_gtz.hpp>
#include <eve/function/is_imag.hpp>
#include <eve/function/is_infinite.hpp>
#include <eve/function/is_less_equal.hpp>
#include <eve/function/is_lessgreater.hpp>
#include <eve/function/is_less.hpp>
#include <eve/function/is_lez.hpp>
#include <eve/function/is_ltz.hpp>
#include <eve/function/is_nan.hpp>
#include <eve/function/is_negative.hpp>
#include <eve/function/is_nez.hpp>
#include <eve/function/is_ngez.hpp>
#include <eve/function/is_ngtz.hpp>
#include <eve/function/is_nlez.hpp>
#include <eve/function/is_nltz.hpp>
#include <eve/function/is_normal.hpp>
#include <eve/function/is_not_denormal.hpp>
#include <eve/function/is_not_equal.hpp>
#include <eve/function/is_not_finite.hpp>
#include <eve/function/is_not_greater_equal.hp
#include <eve/function/is_not_greater.hpp>
#include <eve/function/is_not_imag.hpp>
#include <eve/function/is_not_infinite.hpp>
#include <eve/function/is_not_less_equal.hpp>
#include <eve/function/is_not_less.hpp>
#include <eve/function/is_not_nan.hpp>
#include <eve/function/is_not_real.hpp>
#include <eve/function/is_odd.hpp>
#include <eve/function/is_ordered.hpp>
#include <eve/function/is_positive.hpp>
#include <eve/function/is_pow2.hpp>
#include <eve/function/is_real.hpp>
#include <eve/function/is_unordered.hpp>
#include <eve/function/itrunc.hpp>
#include <eve/function/lcm.hpp>
#include <eve/function/ldexp.hpp>
#include <eve/function/lerp.hpp>
#include <eve/function/load.hpp>
#include <eve/function/log10.hpp>
#include <eve/function/log1p.hpp>
#include <eve/function/log2.hpp>
#include <eve/function/log.hpp>
#include <eve/function/logical_and.hpp>
#include <eve/function/logical_andnot.hpp>
#include <eve/function/logical.hpp>
#include <eve/function/logical_notand.hpp>
#include <eve/function/logical_not.hpp>
#include <eve/function/logical_notor.hpp>
#include <eve/function/logical_or.hpp>
#include <eve/function/logical_ornot.hpp>
#include <eve/function/logical_xor.hpp>
#include <eve/function/lookup.hpp>
#include <eve/function/manhattan.hpp>
#include <eve/function/mantissa.hpp>
#include <eve/function/max.hpp>
#include <eve/function/maximum.hpp>
#include <eve/function/maxmag.hpp>
#include <eve/function/min.hpp>
#include <eve/function/minimum.hpp>
#include <eve/function/minmag.hpp>
#include <eve/function/minus.hpp>
#include <eve/function/modf.hpp>
#include <eve/function/mul.hpp>
#include <eve/function/musl.hpp>
#include <eve/function/count_true.hpp>
#include <eve/function/nb_values.hpp>
#include <eve/function/nearest.hpp>
#include <eve/function/negate.hpp>
#include <eve/function/negatenz.hpp>
#include <eve/function/nextafter.hpp>
#include <eve/function/next.hpp>
#include <eve/function/none.hpp>
#include <eve/function/numeric.hpp>
#include <eve/function/oneminus.hpp>
#include <eve/function/pedantic.hpp>
#include <eve/function/plain.hpp>
#include <eve/function/plus.hpp>
#include <eve/function/popcount.hpp>
#include <eve/function/pow_abs.hpp>
#include <eve/function/pow.hpp>
#include <eve/function/predicates.hpp>
#include <eve/function/prev.hpp>
#include <eve/function/prod.hpp>
#include <eve/function/quadrant.hpp>
#include <eve/function/radindeg.hpp>
#include <eve/function/radinpi.hpp>
#include <eve/function/raw.hpp>
#include <eve/function/rec.hpp>
#include <eve/function/regular.hpp>
#include <eve/function/reldist.hpp>
#include <eve/function/rem.hpp>
#include <eve/function/rem_pio2.hpp>
#include <eve/function/rempio2.hpp>
#include <eve/function/rotl.hpp>
#include <eve/function/rotr.hpp>
#include <eve/function/round.hpp>
#include <eve/function/roundings.hpp>
#include <eve/function/rshl.hpp>
#include <eve/function/rshr.hpp>
#include <eve/function/rsqrt.hpp>
#include <eve/function/saturated.hpp>
#include <eve/function/saturate.hpp>
#include <eve/function/secd.hpp>
#include <eve/function/sech.hpp>
#include <eve/function/sec.hpp>
#include <eve/function/secpi.hpp>
#include <eve/function/shuffle.hpp>
#include <eve/function/sign.hpp>
#include <eve/function/signnz.hpp>
#include <eve/function/sinc.hpp>
#include <eve/function/sincosd.hpp>
#include <eve/function/sincos.hpp>
#include <eve/function/sincpi.hpp>
#include <eve/function/sindcosd.hpp>
#include <eve/function/sind.hpp>
#include <eve/function/sinhc.hpp>
#include <eve/function/sinhcosh.hpp>
#include <eve/function/sinh.hpp>
#include <eve/function/sin.hpp>
#include <eve/function/sinpicospi.hpp>
#include <eve/function/sinpi.hpp>
#include <eve/function/slice_high.hpp>
#include <eve/function/slice_low.hpp>
#include <eve/function/sqr_abs.hpp>
#include <eve/function/sqr.hpp>
#include <eve/function/sqrt.hpp>
#include <eve/function/stirling.hpp>
#include <eve/function/store.hpp>
#include <eve/function/sub.hpp>
#include <eve/function/sum.hpp>
#include <eve/function/swapbytes.hpp>
#include <eve/function/tand.hpp>
#include <eve/function/tanh.hpp>
#include <eve/function/tan.hpp>
#include <eve/function/tanpi.hpp>
#include <eve/function/tgamma.hpp>
#include <eve/function/trigonometric.hpp>
#include <eve/function/trigo_tags.hpp>
#include <eve/function/trunc.hpp>
#include <eve/function/two_add.hpp>
#include <eve/function/two_prod.hpp>
#include <eve/function/two_split.hpp>
#include <eve/function/ulpdist.hpp>
