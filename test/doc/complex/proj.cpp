#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<8>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<8>>;
using c_t     = eve::complex<float>;

int main()
{
  auto inf = eve::inf(eve::as<float>());
  auto nan = eve::nan(eve::as<float>());
  wide_ft fr = {-1.0f, 2.0f, inf,   inf,  -inf, nan,  nan, nan};
  wide_ft fi = { 1.0f, 0.0f, -3.5f, 32.0f, inf, nan, -inf, 3.0f};
  wide_ct f{fr, fi};

  std::cout << "---- simd" << '\n'
            << "<- fr                   = " << fr << '\n'
            << "<- f                    = " << f << '\n'
            << "-> proj(fr)             = " << eve::proj(fr) << '\n'
            << "-> proj(f)              = " << eve::proj(f)  << '\n';

  float  sfr = -1.0f;
  float  sfi =  3.0f;
  c_t    sf{sfr, sfi};

  std::cout << "---- scalar" << '\n'
            << "<- sf                   = " << sf << '\n'
            << "-> proj(sf)             = " << eve::proj(sf) << '\n';
  return 0;
}
