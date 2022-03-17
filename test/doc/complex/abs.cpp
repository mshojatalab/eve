#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<4>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<4>>;
using c_t     = eve::complex<float>;

int main()
{
  wide_ft fr = {-1.0f, 2.0f, -3.0f, -32768.0f};
  wide_ft fi = {-1.0f, 4.0f, -3.5f, -32.0f};
  wide_ct f{fr, fi};

  std::cout << "---- simd" << '\n'
            << "<- fr                  = " << fr << '\n'
            << "<- f                   = " << f << '\n'
            << "-> abs(fr)             = " << eve::abs(fr) << '\n'
            << "-> abs(f)              = " << eve::abs(f)  << '\n';

  float  sfr = -1.0f;
  float  sfi =  3.0f;
  c_t    sf{sfr, sfi};

  std::cout << "---- scalar" << '\n'
            << "<- sf                  = " << sf << '\n'
            << "-> abs(sf)             = " << eve::abs(sf) << '\n';
  return 0;
}
