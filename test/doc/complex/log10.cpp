#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<4>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<4>>;
using c_t     = eve::complex<float>;

int main()
{
  wide_ft fr = {-1.0f, 2.0f, -3.0f,  0.0f};
  wide_ft fi = { 0.0f, 2.0f,  2.0f, -3.5f};
  wide_ct f{fr, fi};

  std::cout << "---- simd" << '\n'
            << "<- fr                   = " << fr << '\n'
            << "<- f                    = " << f << '\n'
            << "-> log10(fr)            = " << eve::log10(fr) << '\n'
            << "-> cmplx(log10)(fr)     = " << eve::cmplx(eve::log10)(fr) << '\n'
            << "-> log10(f)             = " << eve::log10(f)  << '\n';

  float  sfr = -1.0f;
  float  sfi =  3.0f;
  c_t    sf{sfr, sfi};

  std::cout << "---- scalar" << '\n'
            << "<- sf                   = " << sf << '\n'
            << "-> log10(sf)            = " << eve::log10(sf) << '\n';
  return 0;
}
