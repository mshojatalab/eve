#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<4>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<4>>;
using c_t     = eve::complex<float>;

int main()
{
  wide_ft fr = {-0.0f, 2.0f,  0.0f,  2.0f};
  wide_ft fi = { 0.0f, 0.0f,  2.0f, -3.5f};
  wide_ct f{fr, fi};

  std::cout << "---- simd" << '\n'
            << "<- fr                   = " << fr << '\n'
            << "<- f                    = " << f << '\n'
//            << "-> is_eqz(fr)            = " << eve::is_eqz(fr) << '\n'
            << "-> is_eqz(f)             = " << eve::is_eqz(f)  << '\n';

  float  sfr = -0.0f;
  float  sfi =  3.0f;
  c_t    sf{sfr, sfi};

  std::cout << "---- scalar" << '\n'
            << "<- sfr                  = " << sfr << '\n'
    //          << "-> is_eqz(sfr)           = " << eve::is_eqz(sfr) << '\n'
            << "<- sf                   = " << sf << '\n'
            << "-> is_eqz(sf)            = " << eve::is_eqz(sf) << '\n';
  return 0;
}
