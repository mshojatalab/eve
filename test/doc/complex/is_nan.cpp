#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<4>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<4>>;
using c_t     = eve::complex<float>;

int main()
{
  auto inf = eve::inf(eve::as<float>());
  auto nan = eve::nan(eve::as<float>());
  wide_ft fr = {-0.0f, 2.0f,  nan,  2.0f};
  wide_ft fi = { 0.0f, 0.0f,  2.0f, inf };
  wide_ct f{fr, fi};

  std::cout << "---- simd" << '\n'
            << "<- fr             = " << fr << '\n'
            << "<- f              = " << f << '\n'
            << "-> is_nan(fr)     = " << eve::is_nan(fr) << '\n'
            << "-> is_nan(fi)     = " << eve::is_nan(fi) << '\n'
            << "-> is_nan(f)      = " << eve::is_nan(f)  << '\n';

  float  sfr = -0.0f;
  float  sfi =  inf;
  c_t    sf{sfr, sfi};

  std::cout << "---- scalar" << '\n'
            << "<- sfr            = " << sfr << '\n'
            << "-> is_nan(sfr)    = " << eve::is_nan(sfr) << '\n'
            << "<- sf             = " << sf << '\n'
            << "-> is_nan(sf)     = " << eve::is_nan(sf) << '\n';
  return 0;
}
