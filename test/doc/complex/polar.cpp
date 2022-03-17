#include <eve/module/math.hpp>
#include <eve/module/complex.hpp>
#include <eve/wide.hpp>
#include <iostream>

using wide_ft = eve::wide<float, eve::fixed<4>>;
using wide_ct = eve::wide<eve::complex<float>, eve::fixed<4>>;
using c_t     = eve::complex<float>;

int main()
{
  auto pio_4 = eve::pio_4(eve::as<float>());
  wide_ft rho  = {1.0f, 2.0f, 3.0f,  0.0f};
  wide_ft theta =  wide_ft{ 0.0f, 1.0f,  2.0f, 4.0}*pio_4;

  std::cout << "---- simd" << '\n'
            << "<- rho               = " << rho << '\n'
            << "<- theta             = " << theta << '\n'
            << "-> polar(rho, theta) = " << eve::polar(rho, theta) << '\n';

  float  r = 2.0f;
  float  t = pio_4;

  std::cout << "---- scalar" << '\n'
            << "<- r                = " << r << '\n'
            << "<- t                = " << t << '\n'
            << "-> polar(r, t)      = " << eve::polar(r, t) << '\n';
  return 0;
}
