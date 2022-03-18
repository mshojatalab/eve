//==================================================================================================
/**
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
**/
//==================================================================================================
#include "test.hpp"
#include "measures.hpp"
#include <eve/module/complex.hpp>
#include <eve/module/math.hpp>
#include <eve/module/core.hpp>
#include <complex>

EVE_TEST( "Check behavior of cosh on scalar"
        , eve::test::scalar::ieee_reals
        , eve::test::generate(eve::test::randoms(-10, 10)
                             , eve::test::randoms(-10, 10))
        )
  <typename T>(T const& a0, T const& a1 )
{
  using e_t = typename T::value_type;
  using c_t = std::complex<e_t>;
  using z_t = eve::complex<e_t>;
  for(auto e : a0)
  {
    for(auto f : a1)
    {
      std::cout << "e = " << e << " f =  " <<  f << std::endl;
      TTS_ULP_EQUAL(eve::cosh(eve::complex<e_t>(e, f)),  z_t(std::cosh(c_t(e, f))), 2);
      //  TTS_ULP_EQUAL(eve::pedantic(eve::cosh)(eve::complex<e_t>(e, f)),  std::cosh(c_t(e, f)), 2);
      //     TTS_ULP_EQUAL(eve::cosh(eve::complex<e_t>(e, f)),  z_t(std::cosh(c_t(e, f))), 2.0);
    }
  }
};

// EVE_TEST( "Check behavior of cosh on wide"
//         , eve::test::simd::ieee_reals
//         , eve::test::generate(eve::test::randoms(-10, 10)
//                              , eve::test::randoms(-10, 10))
//         )
//   <typename T>(T const& a0, T const& a1 )
// {
//   using e_t = typename T::value_type;
//   using ce_t = eve::complex<e_t>;
//   using z_t = eve::wide<eve::complex<e_t>, typename T::cardinal_type>;
//   using c_t = std::complex<e_t>;
//   auto a = eve::cosh(z_t{a0,a1});
//   z_t b;
//   for(int i = 0; i !=  eve::cardinal_v<T>; ++i)
//   {
//     ce_t z(std::cosh(c_t(a0.get(i), a1.get(i))));
//     //   auto b(map([](auto e,  auto f){ return z_t(std::cosh(c_t(e, f))); }, a0, a1));
//     b.set(i, z);
//   }
//   TTS_EXPECT( eve::all(eve::ulpdist(a, b) <= 2.0) );


//   z_t x = eve::if_else(a0 < 1, a, b);
//   std::cout << "x  " << x << std::endl;
//   std::cout << "a0 " << a0 << std::endl;
//   std::cout << "a1 " << a1 << std::endl;
//   std::cout << "a  " << a  << std::endl;
//   std::cout << "b  " << b  << std::endl;
// };


// EVE_TEST( "Check behavior of pedantic(cosh) on wide"
//         , eve::test::simd::ieee_reals
//         , eve::test::generate(eve::test::randoms(-10, 10)
//                              , eve::test::randoms(-10, 10))
//         )
//   <typename T>(T const& a0, T const& a1 )
// {
//   using e_t = typename T::value_type;
//   using ce_t = eve::complex<e_t>;
//   using z_t = eve::wide<eve::complex<e_t>, typename T::cardinal_type>;
//   using c_t = std::complex<e_t>;
//   auto a = eve::pedantic(eve::cosh)(z_t{a0,a1});
//   z_t b;
//   for(int i = 0; i !=  eve::cardinal_v<T>; ++i)
//   {
//     ce_t z(std::cosh(c_t(a0.get(i), a1.get(i))));
//     //   auto b(map([](auto e,  auto f){ return z_t(std::cosh(c_t(e, f))); }, a0, a1));
//     b.set(i, z);
//   }
//   TTS_EXPECT( eve::all(eve::ulpdist(a, b) <= 2.0) );


//   z_t x = eve::if_else(a0 < 1, a, b);
//   std::cout << "x  " << x << std::endl;
//   std::cout << "a0 " << a0 << std::endl;
//   std::cout << "a1 " << a1 << std::endl;
//   std::cout << "a  " << a  << std::endl;
//   std::cout << "b  " << b  << std::endl;
// };

// EVE_TEST_TYPES( "Check return types of eve::abs", eve::test::scalar::ieee_reals)
//   <typename T>(eve::as<T>)
// {
//   using e_t = eve::element_type_t<T>;
//   using c_t = eve::complex<e_t>;
//   using eve::as;
//   const int N = 20;
//   e_t pie = eve::pi (as<e_t>())-eve::eps(as<e_t>());
//   std::array<c_t, N> inputs =
//     { c_t(eve::zero(as<e_t>()),eve::zero(as<e_t>())),//0   c_t(eve::one(as<e_t>())(), eve::zero(as<e_t>()()),//0 ok
//       c_t(eve::inf(as<e_t>()),eve::zero(as<e_t>())), //1   c_t(eve::inf(as<e_t>())(),eve::zero(as<e_t>()()), //1 ok
//       c_t(eve::minf(as<e_t>()),eve::zero(as<e_t>())),//2   c_t(eve::inf(as<e_t>()),eve::zero(as<e_t>()()), //2 ok
//       c_t(eve::nan(as<e_t>()),eve::zero(as<e_t>())), //3   c_t(eve::nan(as<e_t>()),eve::zero(as<e_t>()()), //3 ok
//       c_t(eve::zero(as<e_t>()),eve::inf(as<e_t>())), //4   c_t(eve::nan(as<e_t>()), eve::zero(as<e_t>())()),//4 ok
//       c_t(eve::inf(as<e_t>()),eve::inf(as<e_t>())),  //5   c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>()()),  //5 ok
//       c_t(eve::minf(as<e_t>()),eve::inf(as<e_t>())), //6   c_t(eve::inf(as<e_t>(),eve::nan(as<e_t>()),  //6
//       c_t(eve::nan(as<e_t>()),eve::inf(as<e_t>())),  //7   c_t(eve::nan(as<e_t>(),eve::nan(as<e_t>()),  //7
//       c_t(eve::zero(as<e_t>()),eve::minf(as<e_t>())),//8   c_t(eve::nan(as<e_t>(), eve::zero(as<e_t>()),//8
//       c_t(eve::inf(as<e_t>()),eve::minf(as<e_t>())), //9   c_t(eve::inf(as<e_t>(),eve::nan(as<e_t>()),  //9
//       c_t(eve::minf(as<e_t>()),eve::minf(as<e_t>())),//10  c_t(eve::inf(as<e_t>(),eve::nan(as<e_t>()),  //10
//       c_t(eve::nan(as<e_t>()),eve::minf(as<e_t>())), //11  c_t(eve::nan(as<e_t>(),eve::nan(as<e_t>()),  //11
//       c_t(eve::zero(as<e_t>()),eve::nan(as<e_t>())), //12  c_t(eve::zero(as<e_t>(),eve::nan(as<e_t>()), //12
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //13  c_t(eve::inf(as<e_t>(),eve::nan(as<e_t>()),  //13
//       c_t(eve::minf(as<e_t>()),eve::nan(as<e_t>())), //14  c_t(eve::inf(as<e_t>(),eve::nan(as<e_t>()),  //14
//       c_t(eve::nan(as<e_t>()),eve::nan(as<e_t>())),  //15  c_t(eve::nan(as<e_t>(),eve::nan(as<e_t>()),  //15
//       c_t(eve::zero(as<e_t>()),pie),                 //16  c_t(eve::mone(as<e_t>(),eve::zero(as<e_t>()),//16
//       c_t(eve::inf(as<e_t>()),pie),                  //17  c_t(eve::minf(as<e_t>(),eve::inf(as<e_t>()), //17
//       c_t(eve::minf(as<e_t>()),pie),                 //18  c_t(eve::minf(as<e_t>(),eve::minf(as<e_t>()),//18
//       c_t(eve::nan(as<e_t>()),pie),                  //19  c_t(eve::nan(as<e_t>(),eve::nan (as<e_t>()), //19
//     };

//   std::array<c_t, N> results =
//     { c_t(eve::one(as<e_t>()), eve::zero(as<e_t>())),//0
//       c_t(eve::inf(as<e_t>()),eve::zero(as<e_t>())), //1
//       c_t(eve::inf(as<e_t>()),eve::zero(as<e_t>())), //2
//       c_t(eve::nan(as<e_t>()),eve::zero(as<e_t>())), //3
//       c_t(eve::nan(as<e_t>()), eve::zero(as<e_t>())),//4
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //5
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //6
//       c_t(eve::nan(as<e_t>()),eve::nan(as<e_t>())),  //7
//       c_t(eve::nan(as<e_t>()), eve::zero(as<e_t>())),//8
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //9
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //10
//       c_t(eve::nan(as<e_t>()),eve::nan(as<e_t>())),  //11
//       c_t(eve::nan(as<e_t>()),eve::zero(as<e_t>())), //12
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //13
//       c_t(eve::inf(as<e_t>()),eve::nan(as<e_t>())),  //14
//       c_t(eve::nan(as<e_t>()),eve::nan(as<e_t>())),  //15
//       c_t(eve::mone(as<e_t>()),eve::zero(as<e_t>())),//16
//       c_t(eve::minf(as<e_t>()),eve::inf(as<e_t>())), //17
//       c_t(eve::minf(as<e_t>()),eve::minf(as<e_t>())),//18
//       c_t(eve::nan(as<e_t>()),eve::nan (as<e_t>())), //19
//     };

//   auto eq = eve::numeric(eve::is_equal);
//   auto ch = eve::pedantic(eve::cosh);
//   for(int i=0; i < N; ++i)
//   {
//     std::cout << std::setprecision(10) << std::endl;
//     std::cout << i << std::endl;
//     std::cout << "input  " << inputs[i]<< std::endl;
//     std::cout << "eval   " << eve::cosh(inputs[i]) << std::endl;
//     std::cout << "expect " << results[i] << std::endl;
//     TTS_EXPECT(eq(ch(inputs[i]), results[i]));
//     TTS_EXPECT(eq(ch(-inputs[i]), ch(inputs[i])));
//     TTS_EXPECT(eq(ch(inputs[i]), eve::pedantic(eve::cos)(eve::mul_i(inputs[i]))));
//   }
// };
