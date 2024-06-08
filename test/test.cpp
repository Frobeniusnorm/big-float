#include <bitset>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../bigfloat.hpp"
#include "doctest.h"

TEST_SUITE("Floating Point semantics") {
  TEST_CASE("Addition") {
    BigFloat a(5.0);
    BigFloat b(6.0);
    BigFloat c = a + b;
    a += b;
    CHECK(*c == *a);
    CHECK_EQ(doctest::Approx(11.).epsilon(0.0000000001), *a);
    CHECK_EQ(doctest::Approx(11.).epsilon(0.0000000001), *c);
    // more complex shifting
    BigFloat d(128.0);
    BigFloat pi(3.141592);
    CHECK_EQ(doctest::Approx(128. + 3.141592).epsilon(0.0000000001), *(d + pi));
    CHECK_EQ(doctest::Approx(128. + 3.141592).epsilon(0.0000000001), *(pi + d));
    // sub decimal
    BigFloat f(0.00009876);
    BigFloat g(12345.6789);
    BigFloat h = f + g;
    CHECK_EQ(doctest::Approx(0.00009876 + 12345.6789), *(g + f));
    CHECK_EQ(doctest::Approx(0.00009876 + 12345.6789), *(f + g));
    BigFloat j(0.0012345);
    BigFloat k(0.0154321);
    double l = 0.0166666;
    CHECK_EQ(doctest::Approx(l), *(j + k));
  }
  TEST_CASE("Subtraction") {
    BigFloat a(5.0);
    BigFloat b(6.0);
    BigFloat c = b - a;
    b -= a;
    CHECK(*c == *b);
    CHECK_EQ(doctest::Approx(1.).epsilon(0.0000000001), *c);
    CHECK_EQ(doctest::Approx(1.).epsilon(0.0000000001), *b);
  }
  TEST_CASE("Wrap & Unwrap") {
    double a = 5.0;
    BigFloat b(a, 16);
    CHECK_EQ(doctest::Approx(a).epsilon(0.0000000001), *b);
    double pi = 3.141592;
    BigFloat c(pi, 19);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
    c = c.to_precision(25);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
    c = c.to_precision(20);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
    c = c.to_precision(16);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
    BigFloat d(0.00012345);
    CHECK_EQ(doctest::Approx(.00012345).epsilon(0.0000000001), *d);
  }
}
// TEST_CASE("Test Loading") {
//   std::cout << "===========" << std::endl;
//   DoubleAndLong conv;
//   conv.val = 5.0;
//   std::cout << std::bitset<64>(conv.binary) << std::endl;
//   BigFloat a(5.0);
//   for (long i = a.data.size() - 1; i >= 0; i--) {
//     std::cout << std::bitset<8>(a.data[i]);
//     if (i == a.size_mantissa)
//       std::cout << "|";
//   }
//   std::cout << std::endl;
// }
