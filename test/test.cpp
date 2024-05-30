#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../bigfloat.hpp"
#include "doctest.h"

TEST_SUITE("Floating Point semantics") {
  TEST_CASE("Wrap & Unwrap") {
    double a = 5.0;
    BigFloat b(a, 16);
    CHECK_EQ(doctest::Approx(a).epsilon(0.0000000001), *b);
	double pi = 3.141592;
	BigFloat c(pi, 19);
	c = c.to_precision(25);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
	c = c.to_precision(20);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
	c = c.to_precision(16);
    CHECK_EQ(doctest::Approx(pi).epsilon(0.0000000001), *c);
  }
  TEST_CASE("Addition") {
    BigFloat a(5.0);
    BigFloat b(6.0);
    BigFloat c = a + b;
    a += b;
    CHECK(c == a);
  }
}
