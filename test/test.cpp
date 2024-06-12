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
    BigFloat d(6.0);
    BigFloat e = a - d;
    CHECK_EQ(doctest::Approx(-1.).epsilon(0.0000000001), *e);
    BigFloat f(0.001);
    BigFloat g(1.5);
    f -= g;
    CHECK_EQ(doctest::Approx(-1.499).epsilon(0.0000000001), *f);
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
    BigFloat e(3.14159265359, 64);
    BigFloat f(e);
    f = f.to_precision(16);
    e = e.to_precision(100);
    CHECK_EQ(doctest::Approx(*e).epsilon(0.0000000001), *f);
  }
  TEST_CASE("Comparison") {
    double a = 5.0;
    double b = 6.0;
    double c = 5.0;
    double d = 3.141592;
    double e = 0.0012345;
    auto bfa = BigFloat(a);
    auto bfb = BigFloat(b);
    auto bfc = BigFloat(c);
    auto bfd = BigFloat(d);
    auto bfe = BigFloat(e);
    CHECK(bfa < bfb);
    CHECK(bfa != bfb);
    CHECK(bfa == bfc);
    CHECK(bfa > bfd);
    CHECK(bfa > bfe);
    CHECK(bfe < bfd);
    CHECK(bfe != bfd);
    CHECK(bfb != bfc);
    CHECK(bfb > bfd);
    CHECK((bfa + BigFloat(1.)) == bfb);
  }
  TEST_CASE("Negative Cases") {
    BigFloat a(100.12345);
    BigFloat b(101.);
    CHECK_EQ(doctest::Approx(100.12345 - 101.).epsilon(0.0000000001), *(a - b));
    BigFloat c(-1.);
	b += c;
    CHECK_EQ(doctest::Approx(100.).epsilon(0.0000000001), *b);
  }
  TEST_CASE("Multiplication") {
    BigFloat a(1.);
    BigFloat b(2.);
    CHECK_EQ(*(a * b), *b);
  }
}
