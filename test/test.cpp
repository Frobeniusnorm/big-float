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
    BigFloat m = BigFloat(0.000014325) + BigFloat(41347123.12312);
    CHECK_EQ(doctest::Approx(0.000014325 + 41347123.12312), *m);
    BigFloat o(21505.2), p(29105.1);
    CHECK_EQ(doctest::Approx(21505.2 + 29105.1), *(o + p));
    BigFloat q = BigFloat(31649.6, 12) + BigFloat(94.8853, 12);
    double r = 31649.6 + 94.8853;
    CHECK_EQ(doctest::Approx(r), *q);
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
    // 39534.3 - 10348.2
    double h = 10348.2 - 39534.3;
    CHECK(BigFloat(39534.3) > BigFloat(10348.2));
    BigFloat i = BigFloat(10348.2);
    i -= BigFloat(39534.3);
    CHECK_EQ(doctest::Approx(h), *i);
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
    CHECK_EQ(doctest::Approx(31649.6), *BigFloat(31649.6));
    CHECK_EQ(doctest::Approx(94.8853), *BigFloat(94.8853));
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
    BigFloat a(0.5);
    BigFloat b(2.);
    CHECK_EQ(*(a * b), 1.);

    BigFloat c(3.141592);
    BigFloat d(0.001234);
    BigFloat e = c * d;
    CHECK_EQ(doctest::Approx(3.141592 * 0.001234), *e);
    BigFloat f = BigFloat(0.00324) * BigFloat(333.);
    CHECK_EQ(doctest::Approx(333. * 0.00324), *f);
    BigFloat g = BigFloat(0.000014325, 12) * BigFloat(41347123.12312, 12);
    CHECK_EQ(doctest::Approx(0.000014325 * 41347123.12312), *g);
    CHECK_EQ(doctest::Approx(12345 * 54321),
             *(BigFloat(12345.) * BigFloat(54321.)));
    CHECK_EQ(doctest::Approx(21505.2 * 29105.1),
             *(BigFloat(21505.2) * BigFloat(29105.1)));
    double h = 31436.5 * 8106.82;
    BigFloat k = BigFloat(31436.5) * BigFloat(8106.82);
    CHECK_EQ(doctest::Approx(h), *k);
  }
  TEST_CASE("Randomized") {
    for (int i = 0; i < 10000; i++) {
      double a = rand() / 50000.0;
      double b = rand() / 50000.0;
      BigFloat x(a);
      BigFloat y(b);
      // wrap and unwrap
      CHECK_EQ(a, *x);
      CHECK_EQ(b, *y);
      // comparison
      if (a < b)
        CHECK(x < y);
      else if (a > b)
        CHECK(x > y);
      else
        CHECK(x == y);
      // arithmetic
      CHECK_EQ(doctest::Approx(a - b), *(x - y));
      CHECK_EQ(doctest::Approx(a + b), *(x + y));
      CHECK_EQ(doctest::Approx(a + (-b)), *(x + (-y)));
      CHECK_EQ(doctest::Approx((-a) + (-b)), *((-x) + (-y)));
      CHECK_EQ(doctest::Approx((-a) + b), *((-x) + y));
      CHECK_EQ(doctest::Approx((-a) - b), *((-x) - y));
      CHECK_EQ(doctest::Approx((-a) - (-b)), *((-x) - (-y)));
      CHECK_EQ(doctest::Approx(b - a), *(y - x));
      CHECK_EQ(doctest::Approx(a * b), *(x * y));
    }
  }
}
