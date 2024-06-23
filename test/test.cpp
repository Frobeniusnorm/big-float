#include <bitset>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../bigfloat.hpp"
#include "doctest.h"
// #include <sycl/sycl.hpp>

// using namespace sycl;
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
    double f = 8;
    double g = (1 << 16);
    auto bfa = BigFloat(a);
    auto bfb = BigFloat(b);
    auto bfc = BigFloat(c);
    auto bfd = BigFloat(d);
    auto bfe = BigFloat(e);
    auto bff = FixedFloat<16>(f);
    auto bfg = FixedFloat<16>(g);
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
    CHECK_EQ(0, *(FixedFloat<16>(0) * FixedFloat<16>(0)));
  }
  TEST_CASE("Mandelbrot") {
    double x_max = 1, x_min = -2, y_max = 1, y_min = -1;
    int res_width = 1920, res_height = 1080;
    for (int xi = 0; xi < 1; xi++) {
      for (int yi = 0; yi < 1; yi++) {
        int iter_normal;
        {
          const double x =
              xi / (double)(res_width - 1); // screen space in [0, 1]
          const double y =
              yi / (double)(res_height -
                            1); // screen space);  // screen space in [0, 1]

          const double x0 = (x * std::fabs(x_max - x_min) +
                             x_min); // complex plane in [x_min, x_max]
          const double y0 = y * std::fabs(y_max - y_min) +
                            y_min; // complex plane in [y_min, y_max]
          double c_x(0.0);
          double c_y(0.0);
          const size_t max_iter = 10;
          // simulate complex conjecture
          size_t iter = 0;
          for (; iter < max_iter; iter++) {
            double cx_squared = c_x * c_x;
            double cy_squared = c_y * c_y;
            std::cout << iter << " iteration: " << cx_squared + cy_squared
                      << std::endl;
            if (cx_squared + cy_squared > (1 << 16))
              break;
            double newx = cx_squared - cy_squared + x0;
            double newy = 2. * c_x * c_y + y0;
            c_x = newx;
            c_y = newy;
          }
          iter_normal = iter;
        }
        {
          const FixedFloat<16> x(
              xi / (double)(res_width - 1)); // screen space in [0, 1]
          const FixedFloat<16> y(
              yi / (double)(res_height -
                            1)); // screen space);  // screen space in [0, 1]

          const FixedFloat<16> x0 =
              (x * FixedFloat<16>(std::fabs(x_max - x_min) +
                                  x_min)); // complex plane in [x_min, x_max]
          const FixedFloat<16> y0 =
              (y * FixedFloat<16>(std::fabs(y_max - y_min) +
                                  y_min)); // complex plane in [y_min, y_max]
          FixedFloat<16> c_x(0.0);
          FixedFloat<16> c_y(0.0);
          const size_t max_iter = 10;
          // simulate complex conjecture
          size_t iter = 0;
          for (; iter < max_iter; iter++) {
            FixedFloat<16> cx_squared = c_x * c_x;
	    std::cout << *c_x << "^2 = " << *cx_squared << std::endl;
            FixedFloat<16> cy_squared = c_y * c_y;
	    std::cout << *c_y << "^2 = " << *cx_squared << std::endl;
            FixedFloat<16> l = cx_squared + cy_squared;
            FixedFloat<16> r((double)(1 << 16));
            std::cout << iter << " iteration: " << *l << std::endl;
            std::cout << *l << " > " << *r << " is " << (l > r) << std::endl;
            if (l > r)
              break;
            FixedFloat<16> newx = cx_squared - cy_squared + x0;
            FixedFloat<16> newy = FixedFloat<16>(2.) * c_x * c_y + y0;
            c_x = newx;
            c_y = newy;
          }
          CHECK_EQ(iter_normal, iter);
        }
      }
    }
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
TEST_SUITE("SyCL compatibility") {
  TEST_CASE("Templated Version") {
    for (int i = 0; i < 10000; i++) {
      double a = rand() / 50000.0;
      double b = rand() / 50000.0;
      FixedFloat<16> x(a);
      FixedFloat<16> y(b);
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
    for (int i = 0; i < 1000; i++) {
      double a = rand() / 500000.0;
      double b = rand() / 5.0;
      FixedFloat<16> x(a);
      FixedFloat<16> y(b);
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
  /* TEST_CASE("In SyCL") {

     default_selector device_selector;
     queue Q(device_selector);
     std::vector<double> res(500);
     buffer image_buffer(res.data(), range<1>{res.size()});
     Q.submit([&](auto &h) {
       accessor img(image_buffer, h, write_only, no_init);
       h.parallel_for(500, [=](item<1> i) {
         FixedFloat<16> a(5.3);
         FixedFloat<16> b((double)(i * 0.5));
         FixedFloat<16> c = a * (b - a);
         c += b * a;
         img[i] = *c;
       });
     });
     Q.wait();
     host_accessor h_acc(image_buffer);
     for (int i = 0; i < res.size(); i++) {
       double a = (5.3);
       double b = (i * 0.5);
       CHECK_EQ(doctest::Approx(a * (b - a) + b * a), h_acc[i]);
     }
   }*/
}
