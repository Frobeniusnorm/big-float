#ifndef BIG_FLOAT
#define BIG_FLOAT
#include <vector>
typedef unsigned long size_t;
struct BigFloat {
  // it should always have more bytes than a double
  BigFloat(double init, size_t bytes = 16) : data(bytes, 0) {
    size_mantissa = (size_t)(bytes / 1.3);
    size_exponent = bytes - size_mantissa;
    // init data
    long init_data = (long)init;
    // first the fraction (6 bytes + 4 bits)
    for (int i = 0; i < 6; i++)
      data[i] = init_data & (0xFFl << (i * 8));
    // now the half of the 7th bit
    data[6] = init_data & (0xFl << 48);
    // exponent (1 byte + 3 bits) 
    data[size_mantissa] = init_data & (0xFFl << 52);
    data[size_mantissa + 1] = init_data & (0x7l << 60);
    sign = init_data & (0x1l << 63);
  }
  BigFloat() : BigFloat(0, 64) {}

  BigFloat operator+(const BigFloat b) const {
    // TODO
    return BigFloat(0);
  }
  BigFloat operator+=(const BigFloat b) const {
    // TODO
    return BigFloat(0);
  }
  double operator*() const {
    // TODO
    return 0;
  }
  bool operator==(const BigFloat b) const {
    // TODO
    return true;
  }
  // TODO sub, mul
protected:
  size_t size_mantissa;
  size_t size_exponent;
  char sign = 1;
  std::vector<char> data; // exponent, mantissa
};

#endif
