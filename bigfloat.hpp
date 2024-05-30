#ifndef BIG_FLOAT
#define BIG_FLOAT
#include <algorithm>
#include <bitset>
#include <iostream>
#include <ostream>
#include <vector>
typedef unsigned long size_t;

typedef union DoubleAndLong {
  double val;
  long binary;
} DoubleAndLong;
struct BigFloat {
  // it should always have more bytes than a double
  BigFloat(double init, size_t bytes = 16) : data(bytes, 0) {
    size_mantissa = (size_t)(bytes / 1.3);
    size_exponent = bytes - size_mantissa;
    // init data
    DoubleAndLong conv;
    conv.val = init;
    const long init_data = conv.binary;
    // first the fraction (6 bytes + 4 bits)
    for (int i = 0; i < 6; i++)
      data[i] = (init_data & (0xFFl << (i * 8))) >> (i * 8);
    // now the half of the 7th bit
    data[6] = (init_data & (0xFl << 48)) >> 48l;
    // exponent (1 byte + 3 bits)
    const long exponent = (init_data & (0x7FFl << 52)) >> 52;
    data[size_mantissa] = exponent & 0xFFl;
    data[size_mantissa + 1] = (exponent & 0x700l) >> 8;
    // now adapt the bias by adding as many binary 1s as the exponential has
    // additional digits
    {
      const size_t double_expo = 11;
      char carry = 0;
      for (size_t i = size_mantissa * 8 + double_expo - 1;
           i < (size_mantissa + size_exponent) * 8 - 1; i++) {
        const size_t byte = i / 8;
        const size_t bit = i % 8;
        const unsigned char old_val = (data[byte] & (1 << bit)) == 0 ? 0 : 1;
        // addition
        const char new_val = old_val + carry + 1;
        if (new_val == 1) {
          data[byte] |= (1 << bit); // set to 1 (carry and old_val were 0)
        } else if (new_val == 2) {
          // set to 0 and set carry (either carry or old_val was set)
          data[byte] &= ~(1 << bit); // the digit was 1, we got no carry
                                     // -> it becomes 10, carry is set
          // else the digit was 0, we got a carry and a 1 -> keep carry
          carry = 1;
        } else { // 3
                 // carry is set, value is set
          data[byte] |= (1 << bit);
          carry = 1;
        }
      }
      if (carry == 1) {
        // highest exponent bit becomes 1
        data[size_mantissa + size_exponent - 1] |= 0x80;
      }
    }
    sign = (init_data & (0x1l << 63)) >> 63;
  }
  BigFloat(size_t bytes = 16) : data(bytes, 0) {
    size_mantissa = (size_t)(bytes / 1.3);
    size_exponent = bytes - size_mantissa;
  }

  BigFloat to_precision(size_t bytes) const {
    BigFloat res(bytes);
    // if the float is lowered back, because of the bias we assume that the
    // exponent has to stay the same
    res.size_exponent = size_exponent;
    for (size_t i = 0; i < std::min(res.size_mantissa, size_mantissa); i++)
      res.data[i] = data[i];
    // copy exponent
    for (size_t i = 0; i < std::min(res.size_exponent, size_exponent); i++) {
      res.data[res.size_mantissa + i] = data[size_mantissa + i];
    }
    // the bias has to be adapted
    // the bias is for k = exponent size in bits: 2^(k-1) - 1 -> binary ones to
    // k - 1 then zeros (k-1 is zero)
    if (res.size_exponent > size_exponent) {
      // if k2 > k1 -> bias is now higher!
      // we have to add bits to the fields (k1 - 1) to (k2 - 1)
      char carry = 0;
      for (size_t i = size_exponent * 8 - 1; i < res.size_exponent * 8 - 1;
           i++) {
        const size_t byte = i / 8;
        const size_t bit = i % 8;
        const char old_val = res.data[byte] & (1 << bit);
        // addition
        const char new_val = old_val + carry + 1;
        if (new_val == 1) {
          res.data[byte] |= (1 << bit); // set to 1 (carry and old_val were 0)
        } else if (new_val == 2) {
          // set to 0 and set carry (either carry or old_val was set)
          res.data[byte] &= ~(1 << bit); // the digit was 1, we got no carry
                                         // -> it becomes 10, carry is set
          // else the digit was 0, we got a carry and a 1 -> keep carry
          carry = 1;
        } else { // 3
                 // carry is set, value is set
          res.data[byte] |= (1 << bit);
          carry = 1;
        }
      }
      // now we have to keep simulating the carry
      for (size_t i = res.size_exponent; i < size_exponent && carry; i++) {
        const size_t byte = i / 8;
        const size_t bit = i % 8;
        const char old_val = res.data[byte] & (1 << bit);
        if (old_val == 0 && carry == 1) {
          carry = 0;
          res.data[byte] |= (1 << bit);
        } else if (old_val == 1 && carry == 1) {
          res.data[byte] &= ~(1 << bit);
        }
      }
    }
    res.sign = sign;
    return res;
  }

  BigFloat operator+(BigFloat b) const {
    const size_t total_bytes = size_exponent + size_mantissa;
    const size_t total_other = b.size_mantissa + b.size_exponent;
    if (total_bytes < total_other)
      return to_precision(total_other) + b;
    else if (total_bytes > total_other)
      b = b.to_precision(total_bytes);
    return BigFloat(0.0);
  }
  BigFloat operator+=(BigFloat b) const {
    const size_t total_bytes = size_exponent + size_mantissa;
    const size_t total_other = b.size_mantissa + b.size_exponent;
    if (total_bytes < total_other)
      return to_precision(total_other) + b;
    else if (total_bytes > total_other)
      b = b.to_precision(total_bytes);

    return BigFloat(0.0);
  }
  double operator*() const {
    // copy the exponent
    DoubleAndLong conv;
    long final = 0;
    // final exponent only has 11 bits, just copy 3 bytes from the original
    for (int i = 0; i < 2; i++)
      final |= ((long)data[size_mantissa + i] << (i * 8 + 52));
    // the old bias was 2^(size_exponent - 1) - 1, the new one is 2^10 -1, i.e.
    // we subtract the bits added in the constructor
    // since only one remaining bit (the 11th) is affected by this no loop is
    // necessary
    if ((data[size_mantissa + 1] & 0x4) !=
        0) {                // the first one that is subtracted
      final &= ~(1l << 62); // set it to 0
    } else
      final |= (1l << 62); // set it to 1 (borrow)
    // just copy as much from the mantissa as possible
    for (int i = 0; i < 52; i++) {
      const size_t byte = i / 8;
      const size_t bit = i % 8;
      final |= (long)(data[byte] & (1 << bit)) << (i - bit);
    }
    conv.binary = final;
    return sign == 0 ? conv.val : -conv.val;
  }
  bool operator==(const BigFloat b) const {
    // TODO
    return true;
  }
  // TODO sub, mul
protected:
  size_t size_mantissa; // in bytes
  size_t size_exponent; // in bytes
  char sign = 1;
  std::vector<char> data; // exponent, mantissa, i.e. 0 starts at the mantissa,
                          // it indexes as for byte operations
};

#endif
