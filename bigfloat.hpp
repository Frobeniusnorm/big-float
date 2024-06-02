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
    BigFloat c = BigFloat(*this);
    c += b;
    return c;
  }
  void operator+=(BigFloat b) {
    const size_t total_bytes = size_exponent + size_mantissa;
    const size_t total_other = b.size_mantissa + b.size_exponent;
    if (total_bytes != total_other)
      b = b.to_precision(total_bytes);
    const long shift = calculate_mantissa_shift(b);
    const long shift_a = shift < 0 ? -shift : 0;
    const long shift_b = shift > 0 ? shift : 0;
    if (shift_a != 0) {
      // copy new exponent of b
      for (size_t i = 0; i < size_exponent; i++)
        data[size_mantissa + i] = b.data[size_mantissa + i];
    }
    // now we can carry out integer addition of mantissa
    perform_mantissa_addition(b, shift_a, shift_b);
    // TODO test on sign
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
      const char bit = i % 8;
      final |= (long)(data[byte] & (1 << bit)) << (i - bit);
    }
    // a 1 before the decimal point might be present -> check for it and correct
    // it accordingly only if the first bit of the exponent is set (else it is
    // negative and the problem is taken care of implictly)
    if (data[size_mantissa + size_exponent - 1] & 0x80) {
      // accumulate shift if exponent
      size_t shift = 0;
      for (size_t i = 0; i < size_exponent * 8 - 1; i++) {
        const size_t byte = i / 8;
        const char bit = i % 8;
        if (data[size_mantissa + byte] & (1 << bit))
          shift += i == 0 ? 1 : 2 << (i - 1);
      }
      // if the bit is set...
      if (shift >= 1 && shift < size_mantissa * 8) {
        size_t shift_byte = shift / 8;
        char shift_bit = shift % 8;
        if (data[size_mantissa - 1 - shift_byte] & (0x80 >> (shift_bit - 1))) {
          // set the according bit in the result
          long shift_res = 0;
          for (size_t i = 0; i < 11; i++) {
            const size_t byte = i / 8;
            const char bit = i % 8;
            if (data[size_mantissa + byte] & (1 << bit))
              shift_res += i == 0 ? 1 : 2l << (i - 1);
          }
          final |= 1l << (52 - shift_res);
        }
      }
    }
    if (sign)
      final |= 1l << 63;
    else
      final &= ~(1l << 63);
    conv.binary = final;
    return conv.val;
  }
  bool operator==(const BigFloat b) const {
    // TODO
    return true;
  }
  // TODO sub, mul
  // protected:
  size_t size_mantissa; // in bytes
  size_t size_exponent; // in bytes
  char sign = 1;
  std::vector<char> data; // exponent, mantissa, i.e. 0 starts at the mantissa,
                          // it indexes as for byte operations
  size_t calculate_mantissa_shift(const BigFloat &b) const {
    // exponent may differ -> the one with the lower exponent has to be right
    // shifted to match the higher
    // count difference between this and b in exponent. If > 0 -> b has to be
    // right shifted, if < 0 this has to be right shifted
    // we manually subtract a - b in shift_a and b - a in shift_b s.t. we don't
    // have to care about negative numbers
    long shift_a = 0;
    long shift_b = 0;
    char carry_b = 0;
    char carry_a = 0;
    for (size_t i = 0; i < size_exponent * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const char data_a = (data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char data_b = (b.data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char sum_a = data_a - data_b - carry_a;
      if (sum_a == 0)
        carry_a = 0;
      else if (sum_a == 1)
        shift_a |= (1l << i);
      else if (sum_a == -1) {
        // borrow 2
        shift_a |= (1l << i);
        carry_a = 1;
      } else if (sum_a == -2) {
        carry_a = 1;
      }
      const char sum_b = data_b - data_a - carry_b;
      if (sum_b == 0)
        carry_b = 0;
      else if (sum_b == 1)
        shift_b |= (1l << i);
      else if (sum_b == -1) {
        // borrow 2
        shift_b |= (1l << i);
        carry_b = 1;
      } else if (sum_b == -2) {
        carry_b = 1;
      }
    }
    return carry_a == 0 ? shift_a : -shift_b;
  }
  void perform_mantissa_addition(BigFloat &b, size_t shift_a, size_t shift_b) {
    const size_t byte_shift_a = shift_a / 8;
    const size_t byte_shift_b = shift_b / 8;
    const size_t bit_shift_a = shift_a % 8;
    const size_t bit_shift_b = shift_b % 8;
    char carry = 0;
    for (size_t i = 0; i < size_mantissa * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      char data_a = 0;
      // since it is 1.mantissa, the 1 appears
      if (byte + byte_shift_a == size_mantissa - 1 && bit + bit_shift_a == 8) {
        data_a = 1;
      } else if (byte + byte_shift_a < size_mantissa &&
                 (byte + byte_shift_a != size_mantissa - 1 ||
                  bit + bit_shift_a < 8)) {
        const size_t byte_a = byte + byte_shift_a + ((bit + bit_shift_a) / 8);
        const char bit_a = (bit + bit_shift_a) % 8;
        data_a = (data[byte_a] & (1 << (bit_a))) >> bit_a;
      }
      char data_b = 0;
      if (byte + byte_shift_b == size_mantissa - 1 && bit + bit_shift_b == 8) {
        data_b = 1;
      } else if (byte + byte_shift_b < size_mantissa &&
                 (byte + byte_shift_b != size_mantissa - 1 ||
                  bit + bit_shift_b < 8)) {
        const size_t byte_b = byte + byte_shift_b + ((bit + bit_shift_b) / 8);
        const char bit_b = (bit + bit_shift_b) % 8;
        data_b = (b.data[byte_b] & (1 << (bit_b))) >> bit_b;
      }
      const char val = data_a + data_b + carry;
      if (val % 2 == 0)
        data[byte] &= ~(1 << bit);
      else
        data[byte] |= (1 << bit);
      if (val > 1)
        carry = 1;
      else
        carry = 0;
    }
    if (shift_a == 0 && shift_b == 0)
      carry = 1;
    // add carry to exponent and shift mantissa by 1 to right
    if (carry) {
      for (size_t i = 1; i < size_mantissa * 8; i++) {
        const size_t byte = i / 8;
        const char bit = i % 8;
        const char d = (data[byte] & (1 << bit));
        const size_t d_byte = bit == 0 ? byte - 1 : byte;
        const char d_bit = bit == 0 ? 7 : bit - 1;
        if (d)
          data[d_byte] |= (1 << d_bit);
        else
          data[d_byte] &= ~(1 << d_bit);
      }
      for (size_t i = 0; i < size_exponent * 8 && carry != 0; i++) {
        const size_t byte = i / 8;
        const char bit = i % 8;
        const char d = (data[size_mantissa + byte] & (1 << bit));
        if (d) {
          data[size_mantissa + byte] &= ~(1 << bit);
        } else {
          carry = 0;
          data[size_mantissa + byte] |= (1 << bit);
        }
      }
    }
  }
};

#endif
