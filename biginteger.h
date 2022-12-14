#pragma once
#include <compare>
#include <complex>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <vector>

class BigInteger {
 private:
  std::vector<int> digits_;
  bool is_positive_;

  static void fft(std::vector<std::complex<double>>&, bool);
  static size_t reverseBits(size_t, size_t);
  static size_t closestPower(size_t);
  static size_t getBinPower(size_t);
  std::vector<std::complex<double>> convertDigits() const;
  void fixDigits();
  BigInteger stupidMult(int) const;

 public:
  static const int BASE = 1'00;
  static const int DEC_POW = 2;

  BigInteger();
  BigInteger(int number);
  BigInteger(const BigInteger&);
  explicit BigInteger(const std::string&);

  BigInteger operator+();
  BigInteger operator-();
  BigInteger& operator++();
  BigInteger operator++(int);
  BigInteger& operator--();
  BigInteger operator--(int);
  friend void swap(BigInteger&, BigInteger&);
  BigInteger& operator=(const BigInteger&);
  friend BigInteger& operator+=(BigInteger&, const BigInteger&);
  friend BigInteger& operator*=(BigInteger&, const BigInteger&);
  friend BigInteger& operator/=(BigInteger&, const BigInteger&);

  friend std::istream& operator>>(std::istream&, BigInteger&);
  friend std::ostream& operator<<(std::ostream&, const BigInteger&);

  std::string toString(bool) const;
  const std::vector<int>& digits() const;
  size_t digitsSize() const;
  void changeSign();
  bool isPositive() const;
  bool isZero() const;
  bool absCmp(const BigInteger&) const;
  static BigInteger gcd(const BigInteger&, const BigInteger&);
  void shiftBase(int);

  explicit operator bool() const;
};

class Rational {
 private:
  BigInteger numerator_, denominator_;

  void fixDividers();
  bool isZero() const;

 public:
  Rational();
  Rational(int, int);
  Rational(const BigInteger&, const BigInteger&);
  Rational(const Rational&);

  void changeSign();
  bool isPositive() const;
  const BigInteger& numerator() const;
  const BigInteger& denominator() const;

  Rational operator+();
  Rational operator-();
  friend void swap(Rational&, Rational&);
  Rational& operator=(const Rational&);
  friend Rational& operator+=(Rational&, const Rational&);
  friend Rational& operator-=(Rational&, const Rational&);
  friend Rational& operator*=(Rational&, const Rational&);
  friend Rational& operator/=(Rational&, const Rational&);

  std::string toString() const;
  std::string asDecimal(size_t) const;
  explicit operator double() const;
};

BigInteger operator""_bi(const char* c_str) {
  std::string s = c_str;
  return BigInteger(s);
}

BigInteger operator""_bi(const char* c_str, size_t) {
  return BigInteger(c_str);
}

size_t BigInteger::getBinPower(size_t number) {
  for (size_t i = 0; i < INT32_WIDTH; ++i) {
    if (static_cast<size_t>(1 << i) >= number) return i;
  }
  return INT32_WIDTH;
}

size_t BigInteger::reverseBits(size_t number, size_t size) {
  for (size_t i = 0; i < size / 2; ++i) {
    size_t get = static_cast<size_t>((1 << i) + (1 << (size - i - 1)));
    get &= number;
    get = ((get << (size - 2 * i - 1)) & (1 << (size - i - 1))) +
          ((get >> (size - 2 * i - 1)) & (1 << i));
    number ^= number & static_cast<size_t>(((1 << i) + (1 << (size - i - 1))));
    number |= get;
  }
  return number;
}

size_t BigInteger::closestPower(size_t number) {
  for (size_t i = 0; i < INT32_WIDTH; ++i) {
    if (static_cast<size_t>(1 << i) >= number) return (1 << i);
  }
  return INT32_WIDTH;
}

std::vector<std::complex<double>> BigInteger::convertDigits() const {
  std::vector<std::complex<double>> complex(digitsSize());
  for (size_t index = 0; index < digitsSize(); ++index) {
    complex[index] = digits_[index];
  }
  size_t closest_power = BigInteger::closestPower(complex.size());
  while (complex.size() != closest_power) {
    complex.push_back(0.0);
  }
  return complex;
}

void BigInteger::fft(std::vector<std::complex<double>>& coeffs,
                     bool reverse = false) {
  size_t size = coeffs.size();

  std::complex<double> quark =
      std::polar(1.0, 2 * M_PI / static_cast<double>(size));
  if (reverse) {
    quark = std::complex<double>(1.0) / quark;
  }
  if (size == 1) {
    return;
  }

  for (size_t i = 0; i < size; ++i) {
    size_t reversed_i = reverseBits(i, getBinPower(size));
    if (i < reversed_i) std::swap(coeffs[i], coeffs[reversed_i]);
  }
  size_t end = size;
  for (size_t step = 2; step <= size; step *= 2) {
    std::complex<double> current_quark = quark;
    for (size_t i = size; i > step; i /= 2) current_quark *= current_quark;
    for (size_t start = 0; start != end; start += step) {
      size_t middle = start + step / 2;
      std::complex<double> quark_degree = 1.0;
      for (size_t left = start, right = middle; left != middle;
           left++, right++) {
        std::complex<double> first = coeffs[left];
        std::complex<double> second = coeffs[right] * quark_degree;
        coeffs[left] = first + second;
        coeffs[right] = first - second;
        quark_degree *= current_quark;
      }
    }
  }
  if (reverse) {
    for (size_t i = 0; i < size; ++i) coeffs[i] /= static_cast<double>(size);
  }
}

BigInteger::BigInteger()
    : digits_(std::vector<int>(1, 0)), is_positive_(true) {}
BigInteger::BigInteger(int number)
    : digits_(std::vector<int>(1, std::abs(number))),
      is_positive_(number >= 0) {
  fixDigits();
}
BigInteger::BigInteger(const BigInteger& other)
    : digits_(std::vector<int>(other.digits_)),
      is_positive_(other.is_positive_) {}

BigInteger::BigInteger(const std::string& inp_str)
    : digits_(std::vector<int>()), is_positive_(true) {
  std::string buf(DEC_POW, '0');
  std::string str = inp_str;
  buf[DEC_POW] = 0;
  size_t index = 0;
  if (str[index] == '-') {
    is_positive_ = false;
    ++index;
  } else if (str[index] == '+') {
    is_positive_ = true;
    ++index;
  }
  std::reverse(str.begin() + static_cast<int>(index), str.end());

  while (index < str.size()) {
    size_t buf_index = 0;
    std::fill(buf.begin(), buf.end(), '0');
    while (index < str.size() && buf_index < buf.size()) {
      buf[buf_index] = str[index];
      ++index;
      ++buf_index;
    }
    buf[buf_index] = 0;
    std::reverse(buf.begin(), buf.begin() + static_cast<int>(buf_index));
    digits_.push_back(std::stoi(buf));
    buf_index = 0;
  }
  fixDigits();
}

const std::vector<int>& BigInteger::digits() const { return digits_; }

size_t BigInteger::digitsSize() const { return digits_.size(); }

void BigInteger::changeSign() {
  is_positive_ = !is_positive_;
  if (isZero()) is_positive_ = true;
}

bool BigInteger::isPositive() const { return is_positive_; }

bool BigInteger::isZero() const {
  return digits_.size() == 1 && digits_[0] == 0;
}

BigInteger abs(const BigInteger& arg) {
  BigInteger tmp = arg;
  if (!arg.isPositive()) tmp.changeSign();
  return tmp;
}

void BigInteger::shiftBase(int shift_size) {
  if (shift_size >= 0) {
    BigInteger shift;
    shift.digits_.resize(static_cast<size_t>(shift_size + 1), 0);
    *(shift.digits_.end() - 1) = 1;
    *this *= shift;
    return;
  }
  if (std::abs(shift_size) >= static_cast<int>(digitsSize())) {
    digits_.resize(1, 1);
  } else {
    std::reverse(digits_.begin(), digits_.end());
    digits_.resize(digitsSize() + static_cast<size_t>(shift_size));
    std::reverse(digits_.begin(), digits_.end());
  }
}

BigInteger::operator bool() const { return !isZero(); }

void BigInteger::fixDigits() {
  // fix leading zeros
  while (digitsSize() > 1 && digits_[digitsSize() - 1] == 0) {
    digits_.pop_back();
  }

  // change sign
  if (*(digits_.rbegin()) < 0) {
    is_positive_ = !is_positive_;
    for (size_t index = 0; index < digitsSize(); ++index) {
      digits_[index] = -digits_[index];
    }
  }

  // fix negatives
  for (size_t index = digitsSize() - 1; index > 0; --index) {
    if (digits_[index - 1] <= 0) {
      int add = std::abs(digits_[index - 1]) / BASE + 1;
      digits_[index] -= add;
      digits_[index - 1] += add * BASE;
    }
  }

  // fix big digits
  int go = 0, stay = 0;
  for (size_t index = 0; index < digitsSize(); ++index) {
    digits_[index] += go;
    stay = digits_[index] % BASE;
    go = digits_[index] / BASE;
    digits_[index] = stay;
  }
  while (go != 0) {
    digits_.push_back(go % BASE);
    go /= BASE;
  }

  // fix leading zeros
  while (digitsSize() > 1 && digits_[digitsSize() - 1] == 0) {
    digits_.pop_back();
  }

  if (digits_.size() == 1 && digits_[0] == 0) is_positive_ = true;
}

bool BigInteger::absCmp(const BigInteger& other) const {
  if (digitsSize() > other.digitsSize()) return false;
  if (digitsSize() < other.digitsSize()) return true;
  for (size_t index = digitsSize(); index > 0; --index) {
    if (digits()[index - 1] > other.digits()[index - 1]) return false;
    if (digits()[index - 1] < other.digits()[index - 1]) return true;
  }
  return false;
}

void swap(BigInteger& lhs, BigInteger& rhs) {
  std::swap(lhs.is_positive_, rhs.is_positive_);
  std::swap(lhs.digits_, rhs.digits_);
}

BigInteger& BigInteger::operator=(const BigInteger& other) {
  if (this == &other) return *this;
  BigInteger tmp = other;
  swap(*this, tmp);
  return *this;
}

bool operator==(const BigInteger& lhs, const BigInteger& rhs) {
  if (&lhs == &rhs) return true;
  if (lhs.isZero() && rhs.isZero()) return true;
  if (lhs.digitsSize() != rhs.digitsSize()) return false;
  for (size_t index = 0; index < lhs.digitsSize(); ++index) {
    if (lhs.digits()[index] != rhs.digits()[index]) return false;
  }
  if (lhs.isPositive() != rhs.isPositive()) return false;
  return true;
}

// bool operator!=(const BigInteger& lhs, const BigInteger& rhs) {
//   return !(lhs == rhs);
// }

auto operator<=>(const BigInteger& lhs, const BigInteger& rhs) {
  if (!lhs.isPositive() && rhs.isPositive()) return std::weak_ordering::less;
  if (lhs.isPositive() && !rhs.isPositive()) return std::weak_ordering::greater;
  if (lhs == rhs) return std::weak_ordering::equivalent;
  bool abs_less = lhs.absCmp(rhs);
  if (lhs.isPositive() && rhs.isPositive())
    return abs_less ? std::weak_ordering::less : std::weak_ordering::greater;
  return abs_less ? std::weak_ordering::greater : std::weak_ordering ::less;
}

// bool operator<(const BigInteger& lhs, const BigInteger& rhs) {
//   if (!lhs.isPositive() && rhs.isPositive()) return true;
//   if (lhs.isPositive() && !rhs.isPositive()) return false;
//   bool abs_less = lhs.absCmp(rhs);
//   if (lhs.isPositive() && rhs.isPositive()) return abs_less ? true : false;
//   return abs_less ? false : true;
// }

// bool operator>(const BigInteger& lhs, const BigInteger& rhs) {
//   return rhs < lhs;
// }

// bool operator<=(const BigInteger& lhs, const BigInteger& rhs) {
//   return !(rhs > lhs);
// }

// bool operator>=(const BigInteger& lhs, const BigInteger& rhs) {
//   return !(lhs < rhs);
// }

BigInteger& operator+=(BigInteger& lhs, const BigInteger& rhs) {
  // // std::cerr << "+=";
  for (size_t i = 0; i < rhs.digitsSize(); ++i) {
    if (lhs.digitsSize() == i) {
      lhs.digits_.push_back(0);
    }
    lhs.digits_[i] += (lhs.is_positive_ == rhs.is_positive_ ? rhs.digits_[i]
                                                            : -rhs.digits_[i]);
  }
  lhs.fixDigits();
  return lhs;
}

BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp += rhs;
}

BigInteger& operator-=(BigInteger& lhs, const BigInteger& rhs) {
  // // std::cerr << "-=";
  lhs.changeSign();
  lhs += rhs;
  lhs.changeSign();
  return lhs;
}

BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp -= rhs;
}

BigInteger BigInteger::operator+() { return *this; }
BigInteger BigInteger::operator-() {
  BigInteger tmp = *this;
  tmp.changeSign();
  return tmp;
}

BigInteger& BigInteger::operator++() {
  digits_[0] += (is_positive_ ? 1 : -1);
  if (digits_[0] >= static_cast<int>(BASE)) fixDigits();
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger tmp = *this;
  ++(*this);
  return tmp;
}

BigInteger& BigInteger::operator--() {
  digits_[0] -= (is_positive_ ? 1 : -1);
  if (digits_[0] < 0) fixDigits();
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger tmp = *this;
  this->operator--();
  return tmp;
}

BigInteger& operator*=(BigInteger& lhs, const BigInteger& rhs) {
  // // std::cerr << "*=";
  std::vector<std::complex<double>> left = lhs.convertDigits(),
                                    right = rhs.convertDigits();
  left.resize(std::max(left.size(), right.size()) * 2, 0.0);
  right.resize(left.size(), 0.0);
  BigInteger::fft(left);
  BigInteger::fft(right);
  for (size_t index = 0; index < left.size(); ++index) {
    left[index] *= right[index];
  }
  BigInteger::fft(left, true);
  lhs.digits_.resize(left.size());
  long long digit = 0;
  for (size_t index = 0; index < left.size(); ++index) {
    digit += static_cast<long long>(std::floor(left[index].real() + 0.5));
    lhs.digits_[index] = static_cast<int>(digit % BigInteger::BASE);
    digit /= BigInteger::BASE;
  }
  if (digit) {
    lhs.digits_.push_back(static_cast<int>(digit));
  }

  lhs.is_positive_ = !(lhs.is_positive_ ^ rhs.is_positive_);
  lhs.fixDigits();
  return lhs;
}

BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  tmp *= rhs;
  return tmp;
}

std::string BigInteger::toString(bool leading_zeros = false) const {
  // // std::cerr << "toString";
  std::stringstream ss;
  if (!is_positive_) ss << "-";
  if (!leading_zeros)
    ss << digits_[digits_.size() - 1];
  else
    ss << std::setw(DEC_POW) << std::setfill('0')
       << digits_[digits_.size() - 1];
  if (digits_.size() > 1)
    for (size_t index = digits_.size() - 1; index > 0; --index) {
      ss << std::setw(DEC_POW) << std::setfill('0') << digits_[index - 1];
    }
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, const BigInteger& rhs) {
  // // std::cerr << "<<";
  os << rhs.toString();
  return os;
}

std::istream& operator>>(std::istream& is, BigInteger& rhs) {
  // // std::cerr << ">>";
  std::string str;
  is >> str;
  rhs = BigInteger(str);
  return is;
}

BigInteger BigInteger::stupidMult(int number) const {
  BigInteger tmp = *this;
  for (size_t i = 0; i < tmp.digitsSize(); ++i) tmp.digits_[i] *= number;
  tmp.fixDigits();
  return tmp;
}

BigInteger& operator/=(BigInteger& lhs, const BigInteger& rhs) {
  // // std::cerr << "/=";
  if (&lhs == &rhs) return lhs = 1;
  if (lhs.absCmp(rhs)) return lhs = 0;
  if (rhs.isZero()) {
    throw "BigInteger: Division by zero.";
  }
  BigInteger calc;
  BigInteger ans;
  BigInteger shift;
  bool l_sign = lhs.is_positive_, r_sign = rhs.is_positive_;
  ans.digits_.resize(0);
  shift.digits_.resize(lhs.digitsSize() - rhs.digitsSize() + 1, 0);

  while (!lhs.absCmp(rhs)) {
    calc.digits_.resize(rhs.digitsSize());
    std::copy(lhs.digits_.begin() + static_cast<int>(lhs.digitsSize()) -
                  static_cast<int>(calc.digitsSize()),
              lhs.digits_.end(), calc.digits_.begin());
    while (calc.absCmp(rhs)) {
      calc.digits_.insert(
          calc.digits_.begin(),
          lhs.digits_[lhs.digitsSize() - calc.digitsSize() - 1]);
    }

    for (size_t i = 0;
         i + lhs.digitsSize() + 2 < shift.digitsSize() + calc.digitsSize(); ++i)
      ans.digits_.push_back(0);

    int left = 0, right = BigInteger::BASE;
    while (left + 1 < right) {
      int mid = (left + right) / 2;
      if (calc.absCmp(abs(rhs.stupidMult(mid))))
        right = mid;
      else
        left = mid;
    }
    int digit = left;

    shift.digits_.resize(lhs.digitsSize() - calc.digitsSize() + 1, 0);
    BigInteger rhs_c;
    rhs_c.digits_.resize(shift.digitsSize() - 1 + rhs.digitsSize(), 0);
    std::copy(rhs.digits_.begin(), rhs.digits_.end(),
              rhs_c.digits_.begin() + static_cast<int>(shift.digitsSize() - 1));
    *(shift.digits_.end() - 1) = 1;

    lhs -=
        (l_sign ? abs(rhs_c.stupidMult(digit)) : -abs(rhs_c.stupidMult(digit)));
    ans.digits_.push_back(digit);
  }
  std::reverse(ans.digits_.begin(), ans.digits_.end());
  ans *= shift;
  lhs.digits_ = ans.digits_;
  lhs.is_positive_ = !(l_sign ^ r_sign);
  lhs.fixDigits();
  return lhs;
}

BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  tmp /= rhs;
  return tmp;
}

BigInteger& operator%=(BigInteger& lhs, const BigInteger& rhs) {
  return lhs = lhs - lhs / rhs * rhs;
}

BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp %= rhs;
}

Rational::Rational() : numerator_(0ll), denominator_(1ll) {}

Rational::Rational(int numerator, int denominator = 1ll)
    : numerator_(numerator), denominator_(denominator) {
  fixDividers();
}

Rational::Rational(const BigInteger& numerator,
                   const BigInteger& denominator = BigInteger(1ll))
    : numerator_(numerator), denominator_(denominator) {
  fixDividers();
}

Rational::Rational(const Rational& other)
    : numerator_(other.numerator_), denominator_(other.denominator_) {}

void Rational::changeSign() { numerator_.changeSign(); }

bool Rational::isPositive() const { return numerator_.isPositive(); }

const BigInteger& Rational::numerator() const { return numerator_; }

const BigInteger& Rational::denominator() const { return denominator_; }

BigInteger BigInteger::gcd(const BigInteger& lhs, const BigInteger& rhs) {
  if (rhs.isZero()) return abs(lhs);
  return gcd(abs(rhs), abs(lhs) % abs(rhs));
}

void Rational::fixDividers() {
  if (denominator_.isZero()) throw "Rational: Division by zero.";
  bool is_positive = true;
  if (!numerator_.isPositive()) {
    is_positive ^= !numerator_.isPositive();
    numerator_.changeSign();
  }
  if (!denominator_.isPositive()) {
    is_positive ^= !denominator_.isPositive();
    denominator_.changeSign();
  }
  BigInteger gcd = BigInteger::gcd(numerator_, denominator_);
  numerator_ /= gcd;
  denominator_ /= gcd;
  if (numerator_.isPositive() != is_positive) numerator_.changeSign();
}

Rational Rational::operator+() { return *this; }
Rational Rational::operator-() {
  Rational tmp = *this;
  tmp.changeSign();
  return tmp;
}

void swap(Rational& lhs, Rational& rhs) {
  swap(lhs.numerator_, rhs.numerator_);
  swap(lhs.denominator_, rhs.denominator_);
}

Rational& Rational::operator=(const Rational& other) {
  // std::cerr << toString() << " = " << other.toString() << "\n";
  Rational tmp = other;
  swap(*this, tmp);
  return *this;
}

Rational& operator+=(Rational& lhs, const Rational& rhs) {
  // std::cerr << lhs.toString() << " += " << rhs.toString() << "\n";
  lhs.numerator_ =
      lhs.numerator_ * rhs.denominator_ + rhs.numerator_ * lhs.denominator_;
  lhs.denominator_ *= rhs.denominator_;
  lhs.fixDividers();
  return lhs;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
  // // std::cerr << lhs.toString() << " + " << rhs.toString() << "\n";
  Rational tmp = lhs;
  tmp += rhs;
  return tmp;
}

Rational& operator-=(Rational& lhs, const Rational& rhs) {
  // // std::cerr << lhs.toString() << " -= " << rhs.toString() << "\n";
  lhs.numerator_ =
      lhs.numerator_ * rhs.denominator_ - rhs.numerator_ * lhs.denominator_;
  lhs.denominator_ *= rhs.denominator_;
  lhs.fixDividers();
  return lhs;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
  // // std::cerr << lhs.toString() << " - " << rhs.toString() << "\n";
  Rational tmp = lhs;
  return tmp -= rhs;
}

Rational& operator*=(Rational& lhs, const Rational& rhs) {
  // // std::cerr << "*=\n";
  lhs.numerator_ *= rhs.numerator_;
  lhs.denominator_ *= rhs.denominator_;
  lhs.fixDividers();
  return lhs;
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
  // // std::cerr << lhs.toString() << " *= " << rhs.toString() << "\n";
  Rational tmp = lhs;
  return tmp *= rhs;
}

Rational& operator/=(Rational& lhs, const Rational& rhs) {
  // std::cerr << lhs.toString() << " /= " << rhs.toString() << "\n";
  if (&lhs == &rhs) return lhs = 1;
  lhs.numerator_ *= rhs.denominator_;
  lhs.denominator_ *= rhs.numerator_;
  lhs.fixDividers();
  return lhs;
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
  // std::cerr << lhs.toString() << " / " << rhs.toString() << "\n";
  Rational tmp = lhs;
  return tmp /= rhs;
}

bool operator==(const Rational& lhs, const Rational& rhs) {
  // std::cerr << lhs.toString() << " == " << rhs.toString() << "\n";

  return lhs.numerator() == rhs.numerator() &&
         lhs.denominator() == rhs.denominator();
}

// bool operator!=(const Rational& lhs, const Rational& rhs) {
//   return !(lhs == rhs);
// }

auto operator<=>(const Rational& lhs, const Rational& rhs) {
  // std::cerr << lhs.toString() << " <=> " << rhs.toString() << "\n";

  return (lhs.numerator() * rhs.denominator()) <=>
         (rhs.numerator() * lhs.denominator());
}

// bool operator<(const Rational& lhs, const Rational& rhs) {
//   return (lhs.numerator() * rhs.denominator()) <
//          (rhs.numerator() * lhs.denominator());
// }

// bool operator>(const Rational& lhs, const Rational& rhs) { return rhs < lhs;
// }

// bool operator<=(const Rational& lhs, const Rational& rhs) {
//   return !(lhs > rhs);
// }

// bool operator>=(const Rational& lhs, const Rational& rhs) {
//   return !(lhs < rhs);
// }

std::string Rational::toString() const {
  return numerator_.toString() +
         (denominator_ != BigInteger(1ll) ? "/" + denominator_.toString() : "");
}

std::string Rational::asDecimal(size_t precision) const {
  // std::cerr << toString() << " asDecimal\n";
  if (numerator_.isZero()) return "0." + std::string(precision, '0');
  std::string ret = "";
  BigInteger part = abs(numerator_) / denominator_;
  if (!numerator_.isPositive()) ret += "-";
  ret += part.toString();

  part = abs(numerator_) % denominator_;
  part.shiftBase(2 * static_cast<int>(precision) / BigInteger::DEC_POW);
  part /= denominator_;
  std::string parts = part.toString(true);
  std::string fract(
      2 * precision / BigInteger::DEC_POW * BigInteger::DEC_POW - parts.size(),
      '0');
  fract += parts;
  fract.resize(precision + 1, '0');
  if (fract[fract.size() - 1] >= '5') fract[fract.size() - 2] += 1;
  fract.pop_back();
  ret += "." + fract;
  return ret;
}

Rational::operator double() const { return std::stod(asDecimal(40)); }
