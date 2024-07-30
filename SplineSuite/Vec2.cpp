#include "Vec2.h"
#include <stdexcept>

template <typename T> Vec2<T> Vec2<T>::operator+(const Vec2<T> &other) const {
  return Vec2(x_ + other.x_, y_ + other.y_);
}

template <typename T> Vec2<T> &Vec2<T>::operator+=(const Vec2<T> &other) {
  x_ += other.x_;
  y_ += other.y_;
  return *this;
}

template <typename T> Vec2<T> Vec2<T>::operator-(const Vec2<T> &other) const {
  return Vec2(x_ - other.x_, y_ - other.y_);
}

template <typename T> T Vec2<T>::dot(const Vec2<T> &other) const {
  return x_ * other.x_ + y_ * other.y_;
}

template <typename T> Vec2<T> Vec2<T>::operator*(const float scalar) const {
  return Vec2<T>(x_ * scalar, y_ * scalar);
}

// Specialization for float scalar multiplication for int type
template <> Vec2<int> Vec2<int>::operator*(float scalar) const {
  return Vec2<int>(static_cast<int>(std::round(x_ * scalar)),
                   static_cast<int>(std::round(y_ * scalar)));
}

// Reverse scalar multiplication
template <typename T> Vec2<T> operator*(float scalar, const Vec2<T> &vector) {
  return Vec2<T>(scalar * vector.x(), scalar * vector.y());
}

// Specialization for float scalar multiplication for int type (reverse)
template <> Vec2<int> operator*(float scalar, const Vec2<int> &vector) {
  return Vec2<int>(static_cast<int>(std::round(vector.x() * scalar)),
                   static_cast<int>(std::round(vector.y() * scalar)));
}

template <typename T> float Vec2<T>::len() const {
  return sqrt(x_ * x_ + y_ * y_);
}

// Conversion operators
template <> Vec2<float>::operator Vec2<int>() const {
  return Vec2<int>(static_cast<int>(std::round(x_)),
                   static_cast<int>(std::round(y_)));
}

template <> Vec2<int>::operator Vec2<float>() const {
  return Vec2<float>(static_cast<float>(x_), static_cast<float>(y_));
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vec2<T> &vec) {
  os << "(" << vec.x_ << ", " << vec.y_ << ")";
  return os;
}

// Explicit template instantiation for int and float
template class Vec2<int>;
template class Vec2<float>;
