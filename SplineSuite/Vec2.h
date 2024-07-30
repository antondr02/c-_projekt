#pragma once
#include <cmath>
#include <iostream>

// A simple 2d vector that supports basic operations and methods
template <typename T> class Vec2 {
private:
  T x_, y_;

public:
  // Default constructor
  Vec2(T x = 0, T y = 0) : x_(x), y_(y){};
  // ~Vec2() = default;

  // Vector addition
  Vec2<T> operator+(const Vec2 &other) const;

  Vec2<T> &operator+=(const Vec2 &other);

  // Vector subtraction
  Vec2<T> operator-(const Vec2 &other) const;

  // Dot product
  T dot(const Vec2 &other) const;

  // Scalar multiplication (member function)
  Vec2<T> operator*(const float scalar) const;

  // Friend function for scalar multiplication with scalar on the left
  friend Vec2<T> operator*(float scalar, const Vec2<T> &point);

  // Get the vector length
  float len() const;

  // Getters
  T x() const { return x_; }
  T y() const { return y_; }

  // Conversion operators
  explicit operator Vec2<float>() const;
  explicit operator Vec2<int>() const;

  // Overload the << operator for output
  friend std::ostream &operator<<(std::ostream &os, const Vec2<T> &vec);
};
