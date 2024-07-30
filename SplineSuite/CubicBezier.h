#pragma once
#include "Vec2.h"

// A simple function for point interpolation. t must be within [0.0, 1.0].
Vec2<int> interpolate(Vec2<int> startPoint, Vec2<int> endPoint, float t);

// A class representing a single cubic bezier curve
class CubicBezier {
private:
  // Defining Points of the cubic bezier. The startAnchor_ is the starting point
  // of the curve (where t=0) the endAnchor_ is the ending point of the curve
  // (where t=1). The startHandle is the point that influences the curve near
  // the startAnchor. The endHandle is the point that influences the curve near
  // the endAnchor.
  Vec2<int> startAnchor_;
  Vec2<int> startHandle_;
  Vec2<int> endHandle_;
  Vec2<int> endAnchor_;

public:
  // Default constructor
  CubicBezier() = default;
  // Constructor from Vectors
  CubicBezier(Vec2<int> startAnchor, Vec2<int> startHandle, Vec2<int> endHandle,
              Vec2<int> endAnchor)
      : startAnchor_(startAnchor), startHandle_(startHandle),
        endHandle_(endHandle), endAnchor_(endAnchor){};

  // Returns the point on the cubic bezier corresponding to the given t value
  // using De Casteljau's algorithm. t must be within [0.0, 1.0].
  Vec2<int> getPoint(float t) const;

  // Moves the entire CubicBezier based on the given vector
  void translate(const Vec2<int> &offset);

  // Getters
  Vec2<int> startAnchor() const { return startAnchor_; }
  Vec2<int> startHandle() const { return startHandle_; }
  Vec2<int> endHandle() const { return endHandle_; }
  Vec2<int> endAnchor() const { return endAnchor_; }

  // Setters
  void setStartAnchor(const Vec2<int> &point) { startAnchor_ = point; }
  void setStartHandle(const Vec2<int> &point) { startHandle_ = point; }
  void setEndHandle(const Vec2<int> &point) { endHandle_ = point; }
  void setEndAnchor(const Vec2<int> &point) { endAnchor_ = point; }
};
