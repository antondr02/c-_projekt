#pragma once
#include "CubicBezier.h"
#include <cstddef>
#include <vector>

class BezierSpline {
public:
  BezierSpline() = default;

  // Add a new cubic bezier to the spline
  void add(const CubicBezier &bezier);

  // Moves the entire BezierSpline based on the given vector
  void translate(const Vec2<int> &offset);

  // Getters
  const std::vector<CubicBezier> &getBeziers() const { return beziers_; }

  // Setter
  void setBezierAt(size_t index, const CubicBezier &bezier);

private:
  std::vector<CubicBezier> beziers_;
};
