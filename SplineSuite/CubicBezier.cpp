#include "CubicBezier.h"
#include <cmath>
#include <stdexcept>

Vec2<int> interpolate(Vec2<int> startPoint, Vec2<int> endPoint, float t) {
  if (t < 0.0f || t > 1.0f) {
    throw std::out_of_range("Parameter t is out of range [0.0, 1.0]");
  }
  int x = static_cast<int>(round((1 - t) * startPoint.x() + t * endPoint.x()));
  int y = static_cast<int>(round((1 - t) * startPoint.y() + t * endPoint.y()));
  return Vec2<int>(x, y);
}

Vec2<int> CubicBezier::getPoint(float t) const {
  if (t < 0.0f || t > 1.0f) {
    throw std::out_of_range("Parameter t is out of range [0.0, 1.0]");
  }
  Vec2<int> a = interpolate(startAnchor_, startHandle_, t);
  Vec2<int> b = interpolate(startHandle_, endHandle_, t);
  Vec2<int> c = interpolate(endHandle_, endAnchor_, t);

  Vec2<int> d = interpolate(a, b, t);
  Vec2<int> e = interpolate(b, c, t);

  return interpolate(d, e, t);
}

void CubicBezier::translate(const Vec2<int> &offset) {
  startAnchor_ += offset;
  startHandle_ += offset;
  endHandle_ += offset;
  endAnchor_ += offset;
}