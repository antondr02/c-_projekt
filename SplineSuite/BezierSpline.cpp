#include "BezierSpline.h"
#include <iostream>

void BezierSpline::add(const CubicBezier &bezier) {
  beziers_.push_back(bezier);
}

void BezierSpline::translate(const Vec2<int> &offset) {
  for (auto &bezier : beziers_) {
    bezier.translate(offset);
  }
}

void BezierSpline::setBezierAt(size_t index, const CubicBezier &bezier) {
  beziers_[index] = bezier;
}