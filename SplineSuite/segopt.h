#pragma once
#include "CubicBezier.h"
#include "GraphicsGems.h"
#include "Vec2.h"
#include <cmath>
#include <cstddef>
#include <nlopt.hpp>
#include <utility>
#include <vector>

// An object that holds additional user-defined data required by the objective
// function 'opdist'
struct OptimizationData {
  std::vector<Vec2<int>> segment_vertices;
  Vec2<float> start_tangent;
  Vec2<float> end_tangent;
};

CubicBezier segopt(const std::vector<Vec2<int>> &segment_vertices,
                   Vec2<float> start_tangent, Vec2<float> end_tangent,
                   std::vector<float> initial_guess);

// Computes the sum of errors between each segment vertex and its closest point
// on the curve
// TODO: Find an existing implementation for the 'closest point on curve
// problem'
double opdist(const std::vector<double> &distances, std::vector<double> &grad,
              void *f_data);
