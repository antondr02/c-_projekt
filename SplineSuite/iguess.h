#pragma once
#include "BezierSpline.h"
#include "CubicBezier.h"
#include "Vec2.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

struct IguessData {
  std::vector<size_t> knot_indices;
  std::vector<Vec2<int>> vertices;
  std::vector<Vec2<float>> unit_tangents;
  std::vector<float> initial_distances;
};

// Computes an initial guess for the BezierSpline approximation of a line
// defined by an ordered set of vertices
// Invariants:
// - when knot_indices is explicitly specified its length must match
//   num_knots,
// - length of vertices >= 1
// - num_knots >= 1
// TODO: Remove the nullptr arguments
std::pair<BezierSpline, IguessData>
iguess(std::vector<Vec2<int>> vertices, size_t num_knots,
       std::vector<size_t> knot_indices = {});

// Computes default knot positions for iguess based on a
// formula to equally disperse the knots throughout the data.
std::vector<size_t> compute_default_knots(size_t num_vertices,
                                          size_t num_knots);

// Computes the initial distances from the knot points to their adjacent
// control points for the initial guess curve. It returns the vector
// of distances
std::vector<float> compute_distances(std::vector<Vec2<int>> vertices,
                                     std::vector<size_t> knot_indices);

// This function computes the unit tangent vectors at the knot points.
// It uses chord length parameterization to fit a parametric quadratic curve to
// five data points. The unit tangent vectors are approximated by the unit
// tangent vectors for these quadratic functions.
std::vector<Vec2<float>>
compute_unit_tangents(std::vector<Vec2<int>> vertices,
                      std::vector<size_t> knot_indices);

// This function takes the vertices of the line, the knot_point indices, unit
// tangent vectors, distances between successive knot points as input. It then
// assembles an appropriate BezierSpline
BezierSpline compute_control_points(std::vector<Vec2<int>> vertices,
                                    std::vector<size_t> knot_indices,
                                    std::vector<float> knot_distances,
                                    std::vector<Vec2<float>> unit_tangents);
