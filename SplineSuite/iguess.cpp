#include "iguess.h"

std::pair<BezierSpline, IguessData> iguess(std::vector<Vec2<int>> vertices,
                                           size_t num_knots,
                                           std::vector<size_t> knot_indices) {
  size_t num_vertices = vertices.size();

  if (vertices.size() < 5) {
    BezierSpline minimal_spline;
    // Construct a degenerated CubicBezier - a line basically
    CubicBezier minimal_bezier(vertices[0], vertices[0],
                               vertices[vertices.size() - 1],
                               vertices[vertices.size() - 1]);
    minimal_spline.add(minimal_bezier);

    // Populate the OptimizationData sStruct
    Vec2<float> connection_vector =
        static_cast<Vec2<float>>(vertices[vertices.size() - 1] - vertices[0]);
    Vec2<float> unit_vector = connection_vector * (1 / connection_vector.len());
    IguessData minimal_opt_data = {{0, vertices.size() - 1},
                                   {vertices[0], vertices[vertices.size() - 1]},
                                   {unit_vector, unit_vector},
                                   {connection_vector.len() / 3}};

    return std::make_pair(minimal_spline, minimal_opt_data);
  }

  if (knot_indices.size() == 0) {
    knot_indices = compute_default_knots(num_vertices, num_knots);
  }

  // Call to compute the distance between successive knot points
  std::vector<float> knot_distances = compute_distances(vertices, knot_indices);

  // Call to compute the unit tangent vectors at the knot points
  std::vector<Vec2<float>> unit_tangents =
      compute_unit_tangents(vertices, knot_indices);

  std::vector<float> divided_distances = knot_distances;

  // Divide each element by 3 to get a good fit
  for (auto &distance : divided_distances) {
    distance /= 3.0f;
  }

  // Populate the OptimizationData struct
  IguessData opt_data = {knot_indices, vertices, unit_tangents,
                         divided_distances};

  BezierSpline spline = compute_control_points(vertices, knot_indices,
                                               knot_distances, unit_tangents);

  return std::make_pair(spline, opt_data);
}

std::vector<size_t> compute_default_knots(size_t num_vertices,
                                          size_t num_knots) {
  float interval = static_cast<float>(num_vertices - 1) /
                   (num_knots - 1); // Ensure floating-point division

  std::vector<size_t> indices;
  indices.reserve(num_knots); // Reserve memory for num_knots elements
  for (size_t i = 0; i < num_knots; ++i) {
    indices.push_back(static_cast<size_t>(std::round(i * interval)));
  }

  return indices;
}

std::vector<float> compute_distances(std::vector<Vec2<int>> vertices,
                                     std::vector<size_t> knot_indices) {
  std::vector<float> distances;
  distances.reserve(knot_indices.size() - 1);

  for (size_t i = 0; i < knot_indices.size() - 1; ++i) {
    Vec2<int> connecting_vector =
        vertices[knot_indices[i + 1]] - vertices[knot_indices[i]];
    distances.push_back(connecting_vector.len());
  }

  return distances;
}

std::vector<Vec2<float>>
compute_unit_tangents(std::vector<Vec2<int>> vertices,
                      std::vector<size_t> knot_indices) {
  size_t num_vertices = vertices.size();
  size_t num_knots = knot_indices.size();

  std::vector<Vec2<float>> unit_tangents;
  unit_tangents.reserve(num_knots);

  for (size_t i = 0; i < num_knots; ++i) {

    // Index in vertices of the first vertex used for approximation
    size_t start_vertex_index;

    // Index of the current knot in vertices relative to start_vertex_index
    size_t local_knot_index;

    // Adjust indices for boundary conditions
    if (i == 0) {
      start_vertex_index = 0;
      local_knot_index = 0;
    } else if (i == num_knots - 1) {
      start_vertex_index = num_vertices - 5;
      local_knot_index = 4;
    } else {
      start_vertex_index = knot_indices[i] - 2;
      local_knot_index = 2;
    }

    // Create an array to store the five approximation points
    std::array<Vec2<int>, 5> approximation_points;
    for (size_t j = 0; j < 5; ++j) {
      approximation_points[j] = vertices[start_vertex_index + j];
    }

    // Compute the distances between consecutive approximation points
    std::array<float, 4> approximation_point_distances;
    for (size_t j = 0; j < 4; ++j) {
      approximation_point_distances[j] = (vertices[start_vertex_index + j + 1] -
                                          vertices[start_vertex_index + j])
                                             .len();
    }

    // Compute the cumulative chord lengths
    std::array<float, 5> chord_length;
    chord_length[0] = 0.0f;
    for (size_t j = 1; j < 5; ++j) {
      chord_length[j] =
          chord_length[j - 1] + approximation_point_distances[j - 1];
    }

    // Construct the matrix A holding the chord lengths
    Eigen::Matrix<float, 5, 3> A;
    for (size_t j = 0; j < 5; ++j) {
      A(j, 0) = 1.0;
      A(j, 1) = chord_length[j];
      A(j, 2) = chord_length[j] * chord_length[j];
    }

    // Construct the matrix B holding the points used for approximation
    Eigen::Matrix<float, 5, 2> B;
    for (size_t j = 0; j < 5; ++j) {
      B(j, 0) = approximation_points[j].x();
      B(j, 1) = approximation_points[j].y();
    }

    // Solve the least squares problem: A * coefficients = B;
    Eigen::Matrix<float, 3, 2> coefficients = A.colPivHouseholderQr().solve(B);

    // Compute the tangent vector at the local knot index using the derivative
    // of the parametric quadratic curve
    Eigen::Vector2f tangent_vector;
    tangent_vector[0] = coefficients(1, 0) +
                        2 * coefficients(2, 0) * chord_length[local_knot_index];
    tangent_vector[1] = coefficients(1, 1) +
                        2 * coefficients(2, 1) * chord_length[local_knot_index];

    // Normalize the tangent vector
    tangent_vector.normalize();

    // Add tangent vector as Vec2<int> to unit_tangents
    Vec2<float> unit_tangent(tangent_vector[0], tangent_vector[1]);
    unit_tangents.push_back(unit_tangent);
  }

  return unit_tangents;
}

BezierSpline compute_control_points(std::vector<Vec2<int>> vertices,
                                    std::vector<size_t> knot_indices,
                                    std::vector<float> knot_distances,
                                    std::vector<Vec2<float>> unit_tangents) {
  size_t num_knots = knot_indices.size();

  BezierSpline iguess_spline;

  for (size_t i = 0; i < num_knots - 1; ++i) {
    Vec2<int> startAnchor = vertices[knot_indices[i]];
    Vec2<int> startHandle =
        startAnchor +
        static_cast<Vec2<int>>(unit_tangents[i] * (knot_distances[i] / 3));
    Vec2<int> endAnchor = vertices[knot_indices[i + 1]];
    Vec2<int> endHandle =
        endAnchor -
        static_cast<Vec2<int>>(unit_tangents[i + 1] * (knot_distances[i] / 3));

    iguess_spline.add(
        CubicBezier(startAnchor, startHandle, endHandle, endAnchor));
  }

  return iguess_spline;
}
