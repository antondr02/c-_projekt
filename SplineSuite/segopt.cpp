#include "segopt.h"

CubicBezier assembleBezier(const Vec2<int> &startAnchor,
                           const Vec2<int> &endAnchor,
                           Vec2<float> &start_tangent, Vec2<float> &end_tangent,
                           float start_distance, float end_distance) {
  Vec2<int> startHandle =
      startAnchor + static_cast<Vec2<int>>(start_tangent * start_distance);
  Vec2<int> endHandle =
      endAnchor - static_cast<Vec2<int>>(end_tangent * end_distance);
  return CubicBezier(startAnchor, startHandle, endHandle, endAnchor);
}

CubicBezier segopt(const std::vector<Vec2<int>> &segment_vertices,
                   Vec2<float> start_tangent, Vec2<float> end_tangent,
                   std::vector<float> initial_guess) {

  // Convert initial_guess from float to double
  std::vector<double> double_initial_guess(initial_guess.begin(),
                                           initial_guess.end());

  // Create the object for data required by opdist
  OptimizationData opt_data = {segment_vertices, start_tangent, end_tangent};

  // Create the nlopt optimizer object with the Subplex algorithm (a variant of
  // the Nelder-Mead simplex algorithm which is said to be more efficient and
  // robust) in a 2 dimensional solution space
  nlopt::opt opt(nlopt::LN_SBPLX, 2);

  // Set the objective function to be minimized by passing a function pointer
  opt.set_min_objective(opdist, &opt_data);

  // Set a relative tolerance as stopping criterium
  opt.set_xtol_rel(1e-4);

  // Set upper bounds
  std::vector<double> upper_bounds = {initial_guess[0] * 3.5,
                                      initial_guess[1] * 3.5};
  opt.set_upper_bounds(upper_bounds);

  // Vector to store the optimized variables
  std::vector<float> optimized_distances(2);

  // Perform the optimization
  double minf; // Variable to store the minimum objective value
  nlopt::result result = opt.optimize(double_initial_guess, minf);
  (void)result; // Mark as unused. Keep for debugging or error handling

  // Cast the optimized distances to float and store them
  for (size_t i = 0; i < double_initial_guess.size(); ++i) {
    optimized_distances[i] = static_cast<float>(double_initial_guess[i]);
  }

  // Compute the control points
  CubicBezier bezier = assembleBezier(
      segment_vertices[0], segment_vertices[segment_vertices.size() - 1],
      start_tangent, end_tangent, optimized_distances[0],
      optimized_distances[1]);

  return bezier;
}

double opdist(const std::vector<double> &double_distances,
              [[maybe_unused]] std::vector<double> &grad, void *f_data) {

  // Access data members
  OptimizationData *data = static_cast<OptimizationData *>(f_data);
  auto &segment_vertices = data->segment_vertices;
  auto &start_tangent = data->start_tangent;
  auto &end_tangent = data->end_tangent;

  // Convert double distances to float distances
  std::vector<float> distances(double_distances.begin(),
                               double_distances.end());

  CubicBezier bezier = assembleBezier(
      segment_vertices[0], segment_vertices[segment_vertices.size() - 1],
      start_tangent, end_tangent, distances[0], distances[1]);

  // Convert the control points to types used by GraphicsGems
  // Explicit conversion to double to avoid narrowing conversion errors
  Point2 GGstartAnchor = {static_cast<double>(bezier.startAnchor().x()),
                          static_cast<double>(bezier.startAnchor().y())};
  Point2 GGstartHandle = {static_cast<double>(bezier.startHandle().x()),
                          static_cast<double>(bezier.startHandle().y())};
  Point2 GGendHandle = {static_cast<double>(bezier.endHandle().x()),
                        static_cast<double>(bezier.endHandle().y())};
  Point2 GGendAnchor = {static_cast<double>(bezier.endAnchor().x()),
                        static_cast<double>(bezier.endAnchor().y())};

  // Define a C-style array to hold these points
  Point2 c_bezier[4] = {
      GGstartAnchor, // Start Anchor
      GGstartHandle, // Start Handle
      GGendHandle,   // End Handle
      GGendAnchor    // End Anchor
  };

  // Accumulative variable that stores the total error between sample points and
  // the curve
  float total_error = 0;

  for (size_t i = 1; i < segment_vertices.size() - 1; ++i) {
    Point2 GGvertex = {static_cast<double>(segment_vertices[i].x()),
                       static_cast<double>(segment_vertices[i].y())};
    Point2 GGnearestPointOnCurve = NearestPointOnCurve(GGvertex, c_bezier);
    Vec2<float> vertex(static_cast<float>(GGvertex.x),
                       static_cast<float>(GGvertex.y));
    Vec2<float> nearestPointOnCurve(
        static_cast<float>(GGnearestPointOnCurve.x),
        static_cast<float>(GGnearestPointOnCurve.y));
    float error = (vertex - nearestPointOnCurve).len();
    total_error += error * error;
  }
  return static_cast<double>(total_error);
}
