#include "Cortado.h"
#include <cmath>
#include <unistd.h>

void drawPoint(Vec2<int> pos, OpenGLTerminalManager *terminalManager,
               int radius) {
  for (int i = -radius; i <= radius; i++) {
    for (int j = -radius; j <= radius; j++) {
      terminalManager->drawPixel(pos.y() + j, pos.x() + i, true, 1.0);
    }
  }
}

void showHandles(const CubicBezier &bezier,
                 OpenGLTerminalManager *terminalManager) {
  drawPoint(Vec2<int>(bezier.startHandle().x(), bezier.startHandle().y()),
            terminalManager);
  drawPoint(Vec2<int>(bezier.endHandle().x(), bezier.endHandle().y()),
            terminalManager);
  terminalManager->drawLine({bezier.startAnchor(), bezier.startHandle()});
  terminalManager->drawLine({bezier.endAnchor(), bezier.endHandle()});
}

void showHandles(BezierSpline &spline, OpenGLTerminalManager *terminalManager) {
  for (const CubicBezier &bezier : spline.getBeziers()) {
    showHandles(bezier, terminalManager);
  }
}

void drawRectangle(const Vec2<int> &corner1, const Vec2<int> &corner2,
                   OpenGLTerminalManager *terminalManager,
                   float intensity = 1.0f, bool bounded = false) {
  // Find out the actual top_left and bottom_right
  int left = std::min(corner1.x(), corner2.x());
  int right = std::max(corner1.x(), corner2.x());
  int top = std::min(corner1.y(), corner2.y());
  int bottom = std::max(corner1.y(), corner2.y());

  // Draw the rectangle
  for (int y = top; y <= bottom; ++y) {
    for (int x = left; x <= right; ++x) {
      if (y >= 0 && y < terminalManager->numRows() && x >= 0 &&
          x < terminalManager->numCols()) {
        terminalManager->drawPixel(y, x, true, intensity);
      }
    }
  }

  // Draw the bounding lines if bounded is true
  if (bounded) {
    terminalManager->drawLine({Vec2<int>(left, top), Vec2<int>(right, top),
                               Vec2<int>(right, bottom),
                               Vec2<int>(left, bottom), Vec2<int>(left, top)});
  }
}

void clearPixels(OpenGLTerminalManager *terminalManager) {
  for (int row = 0; row < terminalManager->numRows(); row++) {
    for (int col = 0; col < terminalManager->numCols(); col++) {
      terminalManager->drawPixel(row, col, false, 1.0);
    }
  }
}

bool isPointInRectangle(const Vec2<int> &point, const Vec2<int> &corner1,
                        const Vec2<int> &corner2) {
  // Find out the actual top_left and bottom_right
  int left = std::min(corner1.x(), corner2.x());
  int right = std::max(corner1.x(), corner2.x());
  int top = std::min(corner1.y(), corner2.y());
  int bottom = std::max(corner1.y(), corner2.y());

  // Check if the point is within the bounds
  return point.x() >= left && point.x() <= right && point.y() >= top &&
         point.y() <= bottom;
}

Cortado::Cortado(OpenGLTerminalManager *terminalManager)
    : terminalManager_(terminalManager), selectedTool_(Tool::DRAW) {}

void Cortado::run() {
  while (true) {
    UserInput2 userInput = terminalManager_->getUserInput();
    if (userInput.keycode_ == 'q') {
      break;
    }
    if (userInput.keycode_ == 't') {
      toggleTool();
    }
    if (userInput.keycode_ == 'z') {
      if (!splines_.empty()) {
        splines_.pop_back();
        selectSplines(startSelect, endSelect);
        draw(T_STEP);
        terminalManager_->refresh();
      }
    }

    // Handle the user input based on the selected tool
    switch (selectedTool_) {
    case Tool::DRAW:
      handleDraw(userInput);
      break;
    case Tool::SELECT:
      handleSelect(userInput);
      break;
    }

    mouseLastFrame_ = userInput.isMousePressed_;
    usleep(100);
  }
}

void Cortado::draw(float t_step) {
  // Create stack and fill with selectedSplineIndices_ (backwards)
  std::stack<size_t> selectedIndices;
  for (size_t i = selectedSplineIndices_.size(); i > 0; --i) {
    selectedIndices.push(selectedSplineIndices_[i - 1]);
  }

  for (size_t i = 0; i < splines_.size(); ++i) {
    std::vector<CubicBezier> beziers = splines_[i].getBeziers();
    for (const auto &bezier : beziers) {
      //  Calculate the number of points
      size_t num_points = static_cast<size_t>(std::ceil(1.0f / t_step)) + 1;
      std::vector<Vec2<int>> vertices;
      vertices.reserve(num_points); // Preallocate space for the vertices

      for (float t = 0.0f; t <= 1.0f; t += t_step) {
        Vec2<int> p = bezier.getPoint(t);
        if (!selectedIndices.empty() && i == selectedIndices.top()) {
          vertices.push_back(p + move_offset);
        } else {
          vertices.push_back(p);
        }
      }
      terminalManager_->drawLine(vertices);
    }
    if (!selectedIndices.empty() && i == selectedIndices.top()) {
      selectedIndices.pop();
    }
  }
}

void Cortado::verticesToSpline() {
  int num_knots = std::max(3, static_cast<int>(currentStroke_.size() / 10));

  // Get the initial guess for the curve (just the optimization data is required
  // here)
  auto [spline, opt_data] = iguess(currentStroke_, num_knots);
  // (void)spline; // Displaying initial guess not necessary

  // Unpack opt_data
  std::vector<size_t> knot_indices = opt_data.knot_indices;
  std::vector<Vec2<int>> vertices = opt_data.vertices;
  std::vector<Vec2<float>> unit_tangents = opt_data.unit_tangents;
  std::vector<float> initial_distances = opt_data.initial_distances;

  // Skip optimization for minimal spline (line)
  if (vertices.size() < 5) {
    splines_.push_back(spline);
    // Clear currentStroke_ to prepare for the next one
    currentStroke_.clear();

    // Draw all the splines
    draw(0.02f);

    // Render the screen content
    terminalManager_->refresh();
    return;
  }

  // Spline object for the new brush stroke
  BezierSpline optimizedSpline;

  // Perform the segment-wise optimization
  for (size_t i = 0; i < knot_indices.size() - 1; ++i) {
    // Create a vector for the sample vertices of this segment
    // Both indices inclusive
    size_t start_index = knot_indices[i];
    size_t end_index = knot_indices[i + 1] + 1; // Make end_index exclusive
    std::vector<Vec2<int>> segment_vertices(vertices.begin() + start_index,
                                            vertices.begin() + end_index);

    // Create a vector holding the initial guess distances
    std::vector<float> initial_guess = {initial_distances[i],
                                        initial_distances[i]};

    // Get the segment-wise optimized cubic bezier
    CubicBezier optimized_bezier = segopt(segment_vertices, unit_tangents[i],
                                          unit_tangents[i + 1], initial_guess);
    optimizedSpline.add(optimized_bezier);
  }

  splines_.push_back(optimizedSpline);

  // Clear currentStroke_ to prepare for the next one
  currentStroke_.clear();

  // Draw all the splines
  draw(T_STEP);

  // Render the screen content
  terminalManager_->refresh();
}

void Cortado::translateSelectedSplines(const Vec2<int> &offset) {
  for (const auto &index : selectedSplineIndices_) {
    splines_[index].translate(offset);
  }
}

void Cortado::toggleTool() {
  if (selectedTool_ == Tool::DRAW) {
    if (!currentStroke_.empty()) {
      verticesToSpline();
    }
    drawRectangle(startSelect, endSelect, terminalManager_, SELECT_INTENSITY);
    draw(T_STEP);
    terminalManager_->refresh();
    selectedTool_ = Tool::SELECT;
    selectSplines(startSelect, endSelect);
  } else {
    // Clear selections and draw splines
    clearPixels(terminalManager_);
    draw(T_STEP);
    terminalManager_->refresh();
    selectedTool_ = Tool::DRAW;
  }
  mouseLastFrame_ = false;
}

void Cortado::selectSplines(const Vec2<int> &topLeft,
                            const Vec2<int> &bottomRight) {
  // Clear previously selected splines
  selectedSplineIndices_.clear();

  // Loop through every spline
  for (size_t i = 0; i < splines_.size(); ++i) {
    const BezierSpline &spline = splines_[i];
    bool splineSelected = false;

    // Loop through the contained CubicBeziers
    for (const auto &bezier : spline.getBeziers()) {
      // Check if bezier anchors are contained within the selection
      if (isPointInRectangle(bezier.startAnchor(), topLeft, bottomRight) ||
          isPointInRectangle(bezier.endAnchor(), topLeft, bottomRight)) {
        splineSelected = true;
        break;
      }
    }
    if (splineSelected) {
      selectedSplineIndices_.push_back(i);
      continue;
    }
    // Loop through contained CubicBeziers to check for actual curve values
    for (const auto &bezier : spline.getBeziers()) {
      for (float t = 0; t <= 1.0f; t += T_STEP_SELECT) {
        if (isPointInRectangle(bezier.getPoint(t), topLeft, bottomRight)) {
          splineSelected = true;
          break;
        }
      }
      if (splineSelected) {
        selectedSplineIndices_.push_back(i);
        break;
      }
    }
  }
}

void Cortado::handleDraw(UserInput2 &userInput) {
  if (userInput.isMousePressed_) {
    Vec2<int> mousePos(userInput.currentMouseX_, userInput.currentMouseY_);
    // Sample only when mouse has moved enough
    if (currentStroke_.size() == 0 ||
        (currentStroke_[currentStroke_.size() - 1] - mousePos).len() >=
            LINE_SAMPLE_RATE) {
      // Create a sample
      currentStroke_.push_back(Vec2<int>(mousePos));

      // Draw the line with the new vertex
      terminalManager_->drawLine(currentStroke_);
      draw(T_STEP);
      terminalManager_->refresh();
    }
  } else if (mouseLastFrame_) {
    // Compute and add the approximation spline
    verticesToSpline();
  }
}

void Cortado::handleSelect(UserInput2 &userInput) {
  if (userInput.isMousePressed_) {
    // Clamp the mouse coordinates to be within the valid screen range
    int clampedX = std::max(
        0, std::min(userInput.currentMouseX_, terminalManager_->numCols() - 1));
    int clampedY = std::max(
        0, std::min(userInput.currentMouseY_, terminalManager_->numRows() - 1));

    // Define the clamped mouse position
    Vec2<int> mousePos(clampedX, clampedY);

    // Check if the mouse is inside the current selection
    bool isMouseInSelection = isPointInRectangle(
        mousePos, startSelect + move_offset, endSelect + move_offset);

    if (!mouseLastFrame_) {
      // Mouse pressed for the first time
      if (!selectionComplete || !isMouseInSelection) {
        // Start a new selection if no selection is complete or if the mouse is
        // outside the existing selection
        startSelect = mousePos;
        selectionComplete = false;
        move_ongoing = false;
      } else {
        // Otherwise, initiate a move operation
        move_ongoing = true;
        move_anchor = mousePos;
      }
    }

    if (move_ongoing) {
      // Handle moving the selection
      move_offset = mousePos - move_anchor;
    } else {
      // Handle resizing the selection
      endSelect = mousePos;
    }

    // Display the selection boundaries AND the splines
    clearPixels(terminalManager_);
    drawRectangle(startSelect + move_offset, endSelect + move_offset,
                  terminalManager_, SELECT_INTENSITY);
    draw(T_STEP);
    terminalManager_->refresh();
  } else if (mouseLastFrame_) {
    // Mouse released
    if (move_ongoing) {
      // Finalize the move operation
      startSelect += move_offset;
      endSelect += move_offset;
      translateSelectedSplines(move_offset);
      move_offset = Vec2<int>(0, 0);
      move_ongoing = false;
    }

    selectionComplete = true;
    selectSplines(startSelect, endSelect);

    // Render splines on the screen
    draw(T_STEP);
    terminalManager_->refresh();
  }
}
