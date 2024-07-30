#pragma once
#include "./OpenGLTerminalManager.h"
#include "BezierSpline.h"
#include "CubicBezier.h"
#include "Vec2.h"
#include "iguess.h"
#include "segopt.h"
#include <iostream>
#include <stack>
#include <vector>

// Macros

#define LINE_SAMPLE_RATE 2
#define T_STEP 0.02f
#define T_STEP_SELECT 0.05f
#define SELECT_INTENSITY 0.01f

// Utility functions
void drawPoint(Vec2<int> pos, OpenGLTerminalManager *terminalManager,
               int radius = 3);

void showHandles(const CubicBezier &bezier,
                 OpenGLTerminalManager *terminalManager);

void showHandles(BezierSpline &spline, OpenGLTerminalManager *terminalManager);

void drawRectangle(const Vec2<int> &corner1, const Vec2<int> &corner2,
                   OpenGLTerminalManager *terminalManager, float intensity,
                   bool bounded);

void clearPixels(OpenGLTerminalManager *terminalManager);

bool isPointInRectangle(const Vec2<int> &point, const Vec2<int> &corner1,
                        const Vec2<int> &corner2);

enum class Tool { DRAW, SELECT };

class Cortado {
public:
  Cortado(OpenGLTerminalManager *terminalManager);
  ~Cortado() = default;

  // Method to start the application, handle user input and trigger re-rendering
  // events
  void run();

  void draw(float t_step);

private:
  std::vector<BezierSpline> splines_;
  OpenGLTerminalManager *terminalManager_;
  Tool selectedTool_;

  // Stores a truth value about whether the mouse was pressed last frame
  bool mouseLastFrame_ = false;

  // Buffer for the ongoing brush stroke (converted to BezierSpline after
  // completion)
  std::vector<Vec2<int>> currentStroke_;

  // Method to convert the collection of individual vertices (here
  // currentStroke_) to a BezierSpline that can be appended to splines
  void verticesToSpline();

  // Select points
  Vec2<int> startSelect;
  Vec2<int> endSelect;

  // Offset of onging move
  Vec2<int> move_offset = Vec2<int>(0, 0);

  // Anchor for the move operation
  Vec2<int> move_anchor;

  // Is the selection complete or ongoing
  bool selectionComplete = false;

  // Is a move occuring at this time
  bool move_ongoing = false;

  // Indices in 'splines_' of selected splines
  std::vector<size_t> selectedSplineIndices_;

  // Moves every selected spline by specified offset
  void translateSelectedSplines(const Vec2<int> &offset);

  // Updates the field 'selectedSplineIndices_' based on the boundaries of the
  // rectangle
  void selectSplines(const Vec2<int> &topLeft, const Vec2<int> &bottomRight);

  // Utilities
  void toggleTool();

  // Subroutines
  void handleDraw(UserInput2 &userInput);

  void handleSelect(UserInput2 &userInput);
};
