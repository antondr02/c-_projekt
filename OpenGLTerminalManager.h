// Copyright 2022, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Axel Lehmann <lehmann@cs.uni-freiburg.de>,
//         Claudius Korzen <korzen@cs.uni-freiburg.de>.
//         Johannes Kalmbach <kalmbach@cs.uni-freiburg.de>.
//
#ifndef OPENGLTERMINALMANAGER_H_
#define OPENGLTERMINALMANAGER_H_

#include <GLFW/glfw3.h>
#include <cstddef>
#include <ft2build.h>
#include <queue>
#include <string>
#include <string_view>
#include FT_FREETYPE_H
#include "./Cell.h"
#include "./Vec2.h"

// Class for the input from the user.
class UserInput2 {
public:
  bool isKeyUp();
  bool isKeyDown();
  bool isKeyLeft();
  bool isKeyRight();
  // The code of the key pressed.
  int keycode_;
  // Was the event a mousecklick.
  bool isMouseclick_;
  // Is the left mouse button currently pressed
  bool isMousePressed_;
  // If the event was a mousecklick, then the coordinates
  // of the mouseclick are stored here.
  int mouseX_ = -1;
  int mouseY_ = -1;
  // Current positoion of the mouse independent of mouse clicks
  int currentMouseX_ = -1;
  int currentMouseY_ = -1;
};

// A 2D cell with size_t coordinates.
class GlCell {
public:
  size_t row_;
  size_t col_;

  GlCell(size_t row, size_t col) : row_{row}, col_{col} {}
  // To store the GlCells in a hashMap we have
  // to define the equality == between GlCells.
  bool operator==(const GlCell &other) const {
    return row_ == other.row_ && col_ == other.col_;
  }
};

// We use GlCell as the key type to a hash map, so we have to define
// a hash function
template <> struct std::hash<GlCell> {
  size_t operator()(const GlCell &p) const {
    auto hash1 = std::hash<size_t>{}(p.row_);
    auto hash2 = std::hash<size_t>{}(p.col_);
    if (hash1 != hash2) {
      return hash1 ^ hash2;
    }
    // If hash1 == hash2, their XOR is zero.
    return hash1;
  }
};

// A class managing the input and output via the terminal, using ncurses.
class OpenGLTerminalManager {
private:
public:
  // Constructor: initialize the terminal for use with ncurses.
  OpenGLTerminalManager(const std::string &windowName = "OpenGL");

  // Destructor: Clean up the screen.
  ~OpenGLTerminalManager();

  // Get input from the user.
  UserInput2 getUserInput();

  // Draw a "pixel" at the given position and with the given intensity between
  // The intensity has to be in [0.0, 1.0]
  void drawPixel(int row, int col, bool inverse, float intensity);

  // Draw a string at the given position and with the given intensity.
  // The intensity has to be in [0.0, 1.0]
  void drawString(int row, int col, const char *output, float intensity = 1.0);

  // Draw a line with the given vertices
  void drawLine(const std::vector<Vec2<int>> &vertices);

  // Refresh the screen.
  void refresh();

  // Get the dimensions of the screen.
  int numRows() const { return size_; }
  int numCols() const { return size_; }

  // Get the window name
  const std::string &windowName() const { return windowName_; }

private:
  // The title of the window created
  const std::string windowName_;

  // The number of "logical" rows and columns of the screen.
  int size_ = 1000;
  GLFWwindow *window;
  std::queue<unsigned int> keypresses;
  std::queue<std::pair<int, int>> mouseClicks;
  // texture buffer with red values
  std::vector<float> pixels;
  GLuint texture;
  FT_Library ft;
  FT_Face face;
  struct Character {
    unsigned int textureID; // ID handle of the glyph texture
    std::tuple<unsigned int, unsigned int> size; // Size of glyph
    std::tuple<unsigned int, unsigned int>
        bearing;    // Offset from baseline to left/top of glyph
    FT_Pos advance; // Offset to advance to next glyph
  };
  std::unordered_map<char, Character> characters;
  std::unordered_map<GlCell, char> charsToDraw;

  // Store lines to be drawn
  std::vector<std::vector<Vec2<int>>> lines[2]; // Two buffers for lines
  int currentBuffer = 0;                        // Index of the current buffer

  void swapBuffers();

  void renderText(std::string_view, int col, int row);
};

#endif // OPENGLTERMINALMANAGER_H_
