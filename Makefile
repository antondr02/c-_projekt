CXX = clang++
# TODO: Move libraries to $(LIBS)
CXXFLAGS = -std=c++17 -g -Wall -Wextra -Wdeprecated -fsanitize=address -I/usr/include/freetype2 -I./SplineSuite -I./eigen -I./nlopt/build

# Directories
SRC_DIR = .
LIB_DIR = ./SplineSuite
OBJ_DIR = ./obj
BIN_DIR = ./bin
LIBRARY_DIR = ./lib

# Libraries
LIBS = -lglfw -lGL -lX11 -lrt -ldl -lfreetype -lnlopt -lm
TESTLIBS = -lgtest -lgtest_main -lpthread

# Source files
MAIN_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
LIB_SRCS = $(wildcard $(LIB_DIR)/*.cpp)
MAIN_OBJS = $(MAIN_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIB_OBJS = $(LIB_SRCS:$(LIB_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Output binary and library
TARGET = $(BIN_DIR)/Cortado
STATIC_LIB = $(LIBRARY_DIR)/libsplinesuite.a

# Rules
all: $(TARGET) checkstyle

$(TARGET): $(MAIN_OBJS) $(STATIC_LIB) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBRARY_DIR) -lsplinesuite $(LIBS)

$(STATIC_LIB): $(LIB_OBJS) | $(LIBRARY_DIR)
	ar cr $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN_DIR) $(OBJ_DIR) $(LIBRARY_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(LIBRARY_DIR)

checkstyle:
	clang-format-14 --dry-run -Werror $(SRC_DIR)/*.h $(SRC_DIR)/*.cpp $(LIB_DIR)/*.h $(LIB_DIR)/*.cpp

format:
	clang-format-14 -i $(SRC_DIR)/*.h $(SRC_DIR)/*.cpp $(LIB_DIR)/*.h $(LIB_DIR)/*.cpp

.PHONY: all clean
