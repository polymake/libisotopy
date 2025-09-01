CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -Iinclude
SRC = src/isotopy_graph.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = libisotopy.so

TEST_SRC = tests/isotopy_graph.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_BIN = isotopy_graph

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(OBJ) $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN)
