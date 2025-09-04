CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -Iinclude -fPIC
SRC = src/isotopy_graph.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = libisotopy.so

TEST_SRC = tests/isotopy_graph.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_BIN = isotopy_graph

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -shared -o $@ $^

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TARGET) $(TEST_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy


BENCHMARK_SRC = tests/benchmark.cpp
BENCHMARK_OBJ = $(BENCHMARK_SRC:.cpp=.o)
BENCHMARK_BIN = isotopy_graph

	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): $(TARGET) $(BENCHMARK_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN)

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)


