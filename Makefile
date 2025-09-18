CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -Iinclude -fPIC
SRC = src/isotopy_graph.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = libisotopy.so

TEST_SRC = tests/isotopy_graph.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_BIN = isotopy_graph_test

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -shared -o $@ $^

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TARGET) $(TEST_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy


TEST_BATCH_SRC = tests/isotopy_graph_batch.cpp
TEST_BATCH_OBJ = $(TEST_BATCH_SRC:.cpp=.o)
TEST_FULL_BIN = isotopy_graph_test_full

test_full: $(TEST_FULL_BIN)
	./$(TEST_FULL_BIN)

$(TEST_FULL_BIN): $(TARGET) $(TEST_OBJ) $(TEST_BATCH_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy

BENCHMARK_SRC = tests/isotopy_graph_benchmark.cpp
BENCHMARK_OBJ = $(BENCHMARK_SRC:.cpp=.o)
BENCHMARK_BIN = isotopy_graph_benchmark

benchmark: $(BENCHMARK_BIN)
	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): $(TARGET) $(BENCHMARK_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy

PROFILING_SRC = tests/isotopy_graph_profiling.cpp
PROFILING_OBJ = $(PROFILING_SRC:.cpp=.o)
PROFILING_BIN = isotopy_graph_profiling

profiling: $(PROFILING_BIN)
	./$(PROFILING_BIN)

$(PROFILING_BIN): $(TARGET) $(PROFILING_OBJ)
	$(CXX) -o $@ $^ -L. -Wl,-rpath,\$$ORIGIN -lisotopy

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) $(TEST_FULL_BIN)

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)
