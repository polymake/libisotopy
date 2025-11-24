CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -Iinclude -fPIC
SRC = src/isotopy_graph.cpp
OBJ = obj/isotopy_graph.o
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
TEST_BATCH_OBJ = obj/isotopy_graph_batch.o
TEST_FULL_BIN = isotopy_graph_full

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
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) \
	$(TEST_FULL_BIN) $(BENCHMARK_OBJ) $(BENCHMARK_BIN) $(PROFILING_OBJ) \
	$(PROFILING_BIN) obj/*.o wasm_obj/*.o libisotopy_wasm.a

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)

# Emscripten build
.PHONY: emscripten

emscripten: clean emscripten_build

emscripten_build: 
	$(MAKE) OBJ="wasm_obj/isotopy_graph.o" TARGET=libisotopy_wasm.a CXX=em++ CXXFLAGS="$(CXXFLAGS)" libisotopy_wasm.a

libisotopy_wasm.a: $(OBJ)
	ar rcs $@ $(OBJ)

obj/:
	mkdir -p obj

wasm_obj/:
	mkdir -p wasm_obj
# Pattern rule for native object files
obj/%.o: src/%.cpp | obj/
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: tests/%.cpp | obj/
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pattern rule for wasm object files
wasm_obj/%.o: src/%.cpp | wasm_obj/
	em++ $(CXXFLAGS) -c $< -o $@
