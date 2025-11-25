CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -Iinclude -fPIC
WITH_ROOTED_TREES ?= 0
BIN_DIR = bin
TEST_DIR = $(BIN_DIR)/tests
BENCHMARK_DIR = $(BIN_DIR)/benchmarks

SRC = src/isotopy_graph.cpp src/viro_format.cpp src/rooted_trees.cpp

OBJ = $(SRC:src/%.cpp=obj/%.o)
TARGET = libisotopy.so

TEST_SRC = tests/isotopy_graph.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_BIN = $(TEST_DIR)/isotopy_graph_test

ROOTED_TREES_TEST_SRC = tests/rooted_trees.cpp
ROOTED_TREES_TEST_OBJ = $(ROOTED_TREES_TEST_SRC:.cpp=.o)
ROOTED_TREES_TEST_BIN = $(TEST_DIR)/rooted_trees_test

VIRO_FORMAT_TEST_SRC = tests/viro_format.cpp
VIRO_FORMAT_TEST_OBJ = $(VIRO_FORMAT_TEST_SRC:.cpp=.o)
VIRO_FORMAT_TEST_BIN = $(TEST_DIR)/viro_format_test

ROOTED_TREES_BENCHMARK_SRC = tests/rooted_trees_benchmark.cpp
ROOTED_TREES_BENCHMARK_OBJ = $(ROOTED_TREES_BENCHMARK_SRC:.cpp=.o)
ROOTED_TREES_BENCHMARK_BIN = $(BENCHMARK_DIR)/rooted_trees_benchmark
ROOTED_TREES_BENCHMARK_DELTA ?= 8

all: $(TARGET) $(BIN_DIR) $(TEST_DIR) $(BENCHMARK_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(TEST_DIR): | $(BIN_DIR)
	mkdir -p $(TEST_DIR)

$(BENCHMARK_DIR): | $(BIN_DIR)
	mkdir -p $(BENCHMARK_DIR)

$(TARGET): $(OBJ)
	$(CXX) -shared -o $@ $^

ifeq ($(WITH_ROOTED_TREES),1)
test: $(TEST_BIN) $(VIRO_FORMAT_TEST_BIN) $(ROOTED_TREES_TEST_BIN)
	./$(TEST_BIN)
	./$(VIRO_FORMAT_TEST_BIN)
	./$(ROOTED_TREES_TEST_BIN)
else
test: $(TEST_BIN) $(VIRO_FORMAT_TEST_BIN)
	./$(TEST_BIN)
	./$(VIRO_FORMAT_TEST_BIN)
endif

$(TEST_BIN): $(TARGET) $(TEST_OBJ) | $(TEST_DIR)
	$(CXX) -o $@ $(TEST_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy

ifeq ($(WITH_ROOTED_TREES),1)
$(ROOTED_TREES_TEST_BIN): $(TARGET) $(ROOTED_TREES_TEST_OBJ) | $(TEST_DIR)
	$(CXX) -o $@ $(ROOTED_TREES_TEST_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy
else
$(ROOTED_TREES_TEST_BIN):
	@echo "rooted_trees_test unavailable (build with WITH_ROOTED_TREES=1)" >&2
	@exit 1
endif

$(VIRO_FORMAT_TEST_BIN): $(TARGET) $(VIRO_FORMAT_TEST_OBJ) | $(TEST_DIR)
	$(CXX) -o $@ $(VIRO_FORMAT_TEST_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy

ifeq ($(WITH_ROOTED_TREES),1)
benchmark_rooted_trees: $(ROOTED_TREES_BENCHMARK_BIN)
	./$(ROOTED_TREES_BENCHMARK_BIN) $(ROOTED_TREES_BENCHMARK_DELTA)
else
benchmark_rooted_trees:
	@echo "benchmark_rooted_trees unavailable (build with WITH_ROOTED_TREES=1)" >&2
	@exit 1
endif

ifeq ($(WITH_ROOTED_TREES),1)
$(ROOTED_TREES_BENCHMARK_BIN): $(TARGET) $(ROOTED_TREES_BENCHMARK_OBJ) | $(BENCHMARK_DIR)
	$(CXX) -o $@ $(ROOTED_TREES_BENCHMARK_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy
else
$(ROOTED_TREES_BENCHMARK_BIN):
	@echo "rooted_trees_benchmark unavailable (build with WITH_ROOTED_TREES=1)" >&2
	@exit 1
endif


TEST_BATCH_SRC = tests/isotopy_graph_batch.cpp
TEST_BATCH_OBJ = obj/isotopy_graph_batch.o
TEST_FULL_BIN = $(BIN_DIR)/isotopy_graph_full

test_full: $(TEST_FULL_BIN)
	./$(TEST_FULL_BIN)

$(TEST_FULL_BIN): $(TARGET) $(TEST_OBJ) $(TEST_BATCH_OBJ) | $(BIN_DIR)
	$(CXX) -o $@ $(TEST_OBJ) $(TEST_BATCH_OBJ) -L. -Wl,-rpath,\$$ORIGIN/.. -lisotopy

BENCHMARK_SRC = tests/isotopy_graph_benchmark.cpp
BENCHMARK_OBJ = $(BENCHMARK_SRC:.cpp=.o)
BENCHMARK_BIN = $(BENCHMARK_DIR)/isotopy_graph_benchmark

benchmark: $(BENCHMARK_BIN)
	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): $(TARGET) $(BENCHMARK_OBJ) | $(BENCHMARK_DIR)
	$(CXX) -o $@ $(BENCHMARK_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy

PROFILING_SRC = tests/isotopy_graph_profiling.cpp
PROFILING_OBJ = $(PROFILING_SRC:.cpp=.o)
PROFILING_BIN = $(BIN_DIR)/isotopy_graph_profiling

profiling: $(PROFILING_BIN)
	./$(PROFILING_BIN)

$(PROFILING_BIN): $(TARGET) $(PROFILING_OBJ) | $(BIN_DIR)
	$(CXX) -o $@ $(PROFILING_OBJ) -L. -Wl,-rpath,\$$ORIGIN/.. -lisotopy

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) \
	$(TEST_FULL_BIN) $(BENCHMARK_OBJ) $(BENCHMARK_BIN) $(PROFILING_OBJ) \
	$(PROFILING_BIN) $(ROOTED_TREES_TEST_OBJ) $(ROOTED_TREES_TEST_BIN) \
	$(ROOTED_TREES_BENCHMARK_OBJ) $(ROOTED_TREES_BENCHMARK_BIN) \
	$(VIRO_FORMAT_TEST_OBJ) $(VIRO_FORMAT_TEST_BIN) \
	obj/*.o wasm_obj/*.o libisotopy_wasm.a
	rm -rf $(BIN_DIR)

 

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
DELTA ?= 8
