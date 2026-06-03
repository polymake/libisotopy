CXX = g++
MAX_DELTA ?= 20
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -Iinclude -fPIC -DMAX_DELTA=$(MAX_DELTA)
BIN_DIR = bin
TEST_DIR = $(BIN_DIR)/tests
BENCHMARK_DIR = $(BIN_DIR)/benchmarks

# Platform detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    TARGET = libisotopy.dylib
    SDK_PATH = /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
    CXXFLAGS += -isysroot $(SDK_PATH) -isystem $(SDK_PATH)/usr/include/c++/v1
    RPATH_FLAG = -Wl,-rpath,@loader_path/../..
else
    TARGET = libisotopy.so
    RPATH_FLAG = -Wl,-rpath,\$$ORIGIN/../..
endif

SRC = src/isotopy_graph.cpp
OBJ = obj/isotopy_graph.o

TEST_SRC = tests/isotopy_graph.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_BIN = $(TEST_DIR)/isotopy_graph_test

all: $(TARGET) $(BIN_DIR) $(TEST_DIR) $(BENCHMARK_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(TEST_DIR): | $(BIN_DIR)
	mkdir -p $(TEST_DIR)

$(BENCHMARK_DIR): | $(BIN_DIR)
	mkdir -p $(BENCHMARK_DIR)

$(TARGET): $(OBJ)
	$(CXX) -shared -o $@ $^

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TARGET) $(TEST_OBJ) | $(TEST_DIR)
	$(CXX) -o $@ $(TEST_OBJ) -L. $(RPATH_FLAG) -lisotopy


TEST_BATCH_SRC = tests/isotopy_graph_batch.cpp
TEST_BATCH_OBJ = obj/isotopy_graph_batch.o
TEST_FULL_BIN = $(BIN_DIR)/isotopy_graph_full
TEST_YAML = tests/isotopy_tests.yaml
TEST_YAML_XZ = tests/isotopy_tests.yaml.xz

$(TEST_YAML): $(TEST_YAML_XZ)
	xz -dc $< > $@

test_full: $(TEST_FULL_BIN) $(TEST_YAML)
	./$(TEST_FULL_BIN)

$(TEST_FULL_BIN): $(TARGET) $(TEST_OBJ) $(TEST_BATCH_OBJ) | $(BIN_DIR)
	$(CXX) -o $@ $(TEST_OBJ) $(TEST_BATCH_OBJ) -L. -Wl,-rpath,\$$ORIGIN/.. -lisotopy

BENCHMARK_SRC = tests/isotopy_graph_benchmark.cpp
BENCHMARK_OBJ = $(BENCHMARK_SRC:.cpp=.o)
BENCHMARK_BIN = $(BENCHMARK_DIR)/isotopy_graph_benchmark

benchmark: $(BENCHMARK_BIN) $(TEST_YAML)
	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): $(TARGET) $(BENCHMARK_OBJ) | $(BENCHMARK_DIR)
	$(CXX) -o $@ $(BENCHMARK_OBJ) -L. -Wl,-rpath,\$$ORIGIN/../.. -lisotopy

PROFILING_SRC = tests/isotopy_graph_profiling.cpp
PROFILING_OBJ = $(PROFILING_SRC:.cpp=.o)
PROFILING_BIN = $(BIN_DIR)/isotopy_graph_profiling

profiling: $(PROFILING_BIN) $(TEST_YAML)
	./$(PROFILING_BIN)

$(PROFILING_BIN): $(TARGET) $(PROFILING_OBJ) | $(BIN_DIR)
	$(CXX) -o $@ $(PROFILING_OBJ) -L. -Wl,-rpath,\$$ORIGIN/.. -lisotopy

docs:
	doxygen Doxyfile

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) \
	$(TEST_FULL_BIN) $(BENCHMARK_OBJ) $(BENCHMARK_BIN) $(PROFILING_OBJ) \
	$(PROFILING_BIN) obj/*.o wasm_obj/*.o libisotopy_wasm.a $(TEST_YAML)
	rm -rf $(BIN_DIR)

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)

# Emscripten build
.PHONY: emscripten docs

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
