CXX = g++
MAX_DELTA ?= 20
# ARCH is overridable so a -march=native binary (safe on the cluster's fixed
# 2080Ti nodes, where we build on the run host) can be dropped for cross-host
# builds that would otherwise SIGILL.  Build with `make ARCH=` to disable.
ARCH ?= -march=native
CXXFLAGS = -std=c++20 -Wall -Wextra -O3 -flto=auto $(ARCH) -Iinclude -fPIC -DMAX_DELTA=$(MAX_DELTA)
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
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

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

GPU_DIR = gpu
GPU_HEADER = $(GPU_DIR)/isotopy_gpu.cuh
GPU_TEST_SRC = tests/isotopy_gpu.cpp
GPU_TEST_OBJ = obj/isotopy_gpu.o
GPU_TEST_BIN = $(TEST_DIR)/isotopy_gpu_test
GPU_TSAN_OBJ = obj/isotopy_gpu_tsan.o
GPU_TSAN_BIN = $(TEST_DIR)/isotopy_gpu_tsan
GPU_CUDA_BIN = $(BIN_DIR)/isotopy_gpu_cuda
NVCC ?= nvcc
CUDA_ARCH ?= sm_75
# MAX_DELTA sizes libisotopy's arrays; GPU_MAX_DELTA sizes the kernel's shared memory.
# They are independent: the kernel caps at 9 because that is the largest degree whose
# tree code fits one uint64 and whose region count fits a uint32 bitset.
GPU_MAX_DELTA ?= 9
GPU_DEFS = -DMAX_DELTA=$(MAX_DELTA) -DISOTOPY_GPU_MAX_DELTA=$(GPU_MAX_DELTA)
# nvcc rejects host compilers newer than its ceiling (CUDA 11.1 caps at GCC 10),
# so the cluster build points it at an older g++ than the one used for the library.
# NVCC_EXTRA carries host-compiler escape hatches, e.g. -allow-unsupported-compiler
# when the only g++ on the node is newer than the toolkit's supported ceiling.
NVCC_CCBIN ?=
NVCC_EXTRA ?=
NVCCFLAGS = -arch=$(CUDA_ARCH) -std=c++17 -O3 -lineinfo -Iinclude -I$(GPU_DIR) $(GPU_DEFS) \
            $(if $(NVCC_CCBIN),-ccbin $(NVCC_CCBIN),) $(NVCC_EXTRA)

gpu_test: $(GPU_TEST_BIN) $(TEST_YAML)
	./$(GPU_TEST_BIN)

$(GPU_TEST_BIN): $(TARGET) $(GPU_TEST_OBJ) | $(TEST_DIR)
	$(CXX) -pthread -o $@ $(GPU_TEST_OBJ) -L. $(RPATH_FLAG) -lisotopy

$(GPU_TEST_OBJ): $(GPU_TEST_SRC) $(GPU_HEADER) | obj/
	$(CXX) $(CXXFLAGS) -pthread -DISOTOPY_GPU_MAX_DELTA=$(GPU_MAX_DELTA) -I$(GPU_DIR) -c $< -o $@

gpu_test_tsan: $(GPU_TSAN_BIN)
	./$(GPU_TSAN_BIN) "[gpu][threaded]"

$(GPU_TSAN_BIN): $(TARGET) $(GPU_TSAN_OBJ) | $(TEST_DIR)
	$(CXX) -fsanitize=thread -pthread -g -o $@ $(GPU_TSAN_OBJ) -L. $(RPATH_FLAG) -lisotopy

$(GPU_TSAN_OBJ): $(GPU_TEST_SRC) $(GPU_HEADER) | obj/
	$(CXX) -std=c++20 -Wall -Wextra -O1 -g -fsanitize=thread -pthread -fPIC \
	  -DISOTOPY_GPU_HOST_WARP $(GPU_DEFS) -Iinclude -I$(GPU_DIR) -c $< -o $@

gpu_ptxas: $(GPU_DIR)/isotopy_gpu.cu $(GPU_HEADER) | obj/
	$(NVCC) $(NVCCFLAGS) -Xptxas -v -c $< -o obj/isotopy_gpu_device.o

GPU_STOP_AFTER ?= 6
gpu_cuda: $(GPU_CUDA_BIN)

# Stage-truncation builds for portable per-stage timing. Each produces a binary that
# returns after stage k; differencing their throughputs gives the per-stage cost.
gpu_stages: $(GPU_DIR)/isotopy_gpu.cu $(GPU_HEADER) $(TARGET) | $(BIN_DIR)
	@for k in 1 2 3 4 5 6; do \
	  echo "building stage-$$k binary"; \
	  $(NVCC) $(NVCCFLAGS) -DISOTOPY_GPU_MAIN -DISOTOPY_GPU_STOP_AFTER=$$k \
	    $(GPU_DIR)/isotopy_gpu.cu -o $(BIN_DIR)/isotopy_gpu_stage$$k \
	    -L. -Xlinker -rpath -Xlinker '$$$$ORIGIN/..' -lisotopy || exit 1; \
	done

gpu_bench: $(GPU_CUDA_BIN)
	./$(GPU_CUDA_BIN) bench 8 65536 1024 7

$(GPU_CUDA_BIN): $(GPU_DIR)/isotopy_gpu.cu $(GPU_HEADER) $(TARGET) | $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) -DISOTOPY_GPU_MAIN $(GPU_DIR)/isotopy_gpu.cu -o $@ \
	  -L. -Xlinker -rpath -Xlinker '$$ORIGIN/..' -lisotopy

docs:
	doxygen Doxyfile

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) \
	$(TEST_FULL_BIN) $(BENCHMARK_OBJ) $(BENCHMARK_BIN) $(PROFILING_OBJ) \
	$(PROFILING_BIN) $(GPU_TEST_OBJ) $(GPU_TEST_BIN) $(GPU_TSAN_OBJ) $(GPU_TSAN_BIN) \
	$(GPU_CUDA_BIN) $(BIN_DIR)/isotopy_gpu_stage* obj/isotopy_gpu_device.o obj/*.o wasm_obj/*.o libisotopy_wasm.a $(TEST_YAML)
	rm -rf $(BIN_DIR)

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)

# Emscripten build
.PHONY: emscripten docs gpu_test gpu_test_tsan gpu_ptxas gpu_cuda gpu_stages gpu_bench

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
