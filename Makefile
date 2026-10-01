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
GPU_HEADERS = $(GPU_HEADER) $(GPU_DIR)/isotopy_viro.h $(GPU_DIR)/isotopy_types.h \
              $(GPU_DIR)/isotopy_batch.h tests/test_env.h
GPU_SUITES = isotopy_gpu isotopy_viro
GPU_TEST_OBJS = $(addprefix obj/,$(addsuffix _test.o,$(GPU_SUITES)))
GPU_TSAN_OBJS = $(addprefix obj/,$(addsuffix _tsan.o,$(GPU_SUITES)))
GPU_TEST_BINS = $(addprefix $(TEST_DIR)/,$(addsuffix _test,$(GPU_SUITES)))
GPU_TSAN_BINS = $(addprefix $(TEST_DIR)/,$(addsuffix _tsan,$(GPU_SUITES)))
NVCC ?= nvcc
CUDA_ARCH ?= sm_75
# MAX_DELTA sizes libisotopy's arrays; GPU_MAX_DELTA sizes the kernel's shared memory.
# They are independent: the kernel caps at 9 because that is the largest degree whose
# tree code fits one uint64 and whose region count fits a uint32 bitset.
GPU_MAX_DELTA ?= 9
GPU_DEFS = -DMAX_DELTA=$(MAX_DELTA) -DISOTOPY_GPU_MAX_DELTA=$(GPU_MAX_DELTA)
GPU_INC = -Iinclude -I$(GPU_DIR) -Itests
# nvcc rejects host compilers newer than its ceiling (CUDA 11.1 caps at GCC 10),
# so the cluster build points it at an older g++ than the one used for the library.
# NVCC_EXTRA carries host-compiler escape hatches, e.g. -allow-unsupported-compiler
# when the only g++ on the node is newer than the toolkit's supported ceiling.
NVCC_CCBIN ?=
NVCC_EXTRA ?=
GPU_DEVICE_TABLE ?= 0
ifeq ($(GPU_DEVICE_TABLE),1)
GPU_DEFS += -DISOTOPY_GPU_DEVICE_TABLE
endif
GPU_FKYAML ?= 0
ifeq ($(GPU_FKYAML),1)
GPU_DEFS += -DISOTOPY_HAVE_FKYAML
endif
ifeq ($(origin GPU_STOP_AFTER),command line)
GPU_CUDA_BIN = $(BIN_DIR)/isotopy_gpu_stage$(GPU_STOP_AFTER)
GPU_DEFS += -DISOTOPY_GPU_STOP_AFTER=$(GPU_STOP_AFTER)
else
GPU_CUDA_BIN = $(BIN_DIR)/isotopy_gpu_cuda
endif
NVCCFLAGS = -arch=$(CUDA_ARCH) -std=c++17 -O3 -lineinfo $(GPU_INC) $(GPU_DEFS) \
            $(if $(NVCC_CCBIN),-ccbin $(NVCC_CCBIN),) $(NVCC_EXTRA)

gpu_test: $(GPU_TEST_BINS) $(TEST_YAML)
	@for t in $(GPU_TEST_BINS); do echo "== $$t"; ./$$t || exit 1; done

# The fast tags only: everything except [corpus], which needs the 42 MB
# decompressed YAML, and [threaded], which gpu_tsan owns.
gpu_test_fast: $(GPU_TEST_BINS)
	@echo "== $(TEST_DIR)/isotopy_gpu_test"
	@./$(TEST_DIR)/isotopy_gpu_test \
	  "[tables],[random],[reentrant],[fallback],[relabel],[limits]" | tail -3
	@echo "== $(TEST_DIR)/isotopy_viro_test"
	@./$(TEST_DIR)/isotopy_viro_test \
	  "[grammar],[code],[lookup],[batch],[fallback],[unknown]" | tail -3

# GPU_STOP_AFTER only ever reached nvcc on the cluster, so no host binary has
# ever contained a truncated pipeline -- a stage boundary that reads scratch the
# earlier stages never wrote would only surface as a bad per-stage timing run.
# The digests are lane-count dependent, so [stages] asserts termination,
# reproducibility and the flag set, not a value.
GPU_STAGES ?= 1 2 3 4 5
GPU_STAGE_BINS = $(addprefix $(TEST_DIR)/isotopy_gpu_stage,$(GPU_STAGES))

gpu_stage_test: $(GPU_STAGE_BINS)
	@for t in $(GPU_STAGE_BINS); do echo "== $$t"; ./$$t "[stages]" || exit 1; done

$(TEST_DIR)/isotopy_gpu_stage%: tests/isotopy_gpu.cpp $(GPU_HEADERS) $(TARGET) | $(TEST_DIR)
	$(CXX) -std=c++20 -Wall -Wextra -O1 -pthread -fPIC \
	  -DISOTOPY_HAVE_FKYAML -DISOTOPY_GPU_STOP_AFTER=$* \
	  $(GPU_DEFS) $(GPU_INC) $< -o $@ -L. $(RPATH_FLAG) -lisotopy

gpu_tsan: $(GPU_TSAN_BINS)
	@for t in $(GPU_TSAN_BINS); do echo "== $$t"; ./$$t "[threaded]" || exit 1; done

$(GPU_TEST_OBJS): obj/%_test.o: tests/%.cpp $(GPU_HEADERS) | obj/
	$(CXX) $(CXXFLAGS) -pthread -DISOTOPY_HAVE_FKYAML \
	  -DISOTOPY_GPU_MAX_DELTA=$(GPU_MAX_DELTA) $(GPU_INC) -c $< -o $@

$(GPU_TSAN_OBJS): obj/%_tsan.o: tests/%.cpp $(GPU_HEADERS) | obj/
	$(CXX) -std=c++20 -Wall -Wextra -O1 -g -fsanitize=thread -pthread -fPIC \
	  -DISOTOPY_GPU_HOST_WARP -DISOTOPY_HAVE_FKYAML $(GPU_DEFS) $(GPU_INC) -c $< -o $@

$(GPU_TEST_BINS): $(TEST_DIR)/%: obj/%.o $(TARGET) | $(TEST_DIR)
	$(CXX) -pthread -o $@ $< -L. $(RPATH_FLAG) -lisotopy

$(GPU_TSAN_BINS): $(TEST_DIR)/%: obj/%.o $(TARGET) | $(TEST_DIR)
	$(CXX) -fsanitize=thread -pthread -g -o $@ $< -L. $(RPATH_FLAG) -lisotopy

gpu_ptxas: $(GPU_DIR)/isotopy_gpu.cu $(GPU_HEADERS) | obj/
	$(NVCC) $(NVCCFLAGS) -Xptxas -v -c $< -o obj/isotopy_gpu_device.o

gpu_cuda: $(GPU_CUDA_BIN)

$(GPU_CUDA_BIN): $(GPU_DIR)/isotopy_gpu.cu $(GPU_HEADERS) $(TARGET) | $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) -DISOTOPY_GPU_MAIN $(GPU_DIR)/isotopy_gpu.cu -o $@ \
	  -L. -Xlinker -rpath -Xlinker '$$ORIGIN/..' -lisotopy

docs:
	doxygen Doxyfile

clean:
	rm -f $(OBJ) $(TARGET) $(TEST_OBJ) $(TEST_BIN) $(TEST_BATCH_OBJ) \
	$(TEST_FULL_BIN) $(BENCHMARK_OBJ) $(BENCHMARK_BIN) $(PROFILING_OBJ) \
	$(PROFILING_BIN) $(GPU_TEST_OBJS) $(GPU_TEST_BINS) $(GPU_TSAN_OBJS) $(GPU_TSAN_BINS) \
	$(BIN_DIR)/isotopy_gpu_cuda $(BIN_DIR)/isotopy_gpu_stage* obj/isotopy_gpu_device.o \
	obj/*.o wasm_obj/*.o libisotopy_wasm.a $(TEST_YAML)
	rm -rf $(BIN_DIR)

 
debug_test: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE
debug_test: LDFLAGS = -g
debug_test: $(TEST_BIN)
	./$(TEST_BIN)

# Emscripten build
.PHONY: emscripten docs gpu_test gpu_test_fast gpu_stage_test gpu_tsan gpu_ptxas gpu_cuda

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
