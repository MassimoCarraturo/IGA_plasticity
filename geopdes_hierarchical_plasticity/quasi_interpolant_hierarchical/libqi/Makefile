# Makefile for the THB-spline Quasi-Interpolant library (libqi)
#
# Build targets:
#   make            -> build static lib, shared lib, and test binary
#   make test       -> run the test suite
#   make python     -> rebuild the shared library used by the Python wrapper
#   make mex        -> build the MATLAB MEX wrapper (requires MATLAB)
#   make clean      -> remove build artefacts
#
# Detect platform-specific shared-library extension.
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
    SOEXT  := dylib
    SOFLAG := -dynamiclib
else
    SOEXT  := so
    SOFLAG := -shared
endif

CC      ?= gcc
CFLAGS  ?= -O3 -fPIC -std=c99 -Wall -Wextra -Wno-unused-parameter -Wno-unused-but-set-variable -fopenmp
LDFLAGS ?= -fopenmp -lm

INCLUDE := -Iinclude -Isrc

SRC_DIR := src
BUILD   := build
TEST    := test

LIB_SRCS := $(SRC_DIR)/qi.c $(SRC_DIR)/qi_linalg.c $(SRC_DIR)/qi_bspline.c $(SRC_DIR)/qi_spatial.c
LIB_OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD)/%.o,$(LIB_SRCS))

STATIC_LIB := $(BUILD)/libqi.a
SHARED_LIB := $(BUILD)/libqi.$(SOEXT)
TEST_BIN   := $(BUILD)/qi_test

.PHONY: all clean test python mex info

all: $(STATIC_LIB) $(SHARED_LIB) $(TEST_BIN)

# Object files
$(BUILD)/%.o: $(SRC_DIR)/%.c | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(BUILD):
	@mkdir -p $(BUILD)

# Static library
$(STATIC_LIB): $(LIB_OBJS)
	ar rcs $@ $^

# Shared library
$(SHARED_LIB): $(LIB_OBJS)
	$(CC) $(SOFLAG) -o $@ $^ $(LDFLAGS)

# Test binary (links against the static library)
$(TEST_BIN): $(TEST)/qi_test.c $(STATIC_LIB)
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@ $(STATIC_LIB) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

# Optional micro-benchmark (parallel scaling)
bench: $(BUILD)/qi_bench
	./$(BUILD)/qi_bench

$(BUILD)/qi_bench: $(TEST)/qi_bench.c $(STATIC_LIB)
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@ $(STATIC_LIB) $(LDFLAGS)

python: $(SHARED_LIB)
	@echo "Shared library at $(SHARED_LIB) — usable from python/qi.py"

# MATLAB MEX
# Override MATLAB_ROOT on the command line if mex is not on PATH:
#   make mex MATLAB_ROOT=/Applications/MATLAB_R2024a.app
MATLAB_ROOT ?=
ifeq ($(MATLAB_ROOT),)
    MEX := mex
else
    MEX := $(MATLAB_ROOT)/bin/mex
endif

mex: $(LIB_OBJS)
	$(MEX) -R2018a -Iinclude -Isrc \
	    CFLAGS='$$CFLAGS -fopenmp -O3 -std=c99' \
	    LDFLAGS='$$LDFLAGS -fopenmp' \
	    matlab/qi_mex.c $(LIB_OBJS) -outdir matlab

clean:
	rm -rf $(BUILD)
	rm -f matlab/*.mex* python/*.so python/*.dylib

info:
	@echo "Platform   : $(UNAME)"
	@echo "Compiler   : $(CC)"
	@echo "CFLAGS     : $(CFLAGS)"
	@echo "Shared ext : $(SOEXT)"
	@echo "Sources    : $(LIB_SRCS)"
