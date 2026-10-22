# Build each mode separately so switching optimization cannot reuse stale objects.
# Override Eigen's include root with: make EIGEN_INC=/path/to/eigen
BUILD_MODE ?= debug
EIGEN_INC ?= /usr/include/eigen3

ifeq ($(BUILD_MODE),debug)
MODE_FLAGS := -g -O0
SUFFIX := debug
else ifeq ($(BUILD_MODE),releaseO2)
MODE_FLAGS := -g -O2
SUFFIX := O2
else ifeq ($(BUILD_MODE),releaseO3)
MODE_FLAGS := -g -O3
SUFFIX := O3
else
$(error Unknown BUILD_MODE '$(BUILD_MODE)'; use debug, releaseO2, or releaseO3)
endif

OBJ_DIR := build/$(BUILD_MODE)
TARGET := bin/MolecularDynamics_$(SUFFIX)
SRCS := $(shell find src -name '*.cpp' | sort)
OBJS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
CPPFLAGS += -Isrc/include -I$(EIGEN_INC)
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic $(MODE_FLAGS)
# GCC 7's C++17 filesystem implementation lives in this separate library.
LDLIBS += -lstdc++fs

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS) | bin
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJ_DIR)/%.o: src/%.cpp Makefile
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c -o $@ $<

bin:
	mkdir -p $@

clean:
	rm -rf build
	rm -f bin/MolecularDynamics_debug bin/MolecularDynamics_O2 bin/MolecularDynamics_O3

-include $(OBJS:.o=.d)
