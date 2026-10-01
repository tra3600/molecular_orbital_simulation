CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -fopenmp
SRC      := $(filter-out src/main.cpp,$(wildcard src/*.cpp))
OBJ      := $(SRC:src/%.cpp=build/%.o)

all: molorb

build/%.o: src/%.cpp $(wildcard src/*.h)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

molorb: $(OBJ) build/main.o
	$(CXX) $(CXXFLAGS) $^ -o $@

build/main.o: src/main.cpp $(wildcard src/*.h)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

tests: $(OBJ) tests/test_main.cpp
	$(CXX) $(CXXFLAGS) -Isrc $^ -o build/test_runner && ./build/test_runner

clean:
	rm -rf build molorb

.PHONY: all tests clean
