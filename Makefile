.PHONY: run debug

run: ./build/solver
	./build/solver

debug: ./build/solver
	gdb ./build/solver

./build/solver: ./src/solver.cpp ./src/base.cpp ./include/base.hpp
	cmake --build build


