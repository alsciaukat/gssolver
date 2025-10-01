.PHONY: run debug plot

run: ./build/solver
	./build/solver

plot:
	python ./tools/plot.py

debug: ./build/solver
	gdb ./build/solver

./build/solver: ./src/solver.cpp ./src/base.cpp ./src/dataio.cpp ./include/base.hpp ./include/dataio.hpp
	cmake --build build


