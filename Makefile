.PHONY: default run plot build debug

default: build
	./build/solver
	python ./tools/plot.py

run: build
	./build/solver

plot:
	python ./tools/plot.py

build: ./src/solver.cpp ./src/base.cpp ./src/dataio.cpp ./include/base.hpp ./include/dataio.hpp
	cmake --build build

debug: build
	gdb ./build/solver
