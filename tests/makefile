#Tell make to make one .out file for each .cpp file found in the current directory
CC = clang++
CFLAGS = -O3 -march=native -Wno-narrowing -std=c++11 -I..
#CFLAGS = -O0 -g -Wno-narrowing -std=c++11 -march=native -I..

.PHONY: clean


all: $(patsubst %.cpp, %, $(wildcard *.cpp))
%: %.cpp includes.hpp solver.hpp ../sparsefield.hpp ../hcs.hpp
	$(CC) $(CFLAGS) $< -o $@
clean: 
	rm -f $(patsubst %.cpp, %, $(wildcard *.cpp))
	rm -f *.pgm
