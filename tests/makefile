#Tell make to make one .out file for each .cpp file found in the current directory
CC = g++
#CFLAGS = -O3 -Wno-narrowing -std=c++11 -march=native -I..
CFLAGS = -O0 -g -Wno-narrowing -std=c++11 -march=native -I..
all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm *.out 
