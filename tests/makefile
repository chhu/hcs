#Tell make to make one .out file for each .cpp file found in the current directory
CC = g++
#CFLAGS = -O2 -std=c++11 -march=native -I..
CFLAGS = -O3 -std=c++11 -fpermissive -march=native -I..
all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp
	$(CC) $(CFLAGS) $< -o $@ -lm
clean:
	rm *.out 
