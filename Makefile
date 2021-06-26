
# -Wall : warn all | It turns on (almost) all the warnings that g++ can tell you about.
# 		  ref : https://stackoverflow.com/questions/2408038/what-does-wall-in-g-wall-test-cpp-o-test-do

# openmp flag of your compiler
OPENFLAG=-fopenmp
#OPENFLAG=-openmp
#OPENFLAG=-qopenmp

# -g : for debugging
CC = g++ -Wall -g -fopenmp

poissonSolver: Algorithms.o Vector.o PoissonMatrix.o PoissonSolver.o
	$(CC) -o $@ $+

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

Vector.o: Vector.cpp classes.h
	$(CC) -c -o $@ $<

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<

# ref : https://mropengate.blogspot.com/2018/01/makefile.html
.PHONY: clean
clean:
	rm *.o
