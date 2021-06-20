
# -Wall : warn all | It turns on (almost) all the warnings that g++ can tell you about.
# 		  ref : https://stackoverflow.com/questions/2408038/what-does-wall-in-g-wall-test-cpp-o-test-do

# -g : for debugging
CC = g++ -Wall -g

poissonSolver: Algorithms.o Vector.o PoissonSolver.o
	$(CC) -o $@ $+

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

Vector.o: Vector.cpp classes.h
	$(CC) -c -o $@ $<

# ref : https://mropengate.blogspot.com/2018/01/makefile.html
.PHONY: clean
clean:
	rm *.o