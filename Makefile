CFLAGS = -O3 -ffast-math -Wall
#CFLAGS = -g 
FLAGS = CFLAGS -std=c++14
G = mpicxx
C = mpicc
M = Makefile header.h

all: a.out

a.out: main.o matrix_io.o
	$G $(LFLAGS) $^ -o $@

main.o: main.cpp $M
	$C $(LFLAGS) $(CFLAGS) -c $< -o $@

matrix_io.o: matrix_io.cpp matrix_io.h $M
	$C $(LFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm *.o a.out

