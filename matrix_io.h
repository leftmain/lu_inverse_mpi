#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include "header.h"
#define MAX_PRINT 15000 //20
#define SHORT_PRINT
//#define PRINT_LINE

//------------------------------------- common functions

void
print_block_line (double *a, int n, int m, int h, FILE *fp = stdout,
                  int max_h = MAX_PRINT);

void // print matrix n x m
print_matrix (double *a, int n, int m, FILE *fp = stdout);

int
print_matrix (double *a, int n, int m, const char *file_name);


//------------------------------------- MPI functions
// MPI function ensure that all processes have the same return value

int
mpi_read_matrix (double *a, double *line, int n, int m, int p,
                 int my_rank, const char *file_name);

void
mpi_init_matrix (double *a, int n, int m, int p, int my_rank,
                 double (*f)(int, int));

void // print block matrix n x n (m - block size)
mpi_print_block_matrix_simple (double *a, double *lines, int n, int m,
                               int p, int my_rank);

void // print block matrix n x n (m - block size)
mpi_print_block_matrix (double *a, double *line, int n, int m, int p,
                        int my_rank, FILE *fp = stdout);

int
mpi_print_block_matrix (double *a, double *line, int n, int m, int p,
                        int my_rank, const char *file_name);

#endif
