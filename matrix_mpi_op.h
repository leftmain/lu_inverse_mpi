#ifndef MATRIX_MPI_H
#define MATRIX_MPI_H

#include "header.h"
#include "matrix_io.h"

double
mpi_norma (const double *a, double *block, int n, int m,
           int p, int my_rank);

double
mpi_residual (const double *a, const double *b, double *line,
              int n, int m, int p, int my_rank);


#endif
