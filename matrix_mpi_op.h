#ifndef MATRIX_MPI_H
#define MATRIX_MPI_H

#include "header.h"
#include "matrix_io.h"
#include "matrix_op.h"

double
mpi_norma (const double *a, double *block, const int n, const int m,
           const int p, const int my_rank);

double
mpi_residual (const double *a, const double *b, double *workspace, const int n,
              const int m, const int p, const int my_rank, Time &time);

int
mpi_lu_matrix_inverse (double *a, double *b, double *workspace,
                       const int n, const int m, const int p,
                       const int my_rank, Time &time);

#endif
