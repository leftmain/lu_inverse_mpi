#ifndef MATRIX_MPI_H
#define MATRIX_MPI_H

#include "header.h"
#include "matrix_io.h"
#define A(a) MPI_Barrier(MPI_COMM_WORLD);if(my_rank==a){
#define B }MPI_Barrier(MPI_COMM_WORLD);

double
mpi_norma (const double *a, double *block, int n, int m,
           int p, int my_rank);

double
mpi_residual (const double *a, const double *b, double *work_space,
              int n, int m, int p, int my_rank);


#endif
