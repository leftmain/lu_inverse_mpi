#include "matrix_mpi_op.h"
#include "matrix_op.h"

//-----------------------------------------------------------------------------

double
mpi_norma (const double *a, double *block, int n, int m, int p, int my_rank)
{
  int k = n / m;
  int l = n % m;
  int blocks = (l == 0) ? k : k + 1;

  double max = 0.;
  double sum = 0.;
  for (int i = my_rank; i < blocks; i += p)
    {
      int h = (i == k) ? l : m;
      sum = line_norma (a + (i / p) * n * m, block, n, m, h);

      if (sum > max)
        max = sum;
    }
  double norma = 0;
  MPI_Allreduce (&max, &norma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return norma;
}


//-----------------------------------------------------------------------------

double
mpi_residual (const double *a, const double *b, double *lines,
              int n, int m, int p, int my_rank)
{
  for (int step = 0; step < p; ++step)
    {
    }
  return 0;
}


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

