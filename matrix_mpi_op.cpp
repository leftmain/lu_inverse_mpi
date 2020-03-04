#include "matrix_mpi_op.h"
#include "matrix_op.h"


double
mpi_norma (const double *a, double *block, int n, int m, int p, int my_rank)
{
  int k = n / m;
  int l = n % m;
  int blocks = number_of_blocks_kl (k, l);

  double max = 0.;
  double sum = 0.;
  for (int i = my_rank; i < blocks; i += p)
    {
      int h = block_height_iklm (i, k, l, m);
      sum = line_norma (a + (i / p) * n * m, block, n, m, h);

      if (sum > max)
        max = sum;
    }
  double norma = 0;
// BAD: fix true order
  MPI_Allreduce (&max, &norma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return norma;
}


double
mpi_residual (const double *a, double *b, double *lines,
              int n, int m, int p, int my_rank)
{
  double *line = lines;
  double *block_a = line + n * m;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double norm = 0.;

  int k = n / m;
  int l = n % m;
  int blocks = number_of_blocks_kl (k, l);
  int p_blocks = (blocks % p == 0) ? blocks / p : blocks / p + 1;

  MPI_Status status;
  int tag = 0;

  // can be done with p steps using 3n^2 memory
  for (int step = 0; step < blocks; ++step)
    {
      int i = step / p;
      int a_h = block_height_iklm (i, k, l, m);
      if (step % p == 0)
        nullify (line, n * m);

      int diff = (my_rank + step) % p;
      for (int j = diff; j < blocks; j += p)
        {
          int a_w = block_width_jklm (j, k, l, m);
          get_block (a, block_a, n, m, p, my_rank, i * p + my_rank, j);

          for (int t = 0; t < blocks; ++t)
            {
              int t_w = block_width_jklm (t, k, l, m);
              get_block (b, block_b, n, m, p, my_rank, j - diff, t);
              matrix_plus_multiply (block_a, block_b, line + t * a_w * m,
                                    a_h, a_w, t_w);
            }
        }

      if (step % p == p - 1)
        {
          // line - diag E; norm of line
          int diag_ind = i + my_rank;
          // # change a to line
          matrix_minus_E (line + a_h * m * diag_ind,
                          block_height_iklm (diag_ind, k, l, m));
          double max = line_norma (line, block_c, n, m, a_h);
          if (max > norm)
            norm = max;
        }

      MPI_Sendrecv_replace (b, p_blocks * n * m, MPI_DOUBLE,
                            (my_rank - 1 + p) % p, tag, (my_rank + 1) % p,
                            tag, MPI_COMM_WORLD, &status);
    }
  double real_norma = 0.;
// BAD: fix true order
  MPI_Allreduce (&norm, &real_norma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return real_norma;
}

