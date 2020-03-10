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
mpi_residual (const double *a, const double *b, double *work_space,
              int n, int m, int p, int my_rank)
{
  int k = n / m;
  int l = n % m;
  int blocks = number_of_blocks_kl (k, l);
  int p_blocks = p_blocks_blocksp (blocks, p);
  int send_size = p_blocks * m * m;
  int recv_size = p * send_size;

  double *sum_line = work_space;
  double *send_line = sum_line + n;
  double *recv_line = send_line + recv_size;
  double *block_a = recv_line + recv_size;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double *p_buf = block_c + m * m;

  auto get_from_column = [=] (const double *line, double *block, int w, int h, int i)
    {
      int rank = i % p;
      int local_i = i / p;
      int offset = rank * send_size + local_i * w * m;
      memcpy (block, line + offset, w * h * sizeof (double));
    };
  auto add_to_sum = [=] (double *line, const double *block, int i, int h, int w)
    {
      for (int loc_i = 0; loc_i < h; ++loc_i)
        {
          for (int loc_j = 0; loc_j < w; ++loc_j)
            {
              line[i * m + loc_i] += fabs (block[loc_i * w + loc_j]);
            }
        }
    };

  nullify (sum_line, n);
  for (int j = 0; j < blocks; ++j)
    {
      // mpi_gather column j to recv_line
      get_column (b, send_line, n, m, p, my_rank, j);
      MPI_Allgather (send_line, send_size, MPI_DOUBLE,
                     recv_line, send_size, MPI_DOUBLE,
                     MPI_COMM_WORLD);
      // compose good line in send_line
      int w = block_width_jklm (j, k, l, m);
      for (int i = 0; i < blocks; ++i)
        {
          get_from_column (recv_line, send_line + i * m * w,
                           w, block_height_iklm (i, k, l, m), i);
        }
      // scalar product
      for (int i = my_rank; i < blocks; i += p)
        {
          int loc_i = i / p;
          int h = block_height_iklm (i, k, l, m);
          scalar_product (a + loc_i * n * m, send_line, block_c, n, m, h, w);
          if (i == j)
            {
              matrix_minus_E (block_c, w);
            }
          add_to_sum (sum_line, block_c, loc_i, h, w);
        }
    }
  double max = 0;
  for (int i = 0; i < blocks; i += p)
    {
      int loc_i = i / p;
      if (sum_line[loc_i] > max) max = sum_line[loc_i];
    }

  MPI_Allgather (&max, 1, MPI_DOUBLE, p_buf, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  double residual = p_buf[0];
  for (int i = 1; i < p; ++i)
    {
      if (p_buf[i] > residual)
        residual = p_buf[i];
    }

  return residual;
}


