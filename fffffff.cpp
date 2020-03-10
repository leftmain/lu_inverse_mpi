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
  int recv_size = p * p_blocks * m * m;

  double *sum_line = work_space;
  double *send_line = sum_line + n;
  double *recv_line = send_line + p_blocks * m;
  double *block_a = recv_line + recv_size;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double *p_buf = block_c + m * n;

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
      A(0) printf ("0j = %d\n", j); B
      A(1) printf ("1j = %d\n", j); B
      A(2) printf ("2j = %d\n", j); B
      // mpi_gather column j to recv_line
      get_column (b, send_line, n, m, p, my_rank, j);
      A(0)
      printf ("0>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
      print_matrix (send_line, 1, send_size);
      B
      A(1)
      printf ("1>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
      print_matrix (send_line, 1, send_size);
      B
      A(0) printf ("0before\n"); B
      A(1) printf ("1before\n"); B
      A(2) printf ("2before\n"); B
      MPI_Allgather (send_line, send_size, MPI_DOUBLE,
                     recv_line, send_size, MPI_DOUBLE,
                     MPI_COMM_WORLD);
      A(0) printf ("0after\n"); B
      A(1) printf ("1after\n"); B
      A(2) printf ("2after\n"); B
      A(0) printf ("0>\n"); print_matrix (recv_line, 1, recv_size); B
      A(1) printf ("1>\n"); print_matrix (recv_line, 1, recv_size); B
      // compose good line in send_line
      int w = block_width_jklm (j, k, l, m);
      for (int i = 0; i < blocks; ++i)
        {
          get_from_column (recv_line, send_line + i * m * w,
                           w, block_height_iklm (i, k, l, m), i);
        }
      A(0) printf ("#0\n"); print_matrix (send_line, 1, w * n); B
      A(1) printf ("#1\n"); print_matrix (send_line, 1, w * n); B
      // scalar product
      for (int i = my_rank; i < blocks; i += p)
        {
          A(0) printf ("0i = %d\n", i); B
          A(1) printf ("1i = %d\n", i); B
          A(2) printf ("2i = %d\n", i); B
          int loc_i = i / p;
          int h = block_height_iklm (i, k, l, m);
          A(0) printf ("0 i k blocks %d %d %d\n", i, k, blocks); B
          A(1) printf ("1 i k blocks %d %d %d\n", i, k, blocks); B
          A(0) printf ("-0 hl %d,%d\n",h,l); print_matrix (a+loc_i*n*m, h, n); B
          A(1) printf ("-1 hl %d,%d\n",h,l); print_matrix (a+loc_i*n*m, h, n); B
          scalar_product (a + loc_i * n * m, send_line, block_c, n, m, h, w);
          A(0) printf ("##0\n"); print_matrix (block_c, h, w); B
          A(1) printf ("##1\n"); print_matrix (block_c, h, w); B
          if (i == j)
            {
              matrix_minus_E (block_c, w);
            }
          A(0) printf ("###0\n"); print_matrix (block_c, h, w); B
          A(1) printf ("###1\n"); print_matrix (block_c, h, w); B
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


double
mpi_residual_1 (const double *a, double *b, double *lines,
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

