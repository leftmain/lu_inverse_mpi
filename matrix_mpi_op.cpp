#include "matrix_mpi_op.h"
#include "matrix_op.h"


double
mpi_norma (const double *a, double *block, const int n, const int m,
           const int p, const int my_rank)
{
  int k = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k, l);

  double max = 0.;
  double sum = 0.;
  for (int i = my_rank; i < blocks; i += p)
    {
      int h = get_block_size (i, k, l, m);
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
mpi_residual (const double *a, const double *b, double *workspace, const int n,
              const int m, const int p, const int my_rank, Time &time)
{
  GlobalTimer global_timer (time.residual_global_time);
  ProcessTimer process_timer (time.residual_process_time);

  int k = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k, l);
  int p_blocks = get_p_blocks (blocks, p);
  int send_size = p_blocks * m * m;
  int recv_size = p * send_size;

  double *sum_line = workspace;
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
      copy_array (line + offset, block, w * h);
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
      int w = get_block_size (j, k, l, m);
      for (int i = 0; i < blocks; ++i)
        {
          get_from_column (recv_line, send_line + i * m * w,
                           w, get_block_size (i, k, l, m), i);
        }

      // scalar product
      for (int i = my_rank; i < blocks; i += p)
        {
          int loc_i = i / p;
          int h = get_block_size (i, k, l, m);
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


// todo optimize: send not all line
// (while send: line -> line + h * m * k, change size)
// KIJ algorighm
static int
mpi_lu_decomposition (double *a, double *workspace,
                      const int n, const int m, const int p, const int my_rank)
{
  double *block_a = workspace;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double *line = block_c + m * m;

  int k_blocks = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k_blocks, l);

  for (int k = 0; k < blocks - 1; ++k)
    {
      int line_size = m * n; // make m*n-k*m^2
      int root = k % p;

      // set line k from a to line (in process k % p)
      if (my_rank == root)
        {
          get_line (a, line, n, m, p, my_rank, k);

          //  block_a = A_kk (from line)
          get_block_from_line (line, block_a, n, m, m, k); // (h = m for all k)

          //  block_b = (A_kk)^(-1)
          if (lu_matrix_inverse (block_a, block_b, m))
            {
              line[line_size] = -1;
            }
          else
            {
              line[line_size] = 1;
            }

          set_block_to_line (line, block_b, n, m, m, k);
        }

      // send line k to others
      MPI_Bcast (line, line_size + 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      if (line[line_size] < 0)
        {
          return METHOD_ERROR;
        }

      // apply line
      int from = k - root + my_rank;
      if (from <= k)
        from += p;
      for (int i = from; i < blocks; i += p)
        {
          int h = get_block_size (i, k_blocks, l, m);

          //  A_ik = A_ik * (A_kk)^(-1)
          //    block_a = A_ik
          get_block (a, block_a, n, m, p, my_rank, i, k);

          //    block_b = (A_kk)^(-1) (from line)
          get_block_from_line (line, block_b, n, m, m, k);

          //    block_c = A_ik * (A_kk)^(-1) = block_a * block_b
          matrix_multiply (block_a, block_b, block_c, h, m, m);

          //    A_ik = block_c
          set_block (a, block_c, n, m, p, my_rank, i, k);

          //  A_ij -= A_ik * A_kj (block_a = block_c * block_b)
          //  A_kj - from line
          for (int j = k + 1; j < blocks; ++j)
            {
              int w = get_block_size (j, k_blocks, l, m);

              //  block_a = A_ij
              get_block (a, block_a, n, m, p, my_rank, i, j);

              //  block_b = A_kj (from line)
              get_block_from_line (line, block_b, n, m, m, j);

              // block_a -= block_c * block_b
              matrix_minus_multiply (block_c, block_b, block_a, h, m, w);

              // A_ij = block_a
              set_block (a, block_a, n, m, p, my_rank, i, j);
            }
        }
    }
  return ALL_RIGHT;
}


static int
mpi_lu_invert (double *a, const double *b, double *workspace,
               const int n, const int m, const int p, const int my_rank)
{
  double *block_a = workspace;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double *line = block_c + m * m;

  int k_blocks = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k_blocks, l);
  int p_blocks = get_p_blocks (blocks, p);

  double *column = line + p * p_blocks * m * m;

  //  U = U^(-1)
  //    U_ii = (U_ii)^(-1)
  int error = ALL_RIGHT;
  for (int i = my_rank; i < blocks; i += p)
    {
      int h = get_block_size (i, k_blocks, l, m);

      // block_a = U_ii
      get_block (a, block_a, n, m, p, my_rank, i, i);

      // block_b = (U_ii)^(-1) = 1 / block_a
      if (lu_matrix_inverse (block_a, block_b, h))
        {
          error = 1;
          break;
        }

      // U'_ii = block_b
      set_block (a, block_b, n, m, p, my_rank, i, i);
    }
  int sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (sum != ALL_RIGHT)
    return METHOD_ERROR;

  auto get_from_column = [=] (const double *line, double *block, int send_size,
                              int w, int h, int i)
    {
      int rank = i % p;
      int local_i = i / p;
      int offset = rank * send_size + local_i * w * m;
      copy_array (line + offset, block, w * h);
    };

  //    U_ij = -sum
  for (int k = 1; k < blocks; ++k)
    {
      int w = get_block_size (k, k_blocks, l, m);

      //  get column k and share it with others
      get_column (a, column, n, m, p, my_rank, k, 0, k);

      //  share column with others
      int send_size = ((k + 1) / p) * m * w;
      if ((k + 1) % p > 0)
        send_size += w * m;
      MPI_Allgather (column, send_size, MPI_DOUBLE,
                     line, send_size, MPI_DOUBLE,
                     MPI_COMM_WORLD);

      //  put line in column in right order
      for (int i = 0; i <= k; ++i)
        {
          int h = get_block_size (i, k_blocks, l, m);
          get_from_column (line, column + i * m * w, send_size, w, h, i);
        }

      for (int i = my_rank; i < k; i += p)
        {
          int h = get_block_size (i, k_blocks, l, m);

          //  block_c = 0
          nullify (block_c, m * m);

          for (int j = i; j < k; ++j)
            {
              //  block_c -= U'_ij * U_jk

              //    block_a = U'_ij
              get_block (a, block_a, n, m, p, my_rank, i, j);

              //    block_b = U_jk (from column)
              get_block_from_line (column, block_b, n, m, w, j);

              //    block_c -= block_a * block_b
              matrix_minus_multiply (block_a, block_b, block_c, h, m, w);
//              if(my_rank==0){print_matrix(block_a,h,m);print_matrix(block_b,m,w);
//              printf("j=%d\n",j);}
            }

          //  block_a = U'_ij = block_c * U'_kk
          //    block_b = U'_kk (from column)
          get_block_from_line (column, block_b, n, m, w, k);
          
          matrix_multiply (block_c, block_b, block_a, h, w, w);

          //  U'_ij = block_a
          set_block (a, block_a, n, m, p, my_rank, i, k);
        }
    }

  //  L = L^(-1)
  //    L = -L
  int from = my_rank;
  if (from < 1)
    from += p;
  for (int i = from; i < blocks; i += p)
    {
      int h = get_block_size (i, k_blocks, l, m);
      minus_array (a + (i / p) * n * m, i * h * m);
    }

  //    L'_ij = -sum_{s=j}^{i-1} L_is * L'_sj
  for (int k = 1; k < blocks - 1; ++k)
    {
      int root = k % p;

      //  set line k from a to line (in process k % p)
      if (my_rank == root)
        {
          get_line (a, line, n, m, p, my_rank, k, 0, k - 1);
        }

      //  send line k to others
      MPI_Bcast (line, k * m * m, MPI_DOUBLE, root, MPI_COMM_WORLD);

      //  apply line to all rows of process
      from = k - root + my_rank;
      if (from <= k)
        from += p;
      for (int i = from; i < blocks; i += p)
        {
          int h = get_block_size (i, k_blocks, l, m);
          for (int j = 0; j < k; ++j)
            {
              //  block_c = L'_ij
              get_block (a, block_c, n, m, p, my_rank, i, j);

              //  L'_ij -= L_ik * L'_kj
              //    block_a = L_ik (from matrix b)
              get_block (b, block_a, n, m, p, my_rank, i, k);

              //    block_b = L'_kj (from line)
              get_block_from_line (line, block_b, n, m, m, j);

              //    block_c -= block_a * block_b
              matrix_minus_multiply (block_a, block_b, block_c, h, m, m);

              // L'_ij = block_c
              set_block (a, block_c, n, m, p, my_rank, i, j);
            }
        }
    }

  return ALL_RIGHT;
}


static void
mpi_ul_multiply (const double *a, double *b, double *workspace,
                 const int n, const int m, const int p, const int my_rank)
{
  double *block_a = workspace;
  double *block_b = block_a + m * m;
  double *block_c = block_b + m * m;
  double *line = block_c + m * m;

  int k_blocks = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k_blocks, l);
  int p_blocks = get_p_blocks (blocks, p);

  double *column = line + p * p_blocks * m * m;

  auto get_from_column = [=] (const double *line, double *block, int send_size,
                              int w, int h, int i, int k)
    {
      int rank = i % p;
      int local_i = (i - k) / p;
      int offset = rank * send_size + local_i * w * m;
      copy_array (line + offset, block, w * h);
    };


  for (int k = 0; k < blocks; ++k)
    {
      int w = get_block_size (k, k_blocks, l, m);
      //  get column k and share it with others
      get_column (a, column, n, m, p, my_rank, k, k + 1, blocks - 1);

      //  share column with others
      int send_size = ((blocks - k) / p) * w * m;
      if ((blocks - k) % p > 0)
        send_size += w * m;
      MPI_Allgather (column, send_size, MPI_DOUBLE,
                     line, send_size, MPI_DOUBLE,
                     MPI_COMM_WORLD);

      //  put line in column in right order
      for (int i = k + 1; i <= blocks; ++i)
        {
          get_from_column (line, column + i * m * w, send_size,
                           w, get_block_size (i, k_blocks, l, m), i, k + 1);
        }

      for (int i = my_rank; i < blocks; i += p)
        {
          int h = get_block_size (i, k_blocks, l, m);

          //  A_ik = ...
          int from = std::max (i, k);
          if (from == k)
            {
              //  block_a = A_ik
              get_block (a, block_a, n, m, p, my_rank, i, k);
            }
          else
            {
              //  block_a = A_ii * A_ik
              //    block_b = A_ii
              get_block (a, block_b, n, m, p, my_rank, i, i);

              //    block_c = A_ik
              get_block (a, block_c, n, m, p, my_rank, i, k);

              //    block_a = block_b * block_c
              matrix_multiply (block_b, block_c, block_a, h, h, w);
            }
          for (int j = from + 1; j < blocks; ++j)
            {
              int loc_w = get_block_size (j, k_blocks, l, m);
              //  block_a += A_ij * A_jk
              //    block_b = A_ij
              get_block (a, block_b, n, m, p, my_rank, i, j);

              //    block_c = A_jk (from column)
              get_block_from_line (column, block_c, n, m, w, j);

              //    block_a += block_b * block_c
              matrix_plus_multiply (block_b, block_c, block_a, h, loc_w, w);
            }
          //  A_ik = block_a
          set_block (b, block_a, n, m, p, my_rank, i, k);
        }
    }
}


int
mpi_lu_matrix_inverse (double *a, double *b, double *workspace,
                       const int n, const int m, const int p,
                       const int my_rank, Time &time)
{
  GlobalTimer global_timer (time.global_time);
  ProcessTimer process_timer (time.process_time);

  //  LU decomposition
  ProcessTimer decomposition_timer (time.decomposition_time);
  if (mpi_lu_decomposition (a, workspace, n, m, p, my_rank))
    return METHOD_ERROR;
  decomposition_timer.stop ();

#ifdef DEBUG
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0)
    printf ("decomp:\n");
  mpi_print_block_matrix (a, workspace, n, m, p, my_rank);
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  int k = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k, l);
  int p_blocks = get_p_blocks (blocks, p);
  copy_array (a, b, p_blocks * n * m);

  //  inversoin L and U
  ProcessTimer invert_timer (time.invert_time);
  if (mpi_lu_invert (a, b, workspace, n, m, p, my_rank))
    return METHOD_ERROR;
  invert_timer.stop ();

#ifdef DEBUG
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0)
    printf ("invert:\n");
  mpi_print_block_matrix (a, workspace, n, m, p, my_rank);
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  ProcessTimer mult_timer (time.mult_time);
  mpi_ul_multiply (a, b, workspace, n, m, p, my_rank);
  mult_timer.stop ();

  return ALL_RIGHT;
}


