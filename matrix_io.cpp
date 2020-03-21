#include "matrix_io.h"



int
mpi_read_matrix (double *a, double *line, int n, int m, int p,
                 int my_rank, const char *file_name)
{
  Closer closer;
  FILE *fp = nullptr;
  int error = ALL_RIGHT;

  if (my_rank == 0)
    {
      fp = fopen (file_name, "r");
      closer.add (fp);
      if (!fp)
        error = CANNOT_OPEN;
    }

  MPI_Bcast (&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (error != ALL_RIGHT)
    {
      return error;
    }


  int k = n / m;
  int l = n % m;
  int blocks = (l == 0) ? k : k + 1;
  int tag = 0;
  MPI_Status status;

  for (int I = 0; I < blocks; ++I)
    {
      int h = (I == k) ? l : m;

      // read line
      error = ALL_RIGHT;
      if (my_rank == 0)
        {
          double *curr_line = (I % p == 0) ? a + (I / p) * n * m : line;
          for (int i = 0; i < h; ++i)
            {
              for (int j = 0; j < n; ++j)
                {
                  int w = (j >= k * m) ? l : m;
                  int ind = (j / m) * h * m + i * w + j % m;
                  if (fscanf (fp, "%lf", curr_line + ind) != 1)
                    {
                      error = CANNOT_READ;
                      i = h;
                      break;
                    }
                }
            }
        }
      MPI_Bcast (&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (error != ALL_RIGHT)
        {
          return error;
        }

      // send line
      if (I % p != 0)
        {
          if (my_rank == 0)
            {
              MPI_Send (line, h * n, MPI_DOUBLE, I % p, tag, MPI_COMM_WORLD);
            }
          else if (my_rank == I % p)
            {
              MPI_Recv (a + (I / p) * n * m, h * n, MPI_DOUBLE, 0, tag,
                        MPI_COMM_WORLD, &status);
            }
        }
    }

  return ALL_RIGHT;
}


void
mpi_init_matrix (double *a, int n, int m, int p, int my_rank,
                 double (*f)(int, int))
{
  int k = n / m;
  int l = n % m;
  int blocks = (l == 0) ? k : k + 1;

  for (int i = my_rank; i < blocks; i += p)
    {
      int c_m = (i == k) ? l : m;
      for (int j = 0; j < blocks; ++j)
        {
          int s_m = (j == k) ? l : m;
          int block_index = (i / p) * n * m + j * c_m * m;
          // current block has c_m x s_m size
          for (int loc_i = 0; loc_i < c_m; ++loc_i)
            {
              for (int loc_j = 0; loc_j < s_m; ++loc_j)
                {
                  a[block_index + loc_i * s_m + loc_j] = f (i * m + loc_i,
                                                            j * m + loc_j);
                }
            }
        }
    }
}


void
mpi_print_matrix_simple (const double *a, double *lines,
                               int n, int m, int p, int my_rank)
{
  int k = n / m;
  int l = n % m;
  int s_blocks = (l == 0) ? k : k + 1;
  int c_blocks = (s_blocks % p == 0) ? s_blocks / p : s_blocks / p + 1;
  int tag = 0;
  MPI_Status status;

  if (my_rank == 0)
    {
      for (int t = 0; t < p; ++t)
        {
          const double *curr_array = a;
          int columns = k / p;
          if (k % p > 0 && t < k % p)
            columns++;
          columns *= m;
          if (k % p > 0 && t - 1 < k % p && t >= k % p && l > 0)
            columns += l;
          if (t > 0)
            {
              curr_array = lines;
              MPI_Recv (lines, c_blocks * n * m, MPI_DOUBLE, t, tag,
                        MPI_COMM_WORLD, &status);
            }
          printf ("%d:\n", t);
          print_matrix (curr_array, 1, columns * n);
        }
    }
  else
    {
      MPI_Send (a, c_blocks * n * m, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
}


void
mpi_print_matrix_with_prefix (const double *a, double *line, int n, int m,
                              int p, int my_rank, FILE *fp, const char *prefix)
{
  if (my_rank == 0 && prefix != nullptr)
    {
      fprintf (fp, "%s\n", prefix);
    }
  mpi_print_matrix (a, line, n, m, p, my_rank, fp);
}


void
mpi_print_matrix_with_prefix (const double *a, double *line, int n, int m,
                              int p, int my_rank, const char *prefix)
{
  mpi_print_matrix_with_prefix (a, line, n, m, p, my_rank, stdout, prefix);
}


void
mpi_print_matrix (const double *a, double *line, int n, int m,
                        int p, int my_rank, FILE *fp)
{
  int N = n;
  int tag = 0;
  MPI_Status status;

  if (fp == stdout)
    {
      N = (n > MAX_PRINT) ? MAX_PRINT : n;
    }
  int k = n / m;
  int l = n % m;
  int K = N / m;
  int L = N % m;
  int blocks = (L == 0) ? K : K + 1;

  const double *curr_array = a;
  int i = 0;
  int local_line = -1;
  for (i = 0; i < blocks; ++i)
    {
      int h = (i == k) ? l : m;
      int H = (i == K) ? L : m;
      if (i % p == 0)
        {
          local_line++;
          curr_array = a + local_line * n * m;
        }
      else
        {
          if (my_rank == 0)
            {
              MPI_Recv (line, n * h, MPI_DOUBLE, i % p, tag,
                        MPI_COMM_WORLD, &status);
              curr_array = line;
            }
          else if (i % p == my_rank)
            {
              MPI_Send (a + local_line * n * m, n * h, MPI_DOUBLE, 0, tag,
                        MPI_COMM_WORLD);
            }
        }

      if (my_rank == 0)
        {
          print_block_line (curr_array, n, m, h, fp, H);
        }
    }
  if (my_rank == 0 && fp == stdout)
    {
#ifdef PRINT_LINE
      for (int i = 0; i <= 7 * N; i++)
        printf ("-");
      printf ("\n");
#endif
    }
}


int
mpi_print_matrix (const double *a, double *line, int n, int m,
                        int p, int my_rank, const char *file_name)
{
  Closer closer;
  FILE *fp = stdout;
  int error = ALL_RIGHT;

  if (file_name && my_rank == 0)
    {
      fp = fopen (file_name, "w");
      closer.add (fp);
      if (!fp)
        error = CANNOT_OPEN;
    }

  MPI_Bcast (&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (error != ALL_RIGHT)
    {
      return error;
    }

  mpi_print_matrix (a, line, n, m, p, my_rank, fp);

  return ALL_RIGHT;
}


void
print_block_line (const double *line, int n, int m, int h, FILE *fp, int max_h)
{
  int k = n / m;
  int l = n % m;
  int N = n;
  int H = h;

  if (fp == stdout)
    {
      N = (n > MAX_PRINT) ? MAX_PRINT : n;
      H = (h > max_h) ? max_h : h;
    }

//  printf ("h = %d, max_h = %d, H = %d\n", h, max_h, H);

  for (int i = 0; i < H; ++i)
    {
      for (int j = 0; j < N; ++j)
        {
          int w = (j >= k * m) ? l : m;
          int ind = (j / m) * h * m + i * w + j % m;
#ifdef SHORT_PRINT
          fprintf (fp, "\t%.2lf", line[ind]);
#else
          fprintf (fp, " %le", line[ind]);
#endif
        }
      fprintf (fp, "\n");
    }
}


void
print_matrix (const double *a, int n, int m, FILE *fp)
{
  if (m <= 0 || n <= 0)
    return;

  int matrix_width = m;
  if (fp == stdout)
    {
      n = (n > MAX_PRINT) ? MAX_PRINT : n;
      m = (m > MAX_PRINT) ? MAX_PRINT : m;
    }
  for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < m; ++j)
        {
#ifdef SHORT_PRINT
          fprintf (fp, "\t%.2lf", a[i * matrix_width + j]);
#else
          fprintf (fp, " %le", a[i * matrix_width + j]);
#endif
        }
      fprintf (fp, "\n");
    }
  if (fp == stdout)
    {
#ifdef PRINT_LINE
      for (int i = 0; i <= 7 * m; i++)
        printf ("- ");
      printf ("\n");
#endif
    }
}


int
print_matrix (const double *a, int n, int m, const char *file_name)
{
  Closer closer;
  FILE *fp = stdout;
  if (file_name)
    {
      fp = fopen (file_name, "w");
      closer.add (fp);
      if (!fp)
        return CANNOT_OPEN;
    }
  print_matrix (a, n, m, fp);
  return ALL_RIGHT;
}


