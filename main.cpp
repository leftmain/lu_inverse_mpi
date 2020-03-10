#include "matrix_io.h"
#include "matrix_op.h"
#include "matrix_mpi_op.h"

void
start_algorithm (int argc, char *argv[], int p, int my_rank);

static int N;

double
init_func (int i, int j)
{
  return N - std::max (i, j);
  return (i == j) ? 1 : 0;
  return (double)(i * N + j);
  return (double)(i - j + 1);
  return 1. / (i + j + 1);
}


int
main (int argc, char *argv[])
{
  int my_rank;
  int p;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  if (p < 1)
    MPI_Abort (MPI_COMM_WORLD, 1);

  start_algorithm (argc, argv, p, my_rank);

  MPI_Finalize ();


  return 0;
}



void
start_algorithm (int argc, char *argv[], int p, int my_rank)
{
  Deleter<double *> deleter;
  const char *file_name = nullptr;
  double *a = nullptr;
  double *b = nullptr;
  double *work_space = nullptr;
  int blocks = 0; // blocks in string
  int p_blocks = 0; // blocks in column (for each process)
  int n = 0;
  int m = 0;
  int k = 0;
  int l = 0;
  int error = 0;


  // read args
  if (argc < 3 || argc > 4 || (n = atoi (argv[1])) <= 0
      || (m = atoi (argv[2])) <= 0 || m > n)
    {
      if (my_rank == 0)
        printf ("usage: %s n m [file_name]\n", argv[0]);
      return;
    }
  if (argc == 4)
    file_name = argv[3];
  N = n;


  // memory allocation
  k = n / m;
  l = n % m;
  blocks = (l == 0) ? k : k + 1;
//  p_blocks = (blocks % p == 0) ? blocks / p : blocks / p + 1;
//  p_blocks = (blocks + p - 1) / p;
  p_blocks = p_blocks_blocksp (blocks, p);

  int sum = 0;
  error = ALL_RIGHT;
  a = new double [p_blocks * n * m];
  b = new double [p_blocks * n * m];
  work_space = new double [p * p_blocks * m * m + n * m
                           + p_blocks * m + 3 * m * m + p];
  deleter.add (a);
  deleter.add (b);
  deleter.add (work_space);
  if (!a || !b)
    error = MEMORY_ERROR;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (sum != ALL_RIGHT)
    return;


#ifdef TEST
  file_name = "1_matr.txt";
  double *c = nullptr;
  double *bb = nullptr;
  double *v = nullptr;
  double c_norm = 0.;


  auto residual = [](const double *a, const double *b, double *vector, int n)
  {
    double max_string_sum = 0;
    for (int i = 0; i < n; ++i)
      {
        double string_sum = 0.;
        for (int j = 0; j < n; ++j)
          {
            double sum = 0.;
            for (int t = 0; t < n; ++t)
              sum += a[i * n + t] * b[t * n + j];
            string_sum += fabs (sum);
          }
        string_sum -= 1.;
        if (i == 0)
          max_string_sum = string_sum;
        else if (string_sum > max_string_sum)
          max_string_sum = string_sum;
      }

    return max_string_sum;
  };

  if (my_rank == 0)
    {
      c = new double [n * n];
      bb = new double [n * n];
      v = new double [n];
      if (!c || !bb)
        MPI_Abort (MPI_COMM_WORLD, 0);

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          {
            c[i * n + j] = init_func (i, j);
            bb[i * n + j] = init_func (i, j);
          }

      c_norm = residual (c, bb, v, n);
      print_matrix (c, n, n, file_name);
      delete [] c;
      delete [] bb;
      delete [] v;
    }
  MPI_Barrier (MPI_COMM_WORLD);
#endif


  // fill matrix
  error = ALL_RIGHT;
  if (file_name)
    {
      error = mpi_read_matrix (a, b, n, m, p, my_rank, file_name);
    }
  else
    {
      mpi_init_matrix (a, n, m, p, my_rank, init_func);
    }
  if (error != ALL_RIGHT)
    {
      if (my_rank == 0)
        {
          switch (error)
            {
            case CANNOT_OPEN:
              fprintf (stderr, "cannot open %s\n", file_name);
              break;
            case CANNOT_READ:
              fprintf (stderr, "cannot read %s\n", file_name);
              break;
            default:
              fprintf (stderr, "unknown error\n");
            }
        }
      return;
    }


  mpi_print_block_matrix (a, work_space, n, m, p, my_rank);
  memcpy (b, a, p_blocks * n * m * sizeof (double));
  mpi_print_block_matrix (b, work_space, n, m, p, my_rank);

//  mpi_print_block_matrix (a, b, n, m, p, my_rank);
//  double norm = mpi_norma (a, work_space, n, m, p, my_rank);
  double res = mpi_residual (a, b, work_space, n, m, p, my_rank);
  if (my_rank == 0)
    {
      printf ("residual = %lf\n", res);
    }

#ifdef TEST
//  double norm = mpi_norma (a, line, n, m, p, my_rank);
  if (my_rank == 0)
    {
      printf ("cnorm = %lf\n", c_norm);
//      print_matrix (&c_norm, 1, 1, "2_norm.txt");
    }
#endif


#if 0 // =====NEEDREMOVE=====
  error = 0;
  sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif //====================
}
