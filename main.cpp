#include "matrix_io.h"

void
start_algorithm (int argc, char *argv[], int p, int my_rank);

static int N;

double
init_func (int i, int j)
{
  return (double)(i * N + j);
  return (double)(i - j + 1);
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
  p_blocks = (blocks % p == 0) ? blocks / p : blocks / p + 1;

  int sum = 0;
  error = ALL_RIGHT;
  a = new double [p_blocks * n * m];
  b = new double [p_blocks * n * m];
  deleter.add (a);
  deleter.add (b);
  if (!a || !b)
    error = MEMORY_ERROR;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (sum != ALL_RIGHT)
    return;


#if 0
//#ifdef TEST
  file_name = "1_matr.txt";
  double *c = nullptr;

  if (my_rank == 0)
    {
      c = new double [n * n];
      deleter.add (c);
      if (!c)
        MPI_Abort (MPI_COMM_WORLD, 0);

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          c[i * n + j] = init_func (i, j);

      print_matrix (c, n, n, file_name);
    }
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


//  mpi_print_block_matrix_simple (a, b, n, m, p, my_rank);
  mpi_print_block_matrix (a, b, n, m, p, my_rank);


#if 0 // =====NEEDREMOVE=====
  error = 0;
  sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif //====================
}
