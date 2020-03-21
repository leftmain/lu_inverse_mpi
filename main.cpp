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
  // read args
  int n = 0;
  int m = 0;
  if (argc < 3 || argc > 4 || (n = atoi (argv[1])) <= 0
      || (m = atoi (argv[2])) <= 0 || m > n)
    {
      if (my_rank == 0)
        printf ("usage: %s n m [file_name]\n", argv[0]);
      return;
    }
  N = n;

  const char *file_name = nullptr;
  if (argc == 4)
    file_name = argv[3];


  // memory allocation
  int k = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k, l);
  int p_blocks = get_p_blocks (blocks, p);

  double *a = nullptr;
  double *b = nullptr;
  double *workspace = nullptr;
  a = new double [p_blocks * n * m];
  b = new double [p_blocks * n * m];
  workspace = new double [2 * p * p_blocks * m * m + n * m
                          + 3 * m * m + p];

  Deleter<double *> deleter;
  deleter.add (a);
  deleter.add (b);
  deleter.add (workspace);

  Time *time = nullptr;
  time = new Time[p];
  Deleter<Time *> time_deleter;
  time_deleter.add (time);

  int error = ALL_RIGHT;
  if (!a || !b || !workspace || !time)
    error = MEMORY_ERROR;

  int sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (sum != ALL_RIGHT)
    return;


  // fill matrix
  error = ALL_RIGHT;
  if (file_name)
    {
      error = mpi_read_matrix (a, workspace, n, m, p, my_rank, file_name);
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


  mpi_print_matrix_with_prefix (a, workspace, n, m, p, my_rank, "A");

  Time my_time;

  //  find A^(-1)
  error = mpi_lu_matrix_inverse (a, b, workspace, n, m, p, my_rank, my_time);
  if (error != ALL_RIGHT)
    {
      if (my_rank == 0)
        printf ("method cannot be applyed\n");
      return;
    }

  mpi_print_matrix_with_prefix (b, workspace, n, m, p, my_rank, "A^(-1)");

  //  read matrix again
  if (file_name)
    error = mpi_read_matrix (a, workspace, n, m, p, my_rank, file_name);
  else
    mpi_init_matrix (a, n, m, p, my_rank, init_func);

  //  residual
  double residual = 0.;
  if (n < 1000 || p > 1)
    residual = mpi_residual (a, b, workspace, n, m, p, my_rank, my_time);

  // print time
  int send_size = sizeof (Time) / sizeof (double);
  MPI_Allgather (&my_time, send_size, MPI_DOUBLE, time,
                 send_size, MPI_DOUBLE, MPI_COMM_WORLD);
  if (my_rank == 0)
    {
      printf ("algorithm global time: %.2lf\n", time[0].global_time);
      printf ("residual global time: %.2lf\n", time[0].residual_global_time);
      for (int i = 0; i < p; ++i)
        printf ("process %d time: %.2lf\n", i, time[i].process_time);

#if 1
      for (int i = 0; i < p; ++i)
        printf ("decomposition %d time: %.2lf\n", i, time[i].decomposition_time);
      for (int i = 0; i < p; ++i)
        printf ("invert %d time: %.2lf\n", i, time[i].invert_time);
      for (int i = 0; i < p; ++i)
        printf ("mult %d time: %.2lf\n", i, time[i].mult_time);
      for (int i = 0; i < p; ++i)
        printf ("residual %d time: %.2lf\n", i, time[i].residual_process_time);
#endif

      printf ("Residual = %le,  n = %d,  m = %d,  p = %d,  elapsed = %.2lf\n",
              residual, n, m, p, time[0].global_time);
    }
}
