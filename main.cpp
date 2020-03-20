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

double
init_func_1 (int i, int j)
{
  return N - std::max (i, j) - 2;
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
  int blocks = get_number_of_blocks (k, l); //(l == 0) ? k : k + 1;
  int p_blocks = get_p_blocks (blocks, p);
//  p_blocks = (blocks % p == 0) ? blocks / p : blocks / p + 1;
//  p_blocks = (blocks + p - 1) / p;

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

  int error = ALL_RIGHT;
  if (!a || !b)
    error = MEMORY_ERROR;

  int sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (sum != ALL_RIGHT)
    return;


#ifdef RES_TEST
  file_name = "1_matr.txt";
  double *c = nullptr;
  double *bb = nullptr;
  double *v = nullptr;
  double c_norm = 0.;


  auto f_residual = [](const double *a, const double *b, double *vector, int n)
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
            bb[i * n + j] = init_func_1 (i, j);
          }

      c_norm = f_residual (c, bb, v, n);
      print_matrix (c, n, n);
//      print_matrix (c, n, n, file_name);
      printf ("c_norm = %le\n", c_norm);
//      print_matrix (&c_norm, 1, 1, "norma.txt");
      delete [] c;
      delete [] bb;
      delete [] v;
    }
  MPI_Barrier (MPI_COMM_WORLD);
#endif // TEST


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


  if (my_rank == 0)
    {
      printf ("A:\n");
    }
  mpi_print_block_matrix (a, workspace, n, m, p, my_rank);


  double global_time = 0.;
  if (my_rank == 0)
    global_time = get_earth_time ();

  double time = 0.;
  error = mpi_lu_matrix_inverse (a, b, workspace, n, m, p, my_rank, &time);

  if (my_rank == 0)
    {
      global_time = get_earth_time () - global_time;
      if (error != ALL_RIGHT)
        {
          printf ("method cannot be applyed\n");
          return;
        }
    }


  if (my_rank == 0)
    {
      printf ("A^(-1):\n");
    }
  mpi_print_block_matrix (b, workspace, n, m, p, my_rank);


  double *p_buf = workspace;
  MPI_Allgather (&time, 1, MPI_DOUBLE, p_buf, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  if (my_rank == 0)
    {
      printf ("global time: %.2lf\n", global_time);
      for (int i = 0; i < p; ++i)
        {
          printf ("proccess %d time: %.2lf\n", i, p_buf[i]);
        }
    }


#ifdef CHECK
  mpi_print_block_matrix (b, workspace, n, m, p, my_rank, "res.txt");
  double *aa = new double[n*n];
  double *bb = new double[n*n];
  double *cc = new double[n*n];
  FILE *fp;
  MPI_Barrier (MPI_COMM_WORLD);
  bool ok = true;
  if (my_rank == 0)
    {
      fp = fopen ("res.txt","r");
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          {
            aa[i * n + j] = init_func (i, j);
          }
      lu_matrix_inverse (aa, bb, n);
  //    print_matrix (bb, n, n);

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          {
            if (fscanf (fp, "%lf", aa + i * n + j) != 1)
              MPI_Abort (MPI_COMM_WORLD, 0);
            bb[i * n + j] = init_func (i, j);
          }
      fclose (fp);
      matrix_multiply (aa, bb, cc, n, n, n);
  //    print_matrix (cc, n, n);
      double norm = norma (cc, n, n);
      for (int i = 0; i < n; ++i)
        {
          for (int j = 0; j < n; ++j)
            {
              double offset = (i == j) ? 1. : 0.;
              if (fabs (cc[i * n + j] - offset) >= 1e-16 * norm)
                {
                  printf ("BAD_MATRIX:\n");
                  print_matrix (cc, n, n);
                  i = j = n;
                  ok = false;
                  break;
                }
            }
        }
    }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] aa;
  delete [] bb;
  delete [] cc;
#endif // CHECK


  //  read matrix again
  if (file_name)
    {
      error = mpi_read_matrix (a, workspace, n, m, p, my_rank, file_name);
    }
  else
    {
      mpi_init_matrix (a, n, m, p, my_rank, init_func);
    }

  double residual = mpi_residual (a, b, workspace, n, m, p, my_rank);

  if (my_rank == 0)
    {
#ifdef CHECK
      if (ok)
        printf ("OK\t");
      else
        printf ("WRONG\t");
#endif // CHECK
      printf ("n = %d\tm = %d\tp = %d\tresidual = %le\n", n, m, p, residual);
    }

#ifdef RES_TEST
  memcpy (b, a, p_blocks * n * m * sizeof (double));
  mpi_init_matrix (b, n, m, p, my_rank, init_func_1);
  double res = mpi_residual (a, b, workspace, n, m, p, my_rank);
  if (my_rank == 0)
    {
      printf ("residual = %le\n", res);
    }
#endif


#if 0 // =====NEEDREMOVE=====
  error = 0;
  sum = 0;
  MPI_Allreduce (&error, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif //====================
}
