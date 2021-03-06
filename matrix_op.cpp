#include "matrix_op.h"

int
get_block_size (int i, int k, int l, int m)
{
  return (i == k) ? l : m;
}

int
get_number_of_blocks (int k, int l)
{
  return (l == 0) ? k : k + 1;
}

int
get_p_blocks (int blocks, int p)
{
  return (blocks + p - 1) / p;
}

int
p_blocks_klp (int k, int l, int p)
{
  return (get_number_of_blocks (k, l) + p - 1) / p;
}

void
nullify (double *a, int size)
{
  memset (a, 0, size * sizeof (double));
}


void
sum_arrays (double *a, const double *b, int len)
{
  for (int i = 0; i < (len & 7); ++i)
    a[i] += b[i];
  for (int i = (len & 7); i < len; i += 8)
    {
      a[i + 0] += b[i + 0];
      a[i + 1] += b[i + 1];
      a[i + 2] += b[i + 2];
      a[i + 3] += b[i + 3];
      a[i + 4] += b[i + 4];
      a[i + 5] += b[i + 5];
      a[i + 6] += b[i + 6];
      a[i + 7] += b[i + 7];
    }
}


void
minus_array (double *a, int len)
{
  for (int i = 0; i < (len & 7); ++i)
    a[i] = -a[i];
  for (int i = (len & 7); i < len; i += 8)
    {
      a[i + 0] = -a[i + 0];
      a[i + 1] = -a[i + 1];
      a[i + 2] = -a[i + 2];
      a[i + 3] = -a[i + 3];
      a[i + 4] = -a[i + 4];
      a[i + 5] = -a[i + 5];
      a[i + 6] = -a[i + 6];
      a[i + 7] = -a[i + 7];
    }
}


void
copy_array (const double *a, double *b, int size)
{
  memcpy (b, a, size * sizeof (double));
}


void
matrix_minus_E (double *a, int n)
{
  for (int i = 0; i < (n & 7); ++i)
    a[i * n + i] -= 1.;
  for (int i = (n & 7); i < n; i += 8)
    {
      a[(i + 0) * n + i + 0] -= 1.;
      a[(i + 1) * n + i + 1] -= 1.;
      a[(i + 2) * n + i + 2] -= 1.;
      a[(i + 3) * n + i + 3] -= 1.;
      a[(i + 4) * n + i + 4] -= 1.;
      a[(i + 5) * n + i + 5] -= 1.;
      a[(i + 6) * n + i + 6] -= 1.;
      a[(i + 7) * n + i + 7] -= 1.;
    }
}


double
norma (const double *a, int n, int m)
{
  if (m == 0)
    m = n;

  double max = 0.;

  double s0 = 0.;
  double s1 = 0.;
  double s2 = 0.;
  double s3 = 0.;
  double s4 = 0.;
  double s5 = 0.;
  double s6 = 0.;
  double s7 = 0.;

  for (int i = 0; i < n; ++i)
    {
      s0 = s1 = s2 = s3 = s4 = s5 = s6 = s7 = 0.;

      for (int j = 0; j < (m & 7); ++j)
        s0 += fabs (a[i * m + j]);

      for (int j = (m & 7); j < m; j += 8)
        {
          s0 += fabs (a[i * m + j + 0]);
          s1 += fabs (a[i * m + j + 1]);
          s2 += fabs (a[i * m + j + 2]);
          s3 += fabs (a[i * m + j + 3]);
          s4 += fabs (a[i * m + j + 4]);
          s5 += fabs (a[i * m + j + 5]);
          s6 += fabs (a[i * m + j + 6]);
          s7 += fabs (a[i * m + j + 7]);
        }

      s0 += s1 + s2 + s3 + s4 + s5 + s6 + s7;

      if (s0 > max)
        max = s0;
    }
  return max;
}


double
line_norma (const double *line, double *block, int n, int m, int h)
{
  if (h == 0)
    h = m;

  int k = n / m;
  int l = n % m;

  nullify (block, h * m);
  for (int i = 0; i < k; ++i)
    {
      sum_arrays (block, line + i * h * m, h * m);
    }
  if (l > 0)
    {
      line += k * h * m;
      for (int i = 0; i < h; ++i)
        {
          for (int j = 0; j < l; ++j)
            {
              block[i * m] += line[i * l + j];
            }
        }
    }

  return norma (block, h, m);
}


void
get_block (const double *a, double *block, int n, int m,
           int p, int my_rank, int i, int j)
{
  if (i % p != my_rank)
    {
      printf ("VERY-VERY-BAD - came to unreachable code\n");
      MPI_Abort (MPI_COMM_WORLD, 0);
    }

  int k = n / m;
  int l = n % m;
//  int h = (i == k) ? l : m;
//  int w = (j == k) ? l : m;
  int h = get_block_size (i, k, l, m);
  int w = get_block_size (j, k, l, m);

  memcpy (block, a + (i / p) * n * m + j * h * m,
          h * w * sizeof (double));
}


void
get_block_from_line (const double *line, double *block, int n, int m, int h, int j)
{
  int k = n / m;
  int l = n % m;
  int w = get_block_size (j, k, l, m);

  memcpy (block, line + j * h * m,
          h * w * sizeof (double));
}


void
set_block_to_line (double *line, const double *block, int n, int m, int h, int j)
{
  int k = n / m;
  int l = n % m;
  int w = get_block_size (j, k, l, m);

  memcpy (line + j * h * m, block,
          h * w * sizeof (double));
}


void
set_block (double *a, const double *block, int n, int m,
           int p, int my_rank, int i, int j)
{
  if (i % p != my_rank)
    MPI_Abort (MPI_COMM_WORLD, 0);

  int k = n / m;
  int l = n % m;
//  int h = (i == k) ? l : m;
//  int w = (j == k) ? l : m;
  int h = get_block_size (i, k, l, m);
  int w = get_block_size (j, k, l, m);

  memcpy (a + (i / p) * n * m + j * h * m, block,
          h * w * sizeof (double));
}


void
get_line (const double *a, double *line, int n, int m,
          int p, int my_rank, int i, int from, int to)
{
  if (i % p != my_rank)
    MPI_Abort (MPI_COMM_WORLD, 0);

  int k = n / m;
  int l = n % m;
//  int h = (i == k) ? l : m;
  int h = get_block_size (i, k, l, m);
  int size = (to - from) * h * m;
  if (to == k)
    size += h * l;
  else
    size += h * m;

  memcpy (line, a + (i / p) * n * m + from * h * m,
          size * sizeof (double));
}


void
set_line (double *a, const double *line, int n, int m,
          int p, int my_rank, int i, int from, int to)
{
  if (i % p != my_rank)
    return;

  int k = n / m;
  int l = n % m;
//  int h = (i == k) ? l : m;
  int h = get_block_size (i, k, l, m);
  int size = (to == k) ? (to - from) * h * m + h * l
                       : (to - from + 1) * h * m;

  memcpy (a + (i / p) * n * m + from * h * m, line,
          size * sizeof (double));
}


void
get_column (const double *a, double *line, int n, int m,
            int p, int my_rank, int j)
{
  int k = n / m;
  int l = n % m;
  int blocks = get_number_of_blocks (k, l);
  int offset = 0;
  int w = get_block_size (j, k, l, m);
  for (int i = my_rank; i < blocks; i += p)
    {
      get_block (a, line + offset, n, m, p, my_rank, i, j);
      offset += w * m;
    }
}


void
get_column (const double *a, double *line, int n, int m,
            int p, int my_rank, int j, int from, int to)
{
  int k = n / m;
  int l = n % m;
  int offset = 0;
  int w = get_block_size (j, k, l, m);
  int f = from - from % p + my_rank;
  if (f < from)
    f += p;
  for (int i = f; i <= to; i += p)
    {
      get_block (a, line + offset, n, m, p, my_rank, i, j);
      offset += w * m;
    }
}


void
get_line (const double *a, double *line, int n, int m,
          int p, int my_rank, int i)
{
  if (i % p != my_rank)
    MPI_Abort (MPI_COMM_WORLD, 0);

  int k = n / m;
  int l = n % m;
  int h = get_block_size (i, k, l, m);
  memcpy (line, a + (i / p) * m * n, h * n * sizeof (double));
}


void
scalar_product (const double *first_line, const double *second_line,
                double *result, int n, int m, int h, int w)
{
  int k = n / m;
  int l = n % m;
  int n_blocks = get_number_of_blocks (k, l);

  nullify (result, m * m);

  for (int i = 0; i < n_blocks; ++i)
    {
      int w1 = (i == k) ? l : m;
      matrix_plus_multiply (first_line + i * h * m, second_line + i * m * w,
                            result, h, w1, w);
    }
}

void
matrix_multiply (const double *a, const double *b, double *c,
                 int m, int n, int k)
{
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] = d;
        }
    }
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = (k & 3); j < k; ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] = d;
        }
    }
  for (int i = (m & 3); i < m; ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] = d;
        }
    }
  for (int i = (m & 3); i < m; i += 4)
    {
      for (int j = (k & 3); j < k; j += 4)
        {
          double d00 = 0.;
          double d01 = 0.;
          double d02 = 0.;
          double d03 = 0.;

          double d10 = 0.;
          double d11 = 0.;
          double d12 = 0.;
          double d13 = 0.;

          double d20 = 0.;
          double d21 = 0.;
          double d22 = 0.;
          double d23 = 0.;

          double d30 = 0.;
          double d31 = 0.;
          double d32 = 0.;
          double d33 = 0.;

          for (int t = 0; t < n; ++t)
            {
              d00 += a[(i + 0) * n + t] * b[t * k + j + 0];
              d01 += a[(i + 0) * n + t] * b[t * k + j + 1];
              d02 += a[(i + 0) * n + t] * b[t * k + j + 2];
              d03 += a[(i + 0) * n + t] * b[t * k + j + 3];

              d10 += a[(i + 1) * n + t] * b[t * k + j + 0];
              d11 += a[(i + 1) * n + t] * b[t * k + j + 1];
              d12 += a[(i + 1) * n + t] * b[t * k + j + 2];
              d13 += a[(i + 1) * n + t] * b[t * k + j + 3];

              d20 += a[(i + 2) * n + t] * b[t * k + j + 0];
              d21 += a[(i + 2) * n + t] * b[t * k + j + 1];
              d22 += a[(i + 2) * n + t] * b[t * k + j + 2];
              d23 += a[(i + 2) * n + t] * b[t * k + j + 3];

              d30 += a[(i + 3) * n + t] * b[t * k + j + 0];
              d31 += a[(i + 3) * n + t] * b[t * k + j + 1];
              d32 += a[(i + 3) * n + t] * b[t * k + j + 2];
              d33 += a[(i + 3) * n + t] * b[t * k + j + 3];
            }

          c[(i + 0) * k + j + 0] = d00;
          c[(i + 0) * k + j + 1] = d01;
          c[(i + 0) * k + j + 2] = d02;
          c[(i + 0) * k + j + 3] = d03;

          c[(i + 1) * k + j + 0] = d10;
          c[(i + 1) * k + j + 1] = d11;
          c[(i + 1) * k + j + 2] = d12;
          c[(i + 1) * k + j + 3] = d13;

          c[(i + 2) * k + j + 0] = d20;
          c[(i + 2) * k + j + 1] = d21;
          c[(i + 2) * k + j + 2] = d22;
          c[(i + 2) * k + j + 3] = d23;

          c[(i + 3) * k + j + 0] = d30;
          c[(i + 3) * k + j + 1] = d31;
          c[(i + 3) * k + j + 2] = d32;
          c[(i + 3) * k + j + 3] = d33;
        }
    }
}

void
matrix_plus_multiply (const double *a, const double *b, double *c,
                      int m, int n, int k)
{
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] += d;
        }
    }
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = (k & 3); j < k; ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] += d;
        }
    }
  for (int i = (m & 3); i < m; ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] += d;
        }
    }
  for (int i = (m & 3); i < m; i += 4)
    {
      for (int j = (k & 3); j < k; j += 4)
        {
          double d00 = 0.;
          double d01 = 0.;
          double d02 = 0.;
          double d03 = 0.;

          double d10 = 0.;
          double d11 = 0.;
          double d12 = 0.;
          double d13 = 0.;

          double d20 = 0.;
          double d21 = 0.;
          double d22 = 0.;
          double d23 = 0.;

          double d30 = 0.;
          double d31 = 0.;
          double d32 = 0.;
          double d33 = 0.;

          for (int t = 0; t < n; ++t)
            {
              d00 += a[(i + 0) * n + t] * b[t * k + j + 0];
              d01 += a[(i + 0) * n + t] * b[t * k + j + 1];
              d02 += a[(i + 0) * n + t] * b[t * k + j + 2];
              d03 += a[(i + 0) * n + t] * b[t * k + j + 3];

              d10 += a[(i + 1) * n + t] * b[t * k + j + 0];
              d11 += a[(i + 1) * n + t] * b[t * k + j + 1];
              d12 += a[(i + 1) * n + t] * b[t * k + j + 2];
              d13 += a[(i + 1) * n + t] * b[t * k + j + 3];

              d20 += a[(i + 2) * n + t] * b[t * k + j + 0];
              d21 += a[(i + 2) * n + t] * b[t * k + j + 1];
              d22 += a[(i + 2) * n + t] * b[t * k + j + 2];
              d23 += a[(i + 2) * n + t] * b[t * k + j + 3];

              d30 += a[(i + 3) * n + t] * b[t * k + j + 0];
              d31 += a[(i + 3) * n + t] * b[t * k + j + 1];
              d32 += a[(i + 3) * n + t] * b[t * k + j + 2];
              d33 += a[(i + 3) * n + t] * b[t * k + j + 3];
            }

          c[(i + 0) * k + j + 0] += d00;
          c[(i + 0) * k + j + 1] += d01;
          c[(i + 0) * k + j + 2] += d02;
          c[(i + 0) * k + j + 3] += d03;

          c[(i + 1) * k + j + 0] += d10;
          c[(i + 1) * k + j + 1] += d11;
          c[(i + 1) * k + j + 2] += d12;
          c[(i + 1) * k + j + 3] += d13;

          c[(i + 2) * k + j + 0] += d20;
          c[(i + 2) * k + j + 1] += d21;
          c[(i + 2) * k + j + 2] += d22;
          c[(i + 2) * k + j + 3] += d23;

          c[(i + 3) * k + j + 0] += d30;
          c[(i + 3) * k + j + 1] += d31;
          c[(i + 3) * k + j + 2] += d32;
          c[(i + 3) * k + j + 3] += d33;
        }
    }
}

void
matrix_minus_multiply (const double *a, const double *b, double *c,
                       int m, int n, int k)
{
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] -= d;
        }
    }
  for (int i = 0; i < (m & 3); ++i)
    {
      for (int j = (k & 3); j < k; ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] -= d;
        }
    }
  for (int i = (m & 3); i < m; ++i)
    {
      for (int j = 0; j < (k & 3); ++j)
        {
          double d = 0.;
          for (int t = 0; t < n; ++t)
            {
              d += a[i * n + t] * b[t * k + j];
            }
          c[i * k + j] -= d;
        }
    }
  for (int i = (m & 3); i < m; i += 4)
    {
      for (int j = (k & 3); j < k; j += 4)
        {
          double d00 = 0.;
          double d01 = 0.;
          double d02 = 0.;
          double d03 = 0.;

          double d10 = 0.;
          double d11 = 0.;
          double d12 = 0.;
          double d13 = 0.;

          double d20 = 0.;
          double d21 = 0.;
          double d22 = 0.;
          double d23 = 0.;

          double d30 = 0.;
          double d31 = 0.;
          double d32 = 0.;
          double d33 = 0.;

          for (int t = 0; t < n; ++t)
            {
              d00 += a[(i + 0) * n + t] * b[t * k + j + 0];
              d01 += a[(i + 0) * n + t] * b[t * k + j + 1];
              d02 += a[(i + 0) * n + t] * b[t * k + j + 2];
              d03 += a[(i + 0) * n + t] * b[t * k + j + 3];

              d10 += a[(i + 1) * n + t] * b[t * k + j + 0];
              d11 += a[(i + 1) * n + t] * b[t * k + j + 1];
              d12 += a[(i + 1) * n + t] * b[t * k + j + 2];
              d13 += a[(i + 1) * n + t] * b[t * k + j + 3];

              d20 += a[(i + 2) * n + t] * b[t * k + j + 0];
              d21 += a[(i + 2) * n + t] * b[t * k + j + 1];
              d22 += a[(i + 2) * n + t] * b[t * k + j + 2];
              d23 += a[(i + 2) * n + t] * b[t * k + j + 3];

              d30 += a[(i + 3) * n + t] * b[t * k + j + 0];
              d31 += a[(i + 3) * n + t] * b[t * k + j + 1];
              d32 += a[(i + 3) * n + t] * b[t * k + j + 2];
              d33 += a[(i + 3) * n + t] * b[t * k + j + 3];
            }

          c[(i + 0) * k + j + 0] -= d00;
          c[(i + 0) * k + j + 1] -= d01;
          c[(i + 0) * k + j + 2] -= d02;
          c[(i + 0) * k + j + 3] -= d03;

          c[(i + 1) * k + j + 0] -= d10;
          c[(i + 1) * k + j + 1] -= d11;
          c[(i + 1) * k + j + 2] -= d12;
          c[(i + 1) * k + j + 3] -= d13;

          c[(i + 2) * k + j + 0] -= d20;
          c[(i + 2) * k + j + 1] -= d21;
          c[(i + 2) * k + j + 2] -= d22;
          c[(i + 2) * k + j + 3] -= d23;

          c[(i + 3) * k + j + 0] -= d30;
          c[(i + 3) * k + j + 1] -= d31;
          c[(i + 3) * k + j + 2] -= d32;
          c[(i + 3) * k + j + 3] -= d33;
        }
    }
}


static bool
is_equal (const double first, const double second, const double n = -1.)
{
  static double norma_ = 0.;
  if (n >= 0)
    {
      if (n < std::numeric_limits<double>::epsilon ())
        norma_ = std::numeric_limits<double>::epsilon ();
      else
        norma_ = n;
      return true;
    }
  if (fabs (first - second) < std::numeric_limits<double>::epsilon () * norma_)
    return true;
  else
    return false;
}


static int
lu_decomposition (double *a, int n)
{
  // U_{1,j} = (L_{1,1})^(-1) * A_{1,j}
  if (is_equal (a[0], 0.))
    return 1;
  for (int k = 1; k < n; ++k)
    a[k] /= a[0];

  for (int k = 1; k < n; ++k)
    {
      // L_{i,k} = A_{i,k} - sum_{j=1}^{k-1}(L_{i,j} * U_{j,k})
      for (int i = k; i < n; ++i)
        {
          double sum = 0.;
          for (int j = 0; j < k; ++j)
            {
              sum += a[i * n + j] * a[j * n + k];
            }
          a[i * n + k] -= sum;
        }

      // U_{k,i} = (A_{k,i} - sum_{j=1}^{k-1}(L_{k,j} * U_{j,i})) / L_{i,i}
      for (int i = k + 1; i < n; ++i)
        {
          double sum = 0.;
          for (int j = 0; j < k; ++j)
            {
              sum += a[k * n + j] * a[j * n + i];
            }

          if (is_equal (a[k * n + k], 0.))
            return 1;

          a[k * n + i] -= sum;
          a[k * n + i] /= a[k * n + k];
        }
    }
  return 0;
}


static int
lu_invert (double *a, int n)
{
  // L^(-1)_{i,i} = (L_{i,i})^(-1)
  // L^(-1)_{i,j} = -L^(-1)_{i,i} * (sum_{s=j}^{i-1}(L_{i,s} * L^(-1)_{s,j}))
  //
  // L_{0,0} = ..
  if (is_equal (a[0 * n + 0], 0.))
    return 1;
  a[0 * n + 0] = 1. / a[0 * n + 0];
  for (int i = 1; i < n; ++i)
    {
      // L_{i,i} = ..
      if (is_equal(a[i * n + i], 0.))
        return 1;
      a[i * n + i] = 1. / a[i * n + i];

      // L_{i,j} = ..
      for (int j = 0; j < i; ++j)
        {
          double sum = 0.;
          for (int s = j; s < i; ++s)
            {
              sum += a[i * n + s] * a[s * n + j];
            }
          a[i * n + j] = -sum * a[i * n + i];
        }
    }

  // U_{i,j} = -sum_{s=i}^{j-1} (U^(-1)_{i,s} * U_{s,j})
  for (int i = 0; i < n - 1; ++i)
    {
      for (int j = i + 1; j < n; ++j)
        {
          double sum = a[i * n + j];
          for (int s = i + 1; s < j; ++s)
            {
              sum += a[i * n + s] * a[s * n + j];
            }
          a[i * n + j] = -sum;
        }
    }
  return 0;
}


static void
ul_multiply (const double *a, double *b, int n)
{
  for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
        {
          double sum = 0.;
          int k = std::max (i, j);
          if (k == i)
            sum += a[i * n + j];
          else
            sum += a[i * n + k] * a[k * n + j];
          for (++k; k < n; ++k)
            {
              sum += a[i * n + k] * a[k * n + j];
            }
          b[i * n + j] = sum;
        }
    }
}


int
lu_matrix_inverse (double *a, double *b, int n)
{
  is_equal (0, 0, norma (a, n));

  if (lu_decomposition (a, n))
    return 1;

  if (lu_invert (a, n))
    return 1;

  ul_multiply (a, b, n);

  return 0;
}


double
get_earth_time ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
}


double
get_process_time ()
{
  struct rusage time;
  getrusage (RUSAGE_SELF, &time);
  return (double)time.ru_utime.tv_sec +
         (double)time.ru_utime.tv_usec * 1e-6;
}
