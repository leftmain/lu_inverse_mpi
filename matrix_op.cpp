#include "matrix_op.h"

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

