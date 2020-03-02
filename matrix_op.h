#ifndef MATRIX_OP
#define MATRIX_OP

#include "header.h"
#include "matrix_io.h"

void
nullify (double *a, int size);

void
sum_arrays (double *a, const double *b, int len);

// n x m matrix norm
double
norma (const double *a, int n, int m = 0);

double
line_norma (const double *line, double *block, int n, int m, int h = 0);

#endif
