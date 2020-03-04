#ifndef MATRIX_OP
#define MATRIX_OP

#include "header.h"
#include "matrix_io.h"

void
nullify (double *a, int size);

void
sum_arrays (double *a, const double *b, int len);

void
matrix_minus_E (double *a, int n);

// n x m matrix norm
double
norma (const double *a, int n, int m = 0);

double
line_norma (const double *line, double *block, int n, int m, int h = 0);

void
get_block (const double *a, double *block, int n, int m,
           int p, int my_rank, int i, int j);


void
set_block (double *a, const double *block, int n, int m,
           int p, int my_rank, int i, int j);

void
get_row (const double *a, double *line, int n, int m,
         int p, int my_rank, int i, int from, int to);

void
set_row (double *a, const double *line, int n, int m,
         int p, int my_rank, int i, int from, int to);

void
matrix_multiply (const double *a, const double *b, double *c,
                 int m, int n, int k);

void
matrix_plus_multiply (const double *a, const double *b, double *c,
                      int m, int n, int k);

void
matrix_minus_multiply (const double *a, const double *b, double *c,
                       int m, int n, int k);

int
block_height_iklm (int i, int k, int l, int m);

int
block_width_jklm (int j, int k, int l, int m);

int
number_of_blocks_kl (int k, int l);

#endif
