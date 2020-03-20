#ifndef MATRIX_OP
#define MATRIX_OP

#include "header.h"
#include "matrix_io.h"

// ------------------------------------ memory operations

void
nullify (double *a, int size);

void
sum_arrays (double *a, const double *b, int len);

void
minus_array (double *a, int len);

void
copy_array (const double *a, double *b, int size);


// ------------------------------------ math operations

void
matrix_minus_E (double *a, int n);

// n x m matrix norm
double
norma (const double *a, int n, int m = 0);

void
matrix_multiply (const double *a, const double *b, double *c,
                 int m, int n, int k);

void
matrix_plus_multiply (const double *a, const double *b, double *c,
                      int m, int n, int k);

void
matrix_minus_multiply (const double *a, const double *b, double *c,
                       int m, int n, int k);

/* @return: 0 - all ok
            1 - cannot apply method
   result in b                      */
int
lu_matrix_inverse (double *a, double *b, int n);


// ------------------------------------ block operations

double
line_norma (const double *line, double *block, int n, int m, int h = 0);

void
scalar_product (const double *first_line, const double *second_line,
                double *blocks, int n, int m, int h, int w);


// ------------------------------------ setters / getters

void
get_block (const double *a, double *block, int n, int m,
           int p, int my_rank, int i, int j);

void
set_block (double *a, const double *block, int n, int m,
           int p, int my_rank, int i, int j);

void
get_line (const double *a, double *line, int n, int m,
          int p, int my_rank, int i, int from, int to);

void
set_line (double *a, const double *line, int n, int m,
          int p, int my_rank, int i, int from, int to);

void
get_column (const double *a, double *line, int n, int m,
            int p, int my_rank, int j);

void
get_column (const double *a, double *line, int n, int m,
            int p, int my_rank, int j, int from, int to);

void
get_line (const double *a, double *line, int n, int m,
          int p, int my_rank, int i);

void
get_block_from_line (const double *line, double *block, int n, int m, int h, int j);

void
set_block_to_line (double *line, const double *block, int n, int m, int h, int j);

int
get_block_size (int i, int k, int l, int m);

int
get_number_of_blocks (int k, int l);

int
get_p_blocks (int blocks, int p);

int
p_blocks_klp (int k, int l, int p);


#endif
