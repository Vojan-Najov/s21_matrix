#include "s21_matrix.h"

double s21_fabs(double n) {
  union {
    double n;
    unsigned int bits[2];
  } un = {n};

  un.bits[1] &= 0x7FFFFFFFU;

  return (un.n);
}

int s21_gauss_method(matrix_t *A) {
  int rows = A->rows;
  int columns = A->columns;
  double **a = A->matrix;
  int i, j, k, l;
  double r;
  double eps = 1.0e-8;

  i = 0;
  j = 0;
  while (i < rows && j < columns) {
    r = 0.0;
    for (k = i; k < rows; ++k) {
      if (s21_fabs(a[k][j]) > r) {
        l = k;
        r = s21_fabs(a[k][j]);
      }
    }
    if (r <= eps) {
      for (k = i; k < rows; ++k) {
        a[k][j] = 0.0;
      }
      ++j;
      continue;
    }
    if (i != l) {
      for (k = j; k < columns; ++k) {
        r = a[i][k];
        a[i][k] = a[l][k];
        a[l][k] = -r;
      }
    }
    r = a[i][j];
    for (k = i + 1; k < rows; ++k) {
      double c = -a[k][j] / r;
      a[k][j] = 0.0;
      for (l = j + 1; l < columns; ++l) {
        a[k][l] += c * a[i][l];
      }
    }
    ++i;
    ++j;
  }

  return (i);
}
