#include "s21_matrix.h"

static int s21_determinant_aux(matrix_t *A, double *result);

int s21_determinant(matrix_t *A, double *result) {
  int k;
  double **a = A->matrix;
  double determinant = 0.0;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A)) {
    status = S21_MATRIX_ERROR;
  } else if (!MATRIX_SQUARE(A)) {
    status = S21_MATRIX_CALCULATION_ERROR;
  } else {
    k = A->rows;
    if (k == 1) {
      determinant = **a;
    } else if (k == 2) {
      determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    } else if (k == 3) {
      determinant = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] +
                    a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0] -
                    a[1][2] * a[2][1] * a[0][0] - a[2][2] * a[0][1] * a[1][0];
    } else {
      status = s21_determinant_aux(A, &determinant);
    }
  }

  *result = determinant;

  return (status);
}

static int s21_determinant_aux(matrix_t *A, double *result) {
  matrix_t m;
  double **a = A->matrix;
  int k = A->rows;
  int rank;
  double determinant = 0.0;
  int status = S21_MATRIX_OK;

  status = s21_create_matrix(k, k, &m);
  if (status == S21_MATRIX_OK) {
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < k; ++j) {
        m.matrix[i][j] = a[i][j];
      }
    }
    rank = s21_gauss_method(&m);
    if (rank < k) {
      determinant = 0.0;
    } else {
      determinant = 1.0;
      for (int i = 0; i < k; ++i) {
        determinant *= m.matrix[i][i];
      }
    }
    s21_remove_matrix(&m);
  }

  *result = determinant;

  return (status);
}
