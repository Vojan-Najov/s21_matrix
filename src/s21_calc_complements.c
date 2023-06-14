#include "s21_matrix.h"

static int s21_minor(int i, int j, matrix_t *A, double *result);
static double s21_minor_aux(matrix_t *M);

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int n;
  double minor;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else if (!MATRIX_SQUARE(A) || A->rows < 2) {
    MATRIX_RESET(result);
    status = S21_MATRIX_CALCULATION_ERROR;
  } else {
    n = A->rows;
    status = s21_create_matrix(n, n, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; status == S21_MATRIX_OK && i < n; ++i) {
        for (int j = 0; status == S21_MATRIX_OK && j < n; ++j) {
          status = s21_minor(i, j, A, &minor);
          minor = (i + j) % 2 ? -minor : minor;
          result->matrix[i][j] = minor;
        }
      }
    }
    if (status) {
      s21_remove_matrix(result);
      status = S21_MATRIX_ERROR;
    }
  }

  return (status);
}

static int s21_minor(int i, int j, matrix_t *A, double *result) {
  int k, l, n;
  matrix_t m;
  double **a = A->matrix;
  double minor = 0.0;
  int status = S21_MATRIX_OK;

  n = A->rows;
  status = s21_create_matrix(A->rows - 1, A->columns - 1, &m);
  if (status == S21_MATRIX_OK) {
    for (k = 0; k < i; ++k) {
      for (l = 0; l < j; ++l) {
        m.matrix[k][l] = a[k][l];
      }
      for (++l; l < n; ++l) {
        m.matrix[k][l - 1] = a[k][l];
      }
    }
    for (++k; k < n; ++k) {
      for (l = 0; l < j; ++l) {
        m.matrix[k - 1][l] = a[k][l];
      }
      for (++l; l < n; ++l) {
        m.matrix[k - 1][l - 1] = a[k][l];
      }
    }
    minor = s21_minor_aux(&m);
    s21_remove_matrix(&m);
  }

  *result = minor;

  return (status);
}

static double s21_minor_aux(matrix_t *M) {
  double minor = 1.0;
  double **m = M->matrix;
  int n = M->rows;
  int rank;

  if (n == 1) {
    minor = **m;
  } else if (n == 2) {
    minor = m[0][0] * m[1][1] - m[0][1] * m[1][0];
  } else if (n == 3) {
    minor = m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] +
            m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0] -
            m[1][2] * m[2][1] * m[0][0] - m[2][2] * m[0][1] * m[1][0];
  } else {
    rank = s21_gauss_method(M);
    if (rank < n) {
      minor = 0.0;
    } else {
      for (int i = 0; i < n; ++i) {
        minor *= m[i][i];
      }
    }
  }

  return (minor);
}
