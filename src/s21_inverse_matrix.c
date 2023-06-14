#include "s21_matrix.h"

static int s21_gauss_jordan_method(int n, double **a, double **b);
static int s21_gj_method_to_triangular_matrix(int n, double **a, double **b);
static void s21_gj_method_to_unit_matrix(int n, double **a, double **b);

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int n;
  matrix_t copy;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else if (!MATRIX_SQUARE(A)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_CALCULATION_ERROR;
  } else {
    n = A->rows;
    status = s21_create_matrix(n, n, &copy);
    status += s21_create_matrix(n, n, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          copy.matrix[i][j] = A->matrix[i][j];
        }
      }
      for (int i = 0; i < n; ++i) {
        result->matrix[i][i] = 1.0;
      }
      status = s21_gauss_jordan_method(n, copy.matrix, result->matrix);
      s21_remove_matrix(&copy);
      if (status == S21_MATRIX_CALCULATION_ERROR) {
        s21_remove_matrix(result);
      }
    }
  }

  return (status);
}

static int s21_gauss_jordan_method(int n, double **a, double **b) {
  int status = S21_MATRIX_OK;

  status = s21_gj_method_to_triangular_matrix(n, a, b);
  if (status == S21_MATRIX_OK) {
    s21_gj_method_to_unit_matrix(n, a, b);
  }

  return (status);
}

static int s21_gj_method_to_triangular_matrix(int n, double **a, double **b) {
  int i, j, k, l;
  double r, c;
  double eps = 1.0e-8;
  int status = S21_MATRIX_OK;

  i = 0;
  j = 0;
  while (i < n && j < n) {
    r = 0.0;
    for (k = i; k < n; ++k) {
      if (s21_fabs(a[k][j]) > r) {
        l = k;
        r = s21_fabs(a[k][j]);
      }
    }
    if (r <= eps) {
      status = S21_MATRIX_CALCULATION_ERROR;
      break;
    }
    if (i != l) {
      for (k = j; k < n; ++k) {
        r = a[i][k];
        a[i][k] = a[l][k];
        a[l][k] = r;
      }
      for (k = 0; k < n; ++k) {
        r = b[i][k];
        b[i][k] = b[l][k];
        b[l][k] = r;
      }
    }
    r = a[i][j];
    for (k = i + 1; k < n; ++k) {
      c = -a[k][j] / r;
      a[k][j] = 0.0;
      for (l = j + 1; l < n; ++l) {
        a[k][l] += c * a[i][l];
      }
      for (l = 0; l < n; ++l) {
        b[k][l] += c * b[i][l];
      }
    }
    ++i;
    ++j;
  }

  return (status);
}

static void s21_gj_method_to_unit_matrix(int n, double **a, double **b) {
  int i, j, k, l;
  double r;

  i = n - 1;
  j = n - 1;
  while (i > 0 && j > 0) {
    r = a[i][j];
    for (k = i - 1; k >= 0; --k) {
      double c = -a[k][j] / r;
      a[k][j] = 0.0;
      for (l = 0; l < n; ++l) {
        b[k][l] += c * b[i][l];
      }
    }
    --i;
    --j;
  }

  i = 0;
  while (i < n) {
    r = a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < n; ++j) {
      b[i][j] /= r;
    }
    ++i;
  }
}
