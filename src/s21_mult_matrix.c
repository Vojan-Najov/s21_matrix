#include "s21_matrix.h"

static double _mult_aux(int i, int j, matrix_t *A, matrix_t *B) {
  double *a_row = A->matrix[i];
  double **b = B->matrix;
  int tmp = A->columns;
  double sum = 0.0;

  for (int k = 0; k < tmp; ++k) {
    sum += a_row[k] * b[k][j];
  }

  return (sum);
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rows, columns;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A) || MATRIX_WRONG(B)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else if (MATRIX_COMPATIBILITY(A, B)) {
    rows = A->rows;
    columns = B->columns;
    status = s21_create_matrix(rows, columns, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
          result->matrix[i][j] = _mult_aux(i, j, A, B);
        }
      }
    }
  } else {
    status = S21_MATRIX_CALCULATION_ERROR;
  }
  return (status);
}
