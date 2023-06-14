#include "s21_matrix.h"

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rows, columns;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A) || MATRIX_WRONG(B)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else if (!MATRIX_ORDER_EQUAL(A, B)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_CALCULATION_ERROR;
  } else {
    rows = A->rows;
    columns = A->columns;
    status = s21_create_matrix(rows, columns, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  }

  return (status);
}
