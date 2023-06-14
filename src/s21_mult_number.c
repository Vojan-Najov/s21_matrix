#include "s21_matrix.h"

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int rows, columns;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else {
    rows = A->rows;
    columns = A->columns;
    status = s21_create_matrix(rows, columns, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
          result->matrix[i][j] = number * A->matrix[i][j];
        }
      }
    }
  }

  return (status);
}
