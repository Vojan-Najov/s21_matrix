#include "s21_matrix.h"

int s21_transpose(matrix_t *A, matrix_t *result) {
  int rows, columns;
  int status = S21_MATRIX_OK;

  if (MATRIX_WRONG(A)) {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  } else {
    rows = A->columns;
    columns = A->rows;
    status = s21_create_matrix(rows, columns, result);
    if (status == S21_MATRIX_OK) {
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
          result->matrix[i][j] = A->matrix[j][i];
        }
      }
    }
  }

  return (status);
}
