#include "s21_matrix.h"

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  double **a, **b;
  int rows, columns;

  if (!MATRIX_WRONG(A) && !MATRIX_WRONG(B) && MATRIX_ORDER_EQUAL(A, B)) {
    a = A->matrix;
    b = B->matrix;
    rows = A->rows;
    columns = A->columns;
    for (int i = 0; status == SUCCESS && i < rows; ++i) {
      for (int j = 0; status == SUCCESS && j < columns; ++j) {
        if (!MATRIX_ELEMENTS_EQUAL(a[i][j], b[i][j])) {
          status = FAILURE;
        }
      }
    }
  } else {
    status = FAILURE;
  }

  return (status);
}
