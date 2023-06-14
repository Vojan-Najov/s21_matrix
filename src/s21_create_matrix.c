#include <stdlib.h>

#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = S21_MATRIX_OK;
  char *tmp;

  if (rows > 0 && columns > 0) {
    result->matrix = (double **)calloc(
        1, sizeof(double *) * rows + sizeof(double) * columns * rows);
    if (result->matrix != NULL) {
      result->rows = rows;
      result->columns = columns;
      tmp = (char *)result->matrix;
      tmp += sizeof(double *) * rows;
      for (int i = 0; i < rows; ++i) {
        result->matrix[i] = (double *)tmp;
        tmp += sizeof(double) * columns;
      }
    } else {
      result->rows = 0;
      result->columns = 0;
    }
  } else {
    MATRIX_RESET(result);
    status = S21_MATRIX_ERROR;
  }

  return (status);
}
