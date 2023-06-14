#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <stddef.h>

#define S21_MATRIX_OK 0
#define S21_MATRIX_ERROR 1
#define S21_MATRIX_CALCULATION_ERROR 2
#define SUCCESS 1
#define FAILURE 0

#define MATRIX_WRONG(m_ptr) \
  ((m_ptr)->matrix == NULL || (m_ptr)->rows < 1 || (m_ptr)->columns < 1)

#define MATRIX_ORDER_EQUAL(m1_ptr, m2_ptr) \
  ((m1_ptr)->rows == (m2_ptr)->rows && (m1_ptr)->columns == (m2_ptr)->columns)

#define MATRIX_COMPATIBILITY(ml_ptr, mr_ptr) \
  ((ml_ptr)->columns == (mr_ptr)->rows)

#define MATRIX_ELEMENTS_EQUAL(a, b) (-1.0e-7 < (a) - (b) && (a) - (b) < 1.0e-7)

#define MATRIX_RESET(m_ptr) \
  (m_ptr)->matrix = NULL;   \
  (m_ptr)->rows = 0;        \
  (m_ptr)->columns = 0;

#define MATRIX_SQUARE(m_ptr) ((m_ptr)->rows == (m_ptr)->columns)

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

/* auxiliary functions */
double s21_fabs(double x);
int s21_gauss_method(matrix_t *A);
void s21_print_matrix(int r, int c, double **a);

#endif
