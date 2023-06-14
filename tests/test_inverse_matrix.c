#include <check.h>
#include <stdio.h>

#include "s21_matrix.h"

START_TEST(test_01) {
  matrix_t A, B, C;
  double a[3][3] = {{2., 5., 7.}, {6., 3., 4.}, {5., -2., -3.}};
  double b[3][3] = {{1., -1., 1.}, {-38., 41., -34.}, {27., -29., 24.}};
  int ret;

  ret = s21_create_matrix(3, 3, &A);
  ret += s21_create_matrix(3, 3, &B);
  if (ret == 0) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        A.matrix[i][j] = a[i][j];
        B.matrix[i][j] = b[i][j];
      }
    }
    ret = s21_inverse_matrix(&A, &C);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(C.rows, 3);
    ck_assert_int_eq(C.columns, 3);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        ck_assert_double_eq_tol(C.matrix[i][j], B.matrix[i][j], 1.0e-7);
      }
    }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
  }
}
END_TEST

START_TEST(test_02) {
  matrix_t A, B;
  double a[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {5.0, 7.0, 9.0}};
  int ret;

  ret = s21_create_matrix(3, 3, &A);
  if (ret == 0) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        A.matrix[i][j] = a[i][j];
      }
    }
    ret = s21_inverse_matrix(&A, &B);
    ck_assert_int_eq(ret, 2);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
  }
}
END_TEST

START_TEST(test_03) {
  matrix_t A, B;
  int ret;

  ret = s21_create_matrix(1, 1, &A);
  if (ret == 0) {
    A.matrix[0][0] = 125480.4;
    ret = s21_inverse_matrix(&A, &B);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(B.rows, 1);
    ck_assert_int_eq(B.columns, 1);
    ck_assert_double_eq_tol(5.0 / 627402.0, B.matrix[0][0], 1.0e-7);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
  }
}
END_TEST

START_TEST(test_04) {
  double a[6][6] = {
      {12., 47., 58., 47., 42., 14.},   {47., 59., 63., 54., 89., 12.},
      {15., 56., -65., -97., 32., 16.}, {58., 75., -24., 45., -16., 8.},
      {78., 93., 17., 13., -11., 7.},   {84., 65., 19., -35., 54., 18.}};
  double b[6][6] = {{-0.01, -0.001, -0.009, 0.006, -0.005, 0.016},
                    {0.0, 0.008, 0.01, -0.007, 0.018, -0.018},
                    {0.009, -0.005, -0.004, -0.014, 0.011, 0.002},
                    {0.001, 0.003, -0.004, 0.012, -0.009, -0.001},
                    {-0.011, 0.017, 0.005, -0.001, -0.004, -0.005},
                    {0.075, -0.067, -0.012, 0.04, -0.061, 0.058}};

  matrix_t A, B;
  int ret;

  ret = s21_create_matrix(6, 6, &A);
  if (ret == 0) {
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        A.matrix[i][j] = a[i][j];
      }
    }
    ret = s21_inverse_matrix(&A, &B);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(B.rows, 6);
    ck_assert_int_eq(B.columns, 6);
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        ck_assert_double_eq_tol(B.matrix[i][j], b[i][j], 1.0e-3);
      }
    }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
  }
}
END_TEST

START_TEST(error_test_01) {
  int ret;
  matrix_t A, B;

  A.rows = 20;
  A.columns = 199;
  A.matrix = NULL;
  ret = s21_inverse_matrix(&A, &B);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(error_test_02) {
  int ret;
  matrix_t A, B;
  double tmp = 0.0;
  double *tmpp = &tmp;
  double **tmppp = &tmpp;

  A.rows = 1;
  A.columns = -199;
  A.matrix = tmppp;
  ret = s21_inverse_matrix(&A, &B);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(error_test_03) {
  int ret;
  matrix_t A, B;
  double tmp = 0.0;
  double *tmpp = &tmp;
  double **tmppp = &tmpp;

  A.rows = -5;
  A.columns = 199;
  A.matrix = tmppp;
  ret = s21_inverse_matrix(&A, &B);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(error_test_04) {
  int ret;
  matrix_t A, B;
  double tmp = 0.0;
  double *tmpp = &tmp;
  double **tmppp = &tmpp;

  A.rows = 5;
  A.columns = 199;
  A.matrix = tmppp;
  ret = s21_inverse_matrix(&A, &B);
  ck_assert_int_eq(ret, 2);
}
END_TEST

Suite *test_inverse_matrix(void) {
  Suite *s;
  TCase *tc;

  s = suite_create("s21_inverse_matrix");
  tc = tcase_create("s21_inverse_matrix");

  if (s != NULL && tc != NULL) {
    tcase_add_test(tc, test_01);
    tcase_add_test(tc, test_02);
    tcase_add_test(tc, test_03);
    tcase_add_test(tc, test_04);
    tcase_add_test(tc, error_test_01);
    tcase_add_test(tc, error_test_02);
    tcase_add_test(tc, error_test_03);
    tcase_add_test(tc, error_test_04);
    suite_add_tcase(s, tc);
  }

  return (s);
}
