#include <matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// вычисление скалярного произведения
double scalar_product(double **u, double **v, double h1, double h2, int M,
                      int N) {
  double ans = 0.0;
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
      ans += h1 * h2 * u[j][i] * v[j][i];
    }
  }
  return ans;
}

// вычисление нормы L2
double norm(double **u, double h1, double h2, int M, int N) {
  return sqrt(scalar_product(u, u, h1, h2, M, N));
}

// копирование матрицы
void mat_copy(double **src, double **target, int M, int N) {
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      target[j][i] = src[j][i];
    }
  }
}

// установка во все значения матрицы значения val
void mat_set_value(double **u, int M, int N, int val) {
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      u[j][i] = val;
    }
  }
}

// создние матрицы и инициализация всех ее элементов нулем
double **mat_create(int M, int N) {
  double **mat = (double **)calloc(N, sizeof(double *));
  for (int j = 0; j < N; ++j) {
    mat[j] = (double *)calloc(M, sizeof(double));
  }
  mat_set_value(mat, M, N, 0);
  return mat;
}

// освобождение памяти выделенной под матрицу
void mat_free(double **mat, int M, int N) {
  for (int j = 0; j < N; ++j) {
    free(mat[j]);
  }
  free(mat);
}

// сложени матицы u и матрицы v поэлементно
void mat_plus(double **u, double **v, int M, int N, double **ans) {
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
      ans[j][i] = u[j][i] + v[j][i];
    }
  }
  return;
}

// вычитание из матицы u матрицы v поэлементно
void mat_minus(double **u, double **v, int M, int N, double **ans) {
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
      ans[j][i] = u[j][i] - v[j][i];
    }
  }
  return;
}

// умножение всех элементов матрицы на число val
void mat_mul_number(double **mat, double val, int M, int N) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      mat[j][i] *= val;
    }
  }
  return;
}


// вывд матрицы в точностью 5 знаков после запятой
void mat_print(double **mat, int M, int N) {
  printf("[");
  for (int j = N - 1; j >= 0; --j) {
  printf("[");
    for (int i = 0; i < M; ++i) {
      printf("%.5f, ", mat[j][i]);
    }
    printf("],\n");
  }
  printf("]\n\n\n");
}
