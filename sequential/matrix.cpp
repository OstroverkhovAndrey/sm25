#include <matrix.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// вычисление скалярного произведения
double scalar_product(double *u, double *v, double h1, double h2, int M,
                      int N) {
  double ans = 0.0;
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
    ans += h1 * h2 * u[j * M + i] * v[j * M + i];
    }
  }
  return ans;
}

// вычисление нормы L2
double norm(double *u, double h1, double h2, int M, int N) {
  return sqrt(scalar_product(u, u, h1, h2, M, N));
}

// копирование матрицы
void mat_copy(double *src, int M, int N, double *target) {
  int i, j;
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    target[j * M + i] = src[j * M + i];
  }
}

// обмен матриц указателями
void mat_swap(double **src, int M, int N, double **target) {
  double* tmp = *src;
  *src = *target;
  *target = tmp;
}

// установка во все значения матрицы значения val
void mat_set_value(double *u, int M, int N, int val) {
  int i, j;
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    u[j * M + i] = val;
  }
}

// создние матрицы и инициализация всех ее элементов нулем
double *mat_create(int M, int N) {
  double *mat = (double *)calloc(N * M, sizeof(double));
  mat_set_value(mat, M, N, 0);
  return mat;
}

// освобождение памяти выделенной под матрицу
void mat_free(double *mat, int M, int N) { free(mat); }

// сложени матицы u и матрицы v поэлементно
void mat_plus(double *u, double *v, int M, int N, double *ans) {
  int i, j;
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    ans[j * M + i] = u[j * M + i] + v[j * M + i];
  }
  return;
}

// вычитание из матицы u матрицы v поэлементно
void mat_minus(double *u, double *v, int M, int N, double *ans) {
  int i, j;
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    ans[j * M + i] = u[j * M + i] - v[j * M + i];
  }
  return;
}

// умножение всех элементов матрицы на число val
void mat_mul_number(double *mat, double val, int M, int N) {
  int i, j;
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    mat[j * M + i] *= val;
  }
  return;
}

// вывд матрицы в точностью 5 знаков после запятой
void mat_print(double *mat, int M, int N, Args args) {
  std::ostringstream oss;
  int i, j;
  oss << "[";
  for (int k = 0; k < M * N; ++k) {
    i = k % M;
    j = k / M;
    if (i == 0) {
      oss << "[";
    }
    oss << std::fixed << std::setprecision(args.precision) << mat[j * M + i]
        << ", ";
    if (i == M - 1) {
      oss << "],\n";
    }
  }
  oss << "]\n\n";

  std::string filename = "./txt/w.txt";
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
  }
  file << oss.str();
  file.close();
}
