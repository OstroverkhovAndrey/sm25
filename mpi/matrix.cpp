#include <matrix.hpp>

#include <iomanip>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#include <args.hpp>

// вычисление скалярного произведения
double scalar_product(double *u, double *v, double h1, double h2, int M, int N,
                      Args args) {
  double ans = 0.0;
  int i, j;
  for (int k = 0; k < (M-2) * (N-2); ++k) {
    i = k % (M - 2) + 1;
    j = k / (M - 2) + 1;
    ans += h1 * h2 * u[j * M + i] * v[j * M + i];
  }

  double global_ans = 0.0;
  MPI_Allreduce(&ans, &global_ans, 1, MPI_DOUBLE, MPI_SUM, args.comm2d);
  return global_ans;
}

// вычисление нормы L2
double norm(double *u, double h1, double h2, int M, int N, Args args) {
  return sqrt(scalar_product(u, u, h1, h2, M, N, args));
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
  oss << "coords=(" << args.coords[0] << "," << args.coords[1] << ")\n";

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

  std::string filename = "./txt/" + std::to_string(args.world_size) + "_" + std::to_string(args.rank2d) + ".txt";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
    }
    file << oss.str();
    file.close();
}
