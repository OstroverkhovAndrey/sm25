
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <args.hpp>
#include <matrix.hpp>
#include <variant.hpp>

// оператор A
void A_fun(double *a, double *b, double *w, int M, int N, double h1, double h2,
           double *ans) {
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {

      double s1 = ((a[j * (M + 1) + (i + 1)] *
                    (w[j * (M + 1) + (i + 1)] - w[j * (M + 1) + i]) / h1) -
                   (a[j * (M + 1) + i] *
                    (w[j * (M + 1) + i] - w[j * (M + 1) + (i - 1)]) / h1)) /
                  h1;
      double s2 = ((b[(j + 1) * (M + 1) + i] *
                    (w[(j + 1) * (M + 1) + i] - w[j * (M + 1) + i]) / h2) -
                   (b[j * (M + 1) + i] *
                    (w[j * (M + 1) + i] - w[(j - 1) * (M + 1) + i]) / h2)) /
                  h2;
      ans[j * (M + 1) + i] = -s1 - s2;
    }
  }
  return;
}

// оператор D
void D_fun(double *a, double *b, double *w, int M, int N, double h1, double h2,
           double *ans) {
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      double d = ((a[j * (M + 1) + (i + 1)] + a[j * (M + 1) + i]) / (h1 * h1)) +
                 ((b[(j + 1) * (M + 1) + i] + b[j * (M + 1) + i]) / (h2 * h2));
      ans[j * (M + 1) + i] = w[j * (M + 1) + i] / d;
    }
  }
  return;
}

// вычисление a_ij
double calc_a_ij(double h2, Point p1, Point p2, double eps) {
  double l_ij = calc_l_ij(p1, p2);
  double a_ij = (l_ij / h2) + (1 - l_ij / h2) / eps;
  return a_ij;
}

// вычисление b_ij
double calc_b_ij(double h1, Point p1, Point p2, double eps) {
  double p_ij = calc_p_ij(p1, p2);
  double b_ij = (p_ij / h1) + (1 - p_ij / h1) / eps;
  return b_ij;
}

// вычисление F_ij
double calc_F_ij(double h1, double h2, Point p) {
  double S = calc_S_ij(h1, h2, p);
  double F_ij = S / (h1 * h2); // f(x_i, y_i) == 1;
  return F_ij;
}

// инициализируем и заполняем матрицу a
double *init_a(int M, int N, double h1, double h2, double eps, Args args) {
  double *a = mat_create(M + 1, N + 1);
// тут считаем от нуля, чтобы корректно считалось в mpi версии
// на последовательный и OpenMP код это не влияет
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < N + 1; ++j) {
    for (int i = 0; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 + args.A1, p1.y = j * h2 - 0.5 * h2 + args.A2;
      Point p2;
      p2.x = i * h1 - 0.5 * h1 + args.A1, p2.y = j * h2 + 0.5 * h2 + args.A2;
      a[j * (M + 1) + i] = calc_a_ij(h2, p1, p2, eps);
    }
  }
  return a;
}

// инициализируем и заполняем матрицу b
double *init_b(int M, int N, double h1, double h2, double eps, Args args) {
  double *b = mat_create(M + 1, N + 1);
// тут считаем от нуля, чтобы корректно считалось в mpi версии
// на последовательный и OpenMP код это не влияет
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < N + 1; ++j) {
    for (int i = 0; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 + args.A1, p1.y = j * h2 - 0.5 * h2 + args.A2;
      Point p2;
      p2.x = i * h1 + 0.5 * h1 + args.A1, p2.y = j * h2 - 0.5 * h2 + args.A2;
      b[j * (M + 1) + i] = calc_b_ij(h1, p1, p2, eps);
    }
  }
  return b;
}

// инициализируем и заполняем матрицу F
double *init_F(int M, int N, double h1, double h2, double eps, Args args) {
  double *F = mat_create(M + 1, N + 1);
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      Point p;
      p.x = i * h1 + args.A1, p.y = j * h2 + args.A2;
      F[j * (M + 1) + i] = calc_F_ij(h1, h2, p);
    }
  }
  return F;
}

// условие остановки итерационного процесса
bool stop_iter(int i, double err, Args args) {
  if (args.count_iter == -1) {
    return err > args.delta;
  } else {
    return i < args.count_iter;
  }
  std::cout << "Error! Default return in stop iter func" << std::endl;
  return false;
}

int main(int argc, char *argv[]) {
  const Args args = parse_args(argc, argv);

  if (args.with_test) {
    // проверяем, что функции из variant правильно работают
    test_calc_l_ij();
    test_calc_p_ij();
    test_calc_S_ij();
  }

  double start_time, end_time;
  start_time = omp_get_wtime();

  // алиасы
  const int &N = args.N, M = args.M;
  const double &h1 = args.h1;
  const double &h2 = args.h2;
  const double &eps = args.eps;

  double alpha_k;

  double *a, *b, *F;
  a = init_a(M, N, h1, h2, eps, args);
  b = init_b(M, N, h1, h2, eps, args);
  F = init_F(M, N, h1, h2, eps, args);

  double *w_k = mat_create(M + 1, N + 1);
  double *w_k_plus1 = mat_create(M + 1, N + 1);

  double *temp1 = mat_create(M + 1, N + 1);

  double *z_k_0 = mat_create(M + 1, N + 1);
  double *p_k_1 = mat_create(M + 1, N + 1);
  double *r_k_0 = mat_create(M + 1, N + 1);
  double *z_k_1 = mat_create(M + 1, N + 1);
  double *p_k_2 = mat_create(M + 1, N + 1);
  double *r_k_1 = mat_create(M + 1, N + 1);

  {
    std::cout << "Zero iteration\n\n";
    // r_k_0 = B - A w_k
    A_fun(a, b, w_k, M, N, h1, h2, temp1);
    mat_minus(F, temp1, M + 1, N + 1, r_k_0);

    // z_k = r_k / D;
    D_fun(a, b, r_k_0, M, N, h1, h2, z_k_0);
    // p_k_1 = z_k_0;
    mat_copy(z_k_0, M + 1, N + 1, p_k_1);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    A_fun(a, b, p_k_1, M, N, h1, h2, temp1);
    alpha_k = scalar_product(z_k_0, r_k_0, h1, h2, M + 1, N + 1) /
              scalar_product(temp1, p_k_1, h1, h2, M + 1, N + 1);

    // w_1 == w_0 + alpha_k * p_1
    mat_copy(p_k_1, M + 1, N + 1, temp1);
    mat_mul_number(temp1, alpha_k, M + 1, N + 1);
    mat_plus(w_k, temp1, M + 1, N + 1, w_k_plus1);

    // w_0 = w_1, для зацикливания
    mat_swap(&w_k_plus1, M + 1, N + 1, &w_k);
  }

  int i = 0;
  double err = args.delta + 1;
  std::cout << "Start iteration\n\n";
  for (; stop_iter(i, err, args); ++i) {

    // r_k_1 == r_k_0 - alpha_k A p_k_1
    A_fun(a, b, p_k_1, M, N, h1, h2, temp1);
    mat_mul_number(temp1, alpha_k, M + 1, N + 1);
    mat_minus(r_k_0, temp1, M + 1, N + 1, r_k_1);

    // z_k = r_k / D;
    D_fun(a, b, r_k_1, M, N, h1, h2, z_k_1);

    // betta == (z_k_1, r_k_1) / (z_k_0, z_k_0)
    double betta = scalar_product(z_k_1, r_k_1, h1, h2, M + 1, N + 1) /
                   scalar_product(z_k_0, r_k_0, h1, h2, M + 1, N + 1);

    // p_k_2 = z_k + betta_k * p_k_1; // тут меняется p_k_1 но он больше не
    // используется
    mat_mul_number(p_k_1, betta, M + 1, N + 1);
    mat_plus(z_k_1, p_k_1, M + 1, N + 1, p_k_2);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    A_fun(a, b, p_k_2, M, N, h1, h2, temp1);
    alpha_k = scalar_product(z_k_1, r_k_1, h1, h2, M + 1, N + 1) /
              scalar_product(temp1, p_k_2, h1, h2, M + 1, N + 1);

    // w_1 == w_0 + alpha_k * p_k_2
    mat_copy(p_k_2, M + 1, N + 1, temp1);
    mat_mul_number(temp1, alpha_k, M + 1, N + 1);
    mat_plus(w_k, temp1, M + 1, N + 1, w_k_plus1);

    // считаем ошибку
    mat_minus(w_k_plus1, w_k, M + 1, N + 1, temp1);
    err = norm(temp1, h1, h2, M + 1, N + 1);
    std::cout << std::fixed << std::setprecision(args.precision)
              << "error: " << err << std::endl;

    // меняем местами значения для зацикливания
    mat_swap(&w_k_plus1, M + 1, N + 1, &w_k);
    mat_swap(&z_k_1, M + 1, N + 1, &z_k_0);
    mat_swap(&p_k_2, M + 1, N + 1, &p_k_1);
    mat_swap(&r_k_1, M + 1, N + 1, &r_k_0);
  }

  std::cout << "\nCount iteration: " << i << std::endl;
  if (args.print_result) {
    mat_print(w_k, M + 1, N + 1, args);
  }

  // очищаем память
  mat_free(a, M + 1, N + 1);
  mat_free(b, M + 1, N + 1);
  mat_free(F, M + 1, N + 1);

  mat_free(temp1, M + 1, N + 1);

  mat_free(w_k, M + 1, N + 1);
  mat_free(w_k_plus1, M + 1, N + 1);

  mat_free(z_k_0, M + 1, N + 1);
  mat_free(p_k_1, M + 1, N + 1);
  mat_free(r_k_0, M + 1, N + 1);
  mat_free(z_k_1, M + 1, N + 1);
  mat_free(p_k_2, M + 1, N + 1);
  mat_free(r_k_1, M + 1, N + 1);

  end_time = omp_get_wtime();
  std::cout << std::fixed << std::setprecision(args.precision) << "Work took "
            << end_time - start_time << " seconds\n";
  return 0;
}