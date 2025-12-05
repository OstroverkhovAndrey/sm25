
#include <iomanip>
#include <iostream>
#include <math.h>
#include <mpi.h>
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
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      double s1 = ((a[j * (M + 2) + (i + 1)] *
                    (w[j * (M + 2) + (i + 1)] - w[j * (M + 2) + i]) / h1) -
                   (a[j * (M + 2) + i] *
                    (w[j * (M + 2) + i] - w[j * (M + 2) + (i - 1)]) / h1)) /
                  h1;
      double s2 = ((b[(j + 1) * (M + 2) + i] *
                    (w[(j + 1) * (M + 2) + i] - w[j * (M + 2) + i]) / h2) -
                   (b[j * (M + 2) + i] *
                    (w[j * (M + 2) + i] - w[(j - 1) * (M + 2) + i]) / h2)) /
                  h2;
      ans[j * (M + 2) + i] = -s1 - s2;
    }
  }
  return;
}

// оператор D
void D_fun(double *a, double *b, double *w, int M, int N, double h1, double h2,
           double *ans) {
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      double d = ((a[j * (M + 2) + (i + 1)] + a[j * (M + 2) + i]) / (h1 * h1)) +
                 ((b[(j + 1) * (M + 2) + i] + b[j * (M + 2) + i]) / (h2 * h2));
      ans[j * (M + 2) + i] = w[j * (M + 2) + i] / d;
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
  double *a = mat_create(M + 2, N + 2);
  // тут считаем от нуля, чтобы корректно считалось в mpi версии
  // на последовательный и OpenMP код это не влияет
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < N + 2; ++j) {
    for (int i = 0; i < M + 2; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 + args.A1_field,
      p1.y = j * h2 - 0.5 * h2 + args.A2_field;
      Point p2;
      p2.x = i * h1 - 0.5 * h1 + args.A1_field,
      p2.y = j * h2 + 0.5 * h2 + args.A2_field;
      a[j * (M + 2) + i] = calc_a_ij(h2, p1, p2, eps);
    }
  }
  return a;
}

// инициализируем и заполняем матрицу b
double *init_b(int M, int N, double h1, double h2, double eps, Args args) {
  double *b = mat_create(M + 2, N + 2);
  // тут считаем от нуля, чтобы корректно считалось в mpi версии
  // на последовательный и OpenMP код это не влияет
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < N + 2; ++j) {
    for (int i = 0; i < M + 2; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 + args.A1_field,
      p1.y = j * h2 - 0.5 * h2 + args.A2_field;
      Point p2;
      p2.x = i * h1 + 0.5 * h1 + args.A1_field,
      p2.y = j * h2 - 0.5 * h2 + args.A2_field;
      b[j * (M + 2) + i] = calc_b_ij(h1, p1, p2, eps);
    }
  }
  return b;
}

// инициализируем и заполняем матрицу F
double *init_F(int M, int N, double h1, double h2, double eps, Args args) {
  double *F = mat_create(M + 2, N + 2);
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      Point p;
      p.x = i * h1 + args.A1_field, p.y = j * h2 + args.A2_field;
      F[j * (M + 2) + i] = calc_F_ij(h1, h2, p);
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

void send_border(double *v, Args args) {

  MPI_Status st;
  int n = args.N_field;
  int m = args.M_field;

  {
    // слева напрао
    double *send_v = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; ++i) {
      send_v[i] = v[(i + 1) * (m + 2) + m];
    }
    MPI_Sendrecv_replace(send_v, n, MPI_DOUBLE, args.right, 100, args.left, 100,
                         args.comm2d, &st);
    if (args.coords[0] != 0) {
      for (int i = 0; i < n; ++i) {
        v[(i + 1) * (m + 2) + 0] = send_v[i];
      }
    }
  }
  {
    // справа на лево
    double *send_v = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; ++i) {
      send_v[i] = v[(i + 1) * (m + 2) + 1];
    }
    MPI_Sendrecv_replace(send_v, n, MPI_DOUBLE, args.left, 200, args.right, 200,
                         args.comm2d, &st);
    if (args.coords[0] != args.dims[0] - 1) {
      for (int i = 0; i < n; ++i) {
        v[(i + 1) * (m + 2) + m + 1] = send_v[i];
      }
    }
  }
  {
    // сверху вниз
    double *send_v = (double *)malloc((size_t)m * sizeof(double));
    for (int i = 0; i < m; ++i) {
      send_v[i] = v[(n) * (m + 2) + i + 1];
    }
    int rc = MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE, args.down, 300,
                                  args.up, 300, args.comm2d, &st);

    if (args.coords[1] != 0) {
      for (int i = 0; i < m; ++i) {
        v[(0) * (m + 2) + i + 1] = send_v[i];
      }
    }
  }
  {
    // снизу вверх
    double *send_v = (double *)malloc((size_t)m * sizeof(double));
    for (int i = 0; i < m; ++i) {
      send_v[i] = v[(1) * (m + 2) + i + 1];
    }
    MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE, args.up, 400, args.down, 400,
                         args.comm2d, &st);
    if (args.coords[1] != args.dims[1] - 1) {
      for (int i = 0; i < m; ++i) {
        v[(n + 1) * (m + 2) + i + 1] = send_v[i];
      }
    }
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  Args args = parse_args(argc, argv);
  if (args.with_test) {
    // проверяем, что функции из variant правильно работают
    test_calc_l_ij();
    test_calc_p_ij();
    test_calc_S_ij();
  }

  double start_time, end_time;
  start_time = MPI_Wtime();

  // алиасы
  const int &N = args.N_field, M = args.M_field;
  const double &h1 = args.h1;
  const double &h2 = args.h2;
  const double &eps = args.eps;

  double alpha_k;

  double *a, *b, *F;
  a = init_a(args.M_field, args.N_field, h1, h2, eps, args);
  b = init_b(args.M_field, args.N_field, h1, h2, eps, args);
  F = init_F(args.M_field, args.N_field, h1, h2, eps, args);

  double *w_k = mat_create(args.M_field + 2, args.N_field + 2);
  double *w_k_plus1 = mat_create(args.M_field + 2, args.N_field + 2);

  double *temp1 = mat_create(args.M_field + 2, args.N_field + 2);

  double *z_k_0 = mat_create(args.M_field + 2, args.N_field + 2);
  double *p_k_1 = mat_create(args.M_field + 2, args.N_field + 2);
  double *r_k_0 = mat_create(args.M_field + 2, args.N_field + 2);
  double *z_k_1 = mat_create(args.M_field + 2, args.N_field + 2);
  double *p_k_2 = mat_create(args.M_field + 2, args.N_field + 2);
  double *r_k_1 = mat_create(args.M_field + 2, args.N_field + 2);

  {
    if (args.rank2d == 0) {
      std::cout << "Zero iteration\n\n";
    }
    // r_k_0 = B - A w_k
    send_border(w_k, args);
    A_fun(a, b, w_k, args.M_field, args.N_field, h1, h2, temp1);
    mat_minus(F, temp1, args.M_field + 2, args.N_field + 2, r_k_0);

    // z_k = r_k / D;
    D_fun(a, b, r_k_0, args.M_field, args.N_field, h1, h2, z_k_0);

    // p_k_1 = z_k_0;
    mat_copy(z_k_0, args.M_field + 2, args.N_field + 2, p_k_1);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    send_border(p_k_1, args);
    A_fun(a, b, p_k_1, args.M_field, args.N_field, h1, h2, temp1);
    alpha_k = scalar_product(z_k_0, r_k_0, h1, h2, args.M_field + 2,
                             args.N_field + 2, args) /
              scalar_product(temp1, p_k_1, h1, h2, args.M_field + 2,
                             args.N_field + 2, args);

    // w_1 == w_0 + alpha_k * p_1
    mat_copy(p_k_1, args.M_field + 2, args.N_field + 2, temp1);
    mat_mul_number(temp1, alpha_k, args.M_field + 2, args.N_field + 2);
    mat_plus(w_k, temp1, args.M_field + 2, args.N_field + 2, w_k_plus1);

    // w_0 = w_1, для зацикливания
    mat_swap(&w_k_plus1, args.M_field + 2, args.N_field + 2, &w_k);
  }

  if (args.rank2d == 0) {
    std::cout << "Start iteration\n\n";
  }
  int i = 0;
  double err = args.delta + 1;
  for (; stop_iter(i, err, args); ++i) {

    // r_k_1 == r_k_0 - alpha_k A p_k_1
    send_border(p_k_1, args);
    A_fun(a, b, p_k_1, M, N, h1, h2, temp1);
    mat_mul_number(temp1, alpha_k, M + 2, N + 2);
    mat_minus(r_k_0, temp1, M + 2, N + 2, r_k_1);

    // z_k = r_k / D;
    D_fun(a, b, r_k_1, M, N, h1, h2, z_k_1);

    // betta == (z_k_1, r_k_1) / (z_k_0, r_k_0)
    double betta = scalar_product(z_k_1, r_k_1, h1, h2, M + 2, N + 2, args) /
                   scalar_product(z_k_0, r_k_0, h1, h2, M + 2, N + 2, args);

    // p_k_2 = z_k + betta_k * p_k_1; // тут меняется p_k_1 но он больше не
    // используется
    mat_mul_number(p_k_1, betta, M + 2, N + 2);
    mat_plus(z_k_1, p_k_1, M + 2, N + 2, p_k_2);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    send_border(p_k_2, args);
    A_fun(a, b, p_k_2, M, N, h1, h2, temp1);
    alpha_k = scalar_product(z_k_1, r_k_1, h1, h2, M + 2, N + 2, args) /
              scalar_product(temp1, p_k_2, h1, h2, M + 2, N + 2, args);

    // w_1 == w_0 + alpha_k * p_k_2
    mat_copy(p_k_2, M + 2, N + 2, temp1);
    mat_mul_number(temp1, alpha_k, M + 2, N + 2);
    mat_plus(w_k, temp1, M + 2, N + 2, w_k_plus1);

    // считаем ошибку
    mat_minus(w_k_plus1, w_k, M + 2, N + 2, temp1);
    err = norm(temp1, h1, h2, M + 2, N + 2, args);
    if (args.rank2d == 0) {
      std::cout << std::fixed << std::setprecision(args.precision)
                << "error: " << err << std::endl;
    }

    // копируем значения для зацикливания
    mat_swap(&w_k_plus1, M + 2, N + 2, &w_k);
    mat_swap(&z_k_1, M + 2, N + 2, &z_k_0);
    mat_swap(&p_k_2, M + 2, N + 2, &p_k_1);
    mat_swap(&r_k_1, M + 2, N + 2, &r_k_0);
  }

  if (args.rank2d == 0) {
    std::cout << "\nCount iteration: " << i << "\n";
  }
  if (args.print_result) {
    mat_print(w_k, M + 2, N + 2, args);
  }

  // очищаем память
  mat_free(a, args.M_field + 2, args.N_field + 2);
  mat_free(b, args.M_field + 2, args.N_field + 2);
  mat_free(F, args.M_field + 2, args.N_field + 2);

  mat_free(temp1, args.M_field + 2, args.N_field + 2);

  mat_free(w_k, args.M_field + 2, args.N_field + 2);
  mat_free(w_k_plus1, args.M_field + 2, args.N_field + 2);
  mat_free(z_k_0, args.M_field + 2, args.N_field + 2);
  mat_free(p_k_1, args.M_field + 2, args.N_field + 2);
  mat_free(r_k_0, args.M_field + 2, args.N_field + 2);
  mat_free(z_k_1, args.M_field + 2, args.N_field + 2);
  mat_free(p_k_2, args.M_field + 2, args.N_field + 2);
  mat_free(r_k_1, args.M_field + 2, args.N_field + 2);

  // Замеряем время
  end_time = MPI_Wtime();
  double local_time = end_time - start_time, max_time;
  MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (args.rank2d == 0) {
    std::cout << "Максимальное время среди процессов: " << max_time
              << " секунд\n";
  }
  MPI_Comm_free(&args.comm2d);
  MPI_Finalize();
  return 0;
}