
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
  // int i, j;
  //  for (int k = 0; k < (N - 0) * (M - 0); ++k) {
  //    i = (k % (M - 0)) + 1;
  //    j = (k / (M - 0)) + 1;
  //    double s1 = ((a[j * (M + 0) + (i + 1)] *
  //                  (w[j * (M + 0) + (i + 1)] - w[j * (M + 0) + i]) / h1) -
  //                 (a[j * (M + 0) + i] *
  //                  (w[j * (M + 0) + i] - w[j * (M + 0) + (i - 1)]) / h1)) /
  //                h1;
  //    double s2 = ((b[(j + 1) * (M + 0) + i] *
  //                  (w[(j + 1) * (M + 0) + i] - w[j * (M + 0) + i]) / h2) -
  //                 (b[j * (M + 0) + i] *
  //                  (w[j * (M + 0) + i] - w[(j - 1) * (M + 0) + i]) / h2)) /
  //                h2;
  //    ans[j * (M + 0) + i] = -s1 - s2;
  //  }
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

// оператор B
void B_fun(double *F, int M, int N, double *ans) {
  // int i, j;
  // for (int k = 0; k < (N - 0) * (M - 0); ++k) {
  //   i = (k % (M - 0)) + 1;
  //   j = (k / (M - 0)) + 1;
  //   ans[j * (M + 2) + i] = F[j * (M + 2) + i];
  // }
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
     ans[j * (M + 2) + i] = F[j * (M + 2) + i];
    }
  }
  return;
}

// оператор D
void D_fun(double *a, double *b, double *w, int M, int N, double h1, double h2,
           double *ans) {
  //int i, j;
  //for (int k = 0; k < (N - 0) * (M - 0); ++k) {
  //  i = (k % (M - 0)) + 1;
  //  j = (k / (M - 0)) + 1;
  //  double d = ((a[j * (M + 2) + (i + 1)] + a[j * (M + 2) + i]) / (h1 * h1)) +
  //             ((b[(j + 1) * (M + 2) + i] + b[j * (M + 2) + i]) / (h2 * h2));
  //  ans[j * (M + 2) + i] = w[j * (M + 2) + i] / d;
  //}
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
  int i, j;
  for (int k = 0; k < (N + 2) * (M + 2); ++k) {
    i = (k % (M + 2));
    j = (k / (M + 2));
    Point p1;
    p1.x = i * h1 - 0.5 * h1 + args.A1_field,
    p1.y = j * h2 - 0.5 * h2 + args.A2_field;
    Point p2;
    p2.x = i * h1 - 0.5 * h1 + args.A1_field,
    p2.y = j * h2 + 0.5 * h2 + args.A2_field;
    a[j * (M + 2) + i] = calc_a_ij(h2, p1, p2, eps);
  }
  // if (args.coords[0] == 0 && args.coords[1] == 0) {
  //   for (int i = 1; i < args.M_field+1; ++i) {
  //     Point p1;
  //     int j = N+1;
  //     p1.x = i * h1 - 0.5 * h1 + args.A1_field, p1.y = j * h2 - 0.5 * h2 +
  //     args.A2_field; Point p2; p2.x = i * h1 - 0.5 * h1 + args.A1_field, p2.y
  //     = j * h2 + 0.5 * h2 + args.A2_field; a[(N+1)*(M+2)+i] = calc_a_ij(h2,
  //     p1, p2, eps);
  //   }
  // }
  // if (args.coords[0] == 0 && args.coords[1] == 1) {
  //   for (int i = 1; i < args.M_field+1; ++i) {
  //     Point p1;
  //     int j = 0;
  //     p1.x = i * h1 - 0.5 * h1 + args.A1_field, p1.y = j * h2 - 0.5 * h2 +
  //     args.A2_field; Point p2; p2.x = i * h1 - 0.5 * h1 + args.A1_field, p2.y
  //     = j * h2 + 0.5 * h2 + args.A2_field; a[i] = calc_a_ij(h2, p1, p2, eps);
  //   }
  // }
  return a;
}

// инициализируем и заполняем матрицу b
double *init_b(int M, int N, double h1, double h2, double eps, Args args) {
  double *b = mat_create(M + 2, N + 2);
  int i, j;
  for (int k = 0; k < (N + 2) * (M + 2); ++k) {
    i = (k % (M + 2));
    j = (k / (M + 2));
    Point p1;
    p1.x = i * h1 - 0.5 * h1 + args.A1_field,
    p1.y = j * h2 - 0.5 * h2 + args.A2_field;
    Point p2;
    p2.x = i * h1 + 0.5 * h1 + args.A1_field,
    p2.y = j * h2 - 0.5 * h2 + args.A2_field;
    b[j * (M + 2) + i] = calc_b_ij(h1, p1, p2, eps);
  }
  return b;
}

// инициализируем и заполняем матрицу F
double *init_F(int M, int N, double h1, double h2, double eps, Args args) {
  double *F = mat_create(M + 2, N + 2);
  int i, j;
  for (int k = 0; k < (N + 2) * (M + 2); ++k) {
    i = (k % (M + 2));
    j = (k / (M + 2));
    Point p;
    p.x = i * h1 + args.A1_field, p.y = j * h2 + args.A2_field;
    F[j * (M + 2) + i] = calc_F_ij(h1, h2, p);
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

void send_v_f(double *v, Args args) {

  MPI_Status st;

  int n = args.N_field;
  int m = args.M_field;

  //{
  //    // слева напрао
  //    double* send_v = (double*)malloc((size_t)n * sizeof(double));
  //    for (int i = 0; i < n; ++i) {
  //        send_v[i] = v[(i+1)*(m+2) + n-1 +1];
  //    }
  //    MPI_Sendrecv_replace(send_v, n, MPI_DOUBLE,
  //                         args.right, 100,  // dest, tag
  //                         args.left,  100,  // source, tag
  //                         args.comm2d, &st);
  //    if (args.coords[0] != 0) {
  //        for (int i = 0; i < n; ++i) {
  //            v[(i+1)*(m+2) + 0] = send_v[i];
  //        }
  //    }
  //}
  //{
  //    // справа на лево
  //    double* send_v = (double*)malloc((size_t)n * sizeof(double));
  //    for (int i = 0; i < n; ++i) {
  //        send_v[i] = v[(i+1)*(m+2) + 1];
  //    }
  //    MPI_Sendrecv_replace(send_v, n, MPI_DOUBLE,
  //                         args.left, 200,  // dest, tag
  //                         args.right,  200,  // source, tag
  //                         args.comm2d, &st);
  //    if (args.coords[0] != args.dims[0]-1) {
  //        for (int i = 0; i < n; ++i) {
  //            v[(i+1)*(m+2) + n+1] = send_v[i];
  //        }
  //    }
  //}
  {
    // сверху вниз
    double *send_v = (double *)malloc((size_t)m * sizeof(double));
    for (int i = 0; i < m; ++i) {
      send_v[i] = v[(n) * (m + 2) + i + 1];
    }
    int rc =
        MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE, args.down, 300, // dest, tag
                             args.up, 300, // source, tag
                             args.comm2d, &st);

    //if (rc != MPI_SUCCESS) {
    //  char err_string[MPI_MAX_ERROR_STRING];
    //  int err_length;
    //  MPI_Error_string(rc, err_string, &err_length);
    //  fprintf(stderr, "Процесс %d: Ошибка передачи/приема: %s\n", args.rank2d,
    //          err_string);
    //  MPI_Abort(MPI_COMM_WORLD, rc);
    //} else {
    //  // Можно проанализировать статус, например, узнать ранг отправителя
    //  printf("Процесс %d: получил значение = %f от процесса %d (тег %d)\n",
    //         args.rank2d, send_v[0], st.MPI_SOURCE, st.MPI_TAG);
    //}

    //if (args.coords[1] != 0) {
      for (int i = 0; i < m; ++i) {
        v[(0) * (m + 2) + i + 1] = send_v[i];
      }
    //}
    //for (int r = 0; r < args.world_size; ++r) {
    //  if (r == args.rank2d) {
    //    printf("[rank2d=%2d coords=(%d,%d) dims=(%d,%d)] "
    //           "L=%2d R=%2d U=%2d D=%2d\n",
    //           args.rank2d, args.coords[0], args.coords[1], args.dims[0],
    //           args.dims[1], args.left, args.right, args.up, args.down);
    //    fflush(stdout);
    //    for (int i = 0; i < m; ++i) {
    //      std::cout << send_v[i] << " ";
    //    }
    //    std::cout << std::endl;
    //    std::cout << std::flush;
    //  }
    //  MPI_Barrier(args.comm2d);
    //}
  }
  {
    // снизу вверх
    double *send_v = (double *)malloc((size_t)m * sizeof(double));
    for (int i = 0; i < m; ++i) {
      send_v[i] = v[(1) * (m + 2) + i + 1];
    }
    MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE, args.up, 400, // dest, tag
                         args.down, 400,                      // source, tag
                         args.comm2d, &st);
    //if (args.coords[1] != args.dims[1] - 1) {
      for (int i = 0; i < m; ++i) {
        v[(n + 1) * (m + 2) + i + 1] = send_v[i];
      }
    //}
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  const Args args = parse_args(argc, argv);
  // if (args.with_test) {
  //   // проверяем, что функции из variant правильно работают
  //   test_calc_l_ij();
  //   test_calc_p_ij();
  //   test_calc_S_ij();
  // }

  double start_time, end_time;
  start_time = omp_get_wtime();

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
  double *temp2 = mat_create(args.M_field + 2, args.N_field + 2);
  double *temp3 = mat_create(args.M_field + 2, args.N_field + 2);
  double *temp4 = mat_create(args.M_field + 2, args.N_field + 2);

  double *z_k_0 = mat_create(args.M_field + 2, args.N_field + 2);
  double *p_k_1 = mat_create(args.M_field + 2, args.N_field + 2);
  double *r_k_0 = mat_create(args.M_field + 2, args.N_field + 2);
  double *z_k_1 = mat_create(args.M_field + 2, args.N_field + 2);
  double *p_k_2 = mat_create(args.M_field + 2, args.N_field + 2);
  double *r_k_1 = mat_create(args.M_field + 2, args.N_field + 2);

  {
    // std::cout << "Zero iteration\n\n";
    // r_k_0 = B - A w_k
    A_fun(a, b, w_k, args.M_field, args.N_field, h1, h2, temp1);
    B_fun(F, args.M_field, args.N_field, temp2);
    mat_minus(temp2, temp1, args.M_field + 2, args.N_field + 2, r_k_0);

    // z_k = r_k / D;
    D_fun(a, b, r_k_0, args.M_field, args.N_field, h1, h2, z_k_0);

    // p_k_1 = z_k_0;
    mat_copy(z_k_0, args.M_field + 2, args.N_field + 2, p_k_1);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    A_fun(a, b, p_k_1, args.M_field, args.N_field, h1, h2, temp3);
    alpha_k = scalar_product(z_k_0, r_k_0, h1, h2, args.M_field + 2,
                             args.N_field + 2, args) /
              scalar_product(temp3, p_k_1, h1, h2, args.M_field + 2,
                             args.N_field + 2, args);

    // w_1 == w_0 + alpha_k * p_1
    mat_copy(p_k_1, args.M_field + 2, args.N_field + 2, temp4);
    mat_mul_number(temp4, alpha_k, args.M_field + 2, args.N_field + 2);
    mat_plus(w_k, temp4, args.M_field + 2, args.N_field + 2, w_k_plus1);

    // w_0 = w_1, для зацикливания
    mat_copy(w_k_plus1, args.M_field + 2, args.N_field + 2, w_k);

    // обмениваемся данными с соседними процессами
    send_v_f(w_k, args);
  }
  // std::cout << "Start iteration\n\n";

  int i = 0;
  double err = args.delta + 1;
  for (; stop_iter(i, err, args); ++i) {

    // r_k_1 == r_k_0 - alpha_k A p_k_1
    A_fun(a, b, p_k_1, M, N, h1, h2, temp1);
    mat_mul_number(temp1, alpha_k, M + 2, N + 2);
    mat_minus(r_k_0, temp1, M + 2, N + 2, r_k_1);

    // z_k = r_k / D;
    D_fun(a, b, r_k_1, M, N, h1, h2, z_k_1);

    // betta == (z_k_1, r_k_1) / (z_k_0, z_k_0)
    double betta = scalar_product(z_k_1, r_k_1, h1, h2, M + 2, N + 2, args) /
                   scalar_product(z_k_0, r_k_0, h1, h2, M + 2, N + 2, args);

    // p_k_2 = z_k + betta_k * p_k_1; // тут меняется p_k_1 но он больше не
    // используется
    mat_mul_number(p_k_1, betta, M + 2, N + 2);
    mat_plus(z_k_1, p_k_1, M + 2, N + 2, p_k_2);

    // alpha_k == (z_0, k_0) / (Ap_1, p_1)
    A_fun(a, b, p_k_2, M, N, h1, h2, temp2);
    alpha_k = scalar_product(z_k_1, r_k_1, h1, h2, M + 2, N + 2, args) /
              scalar_product(temp2, p_k_2, h1, h2, M + 2, N + 2, args);

    // w_1 == w_0 + alpha_k * p_k_2
    mat_copy(p_k_2, M + 2, N + 2, temp3);
    mat_mul_number(temp3, alpha_k, M + 2, N + 2);
    mat_plus(w_k, temp3, M + 2, N + 2, w_k_plus1);

    // считаем ошибку
    mat_copy(w_k, M + 2, N + 2, temp1);
    mat_copy(w_k_plus1, M + 2, N + 2, temp2);
    mat_minus(temp2, temp1, M + 2, N + 2, temp3);
    err = norm(temp3, h1, h2, M + 2, N + 2, args);

    // std::cout << std::fixed << std::setprecision(args.precision)
    //           << "error: " << err << std::endl;

    // копируем значения для зацикливания
    mat_copy(w_k_plus1, M + 2, N + 2, w_k);
    mat_copy(z_k_1, M + 2, N + 2, z_k_0);
    mat_copy(p_k_2, M + 2, N + 2, p_k_1);
    mat_copy(r_k_1, M + 2, N + 2, r_k_0);

    // обмениваемся данными с соседними процессами
    send_v_f(w_k, args);
  }

  // std::cout << "\nCount iteration: " << i << "\n" << "Result w:\n";
  fflush(stdout);
  for (int r = 0; r < args.world_size; ++r) {
    if (r == args.rank2d) {
      printf("[rank2d=%2d coords=(%d,%d) dims=(%d,%d)] "
             "L=%2d R=%2d U=%2d D=%2d\n",
             args.rank2d, args.coords[0], args.coords[1], args.dims[0],
             args.dims[1], args.left, args.right, args.up, args.down);

      // for (int j = 0; j < m+2; ++j) {
      //     for (int i = 0; i < n+2; ++i) {
      //         printf("%06.2f ", v[j][i]);
      //     }
      //     printf("\n");
      // }
      //// printf("recv v:\n");
      //// for (int i = 0; i < m; ++i) {
      ////     printf("%06.2f ", send_v[i]);
      //// }
      // printf("\n\n");
      mat_print(w_k, args.M_field + 2, args.N_field + 2, args);
      fflush(stdout);
    }
    MPI_Barrier(args.comm2d);
  }

  // очищаем память
  mat_free(a, args.M_field + 2, args.N_field + 2);
  mat_free(b, args.M_field + 2, args.N_field + 2);
  mat_free(F, args.M_field + 2, args.N_field + 2);

  mat_free(temp1, args.M_field + 2, args.N_field + 2);
  mat_free(temp2, args.M_field + 2, args.N_field + 2);
  mat_free(temp3, args.M_field + 2, args.N_field + 2);
  mat_free(temp4, args.M_field + 2, args.N_field + 2);

  mat_free(w_k, args.M_field + 2, args.N_field + 2);
  mat_free(w_k_plus1, args.M_field + 2, args.N_field + 2);
  mat_free(z_k_0, args.M_field + 2, args.N_field + 2);
  mat_free(p_k_1, args.M_field + 2, args.N_field + 2);
  mat_free(r_k_0, args.M_field + 2, args.N_field + 2);
  mat_free(z_k_1, args.M_field + 2, args.N_field + 2);
  mat_free(p_k_2, args.M_field + 2, args.N_field + 2);
  mat_free(r_k_1, args.M_field + 2, args.N_field + 2);

  MPI_Finalize();
  end_time = omp_get_wtime();
  // std::cout << std::fixed << std::setprecision(args.precision)
  //           << "Work took " << end_time - start_time << " seconds\n";
  return 0;
}