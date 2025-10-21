
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <variant.h>
#include <matrix.h>

// максимум двую чисел
double max(double a, double b) {
  if (a > b)
    return a;
  return b;
}

// оператор A
void A_fun(double **a, double **b, double **w, int M, int N, double h1,
           double h2, double **ans) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      double s1 = ((a[j][i + 1] * (w[j][i + 1] - w[j][i]) / h1) -
                   (a[j][i] * (w[j][i] - w[j][i - 1]) / h1)) /
                  h1;
      double s2 = ((b[j + 1][i] * (w[j + 1][i] - w[j][i]) / h2) -
                   (b[j][i] * (w[j][i] - w[j - 1][i]) / h2)) /
                  h2;
      ans[j][i] = -s1 - s2;
    }
  }
  return;
}

// оператор B
void B_fun(double **F, int M, int N, double **ans) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      ans[j][i] = F[j][i];
    }
  }
  return;
}

// оператор D
void D_fun(double **a, double **b, double **w, int M, int N, double h1,
           double h2, double **ans) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      double d = ( (a[j][i+1]+a[j][i]) / (h1*h1) ) + ( (b[j+1][i]+b[j][i]) / (h2*h2) );
      ans[j][i] = w[j][i] / d;
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
double **init_a(int M, int N, double h1, double h2, double eps) {
  double **a = mat_create(M + 1, N + 1);
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 -4, p1.y = j * h2 - 0.5 * h2 -1;
      Point p2;
      p2.x = i * h1 - 0.5 * h1 -4, p2.y = j * h2 + 0.5 * h2 -1;
      a[j][i] = calc_a_ij(h2, p1, p2, eps);
    }
  }
  return a;
}

// инициализируем и заполняем матрицу b
double **init_b(int M, int N, double h1, double h2, double eps) {
  double **b = mat_create(M + 1, N + 1);
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1 -4, p1.y = j * h2 - 0.5 * h2 -1;
      Point p2;
      p2.x = i * h1 + 0.5 * h1 -4, p2.y = j * h2 - 0.5 * h2 -1;
      b[j][i] = calc_b_ij(h1, p1, p2, eps);
    }
  }
  return b;
}

// инициализируем и заполняем матрицу F
double **init_F(int M, int N, double h1, double h2, double eps) {
  double **F = mat_create(M, N);
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      Point p;
      p.x = i * h1 -4, p.y = j * h2 -1;
      F[j][i] = calc_F_ij(h1, h2, p);
    }
  }
  return F;
}


int main() {
  // запуск тестов
  if (1) {
    test_calc_l_ij();
    test_calc_p_ij();
    test_calc_S_ij();
  }
  double start, end;
  start = omp_get_wtime();

  const double A1 = -4.0, B1 = 4.0, A2 = -1.0, B2 = 3.0;
  const int N = 100, M = 100;
  const double h1 = (B1 - A1) / (1.0 * M);
  const double h2 = (B2 - A2) / (1.0 * N);
  const int count_iter = 100;
  const double delta = 0.00001;
  const double eps = max(h1, h2) * max(h1, h2);

  double **a, **b, **F;
  a = init_a(M, N, h1, h2, eps);
  b = init_b(M, N, h1, h2, eps);
  F = init_F(M, N, h1, h2, eps);
  //mat_print(F, M, N);

  double **w_k = mat_create(M + 1, N + 1);
  double **w_k_plus1 = mat_create(M + 1, N + 1);

  double **temp1 = mat_create(M+1, N+1);
  double **temp2 = mat_create(M+1, N+1);
  double **temp3 = mat_create(M+1, N+1);
  double **temp4 = mat_create(M+1, N+1);
  double **temp5 = mat_create(M+1, N+1);

  double **z_k_0 = mat_create(M + 1, N + 1);
  double **p_k_1 = mat_create(M + 1, N + 1);
  double **r_k_0 = mat_create(M + 1, N + 1);
  double **z_k_1 = mat_create(M + 1, N + 1);
  double **p_k_2 = mat_create(M + 1, N + 1);
  double **r_k_1 = mat_create(M + 1, N + 1);
  double alpha_k;

  {
    printf("Zero iteration\n\n");
    A_fun(a, b, w_k, M, N, h1, h2, temp1);
    B_fun(F, M, N, temp2);
    mat_minus(temp2, temp1, M + 1, N + 1, r_k_0);

    // z_k = r_k / D;
    D_fun(a, b, r_k_0, M, N, h1, h2, z_k_0);

    //p_k_1 = z_k_0;
    mat_copy(z_k_0, M+1, N+1, p_k_1);

    A_fun(a, b, p_k_1, M, N, h1, h2, temp3);
    alpha_k = scalar_product(z_k_0, r_k_0, h1, h2, M, N) /
                   scalar_product(temp3, p_k_1, h1, h2, M, N);
    printf("alpha: %f\n", alpha_k);
    
    mat_copy(p_k_1, M+1, N+1, temp4);
    mat_mul_number(temp4, alpha_k, M, N);
    
    mat_plus(w_k, temp4, M + 1, N + 1, w_k_plus1);
    mat_copy(w_k_plus1, M + 1, N + 1, w_k);
  }
  printf("Start iteration\n\n");

  int i = 0;
  double err = delta + 1;
  //for (; err > delta; ++i) {
  for (; i < count_iter; ++i) {

    A_fun(a, b, p_k_1, M, N, h1, h2, temp1);
    mat_mul_number(temp1, alpha_k, M+1, N+1);
    mat_minus(r_k_0, temp1, M + 1, N + 1, r_k_1);

    // z_k = r_k / D;
    D_fun(a, b, r_k_1, M, N, h1, h2, z_k_1);

    double betta = scalar_product(z_k_1, r_k_1, h1, h2, M+1, N+1) /
                   scalar_product(z_k_0, r_k_0, h1, h2, M+1, N+1);
    //printf("betta: %f\n", betta);
    //p_k = z_k + betta_k * p_k; // тут меняется p_k_1 но он больше не используется
    mat_mul_number(p_k_1, betta, M+1, N+1);
    mat_plus(z_k_1, p_k_1, M+1, N+1, p_k_2);


    A_fun(a, b, p_k_2, M, N, h1, h2, temp2);
    alpha_k = scalar_product(z_k_1, r_k_1, h1, h2, M, N) /
                   scalar_product(temp2, p_k_2, h1, h2, M, N);
    //printf("alpha: %f\n", alpha_k);
    
    mat_copy(p_k_2, M+1, N+1, temp3);
    mat_mul_number(temp3, alpha_k, M, N);
    
    mat_plus(w_k, temp3, M + 1, N + 1, w_k_plus1);

    {
      mat_copy(w_k, M+1, N+1, temp1);
      mat_copy(w_k_plus1, M+1, N+1, temp2);
      mat_minus(temp2, temp1, M+1, N+1, temp3);
      double err = norm(temp3, h1, h2, M+1, N+1);
      printf("error: %f\n", err);
    }

    mat_copy(w_k_plus1, M + 1, N + 1, w_k);
    mat_copy(z_k_1, M + 1, N + 1, z_k_0);
    mat_copy(p_k_2, M + 1, N + 1, p_k_1);
    mat_copy(r_k_1, M + 1, N + 1, r_k_0);
  }

  printf("\n");
  printf("Count iteration: %i\n", i);
  printf("Result w:\n");
  mat_print(w_k, M+1, N+1);

  mat_free(a, M + 1, N + 1);
  mat_free(b, M + 1, N + 1);
  mat_free(F, M, N);
  mat_free(temp1, M, N);
  mat_free(temp2, M, N);
  mat_free(temp3, M, N);
  mat_free(w_k, M + 1, N + 1);
  mat_free(w_k_plus1, M + 1, N + 1);
  //mat_free(r_k, M + 1, N + 1);

  end = omp_get_wtime();
  printf("Work took %f seconds\n", end - start);
  return 0;
}