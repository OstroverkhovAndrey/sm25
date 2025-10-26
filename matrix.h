
double scalar_product(double *u, double *v, double h1, double h2, int M,
                      int N);

double norm(double *u, double h1, double h2, int M, int N);

void mat_copy(double *src, int M, int N, double *target);

void mat_set_value(double *u, int M, int N, int val);

double *mat_create(int M, int N);

void mat_free(double *mat, int M, int N);

void mat_plus(double *u, double *v, int M, int N, double *ans);

void mat_minus(double *u, double *v, int M, int N, double *ans);

void mat_mul_number(double *mat, double val, int M, int N);

void mat_print(double *mat, int M, int N);