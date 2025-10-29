
// функции результат работы которых зависит от варианта задания

// структура реализующая точку плоскости
struct Point {
  double x;
  double y;
} typedef Point;

double calc_l_ij(Point P1, Point P2);

double calc_p_ij(Point P1, Point P2);

double calc_S_ij(double h1, double h2, Point p);

void test_calc_l_ij();
void test_calc_p_ij();
void test_calc_S_ij();