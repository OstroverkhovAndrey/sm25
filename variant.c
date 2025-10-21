
#include <variant.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// вычисление l_ij
double calc_l_ij(Point P1, Point P2) {
  if (P1.x != P2.x) {
    fprintf(stderr, "Error in calc_l_ij, P1.x != P2.x\n");
    return 0;
  }
  if (P1.y >= P2.y) {
    fprintf(stderr, "Error in calc_l_ij, P1.y >= P2.y\n");
    return 0;
  }
  // отрезок P1 P2 лежит вне треугольника
  if (P2.y <= 0 || P1.x <= -3 || P1.x >= 3 ||
      P1.y >= ((2.0 / 3.0) * P1.x + 2) || P1.y >= ((-2.0 / 3.0) * P1.x + 2)) {
    return 0.0;
  }
  // отрезок P1 P2 лежит полностью в треугольнике
  if (P1.y >= 0 && P2.y <= ((2.0 / 3.0) * P2.x + 2) &&
      P2.y <= ((-2.0 / 3.0) * P2.x + 2)) {
    return P2.y - P1.y;
  }
  // отрезов P1 P2 пересекает треугольник в одной точке внизу
  if (P1.y < 0 && P2.y < ((2.0 / 3.0) * P2.x + 2) &&
      P2.y < ((-2.0 / 3.0) * P2.x + 2)) {
    return P2.y - 0.0;
  }
  // отрезов P1 P2 пересекает треугольник в одной точке вверху
  if (P1.y >= 0 && P1.x <= 0 && (P2.y > ((2.0 / 3.0) * P2.x + 2))) {
      return ((2.0 / 3.0) * P1.x + 2) - P1.y;
  }
  if (P1.y >= 0 && P1.x >= 0 &&  P2.y > ((-2.0 / 3.0) * P2.x + 2)) {
      return ((-2.0 / 3.0) * P1.x + 2) - P1.y;
  }

  // отрезов P1 P2 пересекает треугольник в двух точках
  if (P1.y < 0 && P1.x <= 0 && P2.y > ((2.0 / 3.0) * P2.x + 2)) {
      return ((2.0 / 3.0) * P1.x + 2) - 0.0;
  }
  if (P1.y < 0 && P1.x >= 0 && P2.y > ((-2.0 / 3.0) * P2.x + 2)) {
      return ((-2.0 / 3.0) * P1.x + 2) - 0.0;
  }

  fprintf(stderr, "Error in calc_l_ij, default return\n");
  return 0;
}

// вычисление p_ij
double calc_p_ij(Point P1, Point P2) {
  if (P1.x >= P2.x) {
    fprintf(stderr, "Error in calc_p_ij, P1.x >= P2.x\n");
    return 0;
  }
  if (P1.y != P2.y) {
    fprintf(stderr, "Error in calc_p_ij, P1.y != P2.y\n");
    return 0;
  }
  // отрезок P1 P2 лежит вне треугольника
  if (P1.x >= 3 || P2.x <= -3 || P1.y <= 0 || P1.y >= 2 ||
      P1.y >= ((-2.0 / 3.0) * P1.x + 2.0) ||
      P2.y >= ((2.0 / 3.0) * P2.x + 2.0)) {
    return 0.0;
  }
  // отрезок P1 P2 лежит полностью в треугольнике
  if (P1.y > 0 && P1.y <= ((2.0 / 3.0) * P1.x + 2) &&
      P2.y <= ((-2.0 / 3.0) * P2.x + 2)) {
    return P2.x - P1.x;
  }
  // отрезов P1 P2 пересекает треугольник в одной точке справа
  if (P2.y >= ((-2.0 / 3.0) * P2.x + 2) && P1.y <= ((2.0 / 3.0) * P1.x + 2)) {
    return ((-3.0 / 2.0) * P1.y + 3) - P1.x;
  }
  // отрезов P1 P2 пересекает треугольник в одной точке слева
  if (P1.y >= ((2.0 / 3.0) * P1.x + 2) && P2.y <= ((-2.0 / 3.0) * P2.x + 2)) {
    return P2.x - ((3.0 / 2.0) * P1.y - 3);
  }
  // отрезов P1 P2 пересекает треугольник в двух точках
  if (P2.y >= ((-2.0 / 3.0) * P2.x + 2) && P1.y >= ((2.0 / 3.0) * P1.x + 2)) {
    return ((-3.0 / 2.0) * P1.y + 3) - ((3.0 / 2.0) * P1.y - 3);
  }
  fprintf(stderr, "Error in calc_p_ij, default return\n");
  return 0;
}

struct TestCalclpConfig {
  Point P1;
  Point P2;
  double ans;
} typedef TestCalclpConfig;

void test_calc_l_ij() {
  double eps = 0.0001;
  int count_tests = 14;
  TestCalclpConfig tests_config[14] = {
      {{0, -1}, {0, 0}, 0}, // 0 отрезок лежит вне треугольника
      {{-4, -1}, {-4, 1}, 0},    // 1
      {{4, -2}, {4, 1}, 0},      // 2
      {{-2, 1}, {-2, 2}, 0},     // 3
      {{1, 2}, {1, 3}, 0},       // 4
      {{-1.5, 0}, {-1.5, 1}, 1}, // 5 полностью внутри
      {{1, 0.5}, {1, 1}, 0.5},   // 6 полностью внутри
      {{0, -1}, {0, 1}, 1}, // 7 одна точка пересечения внизу
      {{-1, 1}, {-1, 2}, 1.0 / 3.0}, // 8 одна точка пересечения вверху
      {{0, 1}, {0, 2}, 1}, // 9 одна точка пересечения вверху
      {{2, 0}, {2, 1}, 2.0 / 3.0}, // 10 одна точка пересечения вверху
      {{-2, -1}, {-2, 3}, 2.0 / 3.0}, // 11 две точки пересечения слева
      {{0, -1}, {0, 3}, 2}, // 12 две точки пересечения по центру
      {{-1, -1}, {-1, 3}, 4.0 / 3.0}, // 13 две точки пересечения справа
  };
  for (int i = 0; i < count_tests; ++i) {
    TestCalclpConfig config = tests_config[i];
    double ans = calc_l_ij(config.P1, config.P2);
    double err = ans - config.ans;
    if (-eps < err && err < eps) {
      printf("Pass test %i for calc_l_ij\n", i);
    } else {
      printf("Faild test %i for calc_l_ij\n", i);
    }
  }
  printf("\n");
}

void test_calc_p_ij() {
  double eps = 0.0001;
  int count_tests = 12;
  TestCalclpConfig tests_config[12] = {
      {{-1, -1}, {1, -1}, 0}, // 0 отрезок лежит вне треугольника
      {{-4, 1}, {-3, 1}, 0}, // 1
      {{3, 2}, {4, 2}, 0},   // 2
      {{-3, 1}, {-2, 1}, 0}, // 3
      {{1, 2}, {2, 2}, 0},   // 4
      {{0, 0.5}, {0.7, 0.5}, 0.7}, // 5 отрезок лежит внутри треугольника
      {{-1.5, 1}, {-1, 1}, 0.5}, // 6
      {{-3, 1}, {-1, 1}, 0.5}, // 7 отрезок пересекает треугольник в одной точке
      {{1, 0.5}, {3, 0.5}, 5.0 / 4.0}, // 8
      {{-1, 2}, {1, 2}, 0}, // 9 отрезок пересекает треугольник в точке (0, 2)
      {{-3, 1}, {3, 1}, 3}, // 10 отрезок пересекает треугольник в двух точках
      {{-3, 0}, {3, 0}, 0}, // 11 отрезок совпадает с нижней гранью треугольника
  };
  for (int i = 0; i < count_tests; ++i) {
    TestCalclpConfig config = tests_config[i];
    double ans = calc_p_ij(config.P1, config.P2);
    double err = ans - config.ans;
    if (-eps < err && err < eps) {
      printf("Pass test %i for calc_p_ij\n", i);
    } else {
      printf("Faild test %i for calc_p_ij\n", i);
    }
  }
  printf("\n");
}

// считает площадь методом шнурков
double shoelace_formula(Point points[], int n) {
  if (n < 3) {
    return 0.0;
  }
  double s1 = 0.0, s2 = 0.0;
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    s1 += points[i].x * points[j].y;
    s2 += points[i].y * points[j].x;
  }
  double ans = fabs(s1 - s2) / 2;
  return ans;
}

// вычисление площади пересечения справа от оси Y
double calc_S_ij_right(double h1, double h2, Point p) {
  Point p1, p2, p3, p4;
  p1.x = p.x - h1 / 2; // низ лево
  p1.y = p.y - h2 / 2;
  p2.x = p.x - h1 / 2; // верх лево
  p2.y = p.y + h2 / 2;
  p3.x = p.x + h1 / 2; // верх право
  p3.y = p.y + h2 / 2;
  p4.x = p.x + h1 / 2; // низ право
  p4.y = p.y - h2 / 2;
  double x1 = p1.x, y1 = p1.y;
  double x2 = p3.x, y2 = p3.y;

  // базовые проверки на пересечение
  if (x2 <= 0 || y2 <= 0) {
    return 0.0;
  }
  if (y1 >= ((-2.0 / 3.0) * x1 + 2)) {
    return 0.0;
  }
  // исключаем область в которой точно нет треугольника
  if (x1 < 0) {
    x1 = 0;
  }
  if (y1 < 0) {
    y1 = 0;
  }
  // собираем все точки по часовой стрелки
  int n = 0;
  Point points[8];
  if (y1 < ((-2.0 / 3.0) * x1 + 2)) {
    points[n] = Point{x1, y1};
    ++n;
  }
  if (y1 <= ((-2.0 / 3.0) * x1 + 2) && ((-2.0 / 3.0) * x1 + 2) <= y2) {
    points[n] = Point{x1, ((-2.0 / 3.0) * x1 + 2)};
    ++n;
  }
  if (y2 < ((-2.0 / 3.0) * x1 + 2)) {
    points[n] = Point{x1, y2};
    ++n;
  }
  if (x1 <= ((-3.0 / 2.0) * y2 + 3) && ((-3.0 / 2.0) * y2 + 3) <= x2) {
    points[n] = Point{((-3.0 / 2.0) * y2 + 3), y2};
    ++n;
  }
  if (y2 < ((-2.0 / 3.0) * x2 + 2)) {
    points[n] = Point{x2, y2};
    ++n;
  }
  if (y1 <= ((-2.0 / 3.0) * x2 + 2) && ((-2.0 / 3.0) * x2 + 2) <= y2) {
    points[n] = Point{x2, ((-2.0 / 3.0) * x2 + 2)};
    ++n;
  }
  if (y1 < ((-2.0 / 3.0) * x2 + 2)) {
    points[n] = Point{x2, y1};
    ++n;
  }
  if (x1 <= ((-3.0 / 2.0) * y1 + 3) && ((-3.0 / 2.0) * y1 + 3) <= x2) {
    points[n] = Point{((-3.0 / 2.0) * y1 + 3), y1};
    ++n;
  }
  double S1 = shoelace_formula(points, n);
  return S1;
}

// вычисление площади пересечения слева от оси Y
double calc_S_ij_left(double h1, double h2, Point p) {
  Point p1, p2, p3, p4;
  p1.x = p.x - h1 / 2; // низ лево
  p1.y = p.y - h2 / 2;
  p2.x = p.x - h1 / 2; // верх лево
  p2.y = p.y + h2 / 2;
  p3.x = p.x + h1 / 2; // верх право
  p3.y = p.y + h2 / 2;
  p4.x = p.x + h1 / 2; // низ право
  p4.y = p.y - h2 / 2;
  double x1 = p4.x, y1 = p4.y;
  double x2 = p2.x, y2 = p2.y;

  // базовые проверки на пересечение
  if (x2 >= 0 || y2 <= 0) {
    return 0.0;
  }
  if (y1 > ((2.0 / 3.0) * x1 + 2.0)) {
    return 0.0;
  }
  // исключаем область в которой точно нет треугольника
  if (x1 >= 0) {
    x1 = 0.0;
  }
  if (y1 <= 0) {
    y1 = 0.0;
  }
  // собираем все точки по часовой стрелки
  int n = 0;
  Point points[8];
  if (y1 < ((2.0 / 3.0) * x1 + 2.0)) {
    points[n] = Point{x1, y1};
    ++n;
  }
  if (x2 <= ((3.0 / 2.0) * y1 - 3.0) && ((3.0 / 2.0) * y1 - 3.0) <= x1) {
    points[n] = Point{((3.0 / 2.0) * y1 - 3.0), y1};
    ++n;
  }
  if (y1 < ((2.0 / 3.0) * x2 + 2.0)) {
    points[n] = Point{x2, y1};
    ++n;
  }
  if (y1 <= ((2.0 / 3.0) * x2 + 2.0) && ((2.0 / 3.0) * x2 + 2.0) <= y2) {
    points[n] = Point{x2, ((2.0 / 3.0) * x2 + 2.0)};
    ++n;
  }
  if (y2 < ((2.0 / 3.0) * x2 + 2.0)) {
    points[n] = Point{x2, y2};
    ++n;
  }
  if (x2 <= ((3.0 / 2.0) * y2 - 3.0) && ((3.0 / 2.0) * y2 - 3.0) <= x1) {
    points[n] = Point{((3.0 / 2.0) * y2 - 3.0), y2};
    ++n;
  }
  if (y2 < ((2.0 / 3.0) * x1 + 2.0)) {
    points[n] = Point{x1, y2};
    ++n;
  }
  if (y1 <= ((2.0 / 3.0) * x1 + 2.0) && ((2.0 / 3.0) * x1 + 2.0) <= y2) {
    points[n] = Point{x1, ((2.0 / 3.0) * x1 + 2.0)};
    ++n;
  }
  double S2 = shoelace_formula(points, n);
  return S2;
}

// вычисление S_ij
double calc_S_ij(double h1, double h2, Point p) {
  double S1 = calc_S_ij_right(h1, h2, p);
  double S2 = calc_S_ij_left(h1, h2, p);
  return S1 + S2;
}

struct TestCalcSConfig {
  double h1;
  double h2;
  Point P;
  double ans;
} typedef TestCalcSConfig;

void test_calc_S_ij() {
  double eps = 0.0001;
  int count_tests = 15;
  TestCalcSConfig tests_config[15] = {
      {1, 1, {0.5, 0.5}, 1}, // 0 квадраты в углу внутри треугольника
      {1, 1, {-0.5, 0.5}, 1}, // 1
      {1, 1, {0.0, -2.0}, 0}, // 2 вне треугольника
      {1, 1, {4.0, 3.0}, 0}, // 3
      {1, 1, {-4.0, 3.0}, 0}, // 4
      {2, 2, {0.0, 2.0}, 4.0/3.0}, // 5 на вершине треугольника
      {2, 2, {3.0, 0.0}, 1.0/3.0}, // 6 на правой вершине треугольника
      {2, 2, {-3.0, 0.0}, 1.0/3.0}, // 7 на левой вершине треугольника
      {2, 2, {1.5, 1.0}, 2.0}, // 8 на точке (1,5; 1)
      {2, 2-2.0/3.0, {1.5, 1.0}, 4.0/3.0}, // 9 на точке (1,5; 1)
      {2, 1, {1.5, 1.0}, 1.0}, // 10 на точке (1,5; 1)
      {2, 2, {-1.5, 1.0}, 2.0}, // 11 на точке (-1,5; 1)
      {2, 2-2.0/3.0, {-1.5, 1.0}, 4.0/3.0}, // 12 на точке (-1,5; 1)
      {2, 1, {-1.5, 1.0}, 1.0}, // 12 на точке (-1,5; 1)

  };
  for (int i = 0; i < count_tests; ++i) {
    TestCalcSConfig config = tests_config[i];
    double ans = calc_S_ij(config.h1, config.h2, config.P);
    double err = ans - config.ans;
    if (-eps < err && err < eps) {
      printf("Pass test %i for calc_S_ij\n", i);
    } else {
      printf("Faild test %i for calc_S_ij\n", i);
    }
  }
}
