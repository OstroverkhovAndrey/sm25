#pragma once

struct Args {
    int N{10};
    int M{10};
    int count_iter{-1};
    double delta{0.01};
    bool with_test{false};
    double A1{-4.0};
    double B1{4.0};
    double A2{-1.0};
    double B2{3.0};
    double h1{0.0};
    double h2{0.0};
    double eps{0.0};
    int precision{5};
};

Args parse_args(int argc, char *argv[]);