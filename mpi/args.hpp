#pragma once

#include <mpi.h>

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

  int N_field{-1};
  int M_field{-1};
  int dims[2]{0, 0};
  int world_rank{-1};
  int world_size{-1};
  MPI_Comm comm2d;
  int rank2d;
  int coords[2];

  // соседи
  int left;
  int right;
  int up;
  int down;

  // откуда начинается часть решетки для этого mpi процесса
  double A1_field{-4.0};
  double A2_field{-1.0};
};

Args parse_args(int argc, char *argv[]);