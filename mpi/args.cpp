#include <args.hpp>

#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>

double square_of_max(double a, double b) {
  double m = std::max(a, b);
  return m * m;
}

Args parse_args(int argc, char *argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "N") {
      args.N = std::stoi(argv[i + 1]);
      ++i;
      continue;
    }
    if (std::string(argv[i]) == "M") {
      args.M = std::stoi(argv[i + 1]);
      ++i;
      continue;
    }
    if (std::string(argv[i]) == "count_iter") {
      args.count_iter = std::stoi(argv[i + 1]);
      ++i;
      continue;
    }
    if (std::string(argv[i]) == "delta") {
      args.delta = std::stof(argv[i + 1]);
      ++i;
      continue;
    }
    if (std::string(argv[i]) == "--test") {
      args.with_test = true;
      continue;
    }
  }
  args.h1 = (args.B1 - args.A1) / (1.0 * args.M);
  args.h2 = (args.B2 - args.A2) / (1.0 * args.N);
  args.eps = square_of_max(args.h1, args.h2);

  // mpi  args
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  args.world_rank = world_rank;
  args.world_size = world_size;

  // определяем размер решетки и размер задачи в каждом mpi процессе
  if (args.M == 4 && args.N == 6) {
    if (world_size == 2) {
      args.M_field = 4;
      args.N_field = 3;
      args.dims[0] = 1;
      args.dims[1] = 2;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else if (args.M == 8 && args.N == 12) {
    if (world_size == 2) {
      args.M_field = 8;
      args.N_field = 6;
      args.dims[0] = 1;
      args.dims[1] = 2;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else if (args.M == 40 && args.N == 40) {
    if (world_size == 1) {
      args.M_field = 40;
      args.N_field = 40;
      args.dims[0] = 1;
      args.dims[1] = 1;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else if (args.M == 40 && args.N == 60) {
    if (world_size == 1) {
      args.M_field = 40;
      args.N_field = 60;
      args.dims[0] = 1;
      args.dims[1] = 1;
    } else if (world_size == 2) {
      args.M_field = 40;
      args.N_field = 30;
      args.dims[0] = 1;
      args.dims[1] = 2;
    } else if (world_size == 4) {
      args.M_field = 20;
      args.N_field = 30;
      args.dims[0] = 2;
      args.dims[1] = 2;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else if (args.M == 400 && args.N == 600) {
    if (world_size == 2) {
      args.M_field = 400;
      args.N_field = 300;
      args.dims[0] = 1;
      args.dims[1] = 2;
    } else if (world_size == 4) {
      args.M_field = 200;
      args.N_field = 300;
      args.dims[0] = 2;
      args.dims[1] = 2;
    } else if (world_size == 8) {
      args.M_field = 200;
      args.N_field = 150;
      args.dims[0] = 2;
      args.dims[1] = 4;
    } else if (world_size == 16) {
      args.M_field = 100;
      args.N_field = 150;
      args.dims[0] = 4;
      args.dims[1] = 4;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else if (args.M == 800 && args.N == 1200) {
    if (world_size == 4) {
      args.M_field = 400;
      args.N_field = 300;
      args.dims[0] = 2;
      args.dims[1] = 2;
    } else if (world_size == 8) {
      args.M_field = 200;
      args.N_field = 300;
      args.dims[0] = 2;
      args.dims[1] = 4;
    } else if (world_size == 16) {
      args.M_field = 200;
      args.N_field = 150;
      args.dims[0] = 4;
      args.dims[1] = 4;
    } else if (world_size == 32) {
      args.M_field = 100;
      args.N_field = 150;
      args.dims[0] = 4;
      args.dims[1] = 8;
    } else {
      std::cout << "Error parse M N world_size" << std::endl;
    }
  } else {
    std::cout << "Error parse M N world_size" << std::endl;
  }

  int periods[2] = {0, 0};
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, args.dims, periods, reorder, &args.comm2d);

  MPI_Comm_rank(args.comm2d, &args.rank2d);
  MPI_Cart_coords(args.comm2d, args.rank2d, 2, args.coords);

  MPI_Cart_shift(args.comm2d, 0, 1, &args.left, &args.right);
  MPI_Cart_shift(args.comm2d, 1, 1, &args.up, &args.down);

  args.A1_field =
      ((args.B1 - args.A1) / args.dims[0]) * args.coords[0] + args.A1;
  args.A2_field =
      ((args.B2 - args.A2) / args.dims[1]) * args.coords[1] + args.A2;

  return args;
}