#include <args.hpp>

#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>

double square_of_max(double a, double b) {
  double m = std::max(a, b);
  return m * m;
}

// считаем какая степень двойки x
// если не степень двойки возвращаем 1
int from_power_two(int x) {
  int ans = 0;
  while (x % 2 == 0) {
    ans++;
    x /= 2;
  }
  if (x != 1) {
    std::cout << "Error MPI size not power two!" << std::endl;
    return 1;
  }
  return ans;
}

// возводим число в степень двойки 
int to_power_two(int x) {
  int ans = 1;
  for (int i = 0; i < x; ++i) {
    ans *= 2;
  }
  return ans;
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
    if (std::string(argv[i]) == "print_result") {
      args.print_result = true;
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
  if (args.M < 10 || args.N < 10) {
    std::cout << "Error too small M or N" << std::endl;
  }
  args.N -= 1;
  args.M -= 1;

  int k = from_power_two(world_size);
  int k_m = 0, k_n = 0;
  double M_temp = args.M, N_temp = args.N;
  while (k > 0) {
    if (M_temp > N_temp) {
      M_temp /= 2.0;
      k_m++;
    } else {
      N_temp /= 2.0;
      k_n++;
    }
    k--;
  }

  args.dims[0] = to_power_two(k_m);
  args.dims[1] = to_power_two(k_n);

  int periods[2] = {0, 0};
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, args.dims, periods, reorder, &args.comm2d);

  MPI_Comm_rank(args.comm2d, &args.rank2d);
  MPI_Cart_coords(args.comm2d, args.rank2d, 2, args.coords);

  MPI_Cart_shift(args.comm2d, 0, 1, &args.left, &args.right);
  MPI_Cart_shift(args.comm2d, 1, 1, &args.up, &args.down);

  args.M_field = (int)(args.M / to_power_two(k_m));
  args.N_field = (int)(args.N / to_power_two(k_n));

  //  сколько ечеек в предыдущих блоках
  int prev_cell_m = args.M_field * args.coords[0];
  int prev_cell_n = args.N_field * args.coords[1];

  if (args.M % to_power_two(k_m) > args.coords[0]) {
    args.M_field += 1;
    prev_cell_m += args.coords[0];
  } else {
    prev_cell_m += args.M % to_power_two(k_m);
  }
  if (args.N % to_power_two(k_n) > args.coords[1]) {
    args.N_field += 1;
    prev_cell_n += args.coords[1];
  } else {
    prev_cell_n += args.N % to_power_two(k_n);
  }

  args.A1_field =
      args.h1 * prev_cell_m + args.A1;
  args.A2_field =
      args.h2 * prev_cell_n + args.A2;

  args.N += 1;
  args.M += 1;

  return args;
}