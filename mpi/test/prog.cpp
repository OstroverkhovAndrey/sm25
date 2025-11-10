#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Разложение по двум измерениям под размер world_size
    int dims[2] = {0, 0};                // 0 = подобрать автоматически
    MPI_Dims_create(world_size, 2, dims); // заполняет dims[], dims[0]*dims[1]==world_size

    int periods[2] = {1, 1};             // тор (периодические границы)
    int reorder = 1;
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);

    if (comm2d == MPI_COMM_NULL) {       // на всякий случай
        MPI_Finalize();
        return 0;
    }

    int rank2d, coords[2];
    MPI_Comm_rank(comm2d, &rank2d);
    MPI_Cart_coords(comm2d, rank2d, 2, coords);

    // Соседи по обоим измерениям
    int left, right, up, down;
    MPI_Cart_shift(comm2d, 1, 1, &left, &right); // dim=1: влево/вправо
    MPI_Cart_shift(comm2d, 0, 1, &up, &down);    // dim=0: вверх/вниз

    // Буфер — одно число (для демонстрации)
    double val = (double)rank2d;

    MPI_Status st;

    // Обмен с левым/правым соседом: отправляем вправо, получаем слева
    MPI_Sendrecv_replace(&val, 1, MPI_DOUBLE,
                         right, 100,  // dest, tag
                         left,  100,  // source, tag
                         comm2d, &st);

    // Обмен с верхним/нижним соседом: отправляем вниз, получаем сверху
    MPI_Sendrecv_replace(&val, 1, MPI_DOUBLE,
                         down,  200,
                         up,    200,
                         comm2d, &st);

    // Аккуратный вывод по порядку рангов в comm2d
    for (int r = 0; r < world_size; ++r) {
        if (r == rank2d) {
            printf("[rank2d=%2d coords=(%d,%d) dims=(%d,%d)] "
                   "L=%2d R=%2d U=%2d D=%2d | val=%.1f\n",
                   rank2d, coords[0], coords[1], dims[0], dims[1],
                   left, right, up, down, val);
            fflush(stdout);
        }
        MPI_Barrier(comm2d);
    }

    MPI_Comm_free(&comm2d);
    MPI_Finalize();
    return 0;
}
