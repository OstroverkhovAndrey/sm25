#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <stdlib.h>

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

    int x = coords[0], y = coords[1];
    int n = 15, m = 10;
    int x0 = n*coords[0], y0 = m*coords[1];
    std::vector<std::vector<double> > v(m+2, std::vector<double>(n+2, 0.0));
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            v[j+1][i+1] = dims[0]*n*(y0+j) + (x0 + i);
        }
    }



    // Соседи по обоим измерениям
    int left, right, up, down;
    MPI_Cart_shift(comm2d, 0, 1, &left, &right); // dim=1: влево/вправо
    MPI_Cart_shift(comm2d, 1, 1, &up, &down);    // dim=0: вверх/вниз

    MPI_Status st;

    {
        // слева напрао
        double* send_v = (double*)malloc((size_t)m * sizeof(double));
        for (int i = 0; i < m; ++i) {
            send_v[i] = v[i+1][n-1 +1];
        }
        MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE,
                             right, 100,  // dest, tag
                             left,  100,  // source, tag
                             comm2d, &st);
        if (x != 0) {
            for (int i = 0; i < m; ++i) {
                v[i+1][0] = send_v[i];
            }
        }
    }
    {
        // справа на лево
        double* send_v = (double*)malloc((size_t)m * sizeof(double));
        for (int i = 0; i < m; ++i) {
            send_v[i] = v[i+1][1];
        }
        MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE,
                             left, 100,  // dest, tag
                             right,  100,  // source, tag
                             comm2d, &st);
        if (x != dims[0]-1) {
            for (int i = 0; i < m; ++i) {
                v[i+1][n+1] = send_v[i];
            }
        }
    }
    {
        // сверху вниз
        double* send_v = (double*)malloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; ++i) {
            send_v[i] = v[m-1 +1][i+1];
        }
        MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE,
                             down, 100,  // dest, tag
                             up,  100,  // source, tag
                             comm2d, &st);
        if (y != 0) {
            for (int i = 0; i < n; ++i) {
                v[0][i+1] = send_v[i];
            }
        }
    }
    {
        // справа на лево
        double* send_v = (double*)malloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; ++i) {
            send_v[i] = v[1][i+1];
        }
        // Обмен с левым/правым соседом: отправляем вправо, получаем слева
        MPI_Sendrecv_replace(send_v, m, MPI_DOUBLE,
                             up, 100,  // dest, tag
                             down,  100,  // source, tag
                             comm2d, &st);
        if (y != dims[1]-1) {
            for (int i = 0; i < n; ++i) {
                v[m+1][i+1] = send_v[i];
            }
        }
    }

    // Обмен с верхним/нижним соседом: отправляем вниз, получаем сверху
    //MPI_Sendrecv_replace(&val, 1, MPI_DOUBLE,
    //                     down,  200,
    //                     up,    200,
    //                     comm2d, &st);

    // Аккуратный вывод по порядку рангов в comm2d
    for (int r = 0; r < world_size; ++r) {
        if (r == rank2d) {
            printf("[rank2d=%2d coords=(%d,%d) dims=(%d,%d)] "
                   "L=%2d R=%2d U=%2d D=%2d\n",
                   rank2d, coords[0], coords[1], dims[0], dims[1],
                   left, right, up, down);
            for (int j = 0; j < m+2; ++j) {
                for (int i = 0; i < n+2; ++i) {
                    printf("%06.2f ", v[j][i]);
                }
                printf("\n");
            }
            // printf("recv v:\n");
            // for (int i = 0; i < m; ++i) {
            //     printf("%06.2f ", send_v[i]);
            // }
            printf("\n\n");
            fflush(stdout);
        }
        MPI_Barrier(comm2d);
    }

    MPI_Comm_free(&comm2d);
    MPI_Finalize();
    return 0;
}
