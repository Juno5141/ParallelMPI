#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 5000

int mandelbrot(double x, double y) {
    double zx = 0.0, zy = 0.0, zx_new;
    int iter = 0;
    while (zx*zx + zy*zy <= 4.0 && iter < MAX_ITER) {
        zx_new = zx*zx - zy*zy + x;
        zy = 2.0*zx*zy + y;
        zx = zx_new;
        iter++;
    }
    return iter;
}

int main(int argc, char** argv) {
    int rank, size, i, j, iter, tag = 0;
    double x, y, dx = 4.0/WIDTH, dy = 4.0/HEIGHT;
    int start, end, count;
    int* chunk = (int*)malloc(WIDTH*HEIGHT*sizeof(int));
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0) {
            printf("This program should be run with 2 processes.\n");
        }
        MPI_Finalize();
        exit(1);
    }

    if (rank == 0) {
        count = 0;
        for (i = 0; i < WIDTH*HEIGHT; i++) {
            MPI_Isend(&i, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &request);
            count++;
        }
        start = count;
        end = 0;
        MPI_Isend(&end, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &request);
        printf("Rank 0 sent tasks to rank 1.\n");
    }
    else {
        while (1) {
            MPI_Recv(&i, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
            if (i == 0) {
                break;
            }
            j = i % WIDTH;
            x = -2.0 + dx*j;
            y = -2.0 + dy*(i-j)/WIDTH;
            iter = mandelbrot(x, y);
            chunk[i] = iter;
        }
        MPI_Recv(&end, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        printf("Rank 1 received all tasks from rank 0.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        count = 0;
        for (i = 0; i < WIDTH*HEIGHT; i++) {
            if (chunk[i] > 0) {
                printf("*");
            } else {
                printf(".");
            }
            count++;
            if (count == WIDTH) {
                printf("\n");
                count = 0;
            }
        }
        free(chunk);
    }

    MPI_Finalize();
    return 0;
}
