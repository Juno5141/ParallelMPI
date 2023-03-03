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
    int start, end, chunk_size;
    int* chunk = (int*)malloc(WIDTH*HEIGHT*sizeof(int));

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

    chunk_size = WIDTH*HEIGHT/size;
    start = rank*chunk_size;
    end = (rank+1)*chunk_size - 1;

    for (i = start; i <= end; i++) {
        j = i % WIDTH;
        x = -2.0 + dx*j;
        y = -2.0 + dy*(i-j)/WIDTH;
        iter = mandelbrot(x, y);
        chunk[i] = iter;
    }

    MPI_Status status;
    if (rank == 0) {
        MPI_Send(chunk + chunk_size, chunk_size, MPI_INT, 1, tag, MPI_COMM_WORLD);
        printf("Rank 0 sent chunk to rank 1.\n");
    }
    else if (rank == 1) {
        MPI_Recv(chunk + chunk_size, chunk_size, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        printf("Rank 1 received chunk from rank 0.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        FILE* fp = fopen("mandelbrot.pgm", "wb");
        fprintf(fp, "P2\n%d %d\n%d\n", WIDTH, HEIGHT, MAX_ITER);
        for (i = 0; i < WIDTH*HEIGHT; i++) {
            fprintf(fp, "%d ", chunk[i]);
            if (i % WIDTH == WIDTH-1) {
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
    }

    MPI_Finalize();
    free(chunk);
    return 0;
}
