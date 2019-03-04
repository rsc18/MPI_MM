//
// Created by ram on 2/28/19.
//
#include <mpi.h>
#include <stdio.h>
#include<stdlib.h>


double a[4000][4000];
double b[4000][4000];
double c[4000][4000];
int N=4;

MPI_Status status;

int main(void)
{
    int world_size,partial_rows,offset;
    //N = atoi(argv[1]);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);




    int number;
    if (my_rank == 0) {
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                a[i][j]=1.0;
                b[i][j]=1.0;
            }
        }

        printf("\nPrinting from Master\n");

        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                printf("%f\t",b[i][j]);
            }
            printf("\n");
        }


        partial_rows=N/(world_size);
        printf("%d is my partial rows",partial_rows);
        offset=0;

        //divide matrix A and send to all processes
        for(int i=1;i<world_size;i++) {
            MPI_Send(&offset, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&partial_rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&a[offset][0], partial_rows*N, MPI_DOUBLE,i,0, MPI_COMM_WORLD);
            offset+=partial_rows;
            //send matrix B to all
            MPI_Send(&b, N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }
    else if (my_rank > 0) {
        int count;
//        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
//        MPI_Get_count(&status, MPI_INT, &count);
//        int* number_buf = (int*)malloc(sizeof(int) * count);

        MPI_Recv(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&partial_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, partial_rows*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&b, N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);


        //receive matrix A part and Matrix B full
        //Compute matrix mul
//        for (int k=0; k<N; k++) {
//            for (int i = 0; i < partial_rows; i++) {
//                c[i][k] = 0.0;
//                for (int j = 0; j < N; j++)
//                    c[i][k] = c[i][k] + a[i][j] * b[j][k];
//            }
//        }
        printf("\nPrinting from Process: %d Partial A: \n",my_rank);
        for(int i=0;i<partial_rows;i++){
            for(int j=0;j<N;j++){
                printf("%f\t",a[i][j]);
            }
            printf("\n");
        }
        printf("\nPrinting from Process: %d Full B: \n",my_rank);

        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                printf("%f\t",b[i][j]);
            }
            printf("\n");
        }
    }
/*
    int myrank;
    int i, j, k;
    myrank = *(int *)arg;
    for (i=0; i<N; i++)
        for (j=N/Nthreads*myrank; j < N/Nthreads*(myrank+1)-1; j++)
            for (k=0; k<N; k++)
                c[i][j] = c[i][j] + a[i][k]*b[k][j];*/

    MPI_Finalize();
    return 0;
}
