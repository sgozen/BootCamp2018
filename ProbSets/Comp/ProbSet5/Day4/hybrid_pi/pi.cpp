#include <stdio.h>
#include <time.h>
#include <cmath>
#include <omp.h>
#include <mpi.h>

#define _USE_MATH_DEFINES
static long num_steps = 4e9; 

int main (int argc, char *argv[])
{
    double time = -MPI_Wtime();
    double sum = 0., step = 1. / (double) num_steps;
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   
    // compute remainder and declare displs array via offset
    int nump = num_steps/size; 
    int rmd = num_steps % size;
    int displs[size+1];
    int offset=0;

    for (int i=0; i<rmd; i++) {
        displs[i]=offset;
        offset+=nump+1;
    }

    for (int i=rmd; i<size+1; i++) {
        displs[i]=offset;
        offset+=nump;
    }
    
    #pragma omp parallel for reduction(+:sum)
    for (int i=displs[rank]; i<displs[rank+1]; i++)
    {
        double x = (i + 0.5)*step;
        sum += 4.0 / (1.0 + x*x);
    }
    double pi = step * sum;
    double pi_root;
    
    MPI_Reduce(&pi, &pi_root, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    time += MPI_Wtime();
    if (rank==0) {
        printf("Pi is approximately: %1.6f\n", pi_root);
        printf("Relative error: %E\n", std::fabs(M_PI-pi_root)/M_PI);
        printf("Execution time: %1.6fs\n", time);
    }
    MPI_Finalize();
    return 0;
}
