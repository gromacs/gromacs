

#include <stdio.h>
#include <stdlib.h>




#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef HAVE_SCHED_H
#include <sched.h>
#endif

#ifndef MPICC
#include "tmpi.h"
#else
#include <mpi.h>
#endif

#include "thread_mpi/atomic/cycles.h"
/*#include "thread_mpi/wait.h"*/



#include "send_recv.h"
#include "multicast.h"




#define NITER 100000
#ifdef TMPI_CYCLE_COUNT
typedef struct waiter
{
    tMPI_Cycles_t start;
    tMPI_Cycles_t stop;
} waiter;
#endif

#if 0
static void latency_tester(void)
{
#ifdef TMPI_CYCLE_COUNT
    double lat_av   = 0;
    double lat_avsq = 0;
    int    i;

    for (i = 0; i < NITER; i++)
    {
        mpi_cycles_t start;
    }

#else
    printf("No cycle count.\n");
#endif

}
#endif

static void tester(const void *arg)
{
#ifdef TMPI_CYCLE_COUNT
    int myrank;
    /*struct timeval timeout;
       timeout.tv_sec = 0; timeout.tv_usec = 0; */
    tMPI_Cycles_t start;
    tMPI_Cycles_t stop;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    printf("Starting Thread %d:\n", myrank);
    fflush(stdout);

    start = tMPI_Cycles_read();
    send_recv_data_tester();
    bcast_data_tester();
    gather_data_tester();
    /*gatherv_data_tester(); */
    scatter_data_tester();
    scatterv_data_tester();
    alltoallv_data_tester();

    stop = tMPI_Cycles_read();

    printf("Thread %d: %g cycles\n", myrank, ((double)(stop-start)));
#else
    printf("No cycle count.\n");
#endif
}

int main(int argc, char *argv[])
{
    int arg = 10;

#ifdef TMPI
    MPI_Init(&argc, &argv, main);
#else
    MPI_Init(&argc, &argv);
#endif
    /* only the main thread goes here: we've said main_thread_returns=false,
       so we can safely run the tester. */
    tester(&arg);
    MPI_Finalize();
    return 0;
}
