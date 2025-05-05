

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

#include "thread_mpi.h"

#ifdef CYCLE_COUNT
#include "thread_mpi/atomic/cycles.h"
#include "thread_mpi/wait.h"
#endif


typedef struct
{
    unsigned int  Nbins;
    unsigned int *bins;
    TMPI_YIELD_WAIT_DATA
} hist;


void hist_init(hist *h, unsigned int Nbins)
{
    h->bins  = (unsigned int*)malloc(sizeof(unsigned int)*Nbins);
    h->Nbins = Nbins;
    TMPI_YIELD_WAIT_DATA_INIT(h);
}

void hist_sample(hist *h, unsigned int a)
{
    if (a < h->Nbins)
    {
        h->bins[a]++;
    }
    else
    {
        h->bins[h->Nbins-1]++;
    }
}

void hist_add(hist *h, hist *h2)
{
    unsigned int i;
    unsigned int N = h->Nbins;

    if (h->Nbins >  h2->Nbins)
    {
        N = h2->Nbins;
    }

    for (i = 0; i < N; i++)
    {
        h->bins[i] += h2->bins[i];
    }
}

void hist_write(hist *h, const char *filename)
{
    unsigned int i;
    FILE        *out = fopen(filename, "w");

    for (i = 0; i < h->Nbins; i++)
    {
        fprintf(out, "%u %u\n", i, h->bins[i]);
    }

    fclose(out);
}

void hist_destroy(hist *h)
{
    free(h->bins);
}


#define NITER 5000000

/*tMPI_Atomic_t shared_counter = {0};*/
int *shared_counter;
/* we only fetch every 64th int to avoid cache line conflicts */
#define STRIDE 64
int *counter_id;
int  Ncounters = 1;



#ifdef TMPI_CYCLE_COUNT

static inline void atomic_add(int *ptr, int val)
{
    __asm__ __volatile__("lock ; addl %1, %0;" : "=m" (*ptr) : "r" (val));
    /*__asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(val) :"m"(*ptr), "0"(val) : "memory");*/
}


static void tester(const void *arg)
{
    int           i;
    int           myrank;
    int           Nthreads;
    tmpi_cycles_t start;
    tmpi_cycles_t stop;
    const int     Nbins = 100000;
    hist          h;
    hist         *hsend[1];
    hist        **hrecv = NULL;
    int           n     = 0;
    int          *mycounter;

    hist_init(&h, Nbins);

    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    tMPI_Comm_size(TMPI_COMM_WORLD, &Nthreads);
    printf("Starting Thread %d; hist=%p:\n", myrank, &h);
    fflush(stdout);

    mycounter = &(shared_counter[STRIDE*counter_id[myrank]]);

    for (i = 0; i < NITER; i++)
    {
        n++;
        start = tmpi_cycles_read();
        atomic_add(mycounter, 1);
        /*tMPI_Atomic_fetch_add( &shared_counter, 1);*/
        stop = tmpi_cycles_read();

        TMPI_YIELD_WAIT(&h);
        hist_sample(&h, (unsigned int)(stop-start));
    }
    if (myrank == 0)
    {
        hrecv = (hist**)malloc(sizeof(hist*)*Nthreads);
    }
    hsend[0] = &h;
    tMPI_Gather(hsend, 1, TMPI_POINTER, hrecv, 1, TMPI_POINTER, 0,
                TMPI_COMM_WORLD);

    if (myrank == 0)
    {
        /* add up all hists*/
        for (i = 1; i < Nthreads; i++)
        {
            printf("hist[%d]=%p\n", i, hrecv[i]);
            hist_add(&h, hrecv[i]);
        }
        free(hrecv);
        hist_write(&h, "cycles.dat");
    }
    /* simple sync mechanism */
    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        int sum_total = 0;
        for (i = 0; i < Ncounters; i++)
        {
            sum_total += shared_counter[i*STRIDE];
        }
        printf("Total result: %d, expected: %d\n", sum_total, Nthreads*n);
    }

    hist_destroy(&h);
}

#else

static void tester(const void *arg)
{
    int myrank;

    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        printf("No cycle count.\n");
    }
}

#endif


int main(int argc, char *argv[])
{
    int arg = 10;
    int Nthreads;
    int i;

    /* parse the arguments */
    if (argc < 2)
    {
        goto usage;
    }

    Nthreads   = argc-1;
    counter_id = (int*)malloc(sizeof(int)*Nthreads);

    /* now check each arg */
    for (i = 1; i < argc; i++)
    {
        char *end;
        counter_id[i-1] = strtol(argv[i], &end, 10);
        if (end == argv[i])
        {
            goto usage;
        }
        if (counter_id[i-1]+1 > Ncounters)
        {
            Ncounters = counter_id[i-1]+1;
        }
    }
    shared_counter = (int*)malloc(sizeof(int)*Ncounters*STRIDE);

    printf("\nsync_cyclecount.\n\n");
    printf("Number of threads: %d\n\n", Nthreads);
    printf("Number of counters: %d\n\n", Ncounters);
    tMPI_Init_fn(1, Nthreads, TMPI_AFFINITY_ALL_CORES, tester, &arg);
    /* only the main thread goes here: we've said main_thread_returns=false,
       so we can safely run the tester. */
    tester(&arg);
    tMPI_Finalize();
    return 0;
usage:
    fprintf(stderr, "sync_cyclecount thread1_counter thread2_counter ...\n");
    fprintf(stderr, "where threadi_counter is the counter ID for thread i.\n");
    return 1;
}
