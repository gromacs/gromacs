
#include <stdlib.h>



#include "thread_mpi.h"

#ifdef CYCLE_COUNT

#include "thread_mpi/atomic/cycles.h"




#ifdef TMPI_CYCLE_COUNT

tMPI_Atomic_ptr_t counter;


#define HIST_N 10000
#define HIST_D 1
#define HIST_OFF 200
#define STOP_COND 1000
#define NTRIES 1000000
#define SHIFT 32
#define MASK 0xFFFFFFFFL


void cycles_tester(void)
{
    int           myrank;
    tmpi_cycles_t myrankc;
    long int      i;
    /* histogram: */
    long int     *hist;

    hist = (long int*)malloc(sizeof(long int)*HIST_N);
    /* reset hist */
    for (i = 0; i < HIST_N; i++)
    {
        hist[i] = 0;
    }

    if (sizeof(tmpi_cycles_t) != sizeof(void*))
    {
        printf("Size mismatch, can't do cycle_tester\n");
        return;
    }

    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    myrankc = (tmpi_cycles_t)myrank << SHIFT;
    for (i = 0; i < NTRIES; i++)
    {
        tmpi_cycles_t          cyc;
        volatile tmpi_cycles_t other_cyc;
        int                    other_id;

        do
        {
            other_cyc = myrankc | (MASK & tmpi_cycles_read());

            other_cyc = (tmpi_cycles_t)tMPI_Atomic_ptr_swap(&counter,
                                                            (void*)other_cyc);

            cyc      = tmpi_cycles_read();
            other_id = other_cyc >> SHIFT;
        }
        while (other_id == myrank);
        /*printf("%d -> %d (other_cyc=%ld)\n", myrank, other_id, other_cyc); fflush(stdout);*/

        if (other_id == STOP_COND)
        {
            break;
        }

        cyc       = cyc & MASK;
        other_cyc = other_cyc & MASK;
        if (i > 0)
        {
            signed long diff = cyc-other_cyc;
            signed long bin  = ((diff+HIST_OFF)/HIST_D);
            if (bin < 0)
            {
                bin = 0;
            }

            if (bin > HIST_N)
            {
                bin = HIST_N-1;
            }

            hist[bin]++;
        }
    }

    {   /* tell someone else to stop */
        tmpi_cycles_t stop = STOP_COND;
        stop = stop << SHIFT;
        tMPI_Atomic_ptr_swap(&counter, (void*)stop);
    }



    /* output histogram */
    {
        char  filename[256];
        FILE *outf;

        sprintf(filename, "hist_%03d.out", myrank);
        printf("Writing out %s..\n", filename);
        outf = fopen(filename, "w");
        for (i = 0; i < HIST_N; i++)
        {
            fprintf(outf, "%ld %ld\n", HIST_D*i - HIST_OFF, hist[i]);
        }
        fclose(outf);
    }
}
#else
void cycles_tester(void)
{
    printf("No cycle counter implemented on this platform; can't do cycle_tester\n");
}
#endif
#endif
