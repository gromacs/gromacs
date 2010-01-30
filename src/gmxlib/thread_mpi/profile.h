




#ifdef TMPI_PROFILE

#include "thread_mpi/atomic/cycles.h"

enum tmpi_functions 
{
    TMPIFN_Send=0, /* first the point-to-point comm functions */
    TMPIFN_Recv,
    TMPIFN_Sendrecv,
    TMPIFN_Isend,
    TMPIFN_Irecv,
    TMPIFN_Wait,
    TMPIFN_Test,
    TMPIFN_Waitall,

    TMPIFN_Barrier, /* then the barrier */

    TMPIFN_Bcast, /* and now the collective comm functions */
    TMPIFN_Gather,
    TMPIFN_Gatherv,
    TMPIFN_Scatter,
    TMPIFN_Scatterv,
    TMPIFN_Alltoall,
    TMPIFN_Alltoallv,

    TMPIFN_Nfunctions
};

enum tmpi_wait_functions
{
    TMPIWAIT_Send, /* p2p send wait */
    TMPIWAIT_Recv, /* p2p recv wait */
    TMPIWAIT_Waitall, /* p2p waitall */
    TMPIWAIT_Coll_send, /* collective recv wait */
    TMPIWAIT_Coll_recv, /* collective recv wait */
    TMPIWAIT_Barrier, /* collective recv wait */
    TMPIWAIT_Reduce, /* collective (all)reduce wait */

    TMPIWAIT_N 
};


/* thread-specific profiling data structure */
struct tmpi_profile 
{
    unsigned long int mpifn_calls[TMPIFN_Nfunctions]; /* array of counters */

#ifdef TMPI_CYCLE_COUNT
    /* cycle counters */
    tmpi_cycles_t wait_cycles[TMPIWAIT_N]; /* the wait cycles */

    tmpi_cycles_t global_start,global_stop; /* timing start and stop times */
    tmpi_cycles_t start,stop; /* individual timing start and stop times.
                                 This can be here because tmpi_profile is
                                 thread-specific. */
    double totals;            /* totals counter for reporting end results */
#endif
};

/* initialize the profile counter */
void tMPI_Profile_init(struct tmpi_profile *prof);

/* counter functions */
void tMPI_Profile_count(enum tmpi_functions fn);
void tMPI_Profile_count_thread(struct tmpi_thread *th,enum tmpi_functions fn);

/* stop counting */
void tMPI_Profile_stop(struct tmpi_profile *prof);

/* wait functions */

#ifdef TMPI_CYCLE_COUNT
/* start waiting cycle count */
void tMPI_Profile_wait_start();
void tMPI_Profile_wait_start_thread(struct tmpi_thread *th);

/* stop waiting cycle count */
void tMPI_Profile_wait_stop(enum tmpi_wait_functions fn);
void tMPI_Profile_wait_stop_thread(struct tmpi_thread *th, 
                                   enum tmpi_wait_functions fn);
#endif



/* output functions */
void tMPI_Profiles_summarize(int Nthreads, struct tmpi_thread *threads);
/* deallocations */
void tMPI_Profiles_destroy(int Nthreads, struct tmpi_thread *threads);

#endif


