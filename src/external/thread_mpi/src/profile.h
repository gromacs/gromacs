/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you want to redistribute modifications, please consider that
   scientific software is very special. Version control is crucial -
   bugs must be traceable. We will be happy to consider code for
   inclusion in the official distribution, but derived work should not
   be called official thread_mpi. Details are found in the README & COPYING
   files.
 */


/* the profiling functions. Many of these are macros, so they're inlined
   forcibly. Profiling is turned on by defining TMPI_PROFILE, but the most
   useful parts depend on the cycle counter, which currently only works for
   x86, x86_64 and ia64. */
#ifdef TMPI_PROFILE

#include "thread_mpi/atomic/cycles.h"

struct tmpi_thread;

enum tmpi_functions
{
    TMPIFN_Send = 0, /* first the point-to-point comm functions */
    TMPIFN_Recv,
    TMPIFN_Sendrecv,
    TMPIFN_Isend,
    TMPIFN_Irecv,
    TMPIFN_Wait,
    TMPIFN_Test,
    TMPIFN_Waitall,
    TMPIFN_Testall,
    TMPIFN_Waitany,
    TMPIFN_Testany,
    TMPIFN_Waitsome,
    TMPIFN_Testsome,

    TMPIFN_Barrier, /* then the barrier */

    TMPIFN_Bcast,   /* and now the collective comm functions */
    TMPIFN_Gather,
    TMPIFN_Gatherv,
    TMPIFN_Scatter,
    TMPIFN_Scatterv,
    TMPIFN_Alltoall,
    TMPIFN_Alltoallv,

    TMPIFN_Reduce,
    TMPIFN_Allreduce,
    TMPIFN_Scan,

    TMPIFN_Nfunctions
};

enum tmpi_wait_functions
{
    TMPIWAIT_P2p,        /* p2p send wait */
    TMPIWAIT_P2p_signal, /* p2p signaling wait */
    TMPIWAIT_Coll_send,  /* collective recv wait */
    TMPIWAIT_Coll_recv,  /* collective recv wait */
    TMPIWAIT_Barrier,    /* collective recv wait */
    TMPIWAIT_Reduce,     /* collective (all)reduce wait */

    TMPIWAIT_N
};


/* thread-specific profiling data structure */
struct tmpi_profile
{
    unsigned long int mpifn_calls[TMPIFN_Nfunctions]; /* array of counters */

    unsigned long int buffered_p2p_xfers;             /* number of buffered p2p transfers */
    unsigned long int total_p2p_xfers;                /* total number of p2p transfers */

    unsigned long int buffered_coll_xfers;            /* number of buffered collective
                                                         transfers */
    unsigned long int total_coll_xfers;               /* total number of collective
                                                         transfers */

#ifdef TMPI_CYCLE_COUNT
    /* cycle counters */
    tMPI_Cycles_t mpifn_cycles[TMPIFN_Nfunctions]; /* array of cycle counters */
    tMPI_Cycles_t wait_cycles[TMPIWAIT_N];         /* the wait cycles */

    tMPI_Cycles_t global_start, global_stop;       /* timing start and stop times */
    tMPI_Cycles_t mpifn_start;                     /* individual timing start times for profiling
                                                      function call times.  This can be here
                                                      because tmpi_profile is thread-specific. */
    enum tmpi_functions fn;                        /* the function being cycle-counted */


    tMPI_Cycles_t wait_start; /* individual timing start times for profiling
                                 wait times. */

    double totals;            /* totals counter for reporting end results */
#endif
};

extern int tMPI_Profile_started;

/* initialize the profile counter */
int tMPI_Profile_init(struct tmpi_profile *prof);

#if 0
/* deallocations */
void tMPI_Profile_destroy(struct tmpi_profile *prof);
#endif

/* stop counting */
void tMPI_Profile_stop(struct tmpi_profile *prof);



/* counter functions */
/* start */
#ifdef TMPI_CYCLE_COUNT
/*void tMPI_Profile_count_start(struct tmpi_thread *th);*/
#define tMPI_Profile_count_start(th) { th->profile.mpifn_start = tMPI_Cycles_read(); }
#else
#define tMPI_Profile_count_start(th) {}
#endif

/* end. this is where the counting actually happens */
/*void tMPI_Profile_count_stop(struct tmpi_thread *th, enum tmpi_functions fn);*/
#ifdef TMPI_CYCLE_COUNT
#define tMPI_Profile_count_stop(th, fn) \
    { \
        tMPI_Cycles_t stop = tMPI_Cycles_read(); \
        th->profile.mpifn_cycles[fn] += (stop - th->profile.mpifn_start); \
        (th->profile.mpifn_calls[fn])++; \
    }
#else
#define tMPI_Profile_count_stop(th, fn) \
    { \
        (th->profile.mpifn_calls[fn])++; \
    }
#endif







/* wait functions */
#ifdef TMPI_CYCLE_COUNT
/* start waiting cycle count */
/*void tMPI_Profile_wait_start(struct tmpi_thread *th);*/
#define tMPI_Profile_wait_start(th) \
    { \
        th->profile.wait_start = tMPI_Cycles_read(); \
    }

/* stop waiting cycle count */
/*void tMPI_Profile_wait_stop(struct tmpi_thread *th,
                            enum tmpi_wait_functions fn);*/
#define tMPI_Profile_wait_stop(th, fn) \
    { \
        tMPI_Cycles_t wait_stop = tMPI_Cycles_read(); \
        th->profile.wait_cycles[fn] += (wait_stop - th->profile.wait_start); \
    }
#else
#define tMPI_Profile_wait_start(th) {}
#define tMPI_Profile_wait_stop(th, fn) {}
#endif


/* count the number of transfers at the receiving end. */
/*void tMPI_Profile_count_buffered_p2p_xfer(struct tmpi_thread *th);
   void tMPI_Profile_count_p2p_xfer(struct tmpi_thread *th);
   void tMPI_Profile_count_buffered_coll_xfer(struct tmpi_thread *th);
   void tMPI_Profile_count_coll_xfer(struct tmpi_thread *th);*/
#define tMPI_Profile_count_buffered_p2p_xfer(th) \
    { \
        (th->profile.buffered_p2p_xfers)++; \
    }

#define tMPI_Profile_count_p2p_xfer(th) \
    { \
        (th->profile.total_p2p_xfers)++; \
    }

#define tMPI_Profile_count_buffered_coll_xfer(th) \
    { \
        (th->profile.buffered_coll_xfers)++; \
    }

#define tMPI_Profile_count_coll_xfer(th) \
    { \
        (th->profile.total_coll_xfers)++; \
    }



/* output functions */
void tMPI_Profiles_summarize(int Nthreads, struct tmpi_thread *threads);

#endif
