/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <fftw_mpi.h>

#include "sched.h"

/* This file contains routines to compute communications schedules for
   all-to-all communications (complete exchanges) that are performed
   in-place.  (That is, the block that processor x sends to processor
   y gets replaced on processor x by a block received from processor y.)

   A schedule, int **sched, is a two-dimensional array where
   sched[pe][i] is the processor that pe expects to exchange a message
   with on the i-th step of the exchange.  sched[pe][i] == -1 for the
   i after the last exchange scheduled on pe.

   Here, processors (pe's, for processing elements), are numbered from
   0 to npes-1.

   There are a couple of constraints that a schedule should satisfy
   (besides the obvious one that every processor has to communicate
   with every other processor exactly once).
   
   * First, and most importantly, there must be no deadlocks.
   
   * Second, we would like to overlap communications as much as possible,
   so that all exchanges occur in parallel.  It turns out that perfect
   overlap is possible if npes is even, and only a single extra step is
   required if npes is odd.

   It turns out that this scheduling problem is actually well-studied,
   and good solutions are known.  The problem is known as a
   "time-tabling" problem, and is specifically the problem of
   scheduling a sports competition (where n teams must compete exactly
   once with every other team).  The problem is discussed and
   algorithms are presented in:

   [1] J. A. M. Schreuder, "Constructing Timetables for Sport
   Competitions," Mathematical Programming Study 13, pp. 58-67 (1980).

   [2] A. Schaerf, "Scheduling Sport Tournaments using Constraint
   Logic Programming," Proc. of 12th Europ. Conf. on
   Artif. Intell. (ECAI-96), pp. 634-639 (Budapest 1996).
   http://hermes.dis.uniromal.it/~aschaerf/publications.html

   (These people actually impose a lot of additional constraints that
   we don't care about, so they are solving harder problems. [1] gives
   a simple enough algorithm for our purposes, though.)

   However, we have to do more: for a particular processor, the
   communications schedule must be sorted in ascending or descending
   order of processor index.  (This is necessary so that the data
   coming in for the transpose does not overwrite data that will be
   sent later; for that processor the incoming and outgoing blocks are
   of different non-zero sizes.)

   Fortunately, it is possible to reorder the schedule to achieve any
   permutation on a given processor while maintaining the two required
   properties above. ...except, when npes is odd (when the schedule
   already contains a stall), our reordering introduces an extra stall
   due to the motion of the self-communication past a stall.  We could
   fix this if it were really important, but it turns out that the
   extra stall is not introduced in the case that we care about (when
   the sorted processor is the first or last processor). */

/* Create a new communications schedule for a given number of processors.
   The schedule is initialized to a deadlock-free, maximum overlap
   schedule.  Returns NULL on an error (may print a message to
   stderr if there is a program bug detected).  */
int **make_comm_schedule(int npes)
{
     int **sched;
     int i;

     sched = (int **) fftw_malloc(sizeof(int *) * npes);
     if (!sched)
	  return NULL;

     for (i = 0; i < npes; ++i)
	  sched[i] = NULL;

     for (i = 0; i < npes; ++i) {
	  sched[i] = (int *) fftw_malloc(sizeof(int) * 10 * (npes + 1));
	  if (!sched[i]) {
	       free_comm_schedule(sched,npes);
	       return NULL;
	  }
     }
     
     empty_comm_schedule(sched,npes);
     fill_comm_schedule(sched,npes);

     if (!check_comm_schedule(sched,npes)) {
	  free_comm_schedule(sched,npes);
	  return NULL;
     }

     return sched;
}

void free_comm_schedule(int **sched, int npes)
{
     if (sched) {
	  int i;

	  for (i = 0; i < npes; ++i)
	       fftw_free(sched[i]);
	  fftw_free(sched);
     }
}

void empty_comm_schedule(int **sched, int npes)
{
     int i;
     for (i = 0; i < npes; ++i)
	  sched[i][0] = -1;
}

static void add_dest_to_comm_schedule(int **sched, int pe, int dest)
{
     int i;
     
     for (i = 0; sched[pe][i] != -1; ++i)
	  ;

     sched[pe][i] = dest;
     sched[pe][i+1] = -1;
}

static void add_pair_to_comm_schedule(int **sched, int pe1, int pe2)
{
     add_dest_to_comm_schedule(sched, pe1, pe2);
     if (pe1 != pe2)
	  add_dest_to_comm_schedule(sched, pe2, pe1);
}

/* Simplification of algorithm presented in [1] (we have fewer
   constraints).  Produces a perfect schedule if npes is even;
   otherwise contains one unavoidable extra step. */

void fill_comm_schedule(int **sched, int npes)
{
     int pe, i, n;

     for (pe = 0; pe < npes; ++pe)
	  add_pair_to_comm_schedule(sched,pe,pe);

     if (npes % 2 == 0)
	  n = npes;
     else
	  n = npes + 1;

     for (pe = 0; pe < n - 1; ++pe) {
	  if (pe != npes - 1)
	       add_pair_to_comm_schedule(sched,pe,npes - 1);
	  
	  for (i = 1; i < n/2; ++i) {
	       int pe_a, pe_b;

	       pe_a = pe - i;
	       if (pe_a < 0)
		    pe_a += n - 1;

	       pe_b = (pe + i) % (n - 1);

	       if (pe_a != npes - 1 && pe_b != npes - 1)
		    add_pair_to_comm_schedule(sched,pe_a,pe_b);
	  }
     }
}

/* Below, we have various checks in case of bugs: */

/* check for deadlocks by simulating the schedule and looking for
   cycles in the dependency list; returns 0 if there are deadlocks
   (or other errors) */
static int check_schedule_deadlock(int **sched, int npes)
{
     int *step, *depend, *visited, pe, pe2, period, done = 0;
     int counter = 0;

     /* step[pe] is the step in the schedule that a given pe is on */
     step = (int *) fftw_malloc(sizeof(int) * npes);

     /* depend[pe] is the pe' that pe is currently waiting for a message
	from (-1 if none) */
     depend = (int *) fftw_malloc(sizeof(int) * npes);

     /* visited[pe] tells whether we have visited the current pe already
	when we are looking for cycles. */
     visited = (int *) fftw_malloc(sizeof(int) * npes);

     if (!step || !depend || !visited) {
	  fftw_free(step); fftw_free(depend); fftw_free(visited);
	  return 0;
     }

     for (pe = 0; pe < npes; ++pe)
	  step[pe] = 0;

     while (!done) {
	  ++counter;

	  for (pe = 0; pe < npes; ++pe)
	       depend[pe] = sched[pe][step[pe]];
	  
	  /* now look for cycles in the dependencies with period > 2: */
	  for (pe = 0; pe < npes; ++pe)
	       if (depend[pe] != -1) {
		    for (pe2 = 0; pe2 < npes; ++pe2)
			 visited[pe2] = 0;

		    period = 0;
		    pe2 = pe;
		    do {
			 visited[pe2] = period + 1;
			 pe2 = depend[pe2];
			 period++;
		    } while (pe2 != -1 && !visited[pe2]);

		    if (pe2 == -1) {
			 fprintf(stderr,
				 "BUG: unterminated cycle in schedule!\n");
			 fftw_free(step); fftw_free(depend);
			 fftw_free(visited);
			 return 0;
		    }
		    if (period - (visited[pe2] - 1) > 2) {
			 fprintf(stderr,"BUG: deadlock in schedule!\n");
			 fftw_free(step); fftw_free(depend);
			 fftw_free(visited);
			 return 0;
		    }

		    if (pe2 == pe)
			 step[pe]++;
	       }

	  done = 1;
	  for (pe = 0; pe < npes; ++pe)
	       if (sched[pe][step[pe]] != -1) {
		    done = 0;
		    break;
	       }
     }

     fftw_free(step); fftw_free(depend); fftw_free(visited);
     return (counter > 0 ? counter : 1);
}

/* sanity checks; prints message and returns 0 on failure.
   undocumented feature: the return value on success is actually the
   number of steps required for the schedule to complete, counting
   stalls. */
int check_comm_schedule(int **sched, int npes)
{
     int pe, i, comm_pe;
     
     for (pe = 0; pe < npes; ++pe) {
	  for (comm_pe = 0; comm_pe < npes; ++comm_pe) {
	       for (i = 0; sched[pe][i] != -1 && sched[pe][i] != comm_pe; ++i)
		    ;
	       if (sched[pe][i] == -1) {
		    fprintf(stderr,"BUG: schedule never sends message from "
			    "%d to %d.\n",pe,comm_pe);
		    return 0;  /* never send message to comm_pe */
	       }
	  }
	  for (i = 0; sched[pe][i] != -1; ++i)
	       ;
	  if (i != npes) {
	       fprintf(stderr,"BUG: schedule sends too many messages from "
		       "%d\n",pe);
	       return 0;
	  }
     }
     return check_schedule_deadlock(sched,npes);
}

/* invert the order of all the schedules; this has no effect on
   its required properties. */
void invert_comm_schedule(int **sched, int npes)
{
     int pe, i;

     for (pe = 0; pe < npes; ++pe)
	  for (i = 0; i < npes/2; ++i) {
	       int dummy = sched[pe][i];
	       sched[pe][i] = sched[pe][npes-1-i];
	       sched[pe][npes-1-i] = dummy;
	  }
}

/* Relabel pe1 <-> pe2 in all the schedules.  The required schedule
   properties are invariant under this operation.  */
static void comm_schedule_swap(int **sched, int npes, int pe1, int pe2)
{
     int pe, i, *dummy;

     for (pe = 0; pe < npes; ++pe)
	  for (i = 0; sched[pe][i] != -1; ++i) {
	       if (sched[pe][i] == pe1)
		    sched[pe][i] = pe2;
	       else if (sched[pe][i] == pe2)
		    sched[pe][i] = pe1;
	  }

     dummy = sched[pe1];
     sched[pe1] = sched[pe2];
     sched[pe2] = dummy;
}

/* Sort the schedule for sort_pe in ascending order of processor
   index.  Unfortunately, for odd npes (when schedule has a stall
   to begin with) this will introduce an extra stall due to
   the motion of the self-communication past a stall.  We could
   fix this if it were really important.  Actually, we don't
   get an extra stall when sort_pe == 0 or npes-1, which is sufficient
   for our purposes. */
void sort_comm_schedule(int **sched, int npes, int sort_pe)
{
     int i,j,pe;

     /* Note that we can do this sort in O(npes) swaps because we know
	that the numbers we are sorting are just 0...npes-1. */

     /* find self-communication: */
     for (i = 0; i < npes; ++i)
	  if (sched[sort_pe][i] == sort_pe)
	       break;

     if (i == npes) {
	  fprintf(stderr,"BUG: missing self-communication for %d.",sort_pe);
	  exit(1);
     }

     /* Note that, to maintain communications
	overlap, we perform the same shift
	on the self-communication for all pe's.
	
	The self-communication is always at the
	same step for all pe's (this invariant
	is maintained both by this operation
	and by comm_schedule_swap). */

     /* shift self-communications to the correct place: 
        -- this has to be done separately because we cannot perform 
	   simple swaps of the self-communication elements */
     for (pe = 0; pe < npes; ++pe) {
	  if (sched[pe][i] != pe) {
	       fprintf(stderr,
		       "BUG: self-communication isn't at fixed step.");
	       exit(1);
	  }
	  for (j = i; j+1 < npes; ++j)
	       sched[pe][j] = sched[pe][j+1];
          for (j = npes - 2; j >= sort_pe; --j)
	       sched[pe][j+1] = sched[pe][j];
	  sched[pe][sort_pe] = pe;
     }

     /* Move the remaining items to their sorted positions: */     
     for (pe = 0; pe < npes; ++pe)
	  if (pe != sort_pe) {
	       for (j = 0; j < npes; ++j)
		    if (sched[sort_pe][j] == pe) break;
	       if (j == npes) {
		    fprintf(stderr,"BUG: missing %d in %d schedule.",
			    pe, sort_pe);
		    exit(1);
	       }
	       /* move communications with pe to correct position: */
	       comm_schedule_swap(sched,npes,
				  sched[sort_pe][pe],
				  sched[sort_pe][j]);
	  }
}

/* print the schedule (for debugging purposes) */
void print_comm_schedule(int **sched, int npes)
{
     int pe, i, width;

     if (npes < 10)
	  width = 1;
     else if (npes < 100)
	  width = 2;
     else
	  width = 3;

     for (pe = 0; pe < npes; ++pe) {
	  printf("pe %*d schedule:", width, pe);
	  for (i = 0; sched[pe][i] != -1; ++i)
	       printf("  %*d",width,sched[pe][i]);
	  printf("\n");
     }
}
