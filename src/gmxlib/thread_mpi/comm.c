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

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


#include "impl.h"

/* helper function for tMPI_Comm_split. Splits N entities with color and key
   out so that the output contains Ngroups groups each with elements
   of the same color. The group array contains the entities in each group. */
static void tMPI_Split_colors(int N, const int *color, const int *key,
                              int *Ngroups, int *grp_N, int *grp_color,
                              int *group);






/* communicator query&manipulation functions */
int tMPI_Comm_N(tMPI_Comm comm)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_N(%p)", comm);
#endif
    if (!comm)
        return 0;
    return comm->grp.N;
}

int tMPI_Comm_size(tMPI_Comm comm, int *size)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_size(%p, %p)", comm, size);
#endif
    return tMPI_Group_size(&(comm->grp), size);
}

int tMPI_Comm_rank(tMPI_Comm comm, int *rank)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_rank(%p, %p)", comm, rank);
#endif
    return tMPI_Group_rank(&(comm->grp), rank);
}


int tMPI_Comm_compare(tMPI_Comm comm1, tMPI_Comm comm2, int *result)
{
    int i,j;
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_compare(%p, %p, %p)", comm1, comm2, result);
#endif
    if (comm1 == comm2)
    {
        *result=TMPI_IDENT;
        return TMPI_SUCCESS;
    }

    if ( (!comm1) || (!comm2) )
    {
        *result=TMPI_UNEQUAL;
        return TMPI_SUCCESS;
    }

    if (comm1->grp.N != comm2->grp.N)
    {
        *result=TMPI_UNEQUAL;
        return TMPI_SUCCESS;
    }

    *result=TMPI_CONGRUENT;
    /* we assume that there are two identical comm members within a comm */
    for(i=0;i<comm1->grp.N;i++)
    {
        if (comm1->grp.peers[i] != comm2->grp.peers[i])
        {
            tmpi_bool found=FALSE;

            *result=TMPI_SIMILAR;
            for(j=0;j<comm2->grp.N;j++)
            {
                if (comm1->grp.peers[i] == comm2->grp.peers[j])
                {
                    found=TRUE;
                    break;
                }
            }
            if (!found)
            {
                *result=TMPI_UNEQUAL;
                return TMPI_SUCCESS;
            }
        }
    }
    return TMPI_SUCCESS;
}


tMPI_Comm tMPI_Comm_alloc(tMPI_Comm parent, int N)
{
    struct tmpi_comm_ *ret;
    int i;

    ret=(struct tmpi_comm_*)tMPI_Malloc(sizeof(struct tmpi_comm_));
    ret->grp.peers=(struct tmpi_thread**)tMPI_Malloc(
                                sizeof(struct tmpi_thread*)*Nthreads);
    ret->grp.N=N;

    tMPI_Thread_mutex_init( &(ret->comm_create_lock) );
    tMPI_Thread_cond_init( &(ret->comm_create_prep) );
    tMPI_Thread_cond_init( &(ret->comm_create_finish) );

    ret->split = NULL;
    ret->new_comm = NULL;
    /* we have no topology to start out with */
    ret->cart=NULL;
    /*ret->graph=NULL;*/


    /* initialize the main barrier */
    tMPI_Barrier_init(&(ret->barrier), N);

    /* the reduce barriers */
    {
        /* First calculate the number of reduce barriers */
        int Niter=0; /* the iteration number */
        int Nred=N; /* the number of reduce barriers for this iteration */
        while(Nred>1) {
            /* Nred is now Nred/2 + a rest term because solitary 
               process at the end of the list must still be accounter for */
            Nred = Nred/2 + Nred%2;
            Niter+=1;
        } 

        ret->N_reduce_iter=Niter;
        /* allocate the list */
        ret->reduce_barrier=(tMPI_Barrier_t**)
                  tMPI_Malloc(sizeof(tMPI_Barrier_t*)*(Niter+1));
        ret->N_reduce=(int*)tMPI_Malloc(sizeof(int)*(Niter+1));

        /* we re-set Nred to N */
        Nred=N;
        for(i=0;i<Niter;i++)
        {
            int j;

            Nred = Nred/2 + Nred%2;
            ret->N_reduce[i] = Nred;
            /* allocate the sub-list */
            ret->reduce_barrier[i]=(tMPI_Barrier_t*)
                      tMPI_Malloc(sizeof(tMPI_Barrier_t)*(Nred));
            for(j=0;j<Nred;j++)
            {
                tMPI_Barrier_init(&(ret->reduce_barrier[i][j]),2);
            }
        }
    }

    /* the reduce buffers */
#if 0
    ret->sendbuf=(volatile void**)tMPI_Malloc(sizeof(void*)*Nthreads);
    ret->recvbuf=(volatile void**)tMPI_Malloc(sizeof(void*)*Nthreads);
#else
    ret->reduce_sendbuf=(tMPI_Atomic_ptr_t*)
              tMPI_Malloc(sizeof(tMPI_Atomic_ptr_t)*Nthreads);
    ret->reduce_recvbuf=(tMPI_Atomic_ptr_t*)
              tMPI_Malloc(sizeof(tMPI_Atomic_ptr_t)*Nthreads);
#endif


    if (parent)
    {
        ret->erh=parent->erh;
    }
    else
    {
        ret->erh=TMPI_ERRORS_ARE_FATAL;
    }

    /* coll_env objects */
    ret->cev=(struct coll_env*)tMPI_Malloc(sizeof(struct coll_env)*N_COLL_ENV);
    for(i=0;i<N_COLL_ENV;i++)
        tMPI_Coll_env_init( &(ret->cev[i]), N);
    /* multi_sync objects */
    ret->csync=(struct coll_sync*)tMPI_Malloc(sizeof(struct coll_sync)*N);
    for(i=0;i<N;i++)
        tMPI_Coll_sync_init( &(ret->csync[i]), N);

    /* we insert ourselves in the circular list, after TMPI_COMM_WORLD */
    if (TMPI_COMM_WORLD)
    {
        ret->next=TMPI_COMM_WORLD;
        ret->prev=TMPI_COMM_WORLD->prev;

        TMPI_COMM_WORLD->prev->next = ret;
        TMPI_COMM_WORLD->prev = ret;
    }
    else
    {
        ret->prev=ret->next=ret;
    }

    return ret;
}

void tMPI_Comm_destroy(tMPI_Comm comm)
{
    int i;

    free(comm->grp.peers);
#if 0
    free(comm->reduce_barrier);
    free(comm->N_reduce_barrier);
#endif
    for(i=0;i<comm->N_reduce_iter;i++)
        free(comm->reduce_barrier[i]);
    free(comm->reduce_barrier);
    free(comm->N_reduce);

    for(i=0;i<N_COLL_ENV;i++)
        tMPI_Coll_env_destroy( &(comm->cev[i]) );
    for(i=0;i<comm->grp.N;i++)
        tMPI_Coll_sync_destroy( &(comm->csync[i]) );
    free(comm->cev);
    free(comm->csync);

    tMPI_Thread_mutex_destroy( &(comm->comm_create_lock) );
    tMPI_Thread_cond_destroy( &(comm->comm_create_prep) );
    tMPI_Thread_cond_destroy( &(comm->comm_create_finish) );

    free((void*)comm->reduce_sendbuf);
    free((void*)comm->reduce_recvbuf);

    if ( comm->cart )
    {
        tMPI_Cart_destroy( comm->cart );
        free(comm->cart);
    }

    /* remove ourselves from the circular list */
    if (comm->next)
        comm->next->prev=comm->prev;
    if (comm->prev)
        comm->prev->next=comm->next;

    free(comm);
}

int tMPI_Comm_free(tMPI_Comm *comm)
{
    int myrank=tMPI_Comm_seek_rank(*comm, tMPI_Get_current());
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_free(%p)", comm);
#endif
#ifndef TMPI_STRICT
    if (! *comm)
        return TMPI_SUCCESS;

    if ((*comm)->grp.N > 1)
    {
        /* we remove ourselves from the comm. */
        tMPI_Thread_mutex_lock(&((*comm)->comm_create_lock));
        (*comm)->grp.peers[myrank] = (*comm)->grp.peers[(*comm)->grp.N-1];
        (*comm)->grp.N--;
        tMPI_Thread_mutex_unlock(&((*comm)->comm_create_lock));
    }
    else
    {
        /* we're the last one so we can safely destroy it */
        tMPI_Comm_destroy(*comm);
    }
#else
    /* This is correct if programs actually treat Comm_free as a 
       collective call */
    /* we need to barrier because the comm is a shared structure and
       we have to be sure that nobody else is using it 
       (for example, to get its rank, like above) before destroying it*/
    tMPI_Barrier(*comm);
    /* this is a collective call on a shared data structure, so only 
       one process (rank[0] in this case) should do anything */
    if (myrank==0)
    {
        tMPI_Comm_destroy(*comm);
    }
#endif
    return TMPI_SUCCESS;
}

int tMPI_Comm_dup(tMPI_Comm comm, tMPI_Comm *newcomm)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_dup(%p, %p)", comm, newcomm);
#endif
    /* we just call Comm_split because it already contains all the
       neccesary synchronization constructs. */
    return tMPI_Comm_split(comm, 0, tMPI_Comm_seek_rank(comm, 
                            tMPI_Get_current()), newcomm);
}


int tMPI_Comm_create(tMPI_Comm comm, tMPI_Group group, tMPI_Comm *newcomm)
{
    int color=TMPI_UNDEFINED;
    int key=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_create(%p, %p, %p)", comm, group, newcomm);
#endif
    if (tMPI_In_group(group))
    {
        color=1;
    }
    /* the MPI specs specifically say that this is equivalent */
    return tMPI_Comm_split(comm, color, key, newcomm);
}

static void tMPI_Split_colors(int N, const int *color, const int *key, 
                              int *Ngroups, int *grp_N, int *grp_color, 
                              int *group)
{
    int i,j;
    tmpi_bool found;

    /* reset groups */
    for(i=0;i<N;i++)
        grp_N[i]=0;
    for(i=0;i<N;i++)
    {
        if (color[i] != TMPI_UNDEFINED)
        {
            found=FALSE;
            for(j=0;j<(*Ngroups);j++)
            {
                if (grp_color[j] == color[i])
                {
                    /* we insert where we need to, by counting back */
                    int k=grp_N[j];

                    while (k>0 && ( key[group[N*j + k-1]]>key[i]) )
                    {
                        /* shift up */
                        group[N*j + k]=group[N*j + k-1];
                        k--;
                    }
                    group[N*j+k]=i;
                    grp_N[j]++;
                    found=TRUE;
                }
            }
            if (!found)
            {
                /* not found. just add a new color */
                grp_N[(*Ngroups)]=1;
                grp_color[(*Ngroups)]=color[i];
                group[N*(*Ngroups) + 0]=i;
                (*Ngroups)++;
            }
        }
    }
}

/* this is the main comm creation function. All other functions that create
    comms use this*/
int tMPI_Comm_split(tMPI_Comm comm, int color, int key, tMPI_Comm *newcomm)
{
    int i,j;
    int N=tMPI_Comm_N(comm);
    volatile tMPI_Comm *newcomm_list;
    volatile int colors[MAX_PREALLOC_THREADS]; /* array with the colors 
                                                  of each thread */
    volatile int keys[MAX_PREALLOC_THREADS]; /* same for keys (only one of 
                                                the threads actually suplies 
                                                these arrays to the comm 
                                                structure) */
    tmpi_bool i_am_first=FALSE;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    struct tmpi_split *spl;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Comm_split(%p, %d, %d, %p)", comm, color, key, 
                       newcomm);
#endif
    if (!comm)
    {
        *newcomm=NULL;
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    tMPI_Thread_mutex_lock(&(comm->comm_create_lock));    
    /* first get the colors */
    if (!comm->new_comm)
    {
        /* i am apparently  first */
        comm->split=(struct tmpi_split*)tMPI_Malloc(sizeof(struct tmpi_split));
        comm->new_comm=(tMPI_Comm*)tMPI_Malloc(N*sizeof(tMPI_Comm));
        if (N<=MAX_PREALLOC_THREADS)
        {
            comm->split->colors=colors;
            comm->split->keys=keys;
        }
        else
        {
            comm->split->colors=(int*)tMPI_Malloc(N*sizeof(int));
            comm->split->keys=(int*)tMPI_Malloc(N*sizeof(int));
        }
        comm->split->Ncol_init=tMPI_Comm_N(comm); 
        comm->split->can_finish=FALSE;
        i_am_first=TRUE;
        /* the main communicator contains a list the size of grp.N */
    }
    newcomm_list=comm->new_comm; /* we copy it to the local stacks because
                                    we can later erase comm->new_comm safely */
    spl=comm->split; /* we do the same for spl */
    spl->colors[myrank] = color;
    spl->keys[myrank] = key;
    spl->Ncol_init--;

    if (spl->Ncol_init == 0)
        tMPI_Thread_cond_signal(&(comm->comm_create_prep));

    if (!i_am_first)
    {
        /* all other threads can just wait until the creator thread is 
           finished */
        while(! spl->can_finish )
        {
            tMPI_Thread_cond_wait(&(comm->comm_create_finish) ,
                                  &(comm->comm_create_lock) );
        }
    }
    else
    {
        int Ncomms=0;
        int comm_color_[MAX_PREALLOC_THREADS]; 
        int comm_N_[MAX_PREALLOC_THREADS]; 
        int *comm_color=comm_color_; /* there can't be more comms than N*/
        int *comm_N=comm_N_; /* the number of procs in a group */

        int *comm_groups; /* the groups */
        tMPI_Comm *comms; /* the communicators */

        /* wait for the colors to be done */
        /*if (N>1)*/
        while(spl->Ncol_init > 0)
        {
            tMPI_Thread_cond_wait(&(comm->comm_create_prep), 
                                  &(comm->comm_create_lock));
        }

        /* reset the state so that a new comm creating function can run */
        spl->Ncol_destroy=N;
        comm->new_comm=0;
        comm->split=0;

        comm_groups=(int*)tMPI_Malloc(N*N*sizeof(int));
        if (N>MAX_PREALLOC_THREADS)
        {
            comm_color=(int*)tMPI_Malloc(N*sizeof(int));
            comm_N=(int*)tMPI_Malloc(N*sizeof(int));
        }

        /* count colors, allocate and split up communicators */
        tMPI_Split_colors(N, (int*)spl->colors, 
                             (int*)spl->keys, 
                             &Ncomms, 
                             comm_N, comm_color, comm_groups);


        /* allocate a bunch of communicators */
        comms=(tMPI_Comm*)tMPI_Malloc(Ncomms*sizeof(tMPI_Comm));
        for(i=0;i<Ncomms;i++)
            comms[i]=tMPI_Comm_alloc(comm, comm_N[i]);

        /* now distribute the comms */
        for(i=0;i<Ncomms;i++)
        {
            comms[i]->grp.N=comm_N[i];
            for(j=0;j<comm_N[i];j++)
                comms[i]->grp.peers[j]=
                    comm->grp.peers[comm_groups[i*comm->grp.N + j]];
        }
        /* and put them into the newcomm_list */
        for(i=0;i<N;i++)
        {
            newcomm_list[i]=TMPI_COMM_NULL;
            for(j=0;j<Ncomms;j++)
            {
                if (spl->colors[i] == comm_color[j])
                {
                    newcomm_list[i] = comms[j];
                    break;
                }
            }
        }

#ifdef TMPI_DEBUG
        /* output */
        for(i=0;i<Ncomms;i++)
        {
            printf("Group %d (color %d) has %d members: ",
                    i, comm_color[i], comm_N[i]);
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",comm_groups[comm->grp.N*i + j]);

            printf(" rank: ");
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",spl->keys[comm_groups[N*i + j]]);
            printf(" color: ");
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",spl->colors[comm_groups[N*i + j]]);
            printf("\n");
        }
#endif
        if (N>MAX_PREALLOC_THREADS)
        {
            free((int*)spl->colors);
            free((int*)spl->keys);
            free(comm_color);
            free(comm_N);
        }
        free(comm_groups);
        free(comms);
        spl->can_finish=TRUE;

        /* tell the waiting threads that there's a comm ready */
        tMPI_Thread_cond_broadcast(&(comm->comm_create_finish));
    }
    /* here the individual threads get their comm object */
    *newcomm=newcomm_list[myrank];

    /* free when we have assigned them all, so we can reuse the object*/
    spl->Ncol_destroy--;
    if (spl->Ncol_destroy==0)
    {
        free((void*)newcomm_list);
        free(spl);
    }

    tMPI_Thread_mutex_unlock(&(comm->comm_create_lock));

    return TMPI_SUCCESS;    
}

int tMPI_Comm_seek_rank(tMPI_Comm comm, struct tmpi_thread *th)
{
    int i;
    if (!comm)
        return -1;

    for(i=0;i<comm->grp.N;i++)
    {
        if (comm->grp.peers[i] == th)
            return i;
    }
    return -1;
}



