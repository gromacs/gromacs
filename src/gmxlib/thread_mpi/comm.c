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

To help us fund development, we humbly ask that you cite
any papers on the package - you can find them in the top README file.

*/


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "tmpi_impl.h"

/* helper function for tMPI_Comm_split. Splits N entities with color and key
   out so that the output contains Ngroups groups each with elements
   of the same color. The group array contains the entities in each group. */
static void tMPI_Split_colors(int N, const int *color, const int *key,
                              int *Ngroups, int *grp_N, int *grp_color,
                              int *group);



/* Group query & manipulation functions */

bool tMPI_In_group(tMPI_Group group)
{
    int i;
    struct tmpi_thread *cur;

    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
            return TRUE;
    }
    return FALSE;
}

int tMPI_Group_size(tMPI_Group group, int *size)
{
    if (group)
        *size = group->N;
    else
        *size = 0;
    return TMPI_SUCCESS;
}

int tMPI_Group_rank(tMPI_Group group, int *rank)
{    
    int i;
    struct tmpi_thread *cur;

    if (!group)
        return TMPI_UNDEFINED;

    /* search for my id in the list of peers */
    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
        {
            *rank=i;
            return TMPI_SUCCESS;
        }
    }
    return TMPI_UNDEFINED;
}



tMPI_Group tMPI_Group_alloc(void)
{
    struct tmpi_group_ *ret;

    ret=(struct tmpi_group_*)tMPI_Malloc(sizeof(struct tmpi_group_));
    ret->peers=(struct tmpi_thread**)tMPI_Malloc(
                                sizeof(struct tmpi_thread*)*Nthreads);
    ret->N=0;
#if 0
    ret->Nrefs=1;
#endif

    return ret;
}

int tMPI_Group_free(tMPI_Group *group)
{
    if (group)
    {
        free((*group)->peers);
        free(*group);
    }
    return TMPI_SUCCESS;
}

int tMPI_Comm_group(tMPI_Comm comm, tMPI_Group *group)
{
    int i;
    struct tmpi_group_ *ret=tMPI_Group_alloc();

    ret->N=comm->grp.N;
    for(i=0;i<comm->grp.N;i++)
    {
        ret->peers[i]=comm->grp.peers[i];
    }
    *group=ret;
#if 0
    if (comm)
    {
        *group=&(comm->grp);
    }
    else
    {
        *group=NULL;
    }
#endif

    return TMPI_SUCCESS;
}


int tMPI_Group_incl(tMPI_Group group, int n, int *ranks, tMPI_Group *newgroup)
{
    int i;
    tMPI_Group ng;

    /* just allocate and copy */
    ng=tMPI_Group_alloc();
    ng->N=n;
    for(i=0;i<n;i++)
    {
        if (ranks[i] < 0 || !group || ranks[i] >= group->N)
        {
            return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_GROUP_RANK);
        }
        ng->peers[i]=group->peers[ranks[i]];
    }
    *newgroup=ng;
    return TMPI_SUCCESS;
}




/* communicator query&manipulation functions */
int tMPI_Comm_N(tMPI_Comm comm)
{
    if (!comm)
        return 0;
    return comm->grp.N;
}

int tMPI_Comm_size(tMPI_Comm comm, int *size)
{
    return tMPI_Group_size(&(comm->grp), size);
}

int tMPI_Comm_rank(tMPI_Comm comm, int *rank)
{
    return tMPI_Group_rank(&(comm->grp), rank);
}

tMPI_Comm tMPI_Comm_alloc(tMPI_Comm parent, int N)
{
    struct tmpi_comm_ *ret;
    int i,Nbarriers,Nred;

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

    /* calculate the number of multicast barriers */
    Nbarriers=0;
    Nred=N;
    while(Nred>1) {
        Nbarriers+=1;
        Nred = Nred/2 + Nred%2;
    } 

    ret->Nbarriers=Nbarriers;
    ret->multicast_barrier=(tMPI_Spinlock_barrier_t*)tMPI_Malloc(
                      sizeof(tMPI_Spinlock_barrier_t)*(Nbarriers+1));
    ret->N_multicast_barrier=(int*)tMPI_Malloc(sizeof(int)*(Nbarriers+1));
    Nred=N;
    for(i=0;i<Nbarriers;i++)
    {
        tMPI_Spinlock_barrier_init( &(ret->multicast_barrier[i]), Nred);
        ret->N_multicast_barrier[i]=Nred;
        /* Nred is now Nred/2 + a rest term because solitary 
           process at the end of the list must still be accounter for */
        Nred = Nred/2 + Nred%2;
    }
    ret->sendbuf=(volatile void**)tMPI_Malloc(sizeof(void*)*Nthreads);
    ret->recvbuf=(volatile void**)tMPI_Malloc(sizeof(void*)*Nthreads);


    if (parent)
    {
        ret->erh=parent->erh;
    }
    else
    {
        ret->erh=TMPI_ERRORS_ARE_FATAL;
    }

    /* multi_sync objects */
    ret->msc=(struct multi_sync*)tMPI_Malloc(sizeof(struct multi_sync)*N);
    for(i=0;i<N;i++)
        tMPI_Multi_sync_init( &(ret->msc[i]), N);

    return ret;
}

int tMPI_Comm_free(tMPI_Comm *comm)
{
#ifndef TtMPI_STRICT
    int myrank=tMPI_Comm_seek_rank(*comm, tMPI_Get_current());
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
        free((*comm)->grp.peers);
        free((*comm)->multicast_barrier);
        free((void*)(*comm)->sendbuf);
        free((void*)(*comm)->recvbuf);
        if ( (*comm)->cart)
        {
            free((*comm)->cart->dims);
            free((*comm)->cart->periods);
            free((*comm)->cart);
        }
        free(*comm);
    }
#else
    int myrank=tMPI_Comm_seek_rank(*comm, tMPI_Get_current());
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
        free((*comm)->grp.peers);
        free((*comm)->multicast_barrier);
        free((*comm)->sendbuf);
        free((*comm)->recvbuf);
        if ( (*comm)->cart)
        {
            free((*comm)->cart->dims);
            free((*comm)->cart->periods);
            free((*comm)->cart);
        }
        free(*comm);
    }
#endif
    return TMPI_SUCCESS;
}

int tMPI_Comm_dup(tMPI_Comm comm, tMPI_Comm *newcomm)
{
    /* we just call Comm_split because it already contains all the
       neccesary synchronization constructs. */
    return tMPI_Comm_split(comm, 0, tMPI_Comm_seek_rank(comm, 
                            tMPI_Get_current()), newcomm);
}


int tMPI_Comm_create(tMPI_Comm comm, tMPI_Group group, tMPI_Comm *newcomm)
{
    int color=TMPI_UNDEFINED;
    int key=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
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
    bool found;

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
    bool i_am_first=FALSE;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    struct tmpi_split *spl;

    if (!comm)
    {
        *newcomm=NULL;
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    tMPI_Thread_mutex_lock(&(comm->comm_create_lock));    
    /* first get the colors */
    if (!comm->new_comm)
    {
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
            tMPI_Thread_cond_wait( &(comm->comm_create_finish) ,
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
            tMPI_Thread_cond_wait( &(comm->comm_create_prep), 
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



/* topology functions */
int tMPI_Topo_test(tMPI_Comm comm, int *status)
{
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (comm->cart)
        *status=TMPI_CART;
    /*else if (comm->graph)
        status=MPI_GRAPH;*/
    else 
        *status=TMPI_UNDEFINED;

    return TMPI_SUCCESS;
}

int tMPI_Cartdim_get(tMPI_Comm comm, int *ndims)
{
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
    {
        return TMPI_SUCCESS;
    }
    *ndims=comm->cart->ndims;
    return TMPI_SUCCESS;
}


int tMPI_Cart_get(tMPI_Comm comm, int maxdims, int *dims, int *periods,
                 int *coords)
{
    int i;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;

    tMPI_Cart_coords(comm, myrank, maxdims, coords);

    for(i=0;i<comm->cart->ndims;i++)
    {
        if (i>=maxdims)
        {
            return tMPI_Error(comm, TMPI_ERR_DIMS);
        }
        dims[i]=comm->cart->dims[i];
        periods[i]=comm->cart->periods[i];
    }

    return TMPI_SUCCESS;
}

int tMPI_Cart_rank(tMPI_Comm comm, int *coords, int *rank)
{
    int i,mul=1,ret=0;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;

    /* because of row-major ordering, we count the dimensions down */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        int rcoord=coords[i];
        if (comm->cart->periods[i])
        {
            /* apply periodic boundary conditions */
            rcoord = rcoord % comm->cart->dims[i];
            if (rcoord < 0)
                rcoord += comm->cart->dims[i];
        }
        else
        {
            if (rcoord < 0 || rcoord >= comm->cart->dims[i])
            {
                return tMPI_Error(comm, TMPI_ERR_DIMS);
            }
        }
        ret += mul*rcoord;
        mul *= comm->cart->dims[i];
    }
    *rank=ret;
    return TMPI_SUCCESS;
}

int tMPI_Cart_coords(tMPI_Comm comm, int rank, int maxdims, int *coords)
{
    int i;
    int rank_left=rank;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return TMPI_SUCCESS;
    if (maxdims < comm->cart->ndims)
    {
        return tMPI_Error(comm, TMPI_ERR_DIMS);
    }

    /* again, row-major ordering */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        coords[i]=rank_left%comm->cart->dims[i];
        rank_left /= comm->cart->dims[i];
    }   

    return 0;
}



int tMPI_Cart_map(tMPI_Comm comm, int ndims, int *dims, int *periods, 
                 int *newrank)
{
    /* this function doesn't actually do anything beyond returning the current 
       rank (or TMPI_UNDEFINED if it doesn't fit in the new topology */
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int Ntot=1;
    int i;

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!periods)
    {
        return tMPI_Error(comm, TMPI_ERR_DIMS);
    }
 
    /* calculate the total number of procs in cartesian comm */
    for(i=0;i<ndims;i++)
    {
        Ntot *= dims[i];
    }

    if (myrank >= Ntot)
    {
        *newrank=TMPI_UNDEFINED;
    }
    else
    {
        *newrank=myrank;
    }

    return TMPI_SUCCESS;
}



int tMPI_Cart_create(tMPI_Comm comm_old, int ndims, int *dims, int *periods,
                    int reorder, tMPI_Comm *comm_cart)
{
    int myrank=tMPI_Comm_seek_rank(comm_old, tMPI_Get_current());
    int key=myrank;
    int newrank=-1;
    int color=0;
    int Ntot=1;
    int i;
    

    if (!comm_old)
    {
        return tMPI_Error(comm_old, TMPI_ERR_COMM);
    }
    /* calculate the total number of procs in cartesian comm */
    for(i=0;i<ndims;i++)
    {
        Ntot *= dims[i];
    }
    /* refuse to create if there's not enough procs */
    if (comm_old->grp.N < Ntot)
    {
        *comm_cart=TMPI_COMM_NULL;
#if 1
        return tMPI_Error(comm_old, TMPI_ERR_CART_CREATE_NPROCS);
#endif
    }

    if (key >= Ntot)
        key=TMPI_UNDEFINED;

    if (reorder)
    {
        tMPI_Cart_map(comm_old, ndims, dims, periods, &key);
    }

    if (key==TMPI_UNDEFINED)
    {
        color=TMPI_UNDEFINED;
    }

    tMPI_Comm_split(comm_old, color, key, comm_cart);

    if (*comm_cart)
    {
        tMPI_Comm_rank(*comm_cart, &newrank);
    }

    if (newrank==0)
    {
        (*comm_cart)->cart=(struct cart_topol*)tMPI_Malloc(
                                            sizeof(struct cart_topol));
        (*comm_cart)->cart->dims=(int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->periods=(int*)tMPI_Malloc(ndims*sizeof(int));
        (*comm_cart)->cart->ndims=ndims;
        for(i=0;i<ndims;i++)
        {
            (*comm_cart)->cart->dims[i]=dims[i];
            (*comm_cart)->cart->periods[i]=periods[i];
        }
    }

    /* and we add a barrier to make sure the cart object is seen by 
       every thread that is part of the new communicator */
    if (*comm_cart)
    {
        tMPI_Spinlock_barrier_wait( &( (*comm_cart)->multicast_barrier[0]) );
    }


    return TMPI_SUCCESS;
}

