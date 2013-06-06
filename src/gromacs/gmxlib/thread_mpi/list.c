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

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "thread_mpi/list.h"


void tMPI_Stack_init(tMPI_Stack *st)
{
    tMPI_Atomic_ptr_set(&(st->head), NULL);
}

void tMPI_Stack_destroy(tMPI_Stack *st)
{
    tMPI_Atomic_ptr_set(&(st->head), NULL);
}

void tMPI_Stack_push(tMPI_Stack *st, tMPI_Stack_element *el)
{
    tMPI_Stack_element *head;
    do
    {
        head     = (tMPI_Stack_element*)tMPI_Atomic_ptr_get( &(st->head) );
        el->next = head;
    }
    while (!tMPI_Atomic_ptr_cas(&(st->head), head, el));
}

tMPI_Stack_element *tMPI_Stack_pop(tMPI_Stack *st)
{
    tMPI_Stack_element *head, *next;
    do
    {
        head = (tMPI_Stack_element*)tMPI_Atomic_ptr_get( &(st->head) );
        if (head)
        {
            next = head->next;
        }
        else
        {
            next = NULL;
        }
    }
    while (!tMPI_Atomic_ptr_cas(&(st->head), head, next));

    return head;
}

tMPI_Stack_element *tMPI_Stack_detach(tMPI_Stack *st)
{
    tMPI_Stack_element *head;
    do
    {
        head = (tMPI_Stack_element*)tMPI_Atomic_ptr_get( &(st->head) );
    }
    while (!tMPI_Atomic_ptr_cas(&(st->head), head, NULL));

    return head;
}





#if 0
void tMPI_Queue_init(tMPI_Queue *q)
{
    tMPI_Atomic_ptr_set( &(q->head), NULL);
    tMPI_Atomic_ptr_set( &(q->tail), NULL);
}


void tMPI_Queue_destroy(tMPI_Queue *q)
{
    tMPI_Atomic_ptr_set( &(q->head), NULL);
    tMPI_Atomic_ptr_set( &(q->tail), NULL);
}

void tMPI_Queue_enqueue(tMPI_Queue *q, tMPI_Queue_element *qe)
{
    tMPI_Queue_element *head, *next;

    do
    {
    }
    while (!tMPI_Atomic_ptr_cas(&(q->head), head, next));
}
#endif
