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

#ifndef TMPI_LIST_H_
#define TMPI_LIST_H_

#include "atomic.h"


/** \file
 *
 * \brief Lock-free list data structures.
 *
 */


#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif


/**  Lock-free single-ended stack (FIFO)

   Is a list with push, pop and detach operations */
typedef struct
{
    tMPI_Atomic_ptr_t head;          /**< Pointer to the top stack element. */
} tMPI_Stack;

/**  A single element in stack */
typedef struct tMPI_Stack_element
{
    struct tMPI_Stack_element *next; /**< Pointer to the next stack element. */
    void                      *data; /**< Pointer to data. */
} tMPI_Stack_element;


/**  Initialize a stack */
TMPI_EXPORT
void tMPI_Stack_init(tMPI_Stack *st);

/**  Deallocates a stack */
TMPI_EXPORT
void tMPI_Stack_destroy(tMPI_Stack *st);

/**  Pushes a stack element onto a stack */
TMPI_EXPORT
void tMPI_Stack_push(tMPI_Stack *st, tMPI_Stack_element *el);

/**  Pops a stack element from  a stack */
TMPI_EXPORT
tMPI_Stack_element *tMPI_Stack_pop(tMPI_Stack *st);

/**  Detaches entire stack for use by a single thread */
TMPI_EXPORT
tMPI_Stack_element *tMPI_Stack_detach(tMPI_Stack *st);



#if 0
/**  Lock-free double-ended queue (FIFO)

   Is a list with enqueue and dequeue operations */
typedef struct
{
    tMPI_Atomic_ptr_t head, tail;
} tMPI_Queue;

/**  A single element in a queue */
typedef struct tMPI_Queue_element
{
    struct tMPI_Queue_element *next; /**< Pointer to the next queue element. */
    struct tMPI_Queue_element *prev; /**< Pointer to the prev queue element. */
    void                      *data; /**< Pointer to data. */
} tMPI_Queue_element;

/**  Initialize a queue */
void tMPI_Queue_init(tMPI_Queue *q);

/**  Deallocates a queue */
void tMPI_Queue_destroy(tMPI_Queue *q);

/**  Enqueue an element onto the head of a queue */
void tMPI_Queue_enqueue(tMPI_Queue *q, tMPI_Queue_element *qe);

/**  Dequeue an element from the end a queue */
tMPI_Queue_element *tMPI_Queue_dequeue(tMPI_Queue *q);





/**  Lock-free circular doubly linked list */
typedef struct
{
    tMPI_Atomic_ptr_t head;
} tMPI_List;

/**  Lock-free circular doubly linked list */
typedef struct tMPI_List_element
{
    struct tMPI_List_element *next, *prev;
    void                     *data;
} tMPI_List_element;

/**  Initialize a list */
void tMPI_List_init(tMPI_List *l);
/**  Deallocates a list */
void tMPI_List_destroy(tMPI_List *l);

tMPI_List_element* tMPI_List_first(tMPI_List *l);
tMPI_List_element* tMPI_List_next(tMPI_List         *l,
                                  tMPI_List_element *le);
tMPI_List_element* tMPI_List_prev(tMPI_List         *l,
                                  tMPI_List_element *le);

void tMPI_List_add(tMPI_List *l, tMPI_List_element *le);
void tMPI_List_insert(tMPI_List *l, tMPI_List_element *after,
                      tMPI_List_element *le);
void tMPI_List_remove(tMPI_List *l, tMPI_List_element *le);
#endif


#ifdef __cplusplus
} /* closing extern "C" */
#endif

#endif /* TMPI_LIST_H_ */
