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


#ifdef THREAD_MPI_OPS

/* cpp wizardry follows...

   This file is #included directly from thread_mpi.c, and constructs
   MPI_Reduce operators.

   What this does is create the min, max, sum, prod, etc. functions for a given
   datatype (pre-defined as TYPE, with identifier name TYPENM) and puts pointers
   to these functions in an array called oplist_TYPENM.

   gmx_thread_mpi_reduce.c includes this file once for each type used by MPI,
   and thus builds up a set of arrays of function pointers, that then get used
   in the mpi_datatype_ structure. This way, each operation/datatype entry
   that makes sense can be extracted easily. Note that we don't (yet) support
   user-defined ops */

#define FNAMEr(tp, fn) tMPI_ ## tp ## _ ## fn
#define FNAME(tp, fn) FNAMEr(tp, fn)

/* macros to define functions and prototypes based on a name and an operation */
#define FNr(tp, fname, fn) \
    static void tMPI_ ## tp ## _ ## fname  (void *dest, \
                                            const void *src_a, const void *src_b, \
                                            int count) \
    { \
        /*printf("in function %s, count=%d\n", __FUNCTION__, count);*/ \
        const TYPE *a = (const TYPE*)src_a; \
        const TYPE *b = (const TYPE*)src_b; \
        TYPE *d = (TYPE*)dest; \
        int   i; \
        for (i = 0; i < count; i++) { \
            d[i] = (TYPE)(fn(a[i], b[i])); } \
    }

#define FN(tp, fname, fn) FNr(tp, fname, fn)

#define OPFNr(tp, fname, operator)  \
              static void tMPI_ ## tp ## _ ## fname  (void *dest, \
                                                      const void *src_a, const void *src_b, \
                                                      int count) \
              { \
                  /*printf("in function %s, count=%d\n", __FUNCTION__, count);*/ \
                  const TYPE *a = (const TYPE*)src_a; \
                  const TYPE *b = (const TYPE*)src_b; \
                  TYPE *d = (TYPE*)dest; \
                  int i; \
                  for (i = 0; i < count; i++) { \
                      d[i] = (TYPE)(a[i] operator b[i]); } \
              }

#define OPFN(tp, fname, operator) OPFNr(tp, fname, operator)


/* these are the function prototypes + definitions: */
#if TYPECATEGORY!=LOGICALTYPE
#define MAX(a, b)  (( (a) > (b) ) ? (a) : (b))
FN(TYPENM, max, MAX)
#undef MAX
#define MIN(a, b)  (( (a) < (b) ) ? (a) : (b))
FN(TYPENM, min, MIN)
#undef MIN
OPFN(TYPENM, sum, +)
OPFN(TYPENM, prod, *)
#endif
#if TYPECATEGORY!=FLOATTYPE
OPFN(TYPENM, land, &&)
OPFN(TYPENM, lor, ||)
#define XOR(a, b)  ( (!a) ^ (!b) )
FN(TYPENM, lxor, XOR)
#undef XOR
#endif
#if TYPECATEGORY==INTTYPE
OPFN(TYPENM, band, &)
OPFN(TYPENM, bor, |)
OPFN(TYPENM, bxor, ^)
#endif

#define OPARRAYr(tp) oplist_ ## tp
#define OPARRAY(tp) OPARRAYr(tp)

tMPI_Op_fn OPARRAY(TYPENM)[] =
{
#if TYPECATEGORY==LOGICALTYPE
    0,
    0,
    0,
    0,
    FNAME(TYPENM, land),
    0,
    FNAME(TYPENM, lor),
    0,
    FNAME(TYPENM, lxor),
    0,
#else
    FNAME(TYPENM, max),
    FNAME(TYPENM, min),
    FNAME(TYPENM, sum),
    FNAME(TYPENM, prod),
#if TYPECATEGORY==INTTYPE
    FNAME(TYPENM, land),
    FNAME(TYPENM, band),
    FNAME(TYPENM, lor),
    FNAME(TYPENM, bor),
    FNAME(TYPENM, lxor),
    FNAME(TYPENM, bxor)
#else
    0,
    0,
    0,
    0,
    0,
    0
#endif
#endif
};


#undef FNAME
#undef FNAMEr
#undef OPARRAYr
#undef OPARRAY
#undef FN
#undef FNr
#undef OPFN
#undef OPFNr

#undef TYPE
#undef TYPENM
#undef TYPECATEGORY
#endif
