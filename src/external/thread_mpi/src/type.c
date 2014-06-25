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
#include <string.h>
#if !(defined( _WIN32 ) || defined( _WIN64 ) )
#include <sys/time.h>
#endif

#include "impl.h"


/* this is where all the tMPI_Reduce ops are included from tmpi_ops.h */
#define THREAD_MPI_OPS 1

#define TYPE char
#define TYPENM CHAR
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE short
#define TYPENM SHORT
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE int
#define TYPENM INT
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE long
#define TYPENM LONG
#define INTTYPE 1
#include "tmpi_ops.h"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE long long
#define TYPENM L_LONG
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE long long int
#define TYPENM L_L_INT
#define INTTYPE 1
#include "tmpi_ops.h"

#endif

#define TYPE signed char
#define TYPENM S_CHAR
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE unsigned char
#define TYPENM U_CHAR
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE unsigned short
#define TYPENM U_SHORT
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE unsigned
#define TYPENM UNSIGNED
#define INTTYPE 1
#include "tmpi_ops.h"

#define TYPE unsigned long
#define TYPENM U_LONG
#define INTTYPE 1
#include "tmpi_ops.h"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE unsigned long long
#define TYPENM U_L_LONG
#define INTTYPE 1
#include "tmpi_ops.h"

#endif

#define TYPE float
#define TYPENM FLOAT
#define INTTYPE 0
#include "tmpi_ops.h"

#define TYPE double
#define TYPENM DOUBLE
#define INTTYPE 0
#include "tmpi_ops.h"

#define TYPE long double
#define TYPENM L_DOUBLE
#define INTTYPE 0
#include "tmpi_ops.h"

#define TYPE char
#define TYPENM BYTE
#define INTTYPE 1
#include "tmpi_ops.h"

#ifdef _MSC_VER
#define TYPE __int64
#else
#define TYPE int64_t
#endif
#define TYPENM INT64_T
#define INTTYPE 1
#include "tmpi_ops.h"


/* These are the fundamental data types. They exist as global variables */
tmpi_dt tmpi_char    = {sizeof(char),              oplist_CHAR,     0, NULL, TRUE};
tmpi_dt tmpi_short   = {sizeof(short),             oplist_SHORT,    0, NULL, TRUE};
tmpi_dt tmpi_int     = {sizeof(int),               oplist_INT,      0, NULL, TRUE};
tmpi_dt tmpi_long    = {sizeof(long),              oplist_LONG,     0, NULL, TRUE};
#ifdef SIZEOF_LONG_LONG_INT
tmpi_dt tmpi_l_long  = {sizeof(long long),         oplist_L_LONG,   0, NULL, TRUE};
tmpi_dt tmpi_l_l_int = {sizeof(long long int),     oplist_L_L_INT,  0, NULL, TRUE};
#endif
tmpi_dt tmpi_s_char  = {sizeof(signed char),       oplist_S_CHAR,   0, NULL, TRUE};
tmpi_dt tmpi_u_char  = {sizeof(unsigned char),     oplist_U_CHAR,   0, NULL, TRUE};
tmpi_dt tmpi_u_short = {sizeof(unsigned short),    oplist_U_SHORT,  0, NULL, TRUE};
tmpi_dt tmpi_unsigned = {sizeof(unsigned),          oplist_UNSIGNED, 0, NULL, TRUE};
tmpi_dt tmpi_u_long  = {sizeof(unsigned long),     oplist_U_LONG,   0, NULL, TRUE};
#ifdef SIZEOF_LONG_LONG_INT
tmpi_dt tmpi_u_l_long = {sizeof(unsigned long long), oplist_U_L_LONG, 0, NULL, TRUE};
#endif
tmpi_dt tmpi_float   = {sizeof(float),             oplist_FLOAT,    0, NULL, TRUE};
tmpi_dt tmpi_double  = {sizeof(double),            oplist_DOUBLE,   0, NULL, TRUE};
tmpi_dt tmpi_l_double = {sizeof(long double),       oplist_L_DOUBLE, 0, NULL, TRUE};
tmpi_dt tmpi_byte    = {sizeof(char),              oplist_CHAR,     0, NULL, TRUE};
tmpi_dt tmpi_pointer = {sizeof(void*),             NULL,            0, NULL, TRUE};
tmpi_dt tmpi_int64_t = {8,                         oplist_INT64_T,  0, NULL, TRUE};


/* the variable types as they are referred to from MPI */
const tMPI_Datatype TMPI_CHAR               = &tmpi_char;
const tMPI_Datatype TMPI_SHORT              = &tmpi_short;
const tMPI_Datatype TMPI_INT                = &tmpi_int;
const tMPI_Datatype TMPI_LONG               = &tmpi_long;
#ifdef SIZEOF_LONG_LONG_INT
const tMPI_Datatype TMPI_LONG_LONG          = &tmpi_l_long;
const tMPI_Datatype TMPI_LONG_LONG_INT      = &tmpi_l_l_int;
#endif
const tMPI_Datatype TMPI_SIGNED_CHAR        = &tmpi_s_char;
const tMPI_Datatype TMPI_UNSIGNED_CHAR      = &tmpi_u_char;
const tMPI_Datatype TMPI_UNSIGNED_SHORT     = &tmpi_u_short;
const tMPI_Datatype TMPI_UNSIGNED           = &tmpi_unsigned;
const tMPI_Datatype TMPI_UNSIGNED_LONG      = &tmpi_u_long;
#ifdef SIZEOF_LONG_LONG_INT
const tMPI_Datatype TMPI_UNSIGNED_LONG_LONG = &tmpi_u_l_long;
#endif

const tMPI_Datatype TMPI_FLOAT              = &tmpi_float;
const tMPI_Datatype TMPI_DOUBLE             = &tmpi_double;
const tMPI_Datatype TMPI_LONG_DOUBLE        = &tmpi_l_double;

/*extern tMPI_Datatype tMPI_UNSIGNED_WCHAR*/
const tMPI_Datatype TMPI_BYTE               = &tmpi_byte;

const tMPI_Datatype TMPI_POINTER            = &tmpi_pointer;

const tMPI_Datatype TMPI_INT64_T            = &tmpi_int64_t;





int tMPI_Type_contiguous(int count, tMPI_Datatype oldtype,
                         tMPI_Datatype *newtype)
{
    struct tmpi_datatype_ *ntp;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Type_contiguous(%d, %p, %p)", count, oldtype,
                     newtype);
#endif
    ntp               = (struct tmpi_datatype_*)tMPI_Malloc(sizeof(struct tmpi_datatype_));
    ntp->size         = count*oldtype->size;
    ntp->op_functions = NULL;

    /* establish components */
    ntp->N_comp = 1;
    ntp->comps  = (struct tmpi_datatype_component*)tMPI_Malloc(
                sizeof(struct tmpi_datatype_component)*1);
    ntp->comps[0].type  = oldtype;
    ntp->comps[0].count = 1;
    ntp->committed      = FALSE;

    /* now add it to the list.  */
    tMPI_Spinlock_lock(&(tmpi_global->datatype_lock));
    /* check whether there's space */
    if (tmpi_global->N_usertypes + 1 >= tmpi_global->Nalloc_usertypes)
    {
        /* make space */
        tmpi_global->Nalloc_usertypes = Nthreads*(tmpi_global->N_usertypes) + 1;
        tmpi_global->usertypes        = (struct tmpi_datatype_**)
            tMPI_Realloc(tmpi_global->usertypes,
                         (sizeof(struct tmpi_datatype_ *)*
                          tmpi_global->Nalloc_usertypes)
                         );

    }
    /* add to the list */
    tmpi_global->usertypes[tmpi_global->N_usertypes] = ntp;
    tmpi_global->N_usertypes++;
    *newtype = ntp;
    tMPI_Spinlock_unlock(&(tmpi_global->datatype_lock));

    return TMPI_SUCCESS;
}


int tMPI_Type_commit(tMPI_Datatype *datatype)
{
    int                    i, j;
    struct tmpi_datatype_ *dt = *datatype;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Type_commit(%p)", datatype);
#endif
    if (dt->committed)
    {
        return TMPI_SUCCESS;
    }

    /* search the list for a matching committed type, because if there's
       already a committed type that has the same composition, we just
       make the datatype pointer point to it, ensuring we share datatype
       information across threads. */
    tMPI_Spinlock_lock(&(tmpi_global->datatype_lock));
    for (i = 0; i < tmpi_global->N_usertypes; i++)
    {
        struct tmpi_datatype_ *lt = tmpi_global->usertypes[i];
        if (lt->committed && lt->N_comp == dt->N_comp)
        {
            tmpi_bool found = TRUE;
            for (j = 0; j < lt->N_comp; j++)
            {
                if ( (lt->comps[j].type  != dt->comps[j].type) ||
                     (lt->comps[j].count != dt->comps[j].count) )
                {
                    found = FALSE;
                    break;
                }
            }
            if (found)
            {
                dt = lt;
            }
        }
    }
    if (dt != *datatype)
    {
        tmpi_bool found = FALSE;
        /* we remove the old one from the list */
        for (i = 0; i < tmpi_global->N_usertypes; i++)
        {
            if (tmpi_global->usertypes[i] == *datatype)
            {
                found = TRUE;
                break;
            }
        }
        if (found)
        {
            /* we put the last one in the list in our slot */
            tmpi_global->usertypes[i] = tmpi_global->
                    usertypes[tmpi_global->N_usertypes-1];
            tmpi_global->N_usertypes--;
        }
        free( (*datatype)->comps);
        free(  *datatype );

        /* and overwrite the pointer with the new data type */
        *datatype = dt;
    }
    else
    {
        /* it was the first one of its type */
        dt->committed = TRUE;
    }
    tMPI_Spinlock_unlock(&(tmpi_global->datatype_lock));
    return TMPI_SUCCESS;
}
