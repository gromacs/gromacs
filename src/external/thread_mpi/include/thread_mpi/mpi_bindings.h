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

#ifndef TMPI_MPI_BINDINGS_H_
#define TMPI_MPI_BINDINGS_H_

#include "tmpi.h"

/** \file
   \brief MPI bindings for thread_mpi/tmpi.h

   This file contains only macros and redefinitions to expose the standard
   MPI API with thread_mpi calls.

   This is different from the API exposed through thread_mpi/tmpi.h, which
   uses names like \a tMPI_Send() instead of \a MPI_Send()

   \sa thread_mpi/tmpi.h for documentation of the available data types and
      functions
   \sa http://www.mpi-forum.org/docs/docs.html for MPI documentation.
 */

#ifndef DOXYGEN
#ifdef __cplusplus
extern "C"
{
#endif




/* The MPI_Comm structure contains the group of processes to communicate
   with (defines the scope for global operations such as broadcast) */
typedef struct tmpi_comm_ *MPI_Comm;
/* The group part of the MPI-Comm structure */
typedef struct tmpi_group_ *MPI_Group;
/* Request structure for holding data about non-blocking transfers */
typedef struct tmpi_req_ *MPI_Request;
/* status of receives */
typedef struct tmpi_status_ MPI_Status;
/* data types */
typedef struct tmpi_datatype_ *MPI_Datatype;
/* reduce operations */
typedef tMPI_Op MPI_Op;


#define MPI_CHAR                TMPI_CHAR
#define MPI_SHORT               TMPI_SHORT
#define MPI_INT                 TMPI_INT
#define MPI_LONG                TMPI_LONG
#define MPI_LONG_LONG           TMPI_LONG_LONG
#define MPI_LONG_LONG_INT       TMPI_LONG_LONG_INT
#define MPI_SIGNED_CHAR         TMPI_SIGNED_CHAR
#define MPI_UNSIGNED_CHAR       TMPI_UNSIGNED_CHAR
#define MPI_UNSIGNED_SHORT      TMPI_UNSIGNED_SHORT
#define MPI_UNSIGNED            TMPI_UNSIGNED
#define MPI_UNSIGNED_LONG_LONG  TMPI_UNSIGNED_LONG_LONG
#define MPI_FLOAT               TMPI_FLOAT
#define MPI_DOUBLE              TMPI_DOUBLE
#define MPI_LONG_DOUBLE         TMPI_LONG_DOUBLE
#define MPI_WCHAR               TMPI_WCHAR
#define MPI_BYTE                TMPI_BYTE
#define MPI_INT64_T             TMPI_INT64_T


#define MPI_SUCCESS                 TMPI_SUCCESS
#define MPI_ERR_MALLOC              TMPI_ERR_MALLOC
#define MPI_ERR_INIT                TMPI_ERR_INIT
#define MPI_ERR_FINALIZE            TMPI_ERR_FINALIZE
#define MPI_ERR_GROUP               TMPI_ERR_GROUP
#define MPI_ERR_COMM                TMPI_ERR_COMM
#define MPI_ERR_STATUS              TMPI_ERR_STATUS
#define MPI_ERR_GROUP_RANK          TMPI_ERR_GROUP_RANK
#define MPI_ERR_DIMS                TMPI_ERR_DIMS
#define MPI_ERR_COORDS              TMPI_ERR_COORDS
#define MPI_ERR_CART_CREATE_NPROCS  TMPI_ERR_CART_CREATE_NPROCS
#define MPI_ERR_XFER_COUNTERPART    TMPI_ERR_XFER_COUNTERPART
#define MPI_ERR_XFER_BUFSIZE        TMPI_ERR_XFER_BUFSIZE
#define MPI_ERR_XFER_BUF_OVERLAP    TMPI_ERR_XFER_BUF_OVERLAP
#define MPI_ERR_SEND_DEST           TMPI_ERR_SEND_DEST
#define MPI_ERR_RECV_SRC            TMPI_ERR_RECV_SRC
#define MPI_ERR_BUF                 TMPI_ERR_BUF
#define MPI_ERR_MULTI_MISMATCH      TMPI_ERR_MULTI_MISMATCH
#define MPI_ERR_OP_FN               TMPI_ERR_OP_FN
#define MPI_ERR_ENVELOPES           TMPI_ERR_ENVELOPES
#define MPI_ERR_REQUESTS            TMPI_ERR_REQUESTS
#define MPI_ERR_IN_STATUS           TMPI_ERR_IN_STATUS
#define MPI_FAILURE                 TMPI_FAILURE
#define MPI_ERR_UNKNOWN             TMPI_ERR_UNKNOWN
#define N_MPI_ERR                   N_TMPI_ERR

#define MPI_MAX_ERROR_STRING        TMPI_MAX_ERROR_STRING
#define MPI_UNDEFINED               TMPI_UNDEFINED


#define MPI_Errhandler_fn           tMPI_Errhandler_fn
#define MPI_Errhandler              tMPI_Errhandler
#define MPI_ERRORS_ARE_FATAL        TMPI_ERRORS_ARE_FATAL
#define MPI_ERRORS_RETURN           TMPI_ERRORS_RETURN



/* miscelaneous defines */
#define MPI_ANY_SOURCE          TMPI_ANY_SOURCE
#define MPI_ANY_TAG             TMPI_ANY_TAG

/* comm_compare defines */
#define MPI_IDENT               TMPI_IDENT
#define MPI_CONGRUENT           TMPI_CONGRUENT
#define MPI_SIMILAR             TMPI_SIMILAR
#define MPI_UNEQUAL             TMPI_UNEQUAL


/* topology test defines */
#define MPI_CART                TMPI_CART
#define MPI_GRAPH               TMPI_GRAPH


#define MPI_COMM_WORLD          TMPI_COMM_WORLD
#define MPI_COMM_NULL           TMPI_COMM_NULL


#define MPI_GROUP_NULL          TMPI_GROUP_NULL
#define MPI_GROUP_EMPTY         TMPI_GROUP_EMPTY

#define MPI_MAX_PROCESSOR_NAME  TMPI_MAX_PROCESSOR_NAME


/* MPI status */
#define MPI_STATUS_IGNORE       TMPI_STATUS_IGNORE
#define MPI_STATUSES_IGNORE     TMPI_STATUSES_IGNORE

#define MPI_SOURCE              TMPI_SOURCE
#define MPI_TAG                 TMPI_TAG
#define MPI_ERROR               TMPI_ERROR

#define mpi_status_             tmpi_status_

#define MPI_REQUEST_NULL        TMPI_REQUEST_NULL

#define MPI_IN_PLACE            TMPI_IN_PLACE



#define MPI_MAX          TMPI_MAX
#define MPI_MIN          TMPI_MIN
#define MPI_SUM          TMPI_SUM
#define MPI_PROD         TMPI_PROD
#define MPI_LAND         TMPI_LAND
#define MPI_BAND         TMPI_BAND
#define MPI_LOR          TMPI_LOR
#define MPI_BOR          TMPI_BOR
#define MPI_LXOR         TMPI_LXOR
#define MPI_BXOR         TMPI_BXOR

/* the functions: */
#define MPI_Init                    tMPI_Init
#define MPI_Finalize                tMPI_Finalize
#define MPI_Abort                   tMPI_Abort
#define MPI_Initialized             tMPI_Initialized
#define MPI_Finalized               tMPI_Finalized

#define MPI_Create_errhandler       tMPI_Create_errhandler
#define MPI_Errhandler_free         tMPI_Errhandler_free
#define MPI_Comm_set_errhandler     tMPI_Comm_set_errhandler
#define MPI_Comm_get_errhandler     tMPI_Comm_get_errhandler
#define MPI_Error_string            tMPI_Error_string

#define MPI_Get_processor_name      tMPI_Get_processor_name
#define MPI_Wtime                   tMPI_Wtime

#define MPI_Group_size              tMPI_Group_size
#define MPI_Group_rank              tMPI_Group_rank
#define MPI_Group_incl              tMPI_Group_incl
#define MPI_Comm_group              tMPI_Comm_group
#define MPI_Group_free              tMPI_Group_free

#define MPI_Comm_size               tMPI_Comm_size
#define MPI_Comm_rank               tMPI_Comm_rank
#define MPI_Comm_compare            tMPI_Comm_compare
#define MPI_Comm_free               tMPI_Comm_free
#define MPI_Comm_create             tMPI_Comm_create
#define MPI_Comm_split              tMPI_Comm_split
#define MPI_Comm_dup                tMPI_Comm_dup

#define MPI_Topo_test               tMPI_Topo_test
#define MPI_Cartdim_get             tMPI_Cartdim_get
#define MPI_Cart_get                tMPI_Cart_get
#define MPI_Cart_rank               tMPI_Cart_rank
#define MPI_Cart_coords             tMPI_Cart_coords
#define MPI_Cart_map                tMPI_Cart_map
#define MPI_Cart_create             tMPI_Cart_create
#define MPI_Cart_sub                tMPI_Cart_sub

#define MPI_Type_contiguous         tMPI_Type_contiguous
#define MPI_Type_commit             tMPI_Type_commit

#define MPI_Send                    tMPI_Send
#define MPI_Recv                    tMPI_Recv
#define MPI_Sendrecv                tMPI_Sendrecv
#define MPI_Get_count               tMPI_Get_count

#define MPI_Isend                   tMPI_Isend
#define MPI_Irecv                   tMPI_Irecv
#define MPI_Test                    tMPI_Test
#define MPI_Wait                    tMPI_Wait
#define MPI_Waitall                 tMPI_Waitall

#define MPI_Barrier                 tMPI_Barrier

#define MPI_Bcast                   tMPI_Bcast

#define MPI_Gather                  tMPI_Gather
#define MPI_Gatherv                 tMPI_Gatherv
#define MPI_Scatter                 tMPI_Scatter
#define MPI_Scatterv                tMPI_Scatterv
#define MPI_Alltoall                tMPI_Alltoall
#define MPI_Alltoallv               tMPI_Alltoallv


#define MPI_Reduce                  tMPI_Reduce
#define MPI_Allreduce               tMPI_Allreduce
#define MPI_Scan                    tMPI_Scan

#ifdef __cplusplus
} /* closing extern "C" */
#endif

#endif

#endif /* TMPI_MPI_BINDINGS_H_ */
