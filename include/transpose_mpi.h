/*
 * Copyright (c) 1997,1998 Massachusetts Institute of Technology
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

#ifndef TRANSPOSE_MPI_H
#define TRANSPOSE_MPI_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <mpi.h> /* need access to the MPI type constants */
#include "fftw.h"

/*********************** Transpose Element Type ************************/

/* Here, we define the data type of the individual elements in the
   transposed array. */

typedef fftw_real TRANSPOSE_EL_TYPE;

/* We need to set the MPI type constant corresponding to the type
   TRANSPOSE_EL_TYPE.  (MPI needs to know the data type in case it
   has to make translations on heterogeneous networks.)  Possible
   type constants are MPI_DOUBLE, MPI_FLOAT, MPI_INT, MPI_LONG, ...
   see your MPI manual for more. */
#ifdef DOUBLE
#define MPI_TRANSPOSE_EL_TYPE_CONSTANT MPI_DOUBLE
#else
#define MPI_TRANSPOSE_EL_TYPE_CONSTANT MPI_FLOAT
#endif

/***********************************************************************/

typedef struct {
     int block_num, dest_pe, send_size, recv_size;
} transpose_mpi_block_dest;

typedef struct {
     MPI_Comm comm;
     int n_pes, my_pe;
     
     int nx,ny,local_nx,local_ny;

     transpose_mpi_block_dest *block_dest;
     int num_blocks, block_size;
     int block_size_padded;

     int *perm_block_dest;
     int num_perm_blocks, perm_block_size;

     char *move;
     int move_size;
} transpose_mpi_plan_struct;

typedef transpose_mpi_plan_struct *transpose_mpi_plan;

extern void transpose_mpi_get_local_size(int n, int my_pe, int n_pes,
					 int *local_n, int *local_start);
extern int transpose_mpi_get_local_storage_size(int nx, int ny,
						int my_pe, int n_pes);

extern transpose_mpi_plan transpose_mpi_create_plan(int nx, int ny, 
						    MPI_Comm comm);
extern void transpose_mpi_destroy_plan(transpose_mpi_plan p);

extern void transpose_mpi(transpose_mpi_plan p, 
			  TRANSPOSE_EL_TYPE *local_data, int el_size);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* TRANSPOSE_MPI_H */
