/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_transpose_mpi_h = "$Id$";

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
