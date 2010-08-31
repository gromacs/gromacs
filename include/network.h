/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _network_h
#define _network_h


/*
 * This module defines the interface of the actual communication routines.
 */

#include <stdio.h>

#include "types/simple.h"
#include "types/commrec.h"
#include "typedefs.h"
#include "main.h"
#include "gmx_fatal.h"

#ifdef __cplusplus
extern "C" {
#endif

int gmx_setup(int *argc,char **argv,int *nnodes);
/* Initializes the parallel communication, return the ID of the node */

int gmx_node_num(void);
/* return the number of nodes in the ring */

int gmx_node_rank(void);
/* return the rank of the node */

void gmx_setup_nodecomm(FILE *fplog,t_commrec *cr);
/* Sets up fast global communication for clusters with multi-core nodes */

gmx_bool gmx_mpi_initialized(void);
/* return TRUE when MPI_Init has been called.
 * return FALSE when MPI_Init has not been called OR
 * when GROMACS was compiled without MPI support.
 */

void gmx_barrier(const t_commrec *cr);
/* Wait till all processes in cr->mpi_comm_mygroup have reached the barrier */

void gmx_bcast(int nbytes,void *b,const t_commrec *cr);
/* Broadcast nbytes bytes from the master to cr->mpi_comm_mygroup */

void gmx_bcast_sim(int nbytes,void *b,const t_commrec *cr);
/* Broadcast nbytes bytes from the sim master to cr->mpi_comm_mysim */

void gmx_sumi(int nr,int r[],const t_commrec *cr);
/* Calculate the global sum of an array of ints */

void gmx_sumf(int nr,float r[],const t_commrec *cr);
/* Calculate the global sum of an array of floats */

void gmx_sumd(int nr,double r[],const t_commrec *cr);
/* Calculate the global sum of an array of doubles */

void gmx_sumf_comm(int nr,float r[],MPI_Comm mpi_comm);
/* Calculate the global sum of an array of floats */

void gmx_sumd_comm(int nr,double r[],MPI_Comm mpi_comm);
/* Calculate the global sum of an array of doubles */

void gmx_sumi_sim(int nr,int r[],const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of ints */

void gmx_sumf_sim(int nr,float r[],const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of floats */

void gmx_sumd_sim(int nr,double r[],const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of doubles */

void gmx_abort(int nodeid,int nnodes,int errorno);
/* Abort the parallel run */

void gmx_finalize(void);

/* Finish the parallel run in an ordered manner */

#ifdef GMX_DOUBLE
#define gmx_sum_comm  gmx_sumd_comm
#define gmx_sum       gmx_sumd
#define gmx_sum_sim   gmx_sumd_sim
#else
#define gmx_sum_comm  gmx_sumf_comm
#define gmx_sum       gmx_sumf
#define gmx_sum_sim   gmx_sumf_sim
#endif

#ifdef DEBUG_GMX
#define debug_gmx() do { FILE *fp=debug ? debug : stderr;\
if (bDebugMode()) fprintf(fp,"NODEID=%d, %s  %d\n",gmx_mpi_initialized() ? gmx_node_rank() : -1,__FILE__,__LINE__); fflush(fp); } while (0)
#else
#define debug_gmx()
#endif

#ifdef __cplusplus
}
#endif


#endif	/* _network_h */
