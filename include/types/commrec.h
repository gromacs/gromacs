/*
 * $Id$
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_MPI
#include <mpi.h>
#endif

#define DD_MAXCELL  8
#define DD_MAXICELL 4

typedef struct {
  /* The global charge group division */
  int  *ncg;     /* Number of home charge groups for each node */
  int  *index;   /* Index of nnodes+1 into cg */
  int  *cg;      /* Global charge group index */
  int  *nat;     /* Number of home atoms for each node. */
  int  *ibuf;    /* Buffer for communication */
} gmx_domdec_master_t;

typedef struct {
  /* The cell (node) we communicate with */
  int cell;
  int ncg;
  int nat;
  /* Index of size ncg into the global charge groups */
  int *index_gl;
  /* Index of size ncg into the local charge groups */
  int *index;
  int nalloc;
} gmx_domdec_comm_t;

typedef struct {
  int  j0;       /* j-cell start               */
  int  j1;       /* j-cell end                 */
  int  cg1;      /* i-charge-group end         */
  int  jcg0;     /* j-charge-group start       */
  int  jcg1;     /* j-charge-group end         */
  ivec shift0;   /* Minimum shifts to consider */
  ivec shift1;   /* Maximum shifts to consider */
} gmx_domdec_ns_t;

typedef struct {
  ivec nc;
  int  nnodes;
#ifdef GMX_MPI
  MPI_Comm all;
#endif

  /* The communication setup, identical for each cell */
  int  ncell;
  ivec shift[DD_MAXCELL];

  int  nodeid;

  /* Only available on the master node */
  gmx_domdec_master_t ma;
  /* Switch that tells if the master has the charge group distribution */
  bool bMasterHasAllCG;

  /* Global atom number to local atom number, -1 if not local */
  int  *ga2la;

  /* The following arrays will have size gl.ncell */
  /* Nodes we need to send coordinates to and receive forces from */
  gmx_domdec_comm_t comm0[DD_MAXCELL];
  /* Nodes we need to receive coordinates from and send forces to */
  gmx_domdec_comm_t comm1[DD_MAXCELL];

  /* Index into the coordinates on this node */
  int  at_index[DD_MAXCELL];

  /* For neighborsearching */
  int  ncg_tot;
  int  nicell;
  gmx_domdec_ns_t icell[DD_MAXICELL];

  /* Communication buffers */
  rvec *buf_x;
  int  nalloc_x;
  rvec *buf_v;
  int  nalloc_v;
  rvec *buf_sdx;
  int  nalloc_sdx;
} gmx_domdec_t;

typedef struct {
  int nodeid,nnodes,npmenodes;
  int left,right;
  int threadid,nthreads;
#ifdef GMX_MPI
  MPI_Comm mpi_comm_mygroup;
#endif
  gmx_domdec_t *dd;
} t_commrec;

#define MASTERNODE(cr)     ((cr)->nodeid == 0)
#define MASTERTHREAD(cr)   ((cr)->threadid == 0)
#define MASTER(cr)         (MASTERNODE(cr) && MASTERTHREAD(cr))
#define NODEPAR(cr)        ((cr)->nnodes > 1)
#define THREADPAR(cr)      ((cr)->nthreads > 1)
#define PAR(cr)            (NODEPAR(cr) || THREADPAR(cr))

/*
 * with the help of this macro + enum we are able to find out 
 * what type of pme work the local node has to do 
 *
 * pmeduty = 0 if node does pp only
 *         = 1 if node does pme only
 *         = 2 if ALL nodes do both (no PME/PP node splitting)
 */ 
#define pmeduty(cr)        ((!cr->npmenodes)? 2:(cr->nodeid >= cr->nnodes-cr->npmenodes))

enum {
  epmePPONLY, epmePMEONLY, epmePMEANDPP
};

#define DDMASTERNODE     0
#define DDMASTER(dd)     ((dd)->nodeid == DDMASTERNODE)
