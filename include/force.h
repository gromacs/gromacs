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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _force_h
#define _force_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "typedefs.h"
#include "pbc.h"
#include "nsb.h"
#include "network.h"
#include "tgroup.h"

extern void calc_vir(FILE *log,int nxf,rvec x[],rvec f[],tensor vir);
/* Calculate virial for nxf atoms, and add it to vir */

extern void f_calc_vir(FILE *log,int i0,int i1,rvec x[],rvec f[],tensor vir,
		       t_graph *g,rvec shift_vec[]);
/* Calculate virial taking periodicity into account */

extern t_forcerec *mk_forcerec(void);

extern void make_tables(FILE *fp,t_forcerec *fr,bool bVerbose,char *fn);
/* Make tables for inner loops. When bVerbose the tables are printed
 * to .xvg files
 */
 
extern void pr_forcerec(FILE *log,t_forcerec *fr,t_commrec *cr);

extern void init_forcerec(FILE       *log,     
			  t_forcerec *fr,   
			  t_inputrec *ir,   
			  t_topology *top,
			  t_commrec  *cr,
			  t_mdatoms  *mdatoms,
			  t_nsborder *nsb,
			  matrix     box,
			  bool       bMolEpot,
			  char       *tabfn,
			  bool       bNoSolvOpt);
/* The Force rec struct must be created with mk_forcerec 
 * The booleans have the following meaning:
 * bSetQ:    Copy the charges [ only necessary when they change ]
 * bMolEpot: Use the free energy stuff per molecule
 */
 
extern void update_forcerec(FILE *log,t_forcerec *fr,matrix box);
/* Updates parameters in the forcerec that are time dependent */

/* Compute the average C6 and C12 params for LJ corrections */
extern void set_avcsixtwelve(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms,
			     t_block *excl);

extern void ns(FILE *log,
	       t_forcerec *fr,
	       rvec       x[],
	       rvec       f[],
	       matrix     box,
	       t_groups   *grps,
	       t_grpopts  *opts,
	       t_topology *top,
	       t_mdatoms  *md,
	       t_commrec  *cr,
	       t_nrnb     *nrnb,
	       t_nsborder *nsb,
	       int        step,
	       real       lambda,
	       real       *dvdlambda);
/* Call the neighborsearcher */

extern void force(FILE *log,  
		  int          step,
		  t_forcerec   *fr,
		  t_inputrec   *ir,
		  t_idef       *idef,
		  t_nsborder   *nsb,
		  t_commrec    *cr,
		  t_commrec    *mcr,
		  t_nrnb       *nrnb,
		  t_groups     *grps,
		  t_mdatoms    *md,
		  int          ngener,
		  t_grpopts    *opts,
		  rvec         x[],
		  rvec         f[],    
		  real         epot[], 
		  t_fcdata     *fcd,
		  bool         bVerbose,
		  matrix       box,
		  real         lambda,
		  t_graph      *graph,
		  t_block      *excl,
		  bool         bNBonly,
		  matrix       lr_vir,
		  rvec         mu_tot[2],
		  real         qsum[],
		  bool         bGatherOnly);
/* Call all the force routines */

/* Routine from fnbf.m4 */
extern void do_fnbf(FILE *log,t_commrec *cr,t_forcerec *fr,
		    rvec x[],rvec f[],t_mdatoms *md,
		    real egnb[],real egcoul[],rvec box_size,
		    t_nrnb *nrnb,real lambda,real *dvdlambda,
		    bool bLR,int eNL);

#endif	/* _force_h */
