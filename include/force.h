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

static char *sepdvdlformat="  %-30s V %12.5e  dVdl %12.5e\n";

extern void calc_vir(FILE *fplog,int nxf,rvec x[],rvec f[],tensor vir);
/* Calculate virial for nxf atoms, and add it to vir */

extern void f_calc_vir(FILE *fplog,int i0,int i1,rvec x[],rvec f[],tensor vir,
		       t_graph *g,rvec shift_vec[]);
/* Calculate virial taking periodicity into account */

extern real RF_excl_correction(FILE *fplog,const t_nsborder *nsb,
			       const t_forcerec *fr,t_graph *g,
			       const t_mdatoms *mdatoms,const t_block *excl,
			       rvec x[],rvec f[],rvec *fshift,const t_pbc *pbc,
			       real lambda,real *dvdlambda);
/* Calculate the reaction-field energy correction for this node:
 * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
 * and force correction for all excluded pairs, including self pairs.
 */

extern void calc_rffac(FILE *fplog,int eel,real eps_r,real eps_rf,
		       real Rc,real Temp,
		       real zsq,matrix box,
		       real *kappa,real *krf,real *crf);
/* Determine the reaction-field constants */

extern t_forcerec *mk_forcerec(void);

extern t_forcetable make_tables(FILE *fp,const t_forcerec *fr,
				bool bVerbose,const char *fn,
				real rtab,bool b14only);
/* Return tables for inner loops. When bVerbose the tables are printed
 * to .xvg files
 */
 
extern void pr_forcerec(FILE *fplog,t_forcerec *fr,t_commrec *cr);

extern void init_forcerec(FILE       *fplog,     
			  t_forcerec *fr,   
			  const t_inputrec *ir,   
			  const t_topology *top,
			  const t_commrec  *cr,
			  const t_mdatoms  *mdatoms,
			  const t_nsborder *nsb,
			  matrix     box,
			  bool       bMolEpot,
			  const char *tabfn,
			  const char *tabpfn,
			  bool       bNoSolvOpt);
/* The Force rec struct must be created with mk_forcerec 
 * The booleans have the following meaning:
 * bSetQ:    Copy the charges [ only necessary when they change ]
 * bMolEpot: Use the free energy stuff per molecule
 */
 
extern void update_forcerec(FILE *fplog,t_forcerec *fr,matrix box);
/* Updates parameters in the forcerec that are time dependent */

/* Compute the average C6 and C12 params for LJ corrections */
extern void set_avcsixtwelve(FILE *fplog,t_forcerec *fr,
			     const t_mdatoms *mdatoms,const t_block *excl);

extern void ns(FILE       *fplog,
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
	       real       *dvdlambda,
	       bool       bFillGrid,
	       bool       bDoForces);
/* Call the neighborsearcher */

extern void force(FILE         *fplog,  
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
		  bool         bDoForces,
		  rvec         mu_tot[2],
		  bool         bGatherOnly,
		  t_edsamyn      *edyn);
/* Call all the force routines */


#endif	/* _force_h */
