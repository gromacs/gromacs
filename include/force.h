/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _force_h
#define _force_h

static char *SRCID_force_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) force.h 1.53 5/15/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "pbc.h"
#include "nsb.h"
#include "network.h"
#include "tgroup.h"

extern void calc_vir(FILE *log,int nxf,rvec x[],rvec f[],tensor vir,
		     t_commrec *cr);
/* Calculate virial for nxf atoms, and add it to vir */

extern void f_calc_vir(FILE *log,int i0,int i1,rvec x[],rvec f[],tensor vir,
		       t_commrec *cr,t_graph *g,rvec shift_vec[]);
/* Calculate virial taking periodicity into account */

extern t_forcerec *mk_forcerec(void);

extern void make_tables(t_forcerec *fr,bool bVerbose);
/* Make tables for inner loops. When bVerbose the tables are printed
 * to .xvg files
 */
 
extern void pr_forcerec(FILE *log,t_forcerec *fr,t_commrec *cr);

extern void init_forcerec(FILE       *log,     
			  t_forcerec *fr,   
			  t_inputrec *ir,   
			  t_block    *mols,
			  t_commrec  *cr,
			  t_block    *cgs,
			  t_idef     *idef,
			  t_mdatoms  *mdatoms,
			  matrix     box,
			  bool       bMolEpot);
/* The Force rec struct must be created with mk_forcerec 
 * The booleans have the following meaning:
 * bSetQ:    Copy the charges [ only necessary when they change ]
 * bMolEpot: Use the free energy stuff per molecule
 */
 
extern void update_forcerec(FILE *log,t_forcerec *fr,matrix box);
/* Updates parameters in the forcerec that are time dependent */

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
	       int        step);
/* Call the neighborsearcher */

extern void force(FILE *log,  
		  int          step,
		  t_forcerec   *fr,
		  t_idef       *idef,
		  t_nsborder   *nsb,
		  t_commrec    *cr,
		  t_nrnb       *nrnb,
		  t_groups     *grps,
		  t_mdatoms    *md,
		  int          ngener,
		  t_grpopts    *opts,
		  rvec x[],
		  rvec f[],    tensor virial,
		  real epot[], bool bVerbose,
		  matrix       box,
		  real         lambda,
		  t_graph      *graph,
		  t_block      *excl);
/* Call all the force routines */

#endif	/* _force_h */
