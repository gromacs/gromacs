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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _init_sh_h
#define _init_sh_h

static char *SRCID_init_sh_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "mdebin.h"
	
typedef struct {
  int     nnucl;
  atom_id shell;	        /* The shell id				*/
  atom_id nucl1,nucl2,nucl3;	/* The nuclei connected to the shell	*/
  real    k;		        /* force constant		        */
  real    k_1;		        /* 1 over force constant		*/
} t_shell;

extern t_shell *init_shells(FILE *log,int start,int homenr,
			    t_idef *idef,t_mdatoms *md,int *nshell);

extern int relax_shells(FILE *ene,FILE *log,t_commrec *cr,bool bVerbose,
			int mdstep,t_parm *parm,bool bDoNS,bool bStopCM,
			t_topology *top,real ener[],
			rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
			rvec buf[],t_mdatoms *md,t_nsborder *nsb,t_nrnb *nrnb,
			t_graph *graph,t_groups *grps,tensor vir_part,
			int nshell,t_shell shells[],t_forcerec *fr,
			char *traj,real t,real lambda,
			int natoms,matrix box,t_mdebin *mdebin);

extern  int relax_shells2(FILE *ene,FILE *log,t_commrec *cr,
			  bool bVerbose,int mdstep,
			  t_parm *parm,bool bDoNS,bool bStopCM,
			  t_topology *top,real ener[],
			  rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
			  rvec buf[],t_mdatoms *md,
			  t_nsborder *nsb,t_nrnb *nrnb,
			  t_graph *graph,
			  t_groups *grps,tensor vir_part,
			  int nshell,t_shell shells[],
			  t_forcerec *fr,char *traj,
			  real t,real lambda,
			  int natoms,matrix box,t_mdebin *mdebin);
/* Experimental version with line min */

#endif
