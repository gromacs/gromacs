/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef	_stat_h
#define	_stat_h

#ifdef HAVE_IDENT
#ident	"@(#) stat.h 1.39 8/4/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "force.h"
#include "nsb.h"

extern void global_stat(FILE *log,
			t_commrec *cr,real ener[],
			tensor fvir,tensor svir,
			t_grpopts *opts,t_groups *grps,
			t_nrnb *mynrnb,t_nrnb nrnb[],
			rvec vcm);
/* Communicate statistics around the ring */

extern void write_traj(FILE *log,t_commrec *cr,char *traj,t_nsborder *nsb,
		       int step,real t,real lambda,t_nrnb nr_nb[],
		       int natoms,rvec *xx,rvec *vv,rvec *ff,matrix box);
/* Routine to output statusfiles during a run, as specified in
 * in parm->ir. If any of the pointers xx,vv,ff or ener is not NULL
 * it is written to the trajectory file.
 * Also write the energies etc. to the log file.
 */

extern int do_per_step(int step,int nstep);
/* Return TRUE if io should be done */

extern int do_any_io(int step, t_inputrec *ir);

extern void write_xtc_traj(FILE *log,t_commrec *cr,
			   char *xtc_traj,t_nsborder *nsb,t_mdatoms *md,
			   int step,real t,rvec *xx,
			   matrix box,real prec);

extern void close_xtc_traj(void);

#endif	/* _stat_h */
