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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _sim_util_h
#define _sim_util_h

static char *SRCID_sim_util_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) sim_util.h 1.20 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "network.h"
#include "nsb.h"
#include "mshift.h"
#include "force.h"
#include "tgroup.h"
#include "time.h"

extern void print_time(FILE *out,time_t start,int step,t_inputrec *ir);

extern time_t print_date_and_time(FILE *log,int pid,char *title);

extern void finish_run(FILE *log,t_commrec *cr,
		       char *confout,char *traj,
		       char *xtc_traj,
		       t_nsborder *nsb,
		       t_topology *top,
		       t_parm *parm,
		       real t,real lambda,
		       rvec x[],rvec v[],rvec f[],
		       t_nrnb nrnb[],
		       double dt,int step,
		       bool bWriteStat);

extern void nstop_cm(FILE *log,t_commrec *cr,
		     int start,int nr_atoms,real mass[],rvec x[],rvec v[]);

extern void calc_dispcorr(FILE *log,bool bDispCorr,t_forcerec *fr,int natoms,
			  matrix box,tensor pres,tensor virial,real ener[]);
     
#endif	/* _sim_util_h */







