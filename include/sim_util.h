/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * S  C  A  M  O  R  G
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

extern void do_force(FILE *log,t_commrec *cr,
		     t_parm *parm,t_nsborder *nsb,tensor vir_part,
		     int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
		     rvec x[],rvec v[],rvec f[],rvec buf[],
		     t_mdatoms *mdatoms,real ener[],bool bVerbose,
		     t_graph *graph,
		     bool bNS,bool bMolEpot,t_forcerec *fr);
		     
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
     
#endif	/* _sim_util_h */







