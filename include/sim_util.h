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

extern void calc_dispcorr(FILE *log,int eDispCorr,t_forcerec *fr,int natoms,
			  matrix box,tensor pres,tensor virial,real ener[]);
     
#endif	/* _sim_util_h */







