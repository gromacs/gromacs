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
 * GROningen MAchine for Chemical Simulation
 */

#ifndef	_do_md_h
#define	_do_md_h

#ifdef HAVE_IDENT
#ident	"@(#) do_md.h 1.16 2/2/97"
#endif /* HAVE_IDENT */
#include <stdio.h>
#include "typedefs.h"
#include "network.h"
#include "tgroup.h"
#include "stat.h"

extern time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
		    bool bVerbose,bool bCompact,int stepout,
		    t_parm *parm,t_groups *grps,
		    t_topology *top,real ener[],
		    rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		    rvec buf[],t_mdatoms *mdatoms,
		    t_nsborder *nsb,t_nrnb nrnb[],
		    t_graph *graph);

#endif	/* _do_md_h */
