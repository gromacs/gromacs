/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _sim_util_h
#define _sim_util_h

static char *SRCID_sim_util_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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







