/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
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

void finish_run(FILE *log,t_commrec *cr,char *confout, t_nsborder *nsb,
		t_topology *top, t_parm *parm,t_nrnb nrnb[],double nodetime,
		double realtime,int step,bool bWriteStat);


void calc_dispcorr(FILE *log,int eDispCorr,t_forcerec *fr,int natoms,
		   matrix box,tensor pres,tensor virial,real ener[]);
     
#endif	/* _sim_util_h */







