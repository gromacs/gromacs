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

#ifndef _tgroup_h
#define _tgroup_h

static char *SRCID_tgroup_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) tgroup.h 1.12 2/2/97"
#endif /* HAVE_IDENT */
#ifdef HAVE_IDENT
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "network.h"

extern void init_groups(FILE *log,t_mdatoms *md,
			t_grpopts *opts,t_groups *grps);
/* Allocate memory and set the grpnr array. */

extern void done_groups(t_groups *grps);
/* Free the memory */

extern void accumulate_u(t_commrec *cr,t_grpopts *opts,t_groups *grps);

/*extern void accumulate_ekin(t_commrec *cr,t_grpopts *opts,t_groups *grps);*/
/* Communicate subsystem - group velocities and subsystem ekin respectively
 * and sum them up. Return them in grps.
 */

extern real sum_ekin(t_grpopts *opts,t_groups *grps,tensor ekin,bool bTYZ);
/* Sum the group ekins into total ekin and calc temp per group,
 * return total temperature.
 */
extern void sum_epot(t_grpopts *opts,t_groups *grps,real epot[]);
/* Sum the epot from the group contributions */

extern void update_grps(int start,int homenr,t_groups *grps,
			t_grpopts *opts,rvec v[],t_mdatoms *md,bool bNEMD);
/* Do the update of group velocities (if bNEMD) and
 * (partial) group ekin.
 */

#endif	/* _tgroup_h */
