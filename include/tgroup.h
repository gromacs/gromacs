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
