/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _vsite_parm_h
#define _vsite_parm_h

#include "typedefs.h"
#include "grompp.h"
#include "gpp_atomtype.h"

extern int set_vsites(gmx_bool bVerbose, t_atoms *atoms,  gpp_atomtype_t atype,
		      t_params plist[]);
/* set parameters for vritual sites, return number of virtual sites */

extern void set_vsites_ptype(gmx_bool bVerbose,  gmx_moltype_t *molt);
/* set ptype to VSite for virtual sites */

extern void clean_vsite_bondeds(t_params *ps, int natoms, gmx_bool bRmVSiteBds);
/* remove all bonded interaction (bonds, angles and diherals) that
   have become obsolete due to virtual site constructions */

#endif	/* _vsite_parm_h */
