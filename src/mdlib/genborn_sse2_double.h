/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
#ifndef _genborn_sse2_double_h
#define _genborn_sse2_double_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "grompp.h"


int 
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_localtop_t *top,
							  const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born);

int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda, double *xd, double *f, 
							  double *fshift, double *shift_vec, int gb_algorithm, 
							  gmx_genborn_t *born, t_mdatoms *md);	

int 
calc_gb_rad_hct_obc_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top, 
							const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md, int gb_algorithm);



#endif /* _genborn_sse2_double_h */
