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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _edsam_h
#define _edsam_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void do_edsam(t_inputrec *ir,gmx_large_int_t step,t_mdatoms *md,
                     t_commrec *cr,rvec xs[],rvec v[],matrix box,gmx_edsam_t ed);
/* Essential dynamics constraints, called from constrain() */

gmx_edsam_t ed_open(int nfile,const t_filenm fnm[],unsigned long Flags,t_commrec *cr);
/* Sets the ED input/output filenames, opens output (.edo) file */

void init_edsam(gmx_mtop_t *mtop,t_inputrec *ir,t_commrec *cr,
                       gmx_edsam_t ed, rvec x[], matrix box);
/* Init routine for ED and flooding. Calls init_edi in a loop for every .edi-cycle 
 * contained in the input file, creates a NULL terminated list of t_edpar structures */

void dd_make_local_ed_indices(gmx_domdec_t *dd, gmx_edsam_t ed);
/* Make a selection of the home atoms for the ED groups. 
 * Should be called at every domain decomposition. */
 
void do_flood(FILE *log, t_commrec *cr, rvec x[],rvec force[], gmx_edsam_t ed,
        matrix box, gmx_large_int_t step);
/* Flooding - called from do_force() */

#ifdef __cplusplus
}
#endif

#endif	/* _edsam_h */
