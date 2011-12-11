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

#ifndef _orires_h
#define _orires_h

#include "sysstuff.h"
#include "typedefs.h"

#ifdef __cplusplus
    extern "C" {
#endif

void init_orires(FILE *fplog,const gmx_mtop_t *mtop,
			rvec x[],
			const t_inputrec *ir,
			const gmx_multisim_t *ms,t_oriresdata *od,
			t_state *state);
/* Initializes all the orientation restraint stuff in *od */

real calc_orires_dev(const gmx_multisim_t *ms,
			    int nfa,const t_iatom fa[],const t_iparams ip[],
			    const t_mdatoms *md,const rvec x[],
			    const t_pbc *pbc,t_fcdata *fcd,history_t *hist);
/* 
 * Calculates the time averaged D matrices, the S matrix for each experiment.
 * Returns the weighted RMS deviation of the orientation restraints.
 */

void diagonalize_orires_tensors(t_oriresdata *od);
/*
 * Diagonalizes the order tensor(s) of the orienation restraints.
 * For each experiment eig containts first 3 eigenvalues and then
 * the 3 eigenvectors. The eigenvalues are ordered on magnitude.
 */

void print_orires_log(FILE *log,t_oriresdata *od);
/* Print order parameter, eigenvalues and eigenvectors to the log file */

t_ifunc orires;
/* Does only the orientation restraint force calculation */

void update_orires_history(t_fcdata *fcd,history_t *hist);
/* Copy the new time averages that have been calculated in calc_orires_dev */

#ifdef __cplusplus
}
#endif

#endif	/* _orires_h */
