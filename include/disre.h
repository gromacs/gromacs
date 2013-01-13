
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

#ifndef _disre_h
#define _disre_h

#include "sysstuff.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_disres(FILE *fplog,const gmx_mtop_t *mtop,
                 t_inputrec *ir,const t_commrec *cr,gmx_bool bPartDecomp,
                 t_fcdata *fcd,t_state *state, gmx_bool bIsREMD);
/* Initiate *fcd data, must be called once, nbonds is the number 
 * of iatoms in the ilist of the idef struct.
 * When time averaging is used, the history is initialized in state,
 * unless it was read before from a checkpoint file.
 * The implementation of distance restraints with -multi
 * must differ according to whether REMD is active.
 */

void calc_disres_R_6(const gmx_multisim_t *ms,
			    int nfa,const t_iatom *fa,const t_iparams ip[],
			    const rvec *x,const t_pbc *pbc,
			    t_fcdata *fcd,history_t *hist);
/* Calculates r and r^-3 (inst. and time averaged) for all pairs
 * and the ensemble averaged r^-6 (inst. and time averaged) for all restraints
 */

t_ifunc ta_disres;
/* Calculate the distance restraint forces, return the potential */

void update_disres_history(t_fcdata *fcd,history_t *hist);
/* Copy the new time averages that have been calculated in calc_disres_R_6 */

#ifdef __cplusplus
}
#endif

#endif	/* _disre_h */
