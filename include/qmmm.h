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

#ifndef _QMMM_h
#define _QMMM_h

#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"

#ifdef __cplusplus
extern "C" {
#endif

void atomic_number(int nr, char ***atomtype, int *nucnum);

t_QMMMrec *mk_QMMMrec(void);
/* allocates memory for QMMMrec */

void init_QMMMrec(t_commrec *cr,
			 matrix box,
			 gmx_mtop_t *mtop,
			 t_inputrec *ir,
			 t_forcerec *fr);

/* init_QMMMrec initializes the QMMM record. From
 * topology->atoms.atomname and topology->atoms.atomtype the atom
 * names and types are read; from inputrec->QMcharge
 * resp. inputrec->QMmult the nelecs and multiplicity are determined
 * and md->cQMMM gives numbers of the MM and QM atoms 
 */

void update_QMMMrec(t_commrec *cr,
			   t_forcerec *fr,
			   rvec x[],
			   t_mdatoms *md,
			   matrix box,
			   gmx_localtop_t *top);

/* update_QMMMrec fills the MM stuff in QMMMrec. The MM atoms are
 * taken froom the neighbourlists of the QM atoms. In a QMMM run this
 * routine should be called at every step, since it updates the MM
 * elements of the t_QMMMrec struct.  
 */

real calculate_QMMM(t_commrec *cr,
			   rvec x[], rvec f[],
			   t_forcerec *fr,
			   t_mdatoms *md);

/* QMMM computes the QM forces. This routine makes either function
 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
 * (ab initio)) or generates input files for an external QM package
 * (listed in QMMMrec.QMpackage). The binary of the QM package is
 * called by system().
 */

#ifdef __cplusplus
}
#endif

#endif	/* _QMMM_h */

