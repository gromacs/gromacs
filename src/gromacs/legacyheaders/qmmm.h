/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _QMMM_h
#define _QMMM_h

#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void atomic_number(int nr, char ***atomtype, int *nucnum);

t_QMMMrec *mk_QMMMrec(void);
/* allocates memory for QMMMrec */

void init_QMMMrec(t_commrec  *cr,
                  gmx_mtop_t *mtop,
                  t_inputrec *ir,
                  t_forcerec *fr);

/* init_QMMMrec initializes the QMMM record. From
 * topology->atoms.atomname and topology->atoms.atomtype the atom
 * names and types are read; from inputrec->QMcharge
 * resp. inputrec->QMmult the nelecs and multiplicity are determined
 * and md->cQMMM gives numbers of the MM and QM atoms
 */

void update_QMMMrec(t_commrec      *cr,
                    t_forcerec     *fr,
                    rvec            x[],
                    t_mdatoms      *md,
                    matrix          box,
                    gmx_localtop_t *top);

/* update_QMMMrec fills the MM stuff in QMMMrec. The MM atoms are
 * taken froom the neighbourlists of the QM atoms. In a QMMM run this
 * routine should be called at every step, since it updates the MM
 * elements of the t_QMMMrec struct.
 */

real calculate_QMMM(t_commrec *cr,
                    rvec x[], rvec f[],
                    t_forcerec *fr);

/* QMMM computes the QM forces. This routine makes either function
 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
 * (ab initio)) or generates input files for an external QM package
 * (listed in QMMMrec.QMpackage). The binary of the QM package is
 * called by system().
 */

#ifdef __cplusplus
}
#endif

#endif  /* _QMMM_h */
