/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#ifndef GMX_MDLIB_SETTLE_H
#define GMX_MDLIB_SETTLE_H

#include "gromacs/topology/idef.h"

struct gmx_cmap_t;
struct gmx_mtop_t;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;

/* Abstract type for SETTLE that is defined only in the file that uses it */
typedef struct gmx_settledata *gmx_settledata_t;

gmx_settledata_t settle_init(const gmx_mtop_t *mtop);
/* Initializes and returns a structure with SETTLE parameters */

void settle_free(gmx_settledata_t settled);

void settle_set_constraints(gmx_settledata_t  settled,
                            const t_ilist    *il_settle,
                            const t_mdatoms  *mdatoms);
/* Set up the indices for the settle constraints */

void csettle(gmx_settledata_t    settled,          /* The SETTLE structure */
             int                 nthread,          /* The number of threads used */
             int                 thread,           /* Our thread index */
             const t_pbc        *pbc,              /* PBC data pointer, can be NULL */
             const real          x[],              /* Reference coordinates */
             real                xprime[],         /* New coords, to be settled */
             real                invdt,            /* 1/delta_t */
             real               *v,                /* Also constrain v if v!=NULL */
             bool                bCalcVirial,      /* Calculate the virial contribution */
             tensor              vir_r_m_dr,       /* sum r x m delta_r */
             bool               *bErrorHasOccurred /* True if a settle error occurred */
             );
/* Constrain coordinates using SETTLE.
 * Can be called on any number of threads.
 */

void settle_proj(gmx_settledata_t settled, int econq,
                 int nsettle, t_iatom iatoms[],
                 const t_pbc *pbc,   /* PBC data pointer, can be NULL  */
                 rvec x[],
                 rvec *der, rvec *derp,
                 int CalcVirAtomEnd, tensor vir_r_m_dder);
/* Analytical algorithm to subtract the components of derivatives
 * of coordinates working on settle type constraint.
 */

#endif
