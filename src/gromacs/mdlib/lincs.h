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

#ifndef GMX_MDLIB_LINCS_H
#define GMX_MDLIB_LINCS_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_multisim_t;
struct t_blocka;
struct t_commrec;
struct t_idef;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct t_pbc;

/* Abstract type for LINCS that is defined only in the file that uses it */
typedef struct gmx_lincsdata *gmx_lincsdata_t;

real *lincs_rmsd_data(gmx_lincsdata_t lincsd);
/* Return the data for determining constraint RMS relative deviations */

real lincs_rmsd(gmx_lincsdata_t lincsd);
/* Return the RMSD of the constraint */

gmx_lincsdata_t init_lincs(FILE *fplog, const gmx_mtop_t *mtop,
                           int nflexcon_global, const t_blocka *at2con,
                           gmx_bool bPLINCS, int nIter, int nProjOrder);
/* Initializes and returns the lincs data struct */

void set_lincs(const t_idef *idef, const t_mdatoms *md,
               gmx_bool bDynamics, const t_commrec *cr,
               gmx_lincsdata_t li);
/* Initialize lincs stuff */

gmx_bool
constrain_lincs(FILE *log, gmx_bool bLog, gmx_bool bEner,
                const t_inputrec *ir,
                gmx_int64_t step,
                gmx_lincsdata_t lincsd, t_mdatoms *md,
                const t_commrec *cr,
                const gmx_multisim_t *ms,
                rvec *x, rvec *xprime, rvec *min_proj,
                matrix box, t_pbc *pbc,
                real lambda, real *dvdlambda,
                real invdt, rvec *v,
                gmx_bool bCalcVir, tensor vir_r_m_dr,
                int econ,
                t_nrnb *nrnb,
                int maxwarn, int *warncount);
/* Returns if the constraining succeeded */

#endif
