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
#ifndef GMX_MDLIB_RF_UTIL_H
#define GMX_MDLIB_RF_UTIL_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_blocka;
struct t_forcerec;
struct t_graph;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;

real RF_excl_correction(const t_forcerec *fr,
                        const t_graph    *g,
                        const t_mdatoms  *mdatoms,
                        const t_blocka   *excl,
                        bool              usingDomainDecomposition,
                        rvec              x[],
                        rvec              f[],
                        rvec             *fshift,
                        const t_pbc      *pbc,
                        real              lambda,
                        real             *dvdlambda);
/* Calculate the reaction-field energy correction for this node:
 * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
 * and force correction for all excluded pairs, including self pairs.
 */

void calc_rffac(FILE *fplog, int eel, real eps_r, real eps_rf,
                real Rc, real Temp,
                real zsq, matrix box,
                real *krf, real *crf);
/* Determine the reaction-field constants */

void init_generalized_rf(FILE *fplog,
                         const gmx_mtop_t *mtop, const t_inputrec *ir,
                         t_forcerec *fr);
/* Initialize the generalized reaction field parameters */


#endif
