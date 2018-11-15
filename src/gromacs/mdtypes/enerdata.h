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
#ifndef GMX_MDTYPES_TYPES_ENERDATA_H
#define GMX_MDTYPES_TYPES_ENERDATA_H

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/real.h"

enum {
    egCOULSR, egLJSR, egBHAMSR,
    egCOUL14, egLJ14, egNR
};

struct gmx_grppairener_t
{
    int   nener;      /* The number of energy group pairs     */
    real *ener[egNR]; /* Energy terms for each pair of groups */
};

struct gmx_enerdata_t
{
    real                     term[F_NRE];         /* The energies for all different interaction types */
    struct gmx_grppairener_t grpp;
    double                   dvdl_lin[efptNR];    /* Contributions to dvdl with linear lam-dependence */
    double                   dvdl_nonlin[efptNR]; /* Idem, but non-linear dependence                  */
    /* The idea is that dvdl terms with linear lambda dependence will be added
     * automatically to enerpart_lambda. Terms with non-linear lambda dependence
     * should explicitly determine the energies at foreign lambda points
     * when n_lambda > 0. */

    int                      n_lambda;
    int                      fep_state;           /*current fep state -- just for printing */
    double                  *enerpart_lambda;     /* Partial Hamiltonian for lambda and flambda[], includes at least all perturbed terms */
    real                     foreign_term[F_NRE]; /* alternate array for storing foreign lambda energies */
    struct gmx_grppairener_t foreign_grpp;        /* alternate array for storing foreign lambda energies */
};

#endif
