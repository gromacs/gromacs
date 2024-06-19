/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_MDLIB_ENERDATA_UTILS_H
#define GMX_MDLIB_ENERDATA_UTILS_H

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;
struct gmx_grppairener_t;
struct t_lambda;

void reset_dvdl_enerdata(gmx_enerdata_t* enerd);
/* Resets only the dvdl energy data */

void reset_enerdata(gmx_enerdata_t* enerd);
/* Resets the energy data */

/*! \brief Sums energy group pair contributions into epot */
void sum_epot(const gmx_grppairener_t& grpp, real* epot);

/*! \brief Accumulates potential energy contributions to obtain final potential energies
 *
 * Accumulates energy group pair contributions into the output energy components
 * and sums all potential energies into the total potential energy term.
 * With free-energy also computes the foreign lambda potential energy differences.
 *
 * \param[in,out] enerd    Energy data with components to sum and to accumulate into
 * \param[in]     lambda   Vector of free-energy lambdas
 * \param[in]     fepvals  Foreign lambda energy differences, only summed with !=nullptr
 */
void accumulatePotentialEnergies(gmx_enerdata_t*           enerd,
                                 gmx::ArrayRef<const real> lambda,
                                 const t_lambda*           fepvals);

/*! \brief Accumulates kinetic and constraint contributions to dH/dlambda and foreign energies */
void accumulateKineticLambdaComponents(gmx_enerdata_t*           enerd,
                                       gmx::ArrayRef<const real> lambda,
                                       const t_lambda&           fepvals);

#endif
