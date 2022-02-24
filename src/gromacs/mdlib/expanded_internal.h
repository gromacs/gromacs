/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \internal \file
 * \brief Declares internal functionality for expanded ensemble
 *
 * This file is only used by expanded.cpp and tests/expanded.cpp.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \author Michael Shirts <michael.shirts@colorado.edu>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_EXPANDEDINTERNAL_H
#define GMX_MDLIB_EXPANDEDINTERNAL_H

#include "gromacs/utility/real.h"

enum class LambdaWeightCalculation : int;

namespace gmx
{
/*! \brief Calculates the acceptance weight for a lambda state transition
 *
 * \param[in] calculationMode  How the lambda weights are calculated
 * \param[in] lambdaEnergyDifference  The difference in energy between the two states
 * \return  The acceptance weight
 */
real calculateAcceptanceWeight(LambdaWeightCalculation calculationMode, real lambdaEnergyDifference);
} // namespace gmx

#endif // GMX_MDLIB_EXPANDEDINTERNAL_H
