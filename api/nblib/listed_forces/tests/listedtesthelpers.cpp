/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief
 * A collection of helper utilities that allow setting up both Nblib and
 * GROMACS fixtures for computing listed interactions given sets of parameters
 * and coordinates
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "listedtesthelpers.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

#include "nblib/listed_forces/calculator.h"

#include "gmxcalculator.h"

namespace nblib
{

void compareNblibAndGmxListedImplementations(const ListedInteractionData&  interactionData,
                                             const std::vector<gmx::RVec>& coordinates,
                                             size_t                        numParticles,
                                             int                           numThreads,
                                             const Box&                    box,
                                             real                          tolerance)
{
    ListedForceCalculator calculator(interactionData, numParticles, numThreads, box);

    std::vector<gmx::RVec>            forces(numParticles, gmx::RVec{ 0, 0, 0 });
    std::vector<gmx::RVec>            shiftForces(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 });
    ListedForceCalculator::EnergyType energies;

    calculator.compute(coordinates, forces, shiftForces, energies, true);

    ListedGmxCalculator gmxCalculator(interactionData, numParticles, numThreads, box);

    std::vector<gmx::RVec>            gmxForces(numParticles, gmx::RVec{ 0, 0, 0 });
    std::vector<gmx::RVec>            gmxShiftForces(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 });
    ListedForceCalculator::EnergyType gmxEnergies;

    gmxCalculator.compute(coordinates, gmxForces, gmxShiftForces, gmxEnergies, true);

    gmx::test::FloatingPointTolerance tolSetting(tolerance, tolerance, 1.0e-5, 1.0e-8, 200, 100, false);

    EXPECT_THAT(forces, Pointwise(gmx::test::RVecEq(tolSetting), gmxForces));
    EXPECT_THAT(shiftForces, Pointwise(gmx::test::RVecEq(tolSetting), gmxShiftForces));
}

} // namespace nblib
