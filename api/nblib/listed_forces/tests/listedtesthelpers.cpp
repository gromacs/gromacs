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

#include "nblib/listed_forces/conversionscommon.h"
#include "nblib/vector.h"

#include "gmxcalculator.h"

namespace nblib
{

void compareNblibAndGmxListedImplementations(const ListedInteractionData& interactionData,
                                             const std::vector<Vec3>&     coordinates,
                                             const std::vector<real>&     charges,
                                             size_t                       numParticles,
                                             int                          numThreads,
                                             const Box&                   box,
                                             Vec3                         com,
                                             real                         tolerance)
{
    ListedForceCalculator calculator(
            interactionData, numParticles, numThreads, box, PbcType::Xyz, RefCoordScaling::No);

    std::vector<Vec3> forces(numParticles, Vec3{ 0, 0, 0 });
    std::vector<Vec3> shiftForces(gmx::c_numShiftVectors, Vec3{ 0, 0, 0 });
    std::vector<real> virial(9, 0);
    ListedEnergies    energies;

    if (!charges.empty())
    {
        calculator.compute(coordinates, charges, forces, shiftForces, virial, energies, com);
    }
    else
    {
        std::vector<real> zeroCharges(numParticles, 0.0);
        calculator.compute(coordinates, zeroCharges, forces, shiftForces, virial, energies, com);
    }

    auto [idef, ffparams] = convertToGmxInteractions(interactionData);
    ListedGmxCalculator gmxCalculator(*idef, *ffparams, numParticles, numThreads, box);

    std::vector<Vec3> gmxForces(numParticles, Vec3{ 0, 0, 0 });
    std::vector<Vec3> gmxShiftForces(gmx::c_numShiftVectors, Vec3{ 0, 0, 0 });
    ListedEnergies    gmxEnergies;
    std::vector<real> gmxVirial(9, 0);

    gmxCalculator.compute(coordinates, charges, gmxForces, gmxShiftForces, gmxVirial, gmxEnergies, com);

    gmx::test::FloatingPointTolerance tolSetting(tolerance, tolerance, 1.0e-5, 1.0e-7, 200, 100, false);

    compareVectors(forces.begin(), forces.end(), gmxForces.begin(), tolSetting);
    compareVectors(shiftForces.begin(), shiftForces.end(), gmxShiftForces.begin(), tolSetting);
    compareEnergies(energies.begin(), energies.end(), gmxEnergies.begin(), tolSetting);

    for (std::size_t i = 0; i < virial.size(); ++i)
    {
        EXPECT_REAL_EQ_TOL(virial[i], gmxVirial[i], tolSetting);
    }
}

} // namespace nblib
