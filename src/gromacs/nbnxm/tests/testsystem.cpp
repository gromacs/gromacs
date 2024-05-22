/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * \brief Implements the TestSystem class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "testsystem.h"

#include <numeric>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/listoflists.h"

#include "spc81_coords.h"

namespace gmx
{

namespace test
{

namespace
{

// A 3-site water model
//! The number of atoms in a molecule
constexpr int numAtomsInMolecule = 3;
//! The atom type of the oxygen atom
constexpr int typeO = 0;
//! The atom type of a hydrogen atom with LJ
constexpr int typeHWithLJ = 1;
//! The atom type of a hydrogen atom without LJ
constexpr int typeHWithoutLJ = 2;
//! The charge of the oxygen atom
constexpr real chargeO = -0.8476;
//! The charge of the hydrogen atom
constexpr real chargeH = 0.4238;
//! The LJ sigma parameter of the Oxygen atom
constexpr real sigmaO = 0.316557;
//! The LJ epsilon parameter of the Oxygen atom
constexpr real epsilonO = 0.650194;
//! The LJ sigma parameter of Hydrogen atoms with LJ
constexpr real sigmaH = 0.04;
//! The LJ epsilon parameter Hydrogen atoms with LJ
constexpr real epsilonH = 0.192464;

//! Generate a C6, C12 pair using the combination rule
std::pair<real, real> combineLJParams(const real              sigma0,
                                      const real              epsilon0,
                                      const real              sigma1,
                                      const real              epsilon1,
                                      const LJCombinationRule ljCombinationRule)
{
    real sigma6;
    if (ljCombinationRule == LJCombinationRule::Geometric)
    {
        sigma6 = std::pow(sigma0 * sigma1, 3);
    }
    else
    {
        sigma6 = std::pow(0.5 * (sigma0 + sigma1), 6);
    }
    real c6  = 4 * sqrt(epsilon0 * epsilon1) * sigma6;
    real c12 = c6 * sigma6;

    return { c6, c12 };
}

} // namespace

real TestSystem::maxCharge()
{
    return std::abs(chargeO);
}

TestSystem::TestSystem(const LJCombinationRule ljCombinationRule)
{
    numAtomTypes = 3;
    nonbondedParameters.resize(numAtomTypes * numAtomTypes * 2, 0);
    std::tie(nonbondedParameters[0], nonbondedParameters[1]) =
            combineLJParams(sigmaO, epsilonO, sigmaO, epsilonO, ljCombinationRule);
    std::tie(nonbondedParameters[8], nonbondedParameters[9]) =
            combineLJParams(sigmaH, epsilonH, sigmaH, epsilonH, ljCombinationRule);
    std::tie(nonbondedParameters[2], nonbondedParameters[3]) =
            combineLJParams(sigmaO, epsilonO, sigmaH, epsilonH, ljCombinationRule);
    nonbondedParameters[6] = nonbondedParameters[2];
    nonbondedParameters[7] = nonbondedParameters[3];

    coordinates = spc81Coordinates;
    copy_mat(spc81Box, box);
    put_atoms_in_box(PbcType::Xyz, box, coordinates);

    const int numAtoms = coordinates.size();
    GMX_RELEASE_ASSERT(numAtoms % (3 * numAtomsInMolecule) == 0,
                       "Coordinates should be a multiple of 3 x whole water molecules");

    atomTypes.resize(numAtoms);
    charges.resize(numAtoms);
    atomInfo.resize(numAtoms);

    for (int a = 0; a < numAtoms; a++)
    {
        // The first third of the atoms has no charge to cover all code paths
        const bool hasCharge = (a >= numAtoms / 3);

        if (a % numAtomsInMolecule == 0)
        {
            // Oxgygen
            atomTypes[a] = typeO;
            charges[a]   = hasCharge ? chargeO : 0;
            atomInfo[a] |= gmx::sc_atomInfo_HasVdw;
        }
        else
        {
            // Hydrogen
            // Make the last third of molecules have LJ on all atoms
            if (a >= numAtoms * 2 / 3)
            {
                atomTypes[a] = typeHWithLJ;
                atomInfo[a] |= gmx::sc_atomInfo_HasVdw;
            }
            else
            {
                atomTypes[a] = typeHWithoutLJ;
            }
            charges[a] = hasCharge ? chargeH : 0;
        }
        if (hasCharge)
        {
            atomInfo[a] |= gmx::sc_atomInfo_HasCharge;
        }

        // Set the energy group from 0 to n-1
        atomInfo[a] |= (a / (numAtoms / sc_numEnergyGroups));

        // Generate the exclusions like for water molecules
        excls.pushBackListOfSize(numAtomsInMolecule);
        gmx::ArrayRef<int> exclusionsForAtom   = excls.back();
        const int          firstAtomInMolecule = a - (a % numAtomsInMolecule);
        std::iota(exclusionsForAtom.begin(), exclusionsForAtom.end(), firstAtomInMolecule);
    }
}

} // namespace test

} // namespace gmx
