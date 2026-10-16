/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for the check of the number of perturbed exclusions
 *
 * \author Jan Walther Perthold <jan@ackergarten.at>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/exclusionchecker.h"

#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/mpicomm.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Description of a molecule type and of the block of molecules of that type
struct MoleculeBlockDescription
{
    //! The number of atoms in the molecule
    int numAtoms;
    //! The indices in the molecule of the atoms that are perturbed
    std::vector<int> perturbedAtoms;
    //! The number of molecules in the block
    int numMolecules;
};

/*! \brief Returns a topology with one molecule type and one molecule block for each entry
 * of \p moleculeBlocks
 *
 * All atoms within a molecule exclude each other and themselves.
 * Perturbed atoms have a different atom type in the B-state.
 *
 * \param[in] moleculeBlocks  The molecule types and blocks, in the order of the global atom indices
 * \param[in] exclusionGroup  The global atom indices of the inter-molecular exclusion group
 */
std::unique_ptr<gmx_mtop_t> makeTopology(const std::vector<MoleculeBlockDescription>& moleculeBlocks,
                                         const std::vector<int>& exclusionGroup)
{
    auto mtop = std::make_unique<gmx_mtop_t>();

    mtop->moltype.resize(moleculeBlocks.size());
    for (int moleculeType = 0; moleculeType < gmx::ssize(moleculeBlocks); moleculeType++)
    {
        const MoleculeBlockDescription& moleculeBlock = moleculeBlocks[moleculeType];

        gmx_moltype_t& moltype = mtop->moltype[moleculeType];
        init_t_atoms(&moltype.atoms, moleculeBlock.numAtoms, false);
        for (const int atom : moleculeBlock.perturbedAtoms)
        {
            moltype.atoms.atom[atom].typeB = 1;
        }
        std::vector<int> allAtomsInMolecule(moleculeBlock.numAtoms);
        std::iota(allAtomsInMolecule.begin(), allAtomsInMolecule.end(), 0);
        for (int atom = 0; atom < moleculeBlock.numAtoms; atom++)
        {
            moltype.excls.pushBack(allAtomsInMolecule);
        }

        gmx_molblock_t& molblock = mtop->molblock.emplace_back();
        molblock.type            = moleculeType;
        molblock.nmol            = moleculeBlock.numMolecules;
        mtop->natoms += moleculeBlock.numMolecules * moleculeBlock.numAtoms;
    }
    mtop->intermolecularExclusionGroup = exclusionGroup;
    mtop->finalize();

    return mtop;
}

//! Parameters of a test topology and its number of perturbed exclusions
struct ExclusionCheckerTestParameters
{
    //! Name of the test case
    std::string name;
    //! The molecule types and blocks of the topology
    std::vector<MoleculeBlockDescription> moleculeBlocks;
    //! The global atom indices of the inter-molecular exclusion group
    std::vector<int> exclusionGroup;
    //! The number of exclusions (one-way, including self-exclusions) involving a perturbed atom
    int numPerturbedExclusions;
};

//! Returns the test topologies with their numbers of perturbed exclusions
std::vector<ExclusionCheckerTestParameters> testTopologies()
{
    return {
        // Each of the two molecules with atom 0 perturbed has 3 perturbed exclusions,
        // the molecule without perturbed atoms has none
        { "NoExclusionGroup", { { 3, { 0 }, 2 }, { 3, {}, 1 } }, {}, 6 },
        // Each molecule has 3 perturbed intra-molecular exclusions,
        // and 5 of the 9 atom pairs between the molecules involve a perturbed atom
        { "ExclusionGroupWithOneMoleculeType", { { 3, { 0 }, 2 } }, { 0, 1, 2, 3, 4, 5 }, 11 },
        // 3 intra-molecular and 3 inter-molecular perturbed exclusions
        { "ExclusionGroupWithPerturbedFirstMoleculeType",
          { { 3, { 0 }, 1 }, { 3, {}, 1 } },
          { 0, 1, 2, 3, 4, 5 },
          6 },
        // As above, but with the perturbed molecule after the unperturbed one
        { "ExclusionGroupWithPerturbedSecondMoleculeType",
          { { 3, {}, 1 }, { 3, { 0 }, 1 } },
          { 0, 1, 2, 3, 4, 5 },
          6 },
        // 3 intra-molecular and 1 inter-molecular perturbed exclusions,
        // with the molecule with fewer atoms before the perturbed molecule
        { "ExclusionGroupWithSmallerFirstMoleculeType", { { 1, {}, 1 }, { 3, { 2 }, 1 } }, { 0, 1, 2, 3 }, 4 },
    };
}

//! Returns the name of the test case for GoogleTest
std::string nameOfTest(const testing::TestParamInfo<ExclusionCheckerTestParameters>& info)
{
    return info.param.name;
}

//! Test fixture for checks with the correct count
using ExclusionCheckerTest = testing::TestWithParam<ExclusionCheckerTestParameters>;
//! Test fixture for checks with an incorrect count, which end the process
using ExclusionCheckerDeathTest = testing::TestWithParam<ExclusionCheckerTestParameters>;

TEST_P(ExclusionCheckerTest, AcceptsCorrectCount)
{
    const ExclusionCheckerTestParameters& parameters = GetParam();

    const auto       mtop = makeTopology(parameters.moleculeBlocks, parameters.exclusionGroup);
    const MpiComm    mpiComm(MpiComm::SingleRank{});
    ExclusionChecker exclusionChecker(mpiComm, *mtop, nullptr);

    // With a single rank the check is performed directly and gives a fatal error at mismatch
    exclusionChecker.scheduleCheckOfExclusions(parameters.numPerturbedExclusions);
}

TEST_P(ExclusionCheckerDeathTest, RejectsIncorrectCount)
{
    const ExclusionCheckerTestParameters& parameters = GetParam();

    const auto       mtop = makeTopology(parameters.moleculeBlocks, parameters.exclusionGroup);
    const MpiComm    mpiComm(MpiComm::SingleRank{});
    ExclusionChecker exclusionChecker(mpiComm, *mtop, nullptr);

    // The fatal error reports the expected minus the found count
    GMX_EXPECT_DEATH_IF_SUPPORTED(
            exclusionChecker.scheduleCheckOfExclusions(parameters.numPerturbedExclusions - 1),
            "There are 1 perturbed, excluded non-bonded pair interactions");
    GMX_EXPECT_DEATH_IF_SUPPORTED(
            exclusionChecker.scheduleCheckOfExclusions(parameters.numPerturbedExclusions + 1),
            "There are -1 perturbed, excluded non-bonded pair interactions");
}

INSTANTIATE_TEST_SUITE_P(Topologies, ExclusionCheckerTest, testing::ValuesIn(testTopologies()), nameOfTest);
INSTANTIATE_TEST_SUITE_P(Topologies, ExclusionCheckerDeathTest, testing::ValuesIn(testTopologies()), nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
