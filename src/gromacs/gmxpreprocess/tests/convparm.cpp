/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for convparm.cpp
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/convparm.h"

#include <cmath>

#include <array>
#include <numeric>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/naming.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(ConvertInteractionsTest, DoingNothingWorks)
{
    const int                             numAtomTypes = 0;
    std::array<InteractionsOfType, F_NRE> nonBondedInteractions;
    std::vector<MoleculeInformation>      moleculesInformation;
    const MoleculeInformation*            intermolecularInteractions = nullptr;
    CombinationRule                       combinationRule            = CombinationRule::Geometric;
    const double                          repulsionPower             = 12.0;
    const real                            fudgeQQ                    = 1.0;
    gmx_mtop_t                            mtop;

    convertInteractionsOfType(numAtomTypes,
                              nonBondedInteractions,
                              moleculesInformation,
                              intermolecularInteractions,
                              combinationRule,
                              repulsionPower,
                              fudgeQQ,
                              &mtop);
}

using testing::Eq;
using testing::Pointwise;

//! Fill a vector with size \c size with integer values increasing from 0
std::vector<int> iotaVector(const std::size_t size)
{
    std::vector<int> v(size);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

/*! \brief Fill a vector with iota-style dummy interaction parameter
 * values identical for state A and B
 *
 * Any unused array entries are filled with NaN (on platforms where
 * this is supported). */
std::array<real, MAXFORCEPARAM> iotaParams(const std::size_t size)
{
    GMX_RELEASE_ASSERT(size <= MAXFORCEPARAM / 2, "Size too big for interaction parameter array");
    std::array<real, MAXFORCEPARAM> a;
    std::fill(a.begin(), a.end(), std::nan("0"));
    std::iota(a.begin(), a.begin() + size, 5._real);
    std::iota(a.begin() + size, a.begin() + size * 2, 5._real);
    return a;
}

//! Define a test fixture class taking an integer paramter for ftype
using ConvertInteractionsTest = ::testing::TestWithParam<std::tuple<int>>;

TEST_P(ConvertInteractionsTest, Works)
{
    const int                             numAtomTypes = 0;
    std::array<InteractionsOfType, F_NRE> nonBondedInteractions;
    std::vector<MoleculeInformation>      moleculesInformation;
    char**                                dummyName = nullptr;

    const int ftype = std::get<0>(GetParam());

    // Ensure this function type is handled by convertInteractionsOfType.
    if (!shouldConvertInteractionType(ftype))
    {
        GTEST_SKIP() << "Skipping interaction type that does not represent a interaction with "
                        "parameters converted in grompp";
    }
    for (const auto unsupportedFunctionType :
         { F_GB12_NOLONGERUSED, F_GB13_NOLONGERUSED, F_GB14_NOLONGERUSED, F_GBPOL_NOLONGERUSED, F_NPSOLVATION_NOLONGERUSED })
    {
        if (ftype == unsupportedFunctionType)
        {
            GTEST_SKIP() << "Skipping bonded function type no longer supported";
        }
    }

    // Define a molecule type with a single interaction of type ftype and add it
    // to the molecules information object
    {
        std::array<InteractionsOfType, F_NRE> moleculeInteractions;
        // For function types with no parameters, assign_param()
        // assumes the parameters are all zero, which leads to not
        // appending all-zero parameter sest to the parameter list
        if (interaction_function[ftype].nrfpA == 0 && interaction_function[ftype].nrfpB == 0)
        {
            std::vector<real> zeroParams(MAXFORCEPARAM, 0);
            moleculeInteractions[ftype].interactionTypes.emplace_back(
                    InteractionOfType{ iotaVector(NRAL(ftype)), zeroParams, "name" });
        }
        else
        {
            // Note force parameters end up defined for both FEP states
            moleculeInteractions[ftype].interactionTypes.emplace_back(
                    InteractionOfType{ iotaVector(NRAL(ftype)), iotaParams(NRFPA(ftype)), "name" });
        }
        moleculesInformation.emplace_back(MoleculeInformation{
                dummyName, 0, false, t_atoms{}, t_block{}, ListOfLists<int>{}, moleculeInteractions });
    }

    const MoleculeInformation* intermolecularInteractions = nullptr;
    const double               repulsionPower             = 12.0;
    const real                 fudgeQQ                    = 1.0;
    gmx_mtop_t                 mtop;
    // Add molecule type with index 0
    mtop.moltype.resize(1);
    // Fill the molecule type
    convertInteractionsOfType(numAtomTypes,
                              nonBondedInteractions,
                              moleculesInformation,
                              intermolecularInteractions,
                              CombinationRule::Geometric,
                              repulsionPower,
                              fudgeQQ,
                              &mtop);
    // Add a molecule block with 1 molecule of type 0 which was just filled
    mtop.molblock.emplace_back(gmx_molblock_t{ 0, 1 });

    if (interaction_function[ftype].flags & IF_BOND)
    {
        EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_BOND), 1)
                << "topology has one bonded interaction";
        ASSERT_EQ(gmx_mtop_ftype_count(mtop, ftype), 1) << "topology has one kind of interaction";
        EXPECT_EQ(mtop.moltype[0].ilist[ftype].iatoms[0], 0)
                << "the first interaction of the first molecule type uses the first "
                   "interaction function parameters";
    }
    std::vector<int> expected = iotaVector(NRAL(ftype));

    EXPECT_EQ(mtop.moltype[0].ilist[ftype].iatoms[0], 0)
            << "first interaction has index zero when there is only one interaction added";
    EXPECT_THAT(makeArrayRef(mtop.moltype[0].ilist[ftype].iatoms).subArray(1, NRAL(ftype)), Pointwise(Eq(), expected))
            << "the first interaction of the first molecule type has the expected atom indices";
    ASSERT_EQ(mtop.ffparams.numTypes(), 1)
            << "topology has one set of interaction function parameters for function types "
               "that have parameters";
    // It would be nice to check that the contents of mtop.ffparams.iparams[0] has the expected
    // relationship with moleculesInformation[0].interactions[ftype].interactionTypes[0].forceParam().subArray(0, NRFP(ftype)),
    // but t_iparams is a union that contains a variety of data types and numbers of logical units *and*
    // assign_param() in convparm.cpp sometimes converts units or precomputes squares for the convenience of mdrun.
    // Our code is not yet flexible enough to make that an easy job.
}

std::string ftypeToName(const int ftype)
{
    GMX_RELEASE_ASSERT(ftype < F_NRE, "Must have valid kind of interaction function");
    return interaction_function[ftype].longname;
}

const NameOfTestFromTuple<std::tuple<int>> sc_testNamer{ std::make_tuple(ftypeToName) };

using testing::Combine;
using testing::Range;
INSTANTIATE_TEST_SUITE_P(InteractionFunctionKind,
                         ConvertInteractionsTest,
                         Combine(Range(0, static_cast<int>(F_NRE))),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
