/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
#include "gmxpre.h"

#include <algorithm>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdlib/lincs_constraint_group_sizes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace test
{
namespace
{

using testing::Eq;
using testing::Pointwise;

struct ConstraintSystem
{
    std::vector<AtomPair> constraintsHost;
    std::vector<int>      referenceData;
};

//! Add ethane molecule to constraints
std::vector<AtomPair> addNextConstrainedEthane(const int firstAtom, std::vector<AtomPair>&& constraints)
{
    // Add two groups of three C-H constraints, one for each methyl group
    constraints.push_back({ firstAtom + 0, firstAtom + 1 }); // C1-H1
    constraints.push_back({ firstAtom + 0, firstAtom + 2 }); // C1-H2
    constraints.push_back({ firstAtom + 0, firstAtom + 3 }); // C1-H3
    constraints.push_back({ firstAtom + 5, firstAtom + 6 }); // C2-H4
    constraints.push_back({ firstAtom + 5, firstAtom + 7 }); // C2-H5
    constraints.push_back({ firstAtom + 5, firstAtom + 8 }); // C2-H6
    return constraints;
}

//! Add methane molecule to constraints
std::vector<AtomPair> addNextConstrainedMethane(const int firstAtom, std::vector<AtomPair>&& constraints)
{
    // Add single group of four C-H constraints
    constraints.push_back({ firstAtom + 0, firstAtom + 1 }); // C1-H1
    constraints.push_back({ firstAtom + 0, firstAtom + 2 }); // C1-H2
    constraints.push_back({ firstAtom + 0, firstAtom + 3 }); // C1-H3
    constraints.push_back({ firstAtom + 0, firstAtom + 4 }); // C1-H4
    return constraints;
}

//! Add ethanol molecule to constraints
std::vector<AtomPair> addNextConstrainedEthanol(const int firstAtom, std::vector<AtomPair>&& constraints)
{
    // Add single group of four C-H constraints
    constraints.push_back({ firstAtom + 0, firstAtom + 1 }); // C1-H1
    constraints.push_back({ firstAtom + 0, firstAtom + 2 }); // C1-H2
    constraints.push_back({ firstAtom + 0, firstAtom + 3 }); // C1-H3
    constraints.push_back({ firstAtom + 5, firstAtom + 6 }); // C2-H4
    constraints.push_back({ firstAtom + 5, firstAtom + 7 }); // C2-H5
    constraints.push_back({ firstAtom + 8, firstAtom + 9 }); // O1-HO
    return constraints;
}

static constexpr int c_numConstraintMethane = 4;
static constexpr int c_numConstraintEthane  = 6;
static constexpr int c_numConstraintEthanol = 6;

std::vector<AtomPair> fillConstraints(int numMethane, int numEthane, int numEthanol, int numConstraintThreads)
{
    std::vector<AtomPair> constraints;
    constraints.reserve(numConstraintThreads);

    const int numConstraints = (numMethane * c_numConstraintMethane) + (numEthane * c_numConstraintEthane)
                               + (numEthanol * c_numConstraintEthanol);

    GMX_RELEASE_ASSERT(numConstraints <= numConstraintThreads, "Invalid test setup");

    int indexOfFirstAtom = 0;
    for (int i = 0; i < numMethane; ++i)
    {
        constraints      = addNextConstrainedMethane(indexOfFirstAtom, std::move(constraints));
        indexOfFirstAtom = constraints.back().j + 1;
    }

    for (int i = 0; i < numEthane; ++i)
    {
        constraints      = addNextConstrainedEthane(indexOfFirstAtom, std::move(constraints));
        indexOfFirstAtom = constraints.back().j + 1;
    }

    for (int i = 0; i < numEthanol; ++i)
    {
        constraints      = addNextConstrainedEthanol(indexOfFirstAtom, std::move(constraints));
        indexOfFirstAtom = constraints.back().j + 1;
    }

    for (int i = numConstraints; i < numConstraintThreads; ++i)
    {
        constraints.emplace_back(AtomPair{ -1, -1 });
    }
    return constraints;
}

std::vector<int> generateReferenceData(int numMethane, int numEthane, int numEthanol, int numConstraintThreads)
{
    std::vector<int> referenceData(numConstraintThreads);
    std::fill(referenceData.begin(), referenceData.end(), -1);
    const int numConstraints = (numMethane * c_numConstraintMethane) + (numEthane * c_numConstraintEthane)
                               + (numEthanol * c_numConstraintEthanol);

    GMX_RELEASE_ASSERT(numConstraints <= numConstraintThreads, "Invalid test setup");
    int indexOfFirstAtom = 0;
    for (int i = 0; i < numMethane; ++i)
    {
        referenceData[indexOfFirstAtom] = 3;
        for (int j = indexOfFirstAtom + 1; j < indexOfFirstAtom + 4; ++j)
        {
            referenceData[j] = 0;
        }
        indexOfFirstAtom += 4;
    }
    for (int i = 0; i < numEthane; ++i)
    {
        referenceData[indexOfFirstAtom] = 2;
        for (int j = indexOfFirstAtom + 1; j < indexOfFirstAtom + 3; ++j)
        {
            referenceData[j] = 0;
        }
        indexOfFirstAtom += 3;
        referenceData[indexOfFirstAtom] = 2;
        for (int j = indexOfFirstAtom + 1; j < indexOfFirstAtom + 3; ++j)
        {
            referenceData[j] = 0;
        }
        indexOfFirstAtom += 3;
    }
    for (int i = 0; i < numEthanol; ++i)
    {
        referenceData[indexOfFirstAtom] = 2;
        for (int j = indexOfFirstAtom + 1; j < indexOfFirstAtom + 3; ++j)
        {
            referenceData[j] = 0;
        }
        indexOfFirstAtom += 3;
        referenceData[indexOfFirstAtom] = 1;
        for (int j = indexOfFirstAtom + 1; j < indexOfFirstAtom + 2; ++j)
        {
            referenceData[j] = 0;
        }
        indexOfFirstAtom += 2;
        indexOfFirstAtom++;
    }
    return referenceData;
}


std::vector<ConstraintSystem> c_constraintReorderSystemList = []
{
    std::vector<ConstraintSystem> constraintReorderSystemList;
    {
        const int        numConstraintThreads = 64;
        ConstraintSystem emptySystem;
        emptySystem.constraintsHost = fillConstraints(0, 0, 0, numConstraintThreads);
        emptySystem.referenceData   = generateReferenceData(0, 0, 0, numConstraintThreads);
        constraintReorderSystemList.emplace_back(emptySystem);
    }
    {
        const int        numConstraintThreads = 64;
        const int        numMethane           = numConstraintThreads / c_numConstraintMethane;
        ConstraintSystem methaneOnly;
        methaneOnly.constraintsHost = fillConstraints(numMethane, 0, 0, numConstraintThreads);
        methaneOnly.referenceData   = generateReferenceData(numMethane, 0, 0, numConstraintThreads);
        constraintReorderSystemList.emplace_back(methaneOnly);
    }
    {
        const int        numConstraintThreads = 64;
        const int        numEthane            = numConstraintThreads / c_numConstraintEthane;
        ConstraintSystem ethaneOnly;
        ethaneOnly.constraintsHost = fillConstraints(0, numEthane, 0, numConstraintThreads);
        ethaneOnly.referenceData   = generateReferenceData(0, numEthane, 0, numConstraintThreads);
        constraintReorderSystemList.emplace_back(ethaneOnly);
    }
    {
        const int        numConstraintThreads = 64;
        const int        numEthanol           = numConstraintThreads / c_numConstraintEthanol;
        ConstraintSystem ethanolOnly;
        ethanolOnly.constraintsHost = fillConstraints(0, 0, numEthanol, numConstraintThreads);
        ethanolOnly.referenceData   = generateReferenceData(0, 0, numEthanol, numConstraintThreads);
        constraintReorderSystemList.emplace_back(ethanolOnly);
    }
    {
        const int numConstraintThreads = 128;
        // divide number of entries by two so we have space over for other molecules
        const int numMethane = numConstraintThreads / c_numConstraintMethane / 2;
        const int numEthane =
                (numConstraintThreads - (numMethane * c_numConstraintMethane)) / c_numConstraintEthane;
        ConstraintSystem methaneEthaneMix;
        methaneEthaneMix.constraintsHost = fillConstraints(numMethane, numEthane, 0, numConstraintThreads);
        methaneEthaneMix.referenceData =
                generateReferenceData(numMethane, numEthane, 0, numConstraintThreads);
        constraintReorderSystemList.emplace_back(methaneEthaneMix);
    }
    {
        const int numConstraintThreads = 256;
        // divide number of entries by two so we have space over for other molecules
        const int numMethane = (numConstraintThreads / c_numConstraintMethane) / 2;
        // divide number of entries by two so we have space over for other molecules
        const int numEthane =
                ((numConstraintThreads - (numMethane * c_numConstraintMethane)) / c_numConstraintEthane) / 2;
        const int numEthanol = (numConstraintThreads - (numMethane * c_numConstraintMethane)
                                - (numEthane * c_numConstraintEthane))
                               / c_numConstraintEthanol;
        ConstraintSystem mix;
        mix.constraintsHost = fillConstraints(numMethane, numEthane, numEthanol, numConstraintThreads);
        mix.referenceData = generateReferenceData(numMethane, numEthane, numEthanol, numConstraintThreads);
        constraintReorderSystemList.emplace_back(mix);
    }
    {
        const int        numConstraintThreads = 150659; // Random prime number
        const int        numEthanol           = numConstraintThreads / c_numConstraintEthanol;
        ConstraintSystem largeSystem;
        largeSystem.constraintsHost = fillConstraints(0, 0, numEthanol, numConstraintThreads);
        largeSystem.referenceData   = generateReferenceData(0, 0, numEthanol, numConstraintThreads);
        constraintReorderSystemList.emplace_back(largeSystem);
    }

    return constraintReorderSystemList;
}();

class ConstraintsReorderingTest : public ::testing::TestWithParam<ConstraintSystem>
{
public:
    static void runTest(int numConstraints, ArrayRef<AtomPair> constraintsHost, ArrayRef<const int> referenceData)
    {
        std::vector<int> constraintGroupSize(constraintsHost.size(), -1);

        findConstraintGroupSizes(numConstraints, constraintsHost, constraintGroupSize);
        EXPECT_THAT(constraintGroupSize, Pointwise(Eq(), referenceData))
                << "Constraint group sizes don't match ";
    }
};

TEST_P(ConstraintsReorderingTest, DifferentSystemsWork)
{
    auto constraintsSystem = GetParam();
    auto constraintsHost   = constraintsSystem.constraintsHost;
    auto referenceData     = constraintsSystem.referenceData;
    auto numConstraints    = constraintsHost.size();

    runTest(numConstraints, constraintsHost, referenceData);
}

INSTANTIATE_TEST_SUITE_P(WithParameters,
                         ConstraintsReorderingTest,
                         ::testing::ValuesIn(c_constraintReorderSystemList));

} // namespace
} // namespace test
} // namespace gmx
