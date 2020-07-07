/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests amplitude lookup for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/densityfitting/densityfittingamplitudelookup.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

class DensityFittingAmplitudeLookupTest : public ::testing::Test
{
public:
    DensityFittingAmplitudeLookupTest()
    {
        atoms_.nr      = numberOfAtoms_;
        atoms_.massT   = masses_.data();
        atoms_.chargeA = charges_.data();
    }

protected:
    int               numberOfAtoms_ = 3;
    std::vector<real> masses_        = { 2, 3, 4 };
    std::vector<real> charges_       = { 20, 30, 40 };
    t_mdatoms         atoms_         = {};
    std::vector<int>  lookupIndices_ = { 1, 2 };
};

TEST_F(DensityFittingAmplitudeLookupTest, Unity)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Unity);
    const auto                    lookupResult = lookup(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 1);
    EXPECT_EQ(lookupResult[1], 1);
}

TEST_F(DensityFittingAmplitudeLookupTest, Charge)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Charge);
    const auto                    lookupResult = lookup(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 30);
    EXPECT_EQ(lookupResult[1], 40);
}

TEST_F(DensityFittingAmplitudeLookupTest, Masses)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Mass);
    const auto                    lookupResult = lookup(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 3);
    EXPECT_EQ(lookupResult[1], 4);
}

TEST_F(DensityFittingAmplitudeLookupTest, CanCopyAssign)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Unity);
    DensityFittingAmplitudeLookup lookupCopied = lookup;
    const auto                    lookupResult = lookupCopied(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 1);
    EXPECT_EQ(lookupResult[1], 1);
}

TEST_F(DensityFittingAmplitudeLookupTest, CanCopyConstruct)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Unity);
    DensityFittingAmplitudeLookup lookupCopied(lookup);
    const auto                    lookupResult = lookupCopied(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 1);
    EXPECT_EQ(lookupResult[1], 1);
}

TEST_F(DensityFittingAmplitudeLookupTest, CanMoveAssign)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Unity);
    DensityFittingAmplitudeLookup lookupCopied = std::move(lookup);
    const auto                    lookupResult = lookupCopied(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 1);
    EXPECT_EQ(lookupResult[1], 1);
}

TEST_F(DensityFittingAmplitudeLookupTest, CanMoveConstruct)
{
    DensityFittingAmplitudeLookup lookup(DensityFittingAmplitudeMethod::Unity);
    DensityFittingAmplitudeLookup lookupCopied(std::move(lookup));
    const auto                    lookupResult = lookupCopied(atoms_, lookupIndices_);
    EXPECT_EQ(lookupResult[0], 1);
    EXPECT_EQ(lookupResult[1], 1);
}

} // namespace gmx
