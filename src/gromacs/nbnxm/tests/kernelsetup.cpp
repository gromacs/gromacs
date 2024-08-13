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
 * Tests for Setup of kernels.
 *
 * \author Joe Jordan <ejjordan@kth.se>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/kernel_common.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

TEST(KernelSetupTest, getCoulombKernelTypeRF)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::NotSet, CoulombInteractionType::RF, false),
              CoulombKernelType::ReactionField);
}

TEST(KernelSetupTest, getCoulombKernelTypeCut)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::NotSet, CoulombInteractionType::Cut, false),
              CoulombKernelType::ReactionField);
}

TEST(KernelSetupTest, getCoulombKernelTypeTable)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::Table, CoulombInteractionType::Count, true),
              CoulombKernelType::Table);
}

TEST(KernelSetupTest, getCoulombKernelTypeTableTwin)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::Table, CoulombInteractionType::Count, false),
              CoulombKernelType::TableTwin);
}

TEST(KernelSetupTest, getCoulombKernelTypeEwald)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::NotSet, CoulombInteractionType::Count, true),
              CoulombKernelType::Ewald);
}

TEST(KernelSetupTest, getCoulombKernelTypeEwaldTwin)
{
    EXPECT_EQ(getCoulombKernelType(gmx::EwaldExclusionType::NotSet, CoulombInteractionType::Count, false),
              CoulombKernelType::EwaldTwin);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombGeomNone)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::Geometric,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::None,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBGEOM);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombGeomPotShift)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::Geometric,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::PotShift,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBGEOM);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombLBNone)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::LorentzBerthelot,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::None,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBLB);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombLBPotShift)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::LorentzBerthelot,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::PotShift,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBLB);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombNoneNone)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::None,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::None,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBNONE);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutCombNonePotShift)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::None,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::PotShift,
                               LongRangeVdW::Count),
              vdwktLJCUT_COMBNONE);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutThrows)
{
    EXPECT_ANY_THROW(getVdwKernelType(NbnxmKernelType::NotSet,
                                      LJCombinationRule::Count,
                                      VanDerWaalsType::Cut,
                                      InteractionModifiers::PotShift,
                                      LongRangeVdW::Count));
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutForceSwitch)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::None,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::ForceSwitch,
                               LongRangeVdW::Count),
              vdwktLJFORCESWITCH);
}

TEST(KernelSetupTest, getVdwKernelTypePmeGeom)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::Cpu4x4_PlainC,
                               LJCombinationRule::None,
                               VanDerWaalsType::Pme,
                               InteractionModifiers::Count,
                               LongRangeVdW::Geom),
              vdwktLJEWALDCOMBGEOM);
}

TEST(KernelSetupTest, getVdwKernelTypePmeNone)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::Cpu4x4_PlainC,
                               LJCombinationRule::None,
                               VanDerWaalsType::Pme,
                               InteractionModifiers::Count,
                               LongRangeVdW::Count),
              vdwktLJEWALDCOMBLB);
}

TEST(KernelSetupTest, getVdwKernelTypeLjCutPotSwitch)
{
    EXPECT_EQ(getVdwKernelType(NbnxmKernelType::NotSet,
                               LJCombinationRule::None,
                               VanDerWaalsType::Cut,
                               InteractionModifiers::PotSwitch,
                               LongRangeVdW::Count),
              vdwktLJPOTSWITCH);
}

TEST(KernelSetupTest, getVdwKernelTypeAllCountThrows)
{
    // Count cannot be used for VanDerWaalsType or InteractionModifiers because of calls to
    // enumValueToString(), which require a valid choice to have been made.
    EXPECT_ANY_THROW(getVdwKernelType(NbnxmKernelType::NotSet,
                                      LJCombinationRule::Count,
                                      VanDerWaalsType::Cut,
                                      InteractionModifiers::None,
                                      LongRangeVdW::Count));
}

} // namespace
} // namespace test
} // namespace gmx
