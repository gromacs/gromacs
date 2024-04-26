/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for the update groups functionality.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/updategroups.h"

#include <gtest/gtest.h>

#include "gromacs/topology/topology.h"

#include "testutils/loggertest.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace
{

/* TODO: Actually initialize moltype.atoms.atom when this is converted to C++ */

/*! \brief Returns a flexible ethane united-atom molecule */
gmx_moltype_t flexibleEthaneUA()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr              = 2;
    moltype.ilist[F_BONDS].iatoms = { 0, 0, 1 };

    return moltype;
}

/*! \brief Returns an ethane united-atom molecule */
gmx_moltype_t ethaneUA()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 2;
    moltype.ilist[F_CONSTR].iatoms = { 0, 0, 1 };

    return moltype;
}

/*! \brief Returns a methane molecule */
gmx_moltype_t methane()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 5;
    moltype.ilist[F_CONSTR].iatoms = { 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4 };

    return moltype;
}

/*! \brief Returns an ethane molecule */
gmx_moltype_t ethane()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 8;
    moltype.ilist[F_CONSTR].iatoms = { 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 4, 5, 0, 4, 6, 0, 4, 7 };
    moltype.ilist[F_ANGLES].iatoms = { 1, 1, 0, 2, 1, 1, 0, 3, 1, 2, 0, 3,
                                       1, 5, 4, 6, 1, 5, 4, 7, 1, 6, 4, 7 };

    return moltype;
}

/*! \brief Returns a butane fully-constrained united-atom molecule */
gmx_moltype_t butaneUA()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 4;
    moltype.ilist[F_CONSTR].iatoms = { 0, 0, 1, 0, 1, 2, 0, 2, 3 };

    return moltype;
}

/*! \brief Returns a three-site water molecule */
gmx_moltype_t waterThreeSite()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 3;
    moltype.ilist[F_SETTLE].iatoms = { 0, 0, 1, 2 };

    return moltype;
}

/*! \brief Returns a four-site water molecule with virtual site */
gmx_moltype_t waterFourSite()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 4;
    moltype.ilist[F_SETTLE].iatoms = { 0, 1, 2, 3 };
    moltype.ilist[F_VSITE3].iatoms = { 1, 0, 1, 2, 3 };

    return moltype;
}

/*! \brief Returns a water molecule with flexible angle */
gmx_moltype_t waterFlexAngle()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 3;
    moltype.ilist[F_CONSTR].iatoms = {
        0, 0, 1, 0, 0, 2,
    };
    moltype.ilist[F_ANGLES].iatoms = {
        1,
        1,
        0,
        2,
    };

    return moltype;
}

//! Test fixture class
class UpdateGroupsTest : public ::testing::Test
{
public:
    //! Global toplogy to use in tests
    gmx_mtop_t mtop_;
    //! Default temperature for tests
    real temperature_ = 298;
    //! Logger to use in tests
    test::LoggerTestHelper logHelper_;
};

TEST_F(UpdateGroupsTest, WithEthaneUA)
{
    mtop_.moltype.emplace_back(ethaneUA());
    {
        t_iparams iparams;
        iparams.constr = { 0.3, 0.3 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.3 / 2);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, WithMethane)
{
    mtop_.moltype.emplace_back(methane());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.14);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, WithEthane)
{
    mtop_.moltype.emplace_back(ethane());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
        iparams.harmonic = { 107.800, 276.144, 107.800, 276.144 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 2);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.094746813);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, CheckRadiusCalculationAtDifferentTemperaturesWithEthane)
{
    mtop_.moltype.emplace_back(ethane());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
        iparams.harmonic = { 107.800, 276.144, 107.800, 276.144 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 2);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.094746813);

    // Observe that the temperature affects the radius only when valid
    temperature_ = 0;
    maxRadius    = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.10310466);

    temperature_ = -1;
    maxRadius    = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.125);
}

TEST_F(UpdateGroupsTest, WithButaneUALogsThatUnsuitableForUpdateGroups)
{
    mtop_.moltype.emplace_back(butaneUA());
    {
        t_iparams iparams;
        iparams.constr = { 0.3, 0.3 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), true);

    EXPECT_EQ(std::get<std::string>(result),
              "there are three or more consecutively coupled constraints");

    std::vector<RangePartitioning> updateGroupingsPerMoleculeType;

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.0);
}

TEST_F(UpdateGroupsTest, WithWaterThreeSite)
{
    mtop_.moltype.emplace_back(waterThreeSite());
    {
        t_iparams iparams;
        iparams.settle = { 0.1, 0.1633 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.083887339);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

// Tests update group with virtual site
TEST_F(UpdateGroupsTest, WithWaterFourSite)
{
    mtop_.moltype.emplace_back(waterFourSite());
    {
        t_iparams iparams[2];
        iparams[0].settle = { 0.1, 0.1633 };
        iparams[1].vsite  = { 0.128, 0.128 };
        mtop_.ffparams.iparams.push_back(iparams[0]);
        mtop_.ffparams.iparams.push_back(iparams[1]);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.083887339);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, WithFourAtomsWithSettle)
{
    mtop_.moltype.emplace_back(waterThreeSite());
    mtop_.moltype.back().atoms.nr = 4;
    {
        t_iparams iparams[2];
        iparams[0].settle = { 0.1, 0.1633 };
        iparams[1].vsite  = { 0.128, 0.128 };
        mtop_.ffparams.iparams.push_back(iparams[0]);
        mtop_.ffparams.iparams.push_back(iparams[1]);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 2);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.083887339);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

// Tests groups with two constraints and an angle potential
TEST_F(UpdateGroupsTest, WithWaterFlexAngle)
{
    mtop_.moltype.emplace_back(waterFlexAngle());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
        iparams.harmonic = { 109.47, 383.0, 109.47, 383.0 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.090824135);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, CheckRadiusCalculationAtDifferentTemperaturesWithWaterFlexAngle)
{
    mtop_.moltype.emplace_back(waterFlexAngle());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
        iparams.harmonic = { 109.47, 383.0, 109.47, 383.0 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.090824135);

    // Observe that the temperature affects the radius only when valid
    temperature_ = 0;
    maxRadius    = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.1);

    temperature_ = -1;
    maxRadius    = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.1);
}

TEST_F(UpdateGroupsTest, WithTwoMoltypes)
{
    mtop_.moltype.emplace_back(methane());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    mtop_.moltype.emplace_back(waterThreeSite());
    // Note: iparams not accessed for SETTLE when not computing radius

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 2);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[1].numBlocks(), 1);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0.14);

    logHelper_.expectNoEntries(MDLogger::LogLevel::Info);
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, LogsWhenSizesAreInvalid)
{
    mtop_.moltype.emplace_back(methane());
    {
        t_iparams iparams;
        iparams.constr = { 0.1, 0.1 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "The combination of rlist and box size prohibits");
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), 1e9_real, false, true, true, 1e6_real);
}

TEST_F(UpdateGroupsTest, LogsWhenUpdateGroupsAreNotUseful)
{
    mtop_.moltype.emplace_back(flexibleEthaneUA());
    {
        t_iparams iparams;
        iparams.harmonic = { 0.1, 10.0, 0.1, 10.0 };
        mtop_.ffparams.iparams.push_back(iparams);
    }

    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop_);

    ASSERT_EQ(std::holds_alternative<std::string>(result), false);

    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);

    ASSERT_EQ(updateGroupingsPerMoleculeType.size(), 1);
    EXPECT_EQ(updateGroupingsPerMoleculeType[0].numBlocks(), 2);

    real maxRadius = computeMaxUpdateGroupRadius(mtop_, updateGroupingsPerMoleculeType, temperature_);
    EXPECT_FLOAT_EQ(maxRadius, 0);

    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "No constraints or virtual sites are in use");
    UpdateGroups updateGroups = makeUpdateGroups(
            logHelper_.logger(), std::move(updateGroupingsPerMoleculeType), maxRadius, false, true, false, 1e6_real);
}


} // namespace

} // namespace gmx
