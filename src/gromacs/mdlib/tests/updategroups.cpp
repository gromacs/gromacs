/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Tests for the update groups functionality.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/updategroups.h"

#include <gtest/gtest.h>

#include "gromacs/topology/topology.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace
{

/* TODO: Actually initialize moltype.atoms.atom when this is converted to C++ */

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

/*! \brief Returns a butane united-atom molecule */
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

TEST(UpdateGroups, ethaneUA)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(ethaneUA());
    t_iparams iparams;
    iparams.constr = { 0.3, 0.3 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);

    real temperature = 298;
    real maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.3 / 2);
}

TEST(UpdateGroups, methane)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(methane());
    t_iparams iparams;
    iparams.constr = { 0.1, 0.1 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);

    real temperature = 298;
    real maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.14);
}
TEST(UpdateGroups, ethane)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(ethane());
    t_iparams iparams;
    iparams.constr = { 0.1, 0.1 };
    mtop.ffparams.iparams.push_back(iparams);
    iparams.harmonic = { 107.800, 276.144, 107.800, 276.144 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 2);

    real temperature = 298;
    real maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.094746813);

    temperature = 0;
    maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.10310466);

    temperature = -1;
    maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.125);
}

TEST(UpdateGroups, butaneUA)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(butaneUA());
    t_iparams iparams;
    iparams.constr = { 0.3, 0.3 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    EXPECT_EQ(updateGroups.size(), 0);
}

TEST(UpdateGroups, waterThreeSite)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(waterThreeSite());
    t_iparams iparams;
    iparams.settle = { 0.1, 0.1633 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);

    real temperature = 298;
    real maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.083887339);
}

// Tests update group with virtual site
TEST(UpdateGroups, waterFourSite)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(waterFourSite());
    t_iparams iparams[2];
    iparams[0].settle = { 0.1, 0.1633 };
    iparams[1].vsite  = { 0.128, 0.128 };
    mtop.ffparams.iparams.push_back(iparams[0]);
    mtop.ffparams.iparams.push_back(iparams[1]);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);
}

TEST(UpdateGroups, fourAtomsWithSettle)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(waterThreeSite());
    mtop.moltype.back().atoms.nr = 4;

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 2);
}

// Tests groups with two constraints and an angle potential
TEST(UpdateGroups, waterFlexAngle)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(waterFlexAngle());
    t_iparams iparams;
    iparams.constr = { 0.1, 0.1 };
    mtop.ffparams.iparams.push_back(iparams);
    iparams.harmonic = { 109.47, 383.0, 109.47, 383.0 };
    mtop.ffparams.iparams.push_back(iparams);

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 1);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);

    real temperature = 298;
    real maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.090824135);

    temperature = 0;
    maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.1);

    temperature = -1;
    maxRadius   = computeMaxUpdateGroupRadius(mtop, updateGroups, temperature);
    EXPECT_FLOAT_EQ(maxRadius, 0.1);
}

TEST(UpdateGroups, twoMoltypes)
{
    gmx_mtop_t mtop;

    mtop.moltype.emplace_back(methane());
    t_iparams iparams;
    iparams.constr = { 0.1, 0.1 };
    mtop.ffparams.iparams.push_back(iparams);

    mtop.moltype.emplace_back(waterThreeSite());
    // Note: iparams not accessed for SETTLE when not computing radius

    auto updateGroups = gmx::makeUpdateGroups(mtop);

    ASSERT_EQ(updateGroups.size(), 2);
    EXPECT_EQ(updateGroups[0].numBlocks(), 1);
    EXPECT_EQ(updateGroups[1].numBlocks(), 1);
}

} // namespace

} // namespace gmx
