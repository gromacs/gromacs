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
 * Tests for virtual sites.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/vsite.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Two vsites with the higher depending on the lower one
gmx_moltype_t moltypeVSite12()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr                                        = 3;
    moltype.ilist[InteractionFunction::VirtualSite1].iatoms = { 0, 1, 0 };
    moltype.ilist[InteractionFunction::VirtualSite2].iatoms = { 1, 2, 0, 1 };

    return moltype;
}

//! Two vsites with the lower depending on the higher one
gmx_moltype_t moltypeVSite21()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr                                        = 3;
    moltype.ilist[InteractionFunction::VirtualSite2].iatoms = { 0, 2, 0, 1 };
    moltype.ilist[InteractionFunction::VirtualSite1].iatoms = { 1, 1, 2 };

    return moltype;
}

//! Two vsites of the same type with one depending on the other
gmx_moltype_t moltypeVSite11()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr                                        = 3;
    moltype.ilist[InteractionFunction::VirtualSite1].iatoms = { 0, 1, 0, 1, 2, 1 };

    return moltype;
}

gmx_molblock_t oneMolblock()
{
    return { 0, 1, {}, {} };
}

//! Test fixture class
class VSiteTest : public ::testing::Test
{
public:
    //! Global toplogy to use in tests
    gmx_mtop_t mtop_;
};

TEST_F(VSiteTest, VSiteLowerConstructingWorks)
{
    mtop_.moltype.emplace_back(moltypeVSite12());
    mtop_.molblock.push_back(oneMolblock());
    EXPECT_NO_THROW(makeVirtualSitesHandler(mtop_, nullptr, PbcType::Xyz, {}));
}

TEST_F(VSiteTest, VSiteHigherConstructingThrows)
{
    mtop_.moltype.emplace_back(moltypeVSite21());
    mtop_.molblock.push_back(oneMolblock());
    EXPECT_THROW_GMX(makeVirtualSitesHandler(mtop_, nullptr, PbcType::Xyz, {}), InvalidInputError);
}

TEST_F(VSiteTest, VSiteEqualConstructingThrows)
{
    mtop_.moltype.emplace_back(moltypeVSite11());
    mtop_.molblock.push_back(oneMolblock());
    EXPECT_THROW_GMX(makeVirtualSitesHandler(mtop_, nullptr, PbcType::Xyz, {}), InvalidInputError);
}

} // namespace
} // namespace test
} // namespace gmx
