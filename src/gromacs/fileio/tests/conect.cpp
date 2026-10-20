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
 * Implements test of PDB CONECT generation routines
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "external/scope_guard/scope_guard.h"

#include "gromacs/fileio/pdbio.h"
#include "gromacs/topology/topology.h"

#include "testutils/posixmemstream.h"
#include "testutils/topologyhelpers.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(ConectTest, NoLeaks)
{
    gmx_mtop_t mtop;
    mtop.finalize();

    gmx_conect conect              = gmx_conect_init();
    auto       conectContentsGuard = sg::make_scope_guard([conect]() { gmx_conect_done(conect); });
}

TEST(ConectTest, EmptyToplogy)
{
    gmx_mtop_t mtop;
    mtop.finalize();

    gmx_conect conect              = gmx_conect_generate(mtop);
    auto       conectContentsGuard = sg::make_scope_guard([conect]() { gmx_conect_done(conect); });

    // On non-POSIX systems we can only check for e.g. segfaults above.
    if (PosixMemstream::canCheckBufferContents())
    {
        PosixMemstream memoryStream;
        gmx_conect_dump(memoryStream.stream(), conect);
        EXPECT_EQ(memoryStream.toString(), "");
    }
}

TEST(ConectTest, EthaneToplogy)
{
    gmx_mtop_t mtop;
    const bool bondsAreFlexible = true;
    addEthane(&mtop, bondsAreFlexible);
    mtop.finalize();
    auto atomsGuard = sg::make_scope_guard([&mtop]() { done_atom(&mtop.moltype[0].atoms); });

    gmx_conect conect              = gmx_conect_generate(mtop);
    auto       conectContentsGuard = sg::make_scope_guard([conect]() { gmx_conect_done(conect); });

    if (PosixMemstream::canCheckBufferContents())
    {
        PosixMemstream memoryStream;
        gmx_conect_dump(memoryStream.stream(), conect);
        EXPECT_EQ(memoryStream.toString(),
                  "CONECT    1    2\nCONECT    1    3\nCONECT    1    4\nCONECT    5    6\nCONECT  "
                  "  5    7\nCONECT    5    8\n");
    }
}

TEST(ConectTest, TwoWatersToplogy)
{
    gmx_mtop_t mtop;
    addNWaterMolecules(&mtop, 2);
    mtop.finalize();

    gmx_conect conect              = gmx_conect_generate(mtop);
    auto       conectContentsGuard = sg::make_scope_guard([conect]() { gmx_conect_done(conect); });

    if (PosixMemstream::canCheckBufferContents())
    {
        PosixMemstream memoryStream;
        gmx_conect_dump(memoryStream.stream(), conect);
        EXPECT_EQ(memoryStream.toString(), "") << "SETTLEs are not chemical bonds";
    }
}

TEST(ConectTest, TwoWatersAndEthaneToplogy)
{
    gmx_mtop_t mtop;
    // Work around broken semantics of gmx_moltype_t by avoiding
    // destruction on resize.
    mtop.moltype.reserve(2);
    addNWaterMolecules(&mtop, 2);
    const bool bondsAreFlexible = true;
    addEthane(&mtop, bondsAreFlexible);
    mtop.finalize();

    gmx_conect conect              = gmx_conect_generate(mtop);
    auto       conectContentsGuard = sg::make_scope_guard([conect]() { gmx_conect_done(conect); });

    if (PosixMemstream::canCheckBufferContents())
    {
        PosixMemstream memoryStream;
        gmx_conect_dump(memoryStream.stream(), conect);
        EXPECT_EQ(memoryStream.toString(),
                  "CONECT    7    8\nCONECT    7    9\nCONECT    7   10\nCONECT   11   12\nCONECT  "
                  " 11   13\nCONECT   11   14\n")
                << "two waters with no bonds, then ethane";
    }

    done_atom(&mtop.moltype[0].atoms);
}

} // namespace

} // namespace test

} // namespace gmx
