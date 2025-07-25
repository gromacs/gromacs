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
/*! \internal \file
 * \brief
 * Tests for PosixMemstream
 *
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/posixmemstream.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

TEST(PosixMemstreamTest, ConstructsAndDestructs)
{
    PosixMemstream memstream;
}

TEST(PosixMemstreamTest, HasAStreamUntilClosed)
{
    PosixMemstream memstream;
    EXPECT_NE(memstream.stream(), nullptr);
    memstream.closeStream();
    EXPECT_EQ(memstream.stream(), nullptr);
    // Also tests that destruction works after closing the stream
}

TEST(PosixMemstreamTest, CanFetchString)
{
    PosixMemstream memstream;
    memstream.toString();
}

TEST(PosixMemstreamTest, HasAStreamUntilStringIsFetched)
{
    PosixMemstream memstream;
    EXPECT_NE(memstream.stream(), nullptr);
    memstream.toString();
    EXPECT_EQ(memstream.stream(), nullptr);
}

TEST(PosixMemstreamTest, CanWriteToStream)
{
    PosixMemstream memstream;
    fprintf(memstream.stream(), "Hello world\n");
}

TEST(PosixMemstreamTest, CanCheckBufferWithContents)
{
    PosixMemstream memstream;
    fprintf(memstream.stream(), "Hello world\n");
    if (memstream.canCheckBufferContents())
    {
        EXPECT_EQ(memstream.toString(), "Hello world\n");
    }
}

} // namespace
} // namespace test
} // namespace gmx
