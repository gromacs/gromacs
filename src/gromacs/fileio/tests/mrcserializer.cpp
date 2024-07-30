/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests for mrc file data structure.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/mrcserializer.h"

#include <cstddef>

#include <array>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/inmemoryserializer.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(MrcSerializer, DefaultHeaderValuesAreSet)
{
    const MrcDensityMapHeader inputHeader = {};

    EXPECT_EQ('M', inputHeader.formatIdentifier_[0]);
    EXPECT_EQ('A', inputHeader.formatIdentifier_[1]);
    EXPECT_EQ('P', inputHeader.formatIdentifier_[2]);
    EXPECT_EQ(' ', inputHeader.formatIdentifier_[3]);
}

TEST(MrcSerializer, DefaultHeaderHasRightSerialSize)
{

    InMemorySerializer        serializer;
    const MrcDensityMapHeader inputHeader = {};

    serializeMrcDensityMapHeader(&serializer, inputHeader);
    const auto serializedHeader = serializer.finishAndGetBuffer();

    constexpr size_t c_defaultMrcHeaderSize = 1024;
    EXPECT_EQ(c_defaultMrcHeaderSize, serializedHeader.size());
}

TEST(MrcSerializer, DefaultHeaderIdenticalAfterRoundTrip)
{
    InMemorySerializer        serializer;
    const MrcDensityMapHeader inputHeader = {};

    serializeMrcDensityMapHeader(&serializer, inputHeader);
    const auto serializedHeader = serializer.finishAndGetBuffer();

    InMemoryDeserializer deserializer(serializedHeader, false);
    const auto           deserializedHeader = deserializeMrcDensityMapHeader(&deserializer);

    // comparing serialized results saves MrcDensityHeaders comparison implementation
    serializeMrcDensityMapHeader(&serializer, deserializedHeader);
    const auto roundTripResult = serializer.finishAndGetBuffer();
    EXPECT_THAT(serializedHeader, testing::Pointwise(testing::Eq(), roundTripResult));
}

} // namespace
} // namespace test
} // namespace gmx
