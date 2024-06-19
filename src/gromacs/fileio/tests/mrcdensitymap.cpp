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
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/mrcdensitymap.h"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/inmemoryserializer.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(MrcDensityMap, RoundTripIsIdempotent)
{
    // write header and data to serializer, store the serialized data

    MrcDensityMapHeader header{};
    header.numColumnRowSection_ = { 1, 1, 1 };

    std::vector<float> data(numberOfExpectedDataItems(header));

    MrcDensityMapOfFloatWriter mrcDensityMapWriter(header, data);
    InMemorySerializer         serializer;
    mrcDensityMapWriter.write(&serializer);

    const auto serializedFile = serializer.finishAndGetBuffer();

    // read and write again
    InMemoryDeserializer       deserializer(serializedFile, false);
    MrcDensityMapOfFloatReader mrcDensityMapReader(&deserializer);

    MrcDensityMapOfFloatWriter writerOfDeserializedOutput(mrcDensityMapReader.header(),
                                                          mrcDensityMapReader.constView());
    writerOfDeserializedOutput.write(&serializer);

    const auto roundTripResult = serializer.finishAndGetBuffer();

    // compare serialized data after reading and writing again to serialized data after one-time serialization
    EXPECT_THAT(serializedFile, testing::Pointwise(testing::Eq(), roundTripResult));
}

TEST(MrcDensityMap, ThrowsFileIOErrorWhenFileNotPresent)
{
    EXPECT_THROW(MrcDensityMapOfFloatFromFileReader(""), FileIOError);
}

TEST(MrcDensityMap, ReadsCoordinateTransformationFromFile)
{
    RVec coordinate1(0, 0, 0);
    RVec coordinate2(0, 0, 1);
    RVec coordinate3(1, 0, 0);
    RVec coordinate4(0, 1, 0);

    MrcDensityMapOfFloatFromFileReader mrcFileReader(
            TestFileManager::getInputFilePath("ellipsoid-density.mrc"));
    TranslateAndScale coordinateTransformation(mrcFileReader.transformationToDensityLattice());

    coordinateTransformation({ &coordinate1, &coordinate1 + 1 });
    coordinateTransformation({ &coordinate2, &coordinate2 + 1 });

    EXPECT_REAL_EQ(0, coordinate1[XX]);
    EXPECT_REAL_EQ(-2, coordinate1[YY]);
    EXPECT_REAL_EQ(0, coordinate1[ZZ]);

    EXPECT_REAL_EQ(0, coordinate2[XX]);
    EXPECT_REAL_EQ(-2, coordinate2[YY]);
    EXPECT_REAL_EQ(1.25, coordinate2[ZZ]);

    EXPECT_REAL_EQ(1, coordinate3[XX]);
    EXPECT_REAL_EQ(0, coordinate3[YY]);
    EXPECT_REAL_EQ(0, coordinate3[ZZ]);

    EXPECT_REAL_EQ(0, coordinate4[XX]);
    EXPECT_REAL_EQ(1, coordinate4[YY]);
    EXPECT_REAL_EQ(0, coordinate4[ZZ]);
}

TEST(MrcDensityMap, ReadsDensityDataFromFile)
{
    MrcDensityMapOfFloatFromFileReader mrcFileReader(
            TestFileManager::getInputFilePath("ellipsoid-density.mrc"));
    const auto densityData = mrcFileReader.densityDataCopy();

    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkSequence(begin(densityData.asConstView()),
                          end(densityData.asConstView()),
                          "data ellipsoid density");
}

} // namespace
} // namespace test
} // namespace gmx
