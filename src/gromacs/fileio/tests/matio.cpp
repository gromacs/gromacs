/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Implements tests for matrix file operations
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/fileio/matio.h"

#include <cstdio>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/rgb.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/stringtest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

std::string generateStdReferenceFile(gmx::ArrayRef<const t_mapping> refMaps)
{

    const std::string format       = "%c%c  %20s  %10g  %10g  %10g\n";
    std::string       fileContents = formatString("%ld\n", refMaps.ssize());
    for (const auto& map : refMaps)
    {
        fileContents.append(formatString(format.c_str(),
                                         map.code.c1 ? map.code.c1 : ' ',
                                         map.code.c2 ? map.code.c2 : ' ',
                                         map.desc.c_str(),
                                         map.rgb.r,
                                         map.rgb.g,
                                         map.rgb.b));
    }
    return fileContents;
}

std::vector<t_mapping> getReferenceMapping()
{
    std::vector<t_mapping> reference;
    reference.emplace_back(t_mapping{ { 'A', 'B' }, "101", { 0.23, 0.42, 0.1337 } });
    reference.emplace_back(t_mapping{ { 'A', 0 }, "111", { 0.32, 0.24, 0.3713 } });
    reference.emplace_back(t_mapping{ { 0, 'B' }, "110", { 0.23, 0.42, 0.1337 } });
    reference.emplace_back(t_mapping{ { 'C', 'B' }, "10", { 0.23, 0.42, 0.1337 } });
    return reference;
}

class ColorMapTest : public ::testing::Test
{
public:
    ColorMapTest()
    {
        referenceFilename_ = fileManager_.getTemporaryFilePath("ref.dat");
        testFilename_      = fileManager_.getTemporaryFilePath("test.dat");
        referenceMap_      = getReferenceMapping();
        useStringAsColorMapFile(generateStdReferenceFile(referenceMap()));
        writeColorMapFile();
    }

    const std::filesystem::path& referenceFilename() const { return referenceFilename_; }

    const std::filesystem::path& testFilename() const { return testFilename_; }

    const std::string& referenceContents() const { return referenceContents_; }

    ArrayRef<const t_mapping> referenceMap() const { return referenceMap_; }

    void useStringAsColorMapFile(const std::string& mapString) { referenceContents_ = mapString; }

    void writeColorMapFile() const
    {
        gmx::TextWriter::writeFileFromString(referenceFilename(), referenceContents());
    }

private:
    gmx::test::TestFileManager fileManager_;
    std::filesystem::path      testFilename_;
    std::filesystem::path      referenceFilename_;
    std::string                referenceContents_;
    std::vector<t_mapping>     referenceMap_;
};

TEST_F(ColorMapTest, CanReadFromFile)
{
    auto mappingData = readcmap(referenceFilename().c_str());
    ASSERT_EQ(mappingData.size(), referenceMap().size());
    for (int i = 0; i < gmx::ssize(mappingData); ++i)
    {
        // Depending on what is in the reference structure,
        // different values are set in the output for some reason
        // so we need to check for all possible combinations.
        if (referenceMap()[i].code.c1 != 0)
        {
            EXPECT_EQ(mappingData[i].code.c1, referenceMap()[i].code.c1);
        }
        else
        {
            EXPECT_EQ(mappingData[i].code.c1, referenceMap()[i].code.c2);
        }
        EXPECT_EQ(mappingData[i].code.c2, 0); // for some reason this is always 0 -.-
        EXPECT_EQ(mappingData[i].desc, referenceMap()[i].desc);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.r, referenceMap()[i].rgb.r);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.g, referenceMap()[i].rgb.g);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.b, referenceMap()[i].rgb.b);
    }
}

TEST_F(ColorMapTest, CanWriteToFile)
{
    FILE* out = fopen(testFilename().string().c_str(), "w");
    printcmap(out, referenceMap().size(), referenceMap().data());
    fclose(out);
    gmx::test::StringTestBase::testFilesEqual(referenceFilename(), testFilename());
}

TEST_F(ColorMapTest, RoundTrip)
{
    FILE* out = fopen(testFilename().string().c_str(), "w");
    printcmap(out, referenceMap().size(), referenceMap().data());
    fclose(out);
    auto mappingData = readcmap(testFilename().c_str());
    for (int i = 0; i < gmx::ssize(mappingData); ++i)
    {
        // Depending on what is in the reference structure,
        // different values are set in the output for some reason
        // so we need to check for all possible combinations.
        if (referenceMap()[i].code.c1 != 0)
        {
            EXPECT_EQ(mappingData[i].code.c1, referenceMap()[i].code.c1);
        }
        else
        {
            EXPECT_EQ(mappingData[i].code.c1, referenceMap()[i].code.c2);
        }
        EXPECT_EQ(mappingData[i].code.c2, 0); // for some reason this is always 0 -.-
        EXPECT_EQ(mappingData[i].desc, referenceMap()[i].desc);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.r, referenceMap()[i].rgb.r);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.g, referenceMap()[i].rgb.g);
        EXPECT_FLOAT_EQ(mappingData[i].rgb.b, referenceMap()[i].rgb.b);
    }
}

TEST_F(ColorMapTest, SearchWorks)
{
    t_xpmelmt lookup{ 'A', 0 };
    EXPECT_EQ(1, searchcmap(referenceMap(), lookup));
    lookup = { 'A', 'B' };
    EXPECT_EQ(0, searchcmap(referenceMap(), lookup));
    lookup = { 0, 'B' };
    EXPECT_EQ(2, searchcmap(referenceMap(), lookup));
    lookup = { 'A', 'A' };
    EXPECT_EQ(-1, searchcmap(referenceMap(), lookup));
}

t_matrix generateReferenceMatrix4x3()
{
    t_matrix matrix;
    matrix.nx      = 4;
    matrix.ny      = 3;
    matrix.title   = "Reference 4x3";
    matrix.legend  = "4x3";
    matrix.label_x = "X-Axis 4x3";
    matrix.label_y = "Y-Axis 4x3";
    matrix.axis_x  = { 1, 2, 3, 4 };
    matrix.axis_y  = { 6, 7, 8 };
    matrix.map     = getReferenceMapping();
    matrix.matrix.resize(matrix.nx, matrix.ny);
    for (int i = 0; i < matrix.nx; ++i)
    {
        for (int j = 0; j < matrix.ny; ++j)
        {
            matrix.matrix(i, j) = searchcmap(matrix.map, matrix.map[i].code);
        }
    }
    return matrix;
}

real** compareRealValues(const t_matrix& input, real** values, basic_mdspan<const float, extents<4, 3>> ref)
{
    values = matrix2real(&input, values);
    for (int i = 0; i < input.nx; ++i)
    {
        for (int j = 0; j < input.ny; ++j)
        {
            EXPECT_FLOAT_EQ(ref(i, j), values[i][j]) << "for i = " << i << " and j = " << j;
        }
    }
    return values;
}


class MatioTest : public ::testing::Test
{
public:
    MatioTest()
    {
        referenceFilename_ = fileManager_.getTemporaryFilePath("ref.xpm");
        refRealMat_        = { { 101, 101, 101, 111, 111, 111, 110, 110, 110, 10, 10, 10 } };
    }

    const std::filesystem::path& referenceFilename() const { return referenceFilename_; }

    basic_mdspan<const float, extents<4, 3>> refRealMat() const { return refRealMat_; }

private:
    t_matrix                                               referenceMatrix_;
    gmx::test::TestFileManager                             fileManager_;
    std::filesystem::path                                  testFilename_;
    std::filesystem::path                                  referenceFilename_;
    std::vector<t_mapping>                                 referenceMap_;
    MultiDimArray<std::array<float, 4 * 3>, extents<4, 3>> refRealMat_;
};

TEST_F(MatioTest, CanWriteToFile)
{
    auto  reference = generateReferenceMatrix4x3();
    FILE* out       = fopen(referenceFilename().string().c_str(), "w");
    EXPECT_NO_THROW(write_xpm_m(out, reference));
    fclose(out);
}

TEST_F(MatioTest, CanConvertToExistingRealMatrix)
{
    auto   reference = generateReferenceMatrix4x3();
    real** mat       = nullptr;
    snew(mat, reference.nx);
    for (int i = 0; i < reference.nx; ++i)
    {
        snew(mat[i], reference.ny);
    }
    compareRealValues(reference, mat, refRealMat());
    for (int i = 0; i < reference.nx; ++i)
    {
        sfree(mat[i]);
    }
    sfree(mat);
}

TEST_F(MatioTest, CanConvertToNewRealMatrix)
{
    auto   reference = generateReferenceMatrix4x3();
    real** newMat    = compareRealValues(reference, nullptr, refRealMat());
    for (int i = 0; i < reference.nx; ++i)
    {
        sfree(newMat[i]);
    }
    sfree(newMat);
}

TEST_F(MatioTest, CanReadSingleMatrixAfterWriting)
{
    auto reference = generateReferenceMatrix4x3();
    {
        FILE* out = fopen(referenceFilename().string().c_str(), "w");
        ASSERT_NO_THROW(write_xpm_m(out, reference));
        fclose(out);
    }

    std::vector<t_matrix> reading;
    ASSERT_NO_THROW(reading = read_xpm_matrix(referenceFilename().string().c_str()));
    ASSERT_EQ(reading.size(), 1);
    for (int i = 0; i < gmx::ssize(reading); ++i)
    {
        EXPECT_EQ(reading[i].nx, reference.nx);
        EXPECT_EQ(reading[i].ny, reference.ny);
        ASSERT_EQ(reading[i].axis_x.size(), reference.axis_x.size());
        for (int j = 0; j < gmx::ssize(reading[i].axis_x); ++j)
        {
            EXPECT_EQ(reading[i].axis_x[j], reference.axis_x[j]);
        }
        ASSERT_EQ(reading[i].axis_y.size(), reference.axis_y.size());
        for (int j = 0; j < gmx::ssize(reading[i].axis_y); ++j)
        {
            EXPECT_EQ(reading[i].axis_y[j], reference.axis_y[j]);
        }
        ASSERT_EQ(reading[i].map.size(), reference.map.size());
        for (int j = 0; j < gmx::ssize(reading[i].map); ++j)
        {
            EXPECT_EQ(reading[i].map[j].code, reference.map[j].code);
            EXPECT_EQ(reading[i].map[j].desc, reference.map[j].desc);
            EXPECT_FLOAT_EQ_TOL(reading[i].map[j].rgb.r,
                                reference.map[j].rgb.r,
                                relativeToleranceAsFloatingPoint(reading[i].map[j].rgb.r, 10e-3));
            EXPECT_FLOAT_EQ_TOL(reading[i].map[j].rgb.g,
                                reference.map[j].rgb.g,
                                relativeToleranceAsFloatingPoint(reading[i].map[j].rgb.r, 10e-3));
            EXPECT_FLOAT_EQ_TOL(reading[i].map[j].rgb.b,
                                reference.map[j].rgb.b,
                                relativeToleranceAsFloatingPoint(reading[i].map[j].rgb.r, 10e-3));
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
