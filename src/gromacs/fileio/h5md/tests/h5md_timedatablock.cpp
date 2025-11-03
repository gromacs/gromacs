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
 * \brief Tests for H5mdTimeDataBlock class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_timedatablock.h"

#include <hdf5.h>

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using H5mdTimeDataBlockTest = H5mdTestBase;

TEST_F(H5mdTimeDataBlockTest, BuilderWorks)
{
    constexpr char groupName[] = "float";

    {
        SCOPED_TRACE("Create all time data block sets");
        const H5mdTimeDataBlock<float> dataBlock =
                H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();
        ASSERT_EQ(dataBlock.numFrames(), 0)
                << "Sanity check failed: data block should be empty after creation";
    }
    {
        SCOPED_TRACE("Open the data sets individually to ensure that they were created");
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), groupName));
        EXPECT_NO_THROW(H5mdFrameDataSet<float>(group, "value"))
                << "Value data set was not created in group with given name";
        EXPECT_NO_THROW(H5mdScalarFrameDataSet<int64_t>(group, "step"))
                << "Step data set was not created in group with given name";
        EXPECT_NO_THROW(H5mdScalarFrameDataSet<double>(group, "time"))
                << "Time data set was not created in group with given name";
    }
}

TEST_F(H5mdTimeDataBlockTest, BuilderSetsFrameDimsToValueDataSet)
{
    constexpr char    groupName[] = "float";
    const DataSetDims frameDims   = { 5, 2 };

    H5mdTimeDataBlockBuilder<float>(fileid(), groupName).withFrameDimension(frameDims).build();

    const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), groupName));
    const H5mdFrameDataSet<float> valueDataSet(group, "value");
    EXPECT_EQ(valueDataSet.frameDims(), frameDims);
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorWorks)
{
    constexpr char groupName[] = "float";
    constexpr int  numFrames   = 2;

    {
        SCOPED_TRACE("Create data sets in group with given name and a few values");
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        const H5mdFrameDataSet<float> valueDataSet =
                H5mdFrameDataSetBuilder<float>(group, "value").withNumFrames(numFrames).build();
        const H5mdFrameDataSet<int64_t> stepDataSet =
                H5mdFrameDataSetBuilder<int64_t>(group, "step").withNumFrames(numFrames).build();
        const H5mdFrameDataSet<double> timeDataSet =
                H5mdFrameDataSetBuilder<double>(group, "time").withNumFrames(numFrames).build();
    }
    {
        SCOPED_TRACE("Open the above data sets as a collective time data block");
        H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);
        EXPECT_EQ(dataBlock.numFrames(), numFrames);
    }
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorThrowsIfValueDataSetIsMissing)
{
    constexpr char groupName[] = "float";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
    H5mdFrameDataSetBuilder<double>(group, "time").build();

    EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError);
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    EXPECT_NO_THROW(H5mdTimeDataBlock<float>(fileid(), groupName));
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorThrowsIfStepDataSetIsMissing)
{
    constexpr char groupName[] = "float";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<double>(group, "time").build();

    EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError);
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
    EXPECT_NO_THROW(H5mdTimeDataBlock<float>(fileid(), groupName));
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorWorksIfTimeDataSetIsMissing)
{
    constexpr char groupName[] = "float";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();

    EXPECT_NO_THROW(H5mdTimeDataBlock<float>(fileid(), groupName))
            << "Must not throw for missing time data set since it is optional";
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorThrowsForIncorrectTypes)
{

    {
        SCOPED_TRACE("Value data type does not match");
        constexpr char groupName[]     = "valueAsUnmatchedType";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        {
            // Create the data sets and close them at the end of scope.
            H5mdFrameDataSetBuilder<int64_t>(group, "value").build();
            H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
            H5mdFrameDataSetBuilder<double>(group, "time").build();
        }
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError)
                << "Must throw when opening a float time data block where the value type is int64";
    }
    {
        SCOPED_TRACE("Step data type does not match");
        constexpr char groupName[]     = "stepNotAsInt64";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        {
            // Create the data sets and close them at the end of scope.
            H5mdFrameDataSetBuilder<float>(group, "value").build();
            H5mdFrameDataSetBuilder<int32_t>(group, "step").build();
            H5mdFrameDataSetBuilder<double>(group, "time").build();
        }
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError)
                << "Must throw when opening a time data block where the step data set is not "
                   "int64_t";
    }
    {
        SCOPED_TRACE("Time data type does not match");
        constexpr char groupName[]     = "timeAsNotDouble";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        {
            // Create the data sets and close them at the end of scope.
            H5mdFrameDataSetBuilder<float>(group, "value").build();
            H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
            H5mdFrameDataSetBuilder<float>(group, "time").build();
        }
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError)
                << "Must throw when opening a time data block where the time data set is not "
                   "double";
    }
}

TEST_F(H5mdTimeDataBlockTest, BuilderThrowsForEmptyAndDuplicateName)
{
    EXPECT_THROW(H5mdTimeDataBlockBuilder<float>(fileid(), "").build(), FileIOError)
            << "Must throw for empty name";

    createGroup(fileid(), "goodName");
    EXPECT_THROW(H5mdTimeDataBlockBuilder<float>(fileid(), "goodName").build(), FileIOError)
            << "Must throw if group with name already exists";
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorThrowsForBadName)
{
    H5mdTimeDataBlockBuilder<float>(fileid(), "goodName").build();
    EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), "aBadName"), FileIOError);
    EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), ""), FileIOError);
    EXPECT_NO_THROW(H5mdTimeDataBlock<float>(fileid(), "goodName"));
}

TEST_F(H5mdTimeDataBlockTest, OpenConstructorThrowsIfDataSetsHaveDifferentNumFrames)
{
    {
        SCOPED_TRACE("Value data set has numFrames = 1, others = 0");
        constexpr char groupName[]     = "valueDataSetIsDifferent";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        H5mdFrameDataSetBuilder<float>(group, "value").withNumFrames(1).build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
        H5mdFrameDataSetBuilder<double>(group, "time").build();
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError);
    }
    {
        SCOPED_TRACE("Step data set has numFrames = 1, others = 0");
        constexpr char groupName[]     = "stepDataSetIsDifferent";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        H5mdFrameDataSetBuilder<float>(group, "value").build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").withNumFrames(1).build();
        H5mdFrameDataSetBuilder<double>(group, "time").build();
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError);
    }
    {
        SCOPED_TRACE("Time data set has numFrames = 1, others = 0");
        constexpr char groupName[]     = "timeDataSetIsDifferent";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        H5mdFrameDataSetBuilder<float>(group, "value").build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
        H5mdFrameDataSetBuilder<double>(group, "time").withNumFrames(1).build();
        EXPECT_THROW(H5mdTimeDataBlock<float>(fileid(), groupName), FileIOError);
    }
    {
        SCOPED_TRACE(
                "Edge case: if no time data set exists only the value and step data sets are "
                "checked");
        constexpr char groupName[]     = "noTimeDataSet";
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
        H5mdFrameDataSetBuilder<float>(group, "value").withNumFrames(1).build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").withNumFrames(1).build();
        EXPECT_NO_THROW(H5mdTimeDataBlock<float>(fileid(), groupName));
    }
}

TEST_F(H5mdTimeDataBlockTest, WriteNextFrameNoTimeOverloadWorks)
{
    constexpr char             groupName[] = "float";
    const std::vector<float>   values      = { 5.0, 7.0 };
    const std::vector<int64_t> steps       = { 40, 5000 };
    const std::vector<double>  times       = { 0.1, 10.1 };
    const hsize_t              numFrames   = values.size();

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
    H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);

    for (hsize_t i = 0; i < numFrames; ++i)
    {
        dataBlock.writeNextFrame(constArrayRefFromArray(&values[i], 1), steps[i], times[i]);
        EXPECT_EQ(dataBlock.numFrames(), i + 1);
    }

    H5mdScalarFrameDataSet<float> valueDataSet(group, "value");
    EXPECT_EQ(valueDataSet.numFrames(), numFrames);
    H5mdScalarFrameDataSet<int64_t> stepDataSet(group, "step");
    EXPECT_EQ(stepDataSet.numFrames(), numFrames);

    for (hsize_t i = 0; i < numFrames; ++i)
    {
        float value;
        valueDataSet.readFrame(i, &value);
        EXPECT_FLOAT_EQ(value, values[i]);
        int64_t step;
        stepDataSet.readFrame(i, &step);
        EXPECT_EQ(step, steps[i]);
    }
}

TEST_F(H5mdTimeDataBlockTest, WriteNextFrameNoTimeOverloadThrowsIfTimeIsManaged)
{
    constexpr char groupName[] = "float";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
    H5mdFrameDataSetBuilder<double>(group, "time").build();
    H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);

    constexpr float valueToWrite = 0.0;
    ASSERT_NO_THROW(dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), 0, 0))
            << "Sanity check: Must not throw when passing a time value";
    EXPECT_THROW(dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), 0), gmx::FileIOError)
            << "Must throw when using no-time writeNextFrame with a managed time data set";
}

TEST_F(H5mdTimeDataBlockTest, WriteNextFrameWithTimeWorks)
{
    constexpr char             groupName[] = "float";
    const std::vector<float>   values      = { 5.0, 7.0 };
    const std::vector<int64_t> steps       = { 40, 5000 };
    const std::vector<double>  times       = { 0.1, 10.1 };
    const hsize_t              numFrames   = values.size();

    H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();
    for (hsize_t i = 0; i < numFrames; ++i)
    {
        dataBlock.writeNextFrame(constArrayRefFromArray(&values[i], 1), steps[i], times[i]);
        EXPECT_EQ(dataBlock.numFrames(), i + 1);
    }

    const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), groupName));
    H5mdScalarFrameDataSet<float> valueDataSet(group, "value");
    EXPECT_EQ(valueDataSet.numFrames(), numFrames);
    H5mdScalarFrameDataSet<int64_t> stepDataSet(group, "step");
    EXPECT_EQ(stepDataSet.numFrames(), numFrames);
    H5mdScalarFrameDataSet<double> timeDataSet(group, "time");
    EXPECT_EQ(timeDataSet.numFrames(), numFrames);

    for (hsize_t i = 0; i < numFrames; ++i)
    {
        float value;
        valueDataSet.readFrame(i, &value);
        EXPECT_FLOAT_EQ(value, values[i]);
        int64_t step;
        stepDataSet.readFrame(i, &step);
        EXPECT_EQ(step, steps[i]);
        double time;
        timeDataSet.readFrame(i, &time);
        EXPECT_FLOAT_EQ(time, times[i]);
    }
}

TEST_F(H5mdTimeDataBlockTest, WriteNextFrameWithTimeWorksIfNoTimeDataSetIsManaged)
{
    constexpr char groupName[]     = "float";
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();

    H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);

    constexpr float frameValue = 10;
    dataBlock.writeNextFrame(constArrayRefFromArray(&frameValue, 1), frameValue, -1.0);

    H5mdScalarFrameDataSet<float> valueDataSet(group, "value");
    EXPECT_EQ(valueDataSet.numFrames(), 1);
    H5mdScalarFrameDataSet<int64_t> stepDataSet(group, "step");
    EXPECT_EQ(stepDataSet.numFrames(), 1);

    float value;
    valueDataSet.readFrame(0, &value);
    EXPECT_FLOAT_EQ(value, frameValue);
    int64_t step;
    stepDataSet.readFrame(0, &step);
    EXPECT_EQ(step, frameValue);
}

TEST_F(H5mdTimeDataBlockTest, WriteNextFrameWorksAfterOpen)
{
    constexpr char groupName[] = "float";

    {
        SCOPED_TRACE("Create time data block and write one frame to it before closing.");

        H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();
        constexpr float frameValue = 10;
        dataBlock.writeNextFrame(constArrayRefFromArray(&frameValue, 1), frameValue, frameValue);
    }

    constexpr float valueToAppend = 500;
    {
        SCOPED_TRACE("Reopen the time data block and append a value");
        H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);
        EXPECT_EQ(dataBlock.numFrames(), 1);
        dataBlock.writeNextFrame(constArrayRefFromArray(&valueToAppend, 1), valueToAppend, valueToAppend);
    }
    {
        SCOPED_TRACE("Open the data sets and assert that the values were appended correctly");

        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), groupName));
        H5mdScalarFrameDataSet<float> valueDataSet(group, "value");
        EXPECT_EQ(valueDataSet.numFrames(), 2);
        H5mdScalarFrameDataSet<int64_t> stepDataSet(group, "step");
        EXPECT_EQ(stepDataSet.numFrames(), 2);
        H5mdScalarFrameDataSet<double> timeDataSet(group, "time");
        EXPECT_EQ(timeDataSet.numFrames(), 2);

        float value;
        valueDataSet.readFrame(1, &value);
        EXPECT_FLOAT_EQ(value, valueToAppend);
        int64_t step;
        stepDataSet.readFrame(1, &step);
        EXPECT_EQ(step, valueToAppend);
        double time;
        timeDataSet.readFrame(1, &time);
        EXPECT_FLOAT_EQ(time, valueToAppend);
    }
}

TEST_F(H5mdTimeDataBlockTest, ReadValueAtIndexWorks)
{
    constexpr char groupName[] = "float";
    H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();

    const std::vector<float> valuesToWrite = { 5.0, 11.0, 23.0 };
    for (const auto v : valuesToWrite)
    {
        dataBlock.writeNextFrame(constArrayRefFromArray(&v, 1), 0, 0);
    }

    float readValueBuffer;

    // Read the values out of order to test non-sequential reading
    for (const int i : { 1, 2, 0 })
    {
        EXPECT_TRUE(dataBlock.readValueAtIndex(i, arrayRefFromArray(&readValueBuffer, 1)));
        EXPECT_FLOAT_EQ(readValueBuffer, valuesToWrite[i]);
    }
}

TEST_F(H5mdTimeDataBlockTest, ReadValueAtIndexReturnsFalseIfOutOfBounds)
{
    constexpr char groupName[] = "float";
    H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();

    // Write a single frame to the data sets
    constexpr float valueToWrite = 5.0;
    dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), 0, 0);

    float readValueBuffer;
    EXPECT_FALSE(dataBlock.readValueAtIndex(-1, arrayRefFromArray(&readValueBuffer, 1)))
            << "Must return false for negative frame index";
    EXPECT_FALSE(dataBlock.readValueAtIndex(1, arrayRefFromArray(&readValueBuffer, 1)))
            << "Must return false for frame index >= numFrames";
    EXPECT_TRUE(dataBlock.readValueAtIndex(0, arrayRefFromArray(&readValueBuffer, 1)))
            << "Must return true for valid frame index";
}

TEST_F(H5mdTimeDataBlockTest, ReadStepAndTimeAtIndexWorks)
{
    constexpr char groupName[] = "float";
    H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();

    constexpr hsize_t          numFrames = 3;
    const std::vector<int64_t> steps     = { 10, 501, 10000 };
    const std::vector<double>  times     = { 0.1, 10.1, 5011.35 };
    for (hsize_t i = 0; i < numFrames; ++i)
    {
        constexpr float valueToWrite = 0.0; // dummy value for write operation
        dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), steps[i], times[i]);
    }

    // Read the values out of order to test non-sequential reading
    for (const int i : { 1, 2, 0 })
    {
        EXPECT_EQ(dataBlock.readStepAtIndex(i).value(), steps[i]);
        EXPECT_FLOAT_EQ(dataBlock.readTimeAtIndex(i).value(), times[i]);
    }
}

TEST_F(H5mdTimeDataBlockTest, ReadStepAndTimeAtIndexReturnsNulloptIfOutOfBounds)
{
    constexpr char groupName[] = "float";
    H5mdTimeDataBlock<float> dataBlock = H5mdTimeDataBlockBuilder<float>(fileid(), groupName).build();

    constexpr int64_t stepToWrite  = 11;
    constexpr double  timeToWrite  = 5.1;
    constexpr float   valueToWrite = 0.0; // dummy value for write operation
    dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), stepToWrite, timeToWrite);

    EXPECT_FALSE(dataBlock.readStepAtIndex(-1).has_value());
    EXPECT_FALSE(dataBlock.readStepAtIndex(1).has_value());
    EXPECT_EQ(dataBlock.readStepAtIndex(0).value(), stepToWrite);
    EXPECT_FALSE(dataBlock.readTimeAtIndex(-1).has_value());
    EXPECT_FALSE(dataBlock.readTimeAtIndex(1).has_value());
    EXPECT_EQ(dataBlock.readTimeAtIndex(0).value(), timeToWrite);
}

TEST_F(H5mdTimeDataBlockTest, ReadTimeAtIndexReturnsNulloptIfNotManaged)
{
    constexpr char groupName[]     = "float";
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName));
    H5mdFrameDataSetBuilder<float>(group, "value").build();
    H5mdFrameDataSetBuilder<int64_t>(group, "step").build();

    H5mdTimeDataBlock<float> dataBlock(fileid(), groupName);

    constexpr float   frameValue = 10;
    constexpr int64_t stepValue  = 5;
    constexpr double  timeValue  = 5.1;
    dataBlock.writeNextFrame(constArrayRefFromArray(&frameValue, 1), stepValue, timeValue);

    ASSERT_EQ(dataBlock.readStepAtIndex(0), stepValue)
            << "Sanity check failed: no values were written";
    EXPECT_FALSE(dataBlock.readTimeAtIndex(0).has_value())
            << "Should not return a value if not managing a time data set";
}

TEST_F(H5mdTimeDataBlockTest, WorksForComplexFrameData)
{
    // Read and write frames to a data set of RVec[numFrames][numAtoms] type (e.g. positions)
    constexpr char                        groupName[] = "RVec";
    constexpr hsize_t                     numFrames   = 2;
    constexpr hsize_t                     numAtoms    = 5;
    const DataSetDims                     frameDims   = { numAtoms };
    H5mdTimeDataBlock<BasicVector<float>> dataBlock =
            H5mdTimeDataBlockBuilder<BasicVector<float>>(fileid(), groupName)
                    .withFrameDimension(frameDims)
                    .build();

    std::array<std::array<BasicVector<float>, numAtoms>, numFrames> framesOfData;

    // Fill all frames with some data
    for (hsize_t i = 0; i < numFrames; ++i)
    {
        for (hsize_t j = 0; j < numAtoms; ++j)
        {
            framesOfData[i][j] = { static_cast<float>((i * 100) + (j * 10)) + XX,
                                   static_cast<float>((i * 100) + (j * 10)) + YY,
                                   static_cast<float>((i * 100) + (j * 10)) + ZZ };
        }
    }

    for (hsize_t i = 0; i < numFrames; ++i)
    {
        dataBlock.writeNextFrame(framesOfData[i], i, i);
    }

    std::array<BasicVector<float>, numAtoms> readBuffer;
    EXPECT_TRUE(dataBlock.readValueAtIndex(0, readBuffer));
    EXPECT_EQ(readBuffer, framesOfData[0]);
    EXPECT_TRUE(dataBlock.readValueAtIndex(1, readBuffer));
    EXPECT_EQ(readBuffer, framesOfData[1]);
}

TEST_F(H5mdTimeDataBlockTest, HasTimeWorks)
{
    using ValueType                = float;
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "block"));
    {
        SCOPED_TRACE("Create value and step data sets, hasTime() should return false after open");
        H5mdFrameDataSetBuilder<ValueType>(group, "value").build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").build();

        H5mdTimeDataBlock<ValueType> dataBlock(fileid(), "block");
        EXPECT_FALSE(dataBlock.hasTime());
    }
    {
        SCOPED_TRACE("Create time data set, hasTime() should return true after open");
        H5mdFrameDataSetBuilder<double>(group, "time").build();

        H5mdTimeDataBlock<ValueType> dataBlock(fileid(), "block");
        EXPECT_TRUE(dataBlock.hasTime());
    }
}

TEST_F(H5mdTimeDataBlockTest, ReadFrameWorks)
{
    using ValueType                        = float;
    constexpr int                numValues = 5;
    H5mdTimeDataBlock<ValueType> dataBlock =
            H5mdTimeDataBlockBuilder<ValueType>(fileid(), "testBlockName")
                    .withFrameDimension({ numValues })
                    .build();

    // Prepare unique data for two frames
    const std::vector<ValueType> valuesFrame1 = [&]()
    {
        std::vector<ValueType> values(numValues, 0.0);
        std::iota(values.begin(), values.end(), 0);
        return values;
    }();
    constexpr int64_t            stepFrame1   = 5050;
    constexpr double             timeFrame1   = -5.0;
    const std::vector<ValueType> valuesFrame2 = [&]()
    {
        std::vector<ValueType> values(numValues, 0.0);
        std::iota(values.begin(), values.end(), numValues);
        return values;
    }();
    constexpr int64_t stepFrame2 = 13;
    constexpr double  timeFrame2 = 1030.0;

    // Write the frames
    dataBlock.writeNextFrame(valuesFrame1, stepFrame1, timeFrame1);
    dataBlock.writeNextFrame(valuesFrame2, stepFrame2, timeFrame2);

    std::vector<ValueType> readValuesBuffer(numValues);
    int64_t                readStepBuffer;
    double                 readTimeBuffer;

    {
        SCOPED_TRACE("Read first frame");
        EXPECT_TRUE(dataBlock.readFrame(0, readValuesBuffer, &readStepBuffer, &readTimeBuffer))
                << "readFrame should return true if a frame is read";
        EXPECT_EQ(readValuesBuffer, valuesFrame1);
        EXPECT_EQ(readStepBuffer, stepFrame1);
        EXPECT_EQ(readTimeBuffer, timeFrame1);
    }
    {
        SCOPED_TRACE("Read second frame");
        EXPECT_TRUE(dataBlock.readFrame(1, readValuesBuffer, &readStepBuffer, &readTimeBuffer))
                << "readFrame should return true if a frame is read";
        EXPECT_EQ(readValuesBuffer, valuesFrame2);
        EXPECT_EQ(readStepBuffer, stepFrame2);
        EXPECT_EQ(readTimeBuffer, timeFrame2);
    }
    {
        SCOPED_TRACE("Read invalid indices");
        EXPECT_FALSE(dataBlock.readFrame(2, readValuesBuffer, &readStepBuffer, &readTimeBuffer))
                << "readFrame should return false if a frame index is invalid";
    }
}

TEST_F(H5mdTimeDataBlockTest, ReadFrameWorksWithoutTimeDataSet)
{
    using ValueType                = float;
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "block"));

    {
        SCOPED_TRACE("Create value and step data sets and close them at the end of scope");
        H5mdFrameDataSetBuilder<ValueType>(group, "value").build();
        H5mdFrameDataSetBuilder<int64_t>(group, "step").build();
    }

    H5mdTimeDataBlock<ValueType> dataBlock(fileid(), "block");

    constexpr ValueType valueToWrite = 5;
    constexpr int64_t   stepToWrite  = 5050;
    dataBlock.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1), stepToWrite, 0.0);

    ValueType        readValueBuffer;
    int64_t          readStepBuffer;
    constexpr double timeSentinelValue = 501.0;
    double readTimeBuffer = timeSentinelValue; // sentinel value: should not change after read

    EXPECT_TRUE(dataBlock.readFrame(
            0, arrayRefFromArray(&readValueBuffer, 1), &readStepBuffer, &readTimeBuffer));
    EXPECT_EQ(readValueBuffer, valueToWrite);
    EXPECT_EQ(readStepBuffer, stepToWrite);
    EXPECT_EQ(readTimeBuffer, timeSentinelValue) << "Input time should not be modified";
}

TEST_F(H5mdTimeDataBlockTest, ReadFrameDoesNotReadIntoNullPtrInput)
{
    using ValueType = float;

    H5mdTimeDataBlock<ValueType> dataBlock =
            H5mdTimeDataBlockBuilder<ValueType>(fileid(), "block").build();

    // Write two frames
    constexpr ValueType valueFrame1 = 5;
    dataBlock.writeNextFrame(constArrayRefFromArray(&valueFrame1, 1), -1, -1.0);
    constexpr ValueType valueFrame2 = 3;
    dataBlock.writeNextFrame(constArrayRefFromArray(&valueFrame2, 1), -1, -1.0);

    ValueType readValueBuffer;
    EXPECT_TRUE(dataBlock.readFrame(0, arrayRefFromArray(&readValueBuffer, 1), nullptr, nullptr));
    EXPECT_EQ(readValueBuffer, valueFrame1);
}

TEST_F(H5mdTimeDataBlockTest, DefaultCompression)
{
    using ValueType = double;
    {
        SCOPED_TRACE("Create data sets and close them at the end of scope");
        H5mdTimeDataBlock<ValueType> dataBlockUncompressed =
                H5mdTimeDataBlockBuilder<ValueType>(fileid(), "block").build();
    }
    {
        SCOPED_TRACE("Check value data set default compression");
        H5mdDataSetBase<ValueType> dataSet(fileid(), "block/value");
        const auto [propertyList, propertyListGuard] =
                makeH5mdPropertyListGuard(H5Dget_create_plist(dataSet.id()));

        EXPECT_EQ(H5Pget_nfilters(propertyList), 0)
                << "No filters should be applied to values by default";
    }
    {
        SCOPED_TRACE("Check step data set default compression");
        H5mdDataSetBase<int64_t> dataSet(fileid(), "block/step");
        const auto [propertyList, propertyListGuard] =
                makeH5mdPropertyListGuard(H5Dget_create_plist(dataSet.id()));

        EXPECT_GT(H5Pget_nfilters(propertyList), 0)
                << "At least one filter (deflate) should be applied to steps by default";
        // H5Pget_filter_by_id2 returns: >=0 if the filter (second argument) is set, else <0
        EXPECT_GE(H5Pget_filter_by_id2(propertyList, H5Z_FILTER_DEFLATE, nullptr, nullptr, nullptr, 0, nullptr, nullptr), 0)
                << "Deflate filter should be applied";
    }
    {
        SCOPED_TRACE("Check time data set default compression");
        H5mdDataSetBase<double> dataSet(fileid(), "block/time");
        const auto [propertyList, propertyListGuard] =
                makeH5mdPropertyListGuard(H5Dget_create_plist(dataSet.id()));

        EXPECT_GT(H5Pget_nfilters(propertyList), 0)
                << "At least one filter (deflate) should be applied to times by default";
        // H5Pget_filter_by_id2 returns: >=0 if the filter (second argument) is set, else <0
        EXPECT_GE(H5Pget_filter_by_id2(propertyList, H5Z_FILTER_DEFLATE, nullptr, nullptr, nullptr, 0, nullptr, nullptr), 0)
                << "Deflate filter should be applied";
    }
}

TEST_F(H5mdTimeDataBlockTest, CompressionWorks)
{
    using ValueType         = double;
    constexpr int numValues = 101;

    {
        SCOPED_TRACE(
                "Create data sets, write some data to them and close them at the end of scope");
        H5mdTimeDataBlock<ValueType> dataBlockUncompressed =
                H5mdTimeDataBlockBuilder<ValueType>(fileid(), "uncompressed")
                        .withFrameDimension({ numValues })
                        .build();

        H5mdTimeDataBlock<ValueType> dataBlockCompressed =
                H5mdTimeDataBlockBuilder<ValueType>(fileid(), "compressed")
                        .withFrameDimension({ numValues })
                        .withCompression(H5mdCompression::LosslessShuffle)
                        .build();

        // Write some data to the data sets
        std::vector<ValueType> valuesToWrite(numValues);
        std::iota(valuesToWrite.begin(), valuesToWrite.end(), 0.0);
        dataBlockUncompressed.writeNextFrame(valuesToWrite, 0, 0.0);
        dataBlockCompressed.writeNextFrame(valuesToWrite, 0, 0.0);
    }

    // Open the value data sets as H5mdDataSetBase to access the hid_t handles for below
    H5mdDataSetBase<ValueType> valueDataSetUncompressed(fileid(), "uncompressed/value");
    H5mdDataSetBase<ValueType> valueDataSetCompressed(fileid(), "compressed/value");

    // H5Dget_storage_size returns the size (in bytes) of the stored data inside the set.
    const size_t uncompressedSize = H5Dget_storage_size(valueDataSetUncompressed.id());
    const size_t compressedSize   = H5Dget_storage_size(valueDataSetCompressed.id());
    EXPECT_LT(compressedSize, 0.5 * uncompressedSize)
            << "gzip compression should yield a saving of much more than 50% for this data";
}

} // namespace
} // namespace test
} // namespace gmx
