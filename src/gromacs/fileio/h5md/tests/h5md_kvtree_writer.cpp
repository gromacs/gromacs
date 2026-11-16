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
 * \brief Tests for KeyValueTree writing to H5MD
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_kvtree_writer.h"

#include <hdf5.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "testutils/testasserts.h"


namespace gmx
{
namespace test
{
namespace
{

std::vector<std::string> readFixedSizeStringData(const hid_t container, const char* name)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(H5Dopen(container, name, H5P_DEFAULT));
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    // Get the maximum string size (including the terminating '\0')
    const size_t maxStringSize = H5Tget_size(dataType);

    // Get the total number of values in data set
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    const int            numDims           = H5Sget_simple_extent_ndims(dataSpace);
    std::vector<hsize_t> dims(numDims, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    int numValues = 1;
    for (const hsize_t d : dims)
    {
        numValues *= d;
    }

    std::vector<char> readBuffer(numValues * maxStringSize);

    GMX_H5MD_THROW_UPON_ERROR(
            H5Dread(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, readBuffer.data()) < 0,
            "Error reading fixed-size string data.");

    std::vector<std::string> stringValues;
    stringValues.reserve(numValues);
    auto startOfNextString = readBuffer.cbegin();
    for (int i = 0; i < numValues; ++i)
    {
        stringValues.emplace_back(startOfNextString,
                                  std::find(startOfNextString, startOfNextString + maxStringSize, '\0'));
        startOfNextString += maxStringSize;
    }
    return stringValues;
}

using H5mdKVTreeWriterTest = H5mdTestBase;

TEST_F(H5mdKVTreeWriterTest, WriteSimpleScalarProperties)
{
    const int64_t     numStepsWritten              = 1000;
    const std::string numStepsDescriptionWritten   = "Number of steps to run";
    const double      timestepWritten              = 0.002;
    const std::string timestepUnitWritten          = "ps";
    const std::string timestepDescriptionWritten   = "Integration time step";
    const std::string integratorWritten            = "velocity-verlet";
    const std::string integratorDescriptionWritten = "Integration algorithm";
    const std::string integratorDoiWritten         = "10.xxxx/example.doi";
    const int32_t     factorWritten                = 3;
    const std::string factorDescriptionWritten     = "A simple integer value";
    const float       toleranceWritten             = 0.5F;
    const std::string toleranceUnitWritten         = "kJ/mol";
    const std::string toleranceDescriptionWritten  = "A simple float value";

    {
        SCOPED_TRACE("Write KeyValueTree to H5MD");

        // Build a key-value tree similar to TOML structure
        KeyValueTreeBuilder builder;
        auto                root       = builder.rootObject();
        auto                runControl = root.addObject("run-control");

        auto nsteps = runControl.addObject("nsteps");
        nsteps.addValue<int64_t>("value", numStepsWritten);
        nsteps.addValue<std::string>("description", numStepsDescriptionWritten);

        auto timestepObject = runControl.addObject("timestep");
        timestepObject.addValue<double>("value", timestepWritten);
        timestepObject.addValue<std::string>("unit", timestepUnitWritten);
        timestepObject.addValue<std::string>("description", timestepDescriptionWritten);

        auto algorithm = root.addObject("integrator");
        algorithm.addValue<std::string>("value", integratorWritten);
        algorithm.addValue<std::string>("description", integratorDescriptionWritten);
        algorithm.addValue<std::string>("doi", integratorDoiWritten);

        auto factor = root.addObject("factor");
        factor.addValue<int32_t>("value", factorWritten);
        factor.addValue<std::string>("description", factorDescriptionWritten);

        auto tolerance = root.addObject("tolerance");
        tolerance.addValue<float>("value", toleranceWritten);
        tolerance.addValue<std::string>("unit", toleranceUnitWritten);
        tolerance.addValue<std::string>("description", toleranceDescriptionWritten);

        writeKeyValueTreeToH5md(fileid(), builder.build());
    }

    {
        SCOPED_TRACE("Read back and verify the written data");

        const auto [group1, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), "run-control"));
        ASSERT_TRUE(objectExists(group1, "nsteps"));
        const auto nstepsDataset = H5mdFixedDataSet<int64_t>(group1, "nsteps");
        int64_t    numStepsRead;
        nstepsDataset.readData(arrayRefFromArray(&numStepsRead, 1));
        EXPECT_EQ(numStepsRead, numStepsWritten);
        const auto [nstepsHandle, nstepsGuard] =
                makeH5mdDataSetGuard(H5Dopen(group1, "nsteps", H5P_DEFAULT));
        EXPECT_THAT(getAttribute<std::string>(nstepsHandle, "description"),
                    ::testing::Optional(numStepsDescriptionWritten));

        const auto timestepDataset = H5mdFixedDataSet<double>(group1, "timestep");
        double     timestepRead;
        timestepDataset.readData(arrayRefFromArray(&timestepRead, 1));
        EXPECT_DOUBLE_EQ(timestepRead, timestepWritten);
        const auto [timestepHandle, timestepGuard] =
                makeH5mdDataSetGuard(H5Dopen(group1, "timestep", H5P_DEFAULT));
        EXPECT_THAT(getAttribute<std::string>(timestepHandle, "unit"),
                    ::testing::Optional(timestepUnitWritten));
        EXPECT_THAT(getAttribute<std::string>(timestepHandle, "description"),
                    ::testing::Optional(timestepDescriptionWritten));

        // integrator: root-level string dataset with description and doi as attributes
        const std::vector<std::string> integratorRead =
                readFixedSizeStringData(fileid(), "integrator");
        ASSERT_EQ(integratorRead.size(), 1U);
        EXPECT_EQ(integratorRead[0], integratorWritten);
        const auto [integratorHandle, integratorGuard] =
                makeH5mdDataSetGuard(H5Dopen(fileid(), "integrator", H5P_DEFAULT));
        EXPECT_THAT(getAttribute<std::string>(integratorHandle, "description"),
                    ::testing::Optional(integratorDescriptionWritten));
        EXPECT_THAT(getAttribute<std::string>(integratorHandle, "doi"),
                    ::testing::Optional(integratorDoiWritten));

        // factor: root-level int32 dataset with description as attribute
        auto    factorDataset = H5mdFixedDataSet<int32_t>(fileid(), "factor");
        int32_t factorRead;
        factorDataset.readData(arrayRefFromArray(&factorRead, 1));
        EXPECT_EQ(factorRead, factorWritten);
        const auto [factorHandle, factorGuard] =
                makeH5mdDataSetGuard(H5Dopen(fileid(), "factor", H5P_DEFAULT));
        EXPECT_THAT(getAttribute<std::string>(factorHandle, "description"),
                    ::testing::Optional(factorDescriptionWritten));

        // tolerance: root-level float dataset with unit and description as attributes
        auto  toleranceDataset = H5mdFixedDataSet<float>(fileid(), "tolerance");
        float toleranceRead;
        toleranceDataset.readData(arrayRefFromArray(&toleranceRead, 1));
        EXPECT_FLOAT_EQ(toleranceRead, toleranceWritten);
        const auto [toleranceHandle, toleranceGuard] =
                makeH5mdDataSetGuard(H5Dopen(fileid(), "tolerance", H5P_DEFAULT));
        EXPECT_THAT(getAttribute<std::string>(toleranceHandle, "unit"),
                    ::testing::Optional(toleranceUnitWritten));
        EXPECT_THAT(getAttribute<std::string>(toleranceHandle, "description"),
                    ::testing::Optional(toleranceDescriptionWritten));
    }
}

//! Add \p values as the uniform array \p key of \p object
template<typename ValueType>
void addUniformArrayFromVector(KeyValueTreeObjectBuilder*    object,
                               const std::string&            key,
                               const std::vector<ValueType>& values)
{
    auto array = object->addUniformArray<ValueType>(key);
    for (const ValueType& value : values)
    {
        array.addValue(value);
    }
}

//! Read the values of the numeric data set \p name in \p container
template<typename ValueType>
std::vector<ValueType> readNumericData(const hid_t container, const char* name)
{
    auto                   dataSet = H5mdFixedDataSet<ValueType>(container, name);
    std::vector<ValueType> values(dataSet.numValues());
    dataSet.readData(values);
    return values;
}

TEST_F(H5mdKVTreeWriterTest, WriteArrayProperties)
{
    const std::vector<std::string> groupsWritten        = { "Protein", "LIG", "SOL", "MEMBRANE" };
    const std::vector<int32_t>     indicesWritten       = { 1, 2, 3, 4, 5 };
    const std::vector<int64_t>     indices64Written     = { 1, 2, 3, 4, 5 };
    const std::vector<float>       scalingFactorWritten = { 0.1F, 0.2F, 0.3F };
    const std::vector<double>      configsWritten       = { 0.1, 0.2, 0.3 };
    // Description and doi of each array property, in the order of the arrays above
    const std::vector<std::pair<std::string, std::string>> attributesWritten = {
        { "A list of group names", "10.xxxx/example1.doi" },
        { "A list of indices", "10.xxxx/example2.doi" },
        { "A list of 64-bit indices", "10.xxxx/example3.doi" },
        { "A list of 32-bit float values", "10.xxxx/example4.doi" },
        { "A list of configuration values", "10.xxxx/example5.doi" }
    };
    const std::vector<std::string> names = {
        "groups", "indices", "indices_i64", "scaling_factor_f32", "configs"
    };

    {
        SCOPED_TRACE("Write array properties to H5MD");

        KeyValueTreeBuilder builder;
        auto                arrayProp   = builder.rootObject().addObject("array-property");
        const auto          addProperty = [&arrayProp, &names, &attributesWritten](const Index i)
        {
            auto property = arrayProp.addObject(names[i]);
            property.addValue<std::string>("description", attributesWritten[i].first);
            property.addValue<std::string>("doi", attributesWritten[i].second);
            return property;
        };
        auto groups = addProperty(0);
        addUniformArrayFromVector(&groups, "value", groupsWritten);
        auto indices = addProperty(1);
        addUniformArrayFromVector(&indices, "value", indicesWritten);
        auto indices64 = addProperty(2);
        addUniformArrayFromVector(&indices64, "value", indices64Written);
        auto scalingFactor = addProperty(3);
        addUniformArrayFromVector(&scalingFactor, "value", scalingFactorWritten);
        auto configs = addProperty(4);
        addUniformArrayFromVector(&configs, "value", configsWritten);
        writeKeyValueTreeToH5md(fileid(), builder.build());
    }

    {
        SCOPED_TRACE("Read back and verify array properties");

        const auto [arrayGroup, arrayGuard] =
                makeH5mdGroupGuard(openGroup(fileid(), "array-property"));

        // Each array is a data set directly in the group, with its attributes on the data set
        EXPECT_EQ(readFixedSizeStringData(arrayGroup, "groups"), groupsWritten);
        EXPECT_EQ(readNumericData<int32_t>(arrayGroup, "indices"), indicesWritten);
        EXPECT_EQ(readNumericData<int64_t>(arrayGroup, "indices_i64"), indices64Written);
        EXPECT_EQ(readNumericData<float>(arrayGroup, "scaling_factor_f32"), scalingFactorWritten);
        EXPECT_EQ(readNumericData<double>(arrayGroup, "configs"), configsWritten);

        for (Index i = 0; i < gmx::ssize(names); ++i)
        {
            SCOPED_TRACE("Attributes of " + names[i]);
            const auto [handle, guard] =
                    makeH5mdDataSetGuard(H5Dopen(arrayGroup, names[i].c_str(), H5P_DEFAULT));
            EXPECT_THAT(getAttribute<std::string>(handle, "description"),
                        ::testing::Optional(attributesWritten[i].first));
            EXPECT_THAT(getAttribute<std::string>(handle, "doi"),
                        ::testing::Optional(attributesWritten[i].second));
        }
    }
}


TEST_F(H5mdKVTreeWriterTest, WriteEmptyArray)
{
    const std::string descriptionWritten = "An empty array";

    KeyValueTreeBuilder builder;
    auto                emptyProp = builder.rootObject().addObject("empty-array-prop");
    emptyProp.addValue<std::string>("description", descriptionWritten);
    emptyProp.addUniformArray<std::string>("value"); // no elements added

    writeKeyValueTreeToH5md(fileid(), builder.build());

    // empty-array-prop: zero-length dataset at root level; description attribute on the dataset
    const auto [dataset, datasetGuard] =
            makeH5mdDataSetGuard(H5Dopen(fileid(), "empty-array-prop", H5P_DEFAULT));
    ASSERT_GE(dataset, 0);
    const auto [dataspace, dataspaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataset));
    hsize_t dims[1]                        = { 0 };
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);
    EXPECT_EQ(dims[0], 0U);

    EXPECT_THAT(getAttribute<std::string>(dataset, "description"), ::testing::Optional(descriptionWritten));
}

//! Return the dimensions of the data set \p name in \p container
std::vector<hsize_t> dataSetDims(const hid_t container, const char* name)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(H5Dopen(container, name, H5P_DEFAULT));
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    std::vector<hsize_t> dims(H5Sget_simple_extent_ndims(dataSpace));
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    return dims;
}

TEST_F(H5mdKVTreeWriterTest, WriteArrayWithShape)
{
    const std::vector<double> compressibilityWritten = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    const std::vector<int64_t> shapeWritten = { 3, 3 };
    const std::string          unitWritten  = "bar-1";

    KeyValueTreeBuilder builder;
    auto                compressibility = builder.rootObject().addObject("compressibility");
    addUniformArrayFromVector(&compressibility, "value", compressibilityWritten);
    addUniformArrayFromVector(&compressibility, "shape", shapeWritten);
    compressibility.addValue<std::string>("unit", unitWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    EXPECT_EQ(dataSetDims(fileid(), "compressibility"),
              std::vector<hsize_t>(shapeWritten.begin(), shapeWritten.end()));
    // The values are written in row-major order
    EXPECT_EQ(readNumericData<double>(fileid(), "compressibility"), compressibilityWritten);

    // The shape is consumed as the dimensions, not written as an attribute
    const auto [dsHandle, dsGuard] =
            makeH5mdDataSetGuard(H5Dopen(fileid(), "compressibility", H5P_DEFAULT));
    EXPECT_EQ(H5Aexists(dsHandle, "shape"), 0);
    EXPECT_THAT(getAttribute<std::string>(dsHandle, "unit"), ::testing::Optional(unitWritten));
}

TEST_F(H5mdKVTreeWriterTest, WriteArrayWithoutShapeInOneDimension)
{
    const std::vector<double> compressibilityWritten = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

    KeyValueTreeBuilder builder;
    auto                compressibility = builder.rootObject().addObject("compressibility");
    addUniformArrayFromVector(&compressibility, "value", compressibilityWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    EXPECT_EQ(dataSetDims(fileid(), "compressibility"),
              std::vector<hsize_t>{ compressibilityWritten.size() });
    EXPECT_EQ(readNumericData<double>(fileid(), "compressibility"), compressibilityWritten);
}

TEST_F(H5mdKVTreeWriterTest, WriteStringArrayWithShape)
{
    const std::vector<std::string> namesWritten = { "a", "b", "c", "d", "e", "f" };
    const std::vector<int>         shapeWritten = { 2, 3 };

    KeyValueTreeBuilder builder;
    auto                names = builder.rootObject().addObject("names");
    addUniformArrayFromVector(&names, "value", namesWritten);
    addUniformArrayFromVector(&names, "shape", shapeWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    EXPECT_EQ(dataSetDims(fileid(), "names"),
              std::vector<hsize_t>(shapeWritten.begin(), shapeWritten.end()));
    EXPECT_EQ(readFixedSizeStringData(fileid(), "names"), namesWritten);
}

TEST_F(H5mdKVTreeWriterTest, WriteBoolAsInt)
{
    const std::string trueDescriptionWritten  = "Boolean true stored as int";
    const std::string falseDescriptionWritten = "Boolean false stored as int";

    KeyValueTreeBuilder builder;
    auto                root = builder.rootObject();

    auto trueProp = root.addObject("flag-true");
    trueProp.addValue<bool>("value", true);
    trueProp.addValue<std::string>("description", trueDescriptionWritten);

    auto falseProp = root.addObject("flag-false");
    falseProp.addValue<bool>("value", false);
    falseProp.addValue<std::string>("description", falseDescriptionWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    // Each flag is an int data set at root level holding 1 or 0, with its description attribute
    EXPECT_EQ(readNumericData<int>(fileid(), "flag-true"), std::vector<int>{ 1 });
    const auto [trueHandle, trueGuard] =
            makeH5mdDataSetGuard(H5Dopen(fileid(), "flag-true", H5P_DEFAULT));
    EXPECT_THAT(getAttribute<std::string>(trueHandle, "description"),
                ::testing::Optional(trueDescriptionWritten));

    EXPECT_EQ(readNumericData<int>(fileid(), "flag-false"), std::vector<int>{ 0 });
    const auto [falseHandle, falseGuard] =
            makeH5mdDataSetGuard(H5Dopen(fileid(), "flag-false", H5P_DEFAULT));
    EXPECT_THAT(getAttribute<std::string>(falseHandle, "description"),
                ::testing::Optional(falseDescriptionWritten));
}

TEST_F(H5mdKVTreeWriterTest, WriteDeepNestedObjects)
{
    const double      timestepWritten    = 0.001;
    const std::string unitWritten        = "ps";
    const std::string descriptionWritten = "A deeply nested timestep";

    KeyValueTreeBuilder builder;
    // 4-level nesting: simulation > algorithm > parameters > timestep
    auto timestep = builder.rootObject()
                            .addObject("simulation")
                            .addObject("algorithm")
                            .addObject("parameters")
                            .addObject("timestep");
    timestep.addValue<double>("value", timestepWritten);
    timestep.addValue<std::string>("unit", unitWritten);
    timestep.addValue<std::string>("description", descriptionWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    const auto [simGroup, simGuard]       = makeH5mdGroupGuard(openGroup(fileid(), "simulation"));
    const auto [algoGroup, algoGuard]     = makeH5mdGroupGuard(openGroup(simGroup, "algorithm"));
    const auto [paramsGroup, paramsGuard] = makeH5mdGroupGuard(openGroup(algoGroup, "parameters"));
    // timestep: double dataset directly in paramsGroup; unit and description as attributes
    const std::vector<double> timestepRead = readNumericData<double>(paramsGroup, "timestep");
    ASSERT_EQ(timestepRead.size(), 1U);
    EXPECT_DOUBLE_EQ(timestepRead[0], timestepWritten);

    const auto [tsHandle, tsGuard] = makeH5mdDataSetGuard(H5Dopen(paramsGroup, "timestep", H5P_DEFAULT));
    EXPECT_THAT(getAttribute<std::string>(tsHandle, "unit"), ::testing::Optional(unitWritten));
    EXPECT_THAT(getAttribute<std::string>(tsHandle, "description"),
                ::testing::Optional(descriptionWritten));
}

TEST_F(H5mdKVTreeWriterTest, WriteEmptyRoot)
{
    // The file already holds the H5MD metadata group, so count what the writer adds to it
    H5G_info_t rootBefore;
    ASSERT_GE(H5Gget_info(fileid(), &rootBefore), 0);

    // Generate an empty key-value tree
    KeyValueTreeBuilder builder;
    writeKeyValueTreeToH5md(fileid(), builder.build());

    // Expect no datasets or groups created at root level
    H5G_info_t rootAfter;
    ASSERT_GE(H5Gget_info(fileid(), &rootAfter), 0);
    EXPECT_EQ(rootAfter.nlinks, rootBefore.nlinks);
}

TEST_F(H5mdKVTreeWriterTest, WriteUrlAndVersionAttributes)
{
    const std::string urlWritten     = "https://example.org";
    const std::string versionWritten = "1.0";

    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addValue<double>("value", 1.0);
    property.addValue<std::string>("url", urlWritten);
    property.addValue<std::string>("version", versionWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(H5Dopen(fileid(), "property", H5P_DEFAULT));
    EXPECT_THAT(getAttribute<std::string>(dataSet, "url"), ::testing::Optional(urlWritten));
    EXPECT_THAT(getAttribute<std::string>(dataSet, "version"), ::testing::Optional(versionWritten));
}

TEST_F(H5mdKVTreeWriterTest, WriteAttributesOfEveryType)
{
    const bool                     boolWritten    = true;
    const int                      intWritten     = 3;
    const int64_t                  int64Written   = 4;
    const float                    floatWritten   = 0.5F;
    const double                   doubleWritten  = 0.25;
    const std::string              stringWritten  = "text";
    const std::vector<std::string> stringsWritten = { "a", "b" };
    const std::vector<int>         intsWritten    = { 1, 2 };
    const std::vector<int64_t>     int64sWritten  = { 3, 4 };
    const std::vector<float>       floatsWritten  = { 0.5F, 1.5F };
    const std::vector<double>      doublesWritten = { 0.25, 0.75 };

    // Each group holds its attribute under the same reserved key, with a value of another type
    KeyValueTreeBuilder builder;
    auto                root = builder.rootObject();
    root.addObject("bool").addValue<bool>("version", boolWritten);
    root.addObject("int").addValue<int>("version", intWritten);
    root.addObject("int64").addValue<int64_t>("version", int64Written);
    root.addObject("float").addValue<float>("version", floatWritten);
    root.addObject("double").addValue<double>("version", doubleWritten);
    root.addObject("string").addValue<std::string>("version", stringWritten);
    root.addObject("emptyArray").addUniformArray<std::string>("version");
    auto strings = root.addObject("strings");
    addUniformArrayFromVector(&strings, "version", stringsWritten);
    auto ints = root.addObject("ints");
    addUniformArrayFromVector(&ints, "version", intsWritten);
    auto int64s = root.addObject("int64s");
    addUniformArrayFromVector(&int64s, "version", int64sWritten);
    auto floats = root.addObject("floats");
    addUniformArrayFromVector(&floats, "version", floatsWritten);
    auto doubles = root.addObject("doubles");
    addUniformArrayFromVector(&doubles, "version", doublesWritten);

    writeKeyValueTreeToH5md(fileid(), builder.build());

    const auto version = [this](const char* groupName, auto read)
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), groupName));
        return read(group);
    };
    // A bool is stored as an int
    EXPECT_THAT(version("bool", [](hid_t g) { return getAttribute<int>(g, "version"); }),
                ::testing::Optional(static_cast<int>(boolWritten)));
    EXPECT_THAT(version("int", [](hid_t g) { return getAttribute<int>(g, "version"); }),
                ::testing::Optional(intWritten));
    EXPECT_THAT(version("int64", [](hid_t g) { return getAttribute<int64_t>(g, "version"); }),
                ::testing::Optional(int64Written));
    EXPECT_THAT(version("float", [](hid_t g) { return getAttribute<float>(g, "version"); }),
                ::testing::Optional(floatWritten));
    EXPECT_THAT(version("double", [](hid_t g) { return getAttribute<double>(g, "version"); }),
                ::testing::Optional(doubleWritten));
    EXPECT_THAT(version("string", [](hid_t g) { return getAttribute<std::string>(g, "version"); }),
                ::testing::Optional(stringWritten));
    EXPECT_EQ(version("emptyArray", [](hid_t g) { return H5Aexists(g, "version"); }), 0);
    EXPECT_THAT(
            version("strings", [](hid_t g) { return getAttributeVector<std::string>(g, "version"); }),
            ::testing::Optional(stringsWritten));
    EXPECT_THAT(version("ints", [](hid_t g) { return getAttributeVector<int>(g, "version"); }),
                ::testing::Optional(intsWritten));
    EXPECT_THAT(version("int64s", [](hid_t g) { return getAttributeVector<int64_t>(g, "version"); }),
                ::testing::Optional(int64sWritten));
    EXPECT_THAT(version("floats", [](hid_t g) { return getAttributeVector<float>(g, "version"); }),
                ::testing::Optional(floatsWritten));
    EXPECT_THAT(version("doubles", [](hid_t g) { return getAttributeVector<double>(g, "version"); }),
                ::testing::Optional(doublesWritten));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnNestedObjectInValueProperty)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addValue<double>("value", 1.0);
    property.addObject("nested").addValue<double>("value", 2.0);
    const KeyValueTreeObject tree = builder.build();

    // HDF5 would also fail to create the nested object, so check what is reported
    EXPECT_THAT(
            [&]() { writeKeyValueTreeToH5md(fileid(), tree); },
            ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("cannot hold the nested")));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnKeyThatIsNotAllowed)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addValue<double>("value", 1.0);
    property.addValue<std::string>("comment", "Not an allowed attribute key");

    EXPECT_THROW_GMX(writeKeyValueTreeToH5md(fileid(), builder.build()), H5mdError);
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeThatDoesNotMatchValues)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("compressibility");
    property.addUniformArray<double>("value", { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 });
    property.addUniformArray<int64_t>("shape", { 3, 3 });
    const KeyValueTreeObject tree = builder.build();

    // HDF5 would also fail on the mismatch, so check that the shape is what is reported
    EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("but the array holds 8")));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeThatIsNotPositive)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addUniformArray<double>("value", { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 });
    // As unsigned dimensions, -3 x -3 wraps around to 9 values
    property.addUniformArray<int64_t>("shape", { -3, -3 });
    const KeyValueTreeObject tree = builder.build();

    EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("positive integers")));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeOfScalar)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addValue<double>("value", 0.1);
    property.addUniformArray<int64_t>("shape", { 1 });

    EXPECT_THROW_GMX(writeKeyValueTreeToH5md(fileid(), builder.build()), H5mdError);
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeWithoutValue)
{
    KeyValueTreeBuilder builder;
    builder.rootObject().addObject("group").addUniformArray<int64_t>("shape", { 3, 3 });

    EXPECT_THROW_GMX(writeKeyValueTreeToH5md(fileid(), builder.build()), H5mdError);
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnValueOfUnsupportedType)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addObject("value").addValue<double>("value", 1.0);

    EXPECT_THROW_GMX(writeKeyValueTreeToH5md(fileid(), builder.build()), H5mdError);
}

TEST_F(H5mdKVTreeWriterTest, WritesNothingForScalarOrArrayAtTopLevel)
{
    H5G_info_t rootBefore;
    ASSERT_GE(H5Gget_info(fileid(), &rootBefore), 0);

    KeyValueTreeBuilder builder;
    builder.rootObject().addValue<int>("scalar", 1);
    builder.rootObject().addUniformArray<int>("array", { 1, 2 });
    writeKeyValueTreeToH5md(fileid(), builder.build());

    H5G_info_t rootAfter;
    ASSERT_GE(H5Gget_info(fileid(), &rootAfter), 0);
    EXPECT_EQ(rootAfter.nlinks, rootBefore.nlinks);
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeThatIsNotAnArray)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addUniformArray<double>("value", { 0.1, 0.2, 0.3 });
    property.addValue<int64_t>("shape", 3);
    const KeyValueTreeObject tree = builder.build();

    EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("must be an array")));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnShapeWithNonIntegerElements)
{
    KeyValueTreeBuilder builder;
    auto                property = builder.rootObject().addObject("property");
    property.addUniformArray<double>("value", { 0.1, 0.2, 0.3, 0.4 });
    property.addUniformArray<double>("shape", { 2.0, 2.0 });
    const KeyValueTreeObject tree = builder.build();

    EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("positive integers")));
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnArrayOfUnsupportedElementType)
{
    KeyValueTreeBuilder builder;
    builder.rootObject().addObject("property").addUniformArray<bool>("value", { true, false });

    EXPECT_THROW_GMX(writeKeyValueTreeToH5md(fileid(), builder.build()), H5mdError);
}

TEST_F(H5mdKVTreeWriterTest, ThrowsOnArrayOfMixedTypes)
{
    {
        SCOPED_TRACE("A value array");
        KeyValueTreeBuilder builder;
        auto                values = builder.rootObject().addObject("property").addArray("value");
        values.addRawValue(Any::create<int>(1));
        values.addRawValue(Any::create<double>(2.5));
        const KeyValueTreeObject tree = builder.build();

        EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                    ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("of one type")));
    }
    {
        SCOPED_TRACE("An attribute array");
        KeyValueTreeBuilder builder;
        auto                versions = builder.rootObject().addObject("group").addArray("version");
        versions.addRawValue(Any::create<int>(1));
        versions.addRawValue(Any::create<std::string>("2"));
        const KeyValueTreeObject tree = builder.build();

        EXPECT_THAT([&]() { writeKeyValueTreeToH5md(fileid(), tree); },
                    ::testing::ThrowsMessage<H5mdError>(::testing::HasSubstr("of one type")));
    }
}

} // namespace
} // namespace test
} // namespace gmx
