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
 * Tests for density fitting module options.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/densityfitting/densityfittingoptions.h"

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/densityfitting/densityfittingamplitudelookup.h"
#include "gromacs/applied_forces/densityfitting/densityfittingparameters.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

class DensityFittingOptionsTest : public ::testing::Test
{
public:
    DensityFittingOptionsTest() { init_blocka(&defaultGroups_); }
    ~DensityFittingOptionsTest() override { done_blocka(&defaultGroups_); }

    void setFromMdpValues(const KeyValueTreeObject& densityFittingMdpValues)
    {
        // set up options
        Options densityFittingModuleOptions;
        densityFittingOptions_.initMdpOptions(&densityFittingModuleOptions);

        // Add rules to transform mdp inputs to densityFittingModule data
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);

        densityFittingOptions_.initMdpTransform(transform.rules());

        // Execute the transform on the mdpValues
        auto transformedMdpValues = transform.transform(densityFittingMdpValues, nullptr);
        assignOptionsFromKeyValueTree(&densityFittingModuleOptions, transformedMdpValues.object(), nullptr);
    }

    static KeyValueTreeObject densityFittingSetActiveAsMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue("density-guided-simulation-active", std::string("yes"));
        return mdpValueBuilder.build();
    }

    static IndexGroupsAndNames genericIndexGroupsAndNames()
    {
        std::vector<IndexGroup> indexGroups;
        indexGroups.push_back({ "A", { 0 } });
        indexGroups.push_back({ "protein", { 1 } });
        indexGroups.push_back({ "C", { 2 } });
        return IndexGroupsAndNames(indexGroups);
    }

    static IndexGroupsAndNames differingIndexGroupsAndNames()
    {
        std::vector<IndexGroup> indexGroups;
        indexGroups.push_back({ "protein", { 0 } });
        indexGroups.push_back({ "C", { 1 } });
        indexGroups.push_back({ "A", { 2 } });
        return IndexGroupsAndNames(indexGroups);
    }

    void mangleInternalParameters()
    {
        densityFittingOptions_.setFitGroupIndices(differingIndexGroupsAndNames());
    }

protected:
    t_blocka              defaultGroups_;
    DensityFittingOptions densityFittingOptions_;
};

TEST_F(DensityFittingOptionsTest, DefaultParameters)
{
    const auto defaultParameters = densityFittingOptions_.buildParameters();
    EXPECT_FALSE(defaultParameters.active_);
    EXPECT_EQ(0, defaultParameters.indices_.size());
    EXPECT_EQ(DensitySimilarityMeasureMethod::innerProduct, defaultParameters.similarityMeasureMethod_);
    EXPECT_EQ(DensityFittingAmplitudeMethod::Unity, defaultParameters.amplitudeLookupMethod_);
    EXPECT_REAL_EQ(1e9, defaultParameters.forceConstant_);
    EXPECT_REAL_EQ(0.2, defaultParameters.gaussianTransformSpreadingWidth_);
    EXPECT_REAL_EQ(4.0, defaultParameters.gaussianTransformSpreadingRangeInMultiplesOfWidth_);
}

TEST_F(DensityFittingOptionsTest, OptionSetsActive)
{
    EXPECT_FALSE(densityFittingOptions_.buildParameters().active_);
    setFromMdpValues(densityFittingSetActiveAsMdpValues());
    EXPECT_TRUE(densityFittingOptions_.buildParameters().active_);
}

TEST_F(DensityFittingOptionsTest, OutputNoDefaultValuesWhenInactive)
{
    // Transform module data into a flat key-value tree for output.

    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    densityFittingOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    EXPECT_EQ(stream.toString(),
              std::string(
                      "\n; Density guided simulation\ndensity-guided-simulation-active = false\n"));
}

TEST_F(DensityFittingOptionsTest, OutputDefaultValuesWhenActive)
{
    setFromMdpValues(densityFittingSetActiveAsMdpValues());
    // Transform module data into a flat key-value tree for output.

    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    densityFittingOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();
    std::string expectedString = {
        "\n"
        "; Density guided simulation\n"
        "density-guided-simulation-active = true\n"
        "density-guided-simulation-group = protein\n"
        "; Similarity measure between densities: inner-product, relative-entropy, or "
        "cross-correlation\n"
        "density-guided-simulation-similarity-measure = inner-product\n"
        "; Atom amplitude for spreading onto grid: unity, mass, or charge\n"
        "density-guided-simulation-atom-spreading-weight = unity\n"
        "density-guided-simulation-force-constant = 1e+09\n"
        "density-guided-simulation-gaussian-transform-spreading-width = 0.2\n"
        "density-guided-simulation-gaussian-transform-spreading-range-in-multiples-of-width = 4\n"
        "; Reference density file location as absolute path or relative to the gmx mdrun calling "
        "location\n"
        "density-guided-simulation-reference-density-filename = reference.mrc\n"
        "density-guided-simulation-nst = 1\n"
        "; Normalize the sum of density voxel values to one\n"
        "density-guided-simulation-normalize-densities = true\n"
        "; Apply adaptive force scaling\n"
        "density-guided-simulation-adaptive-force-scaling = false\n"
        "; Time constant for adaptive force scaling in ps\n"
        "density-guided-simulation-adaptive-force-scaling-time-constant = 4\n"
    };

    EXPECT_EQ(expectedString, stream.toString());
}


TEST_F(DensityFittingOptionsTest, CanConvertGroupStringToIndexGroup)
{
    setFromMdpValues(densityFittingSetActiveAsMdpValues());

    const auto indexGroupAndNames = genericIndexGroupsAndNames();
    densityFittingOptions_.setFitGroupIndices(indexGroupAndNames);

    EXPECT_EQ(1, densityFittingOptions_.buildParameters().indices_.size());
    EXPECT_EQ(1, densityFittingOptions_.buildParameters().indices_[0]);
}

TEST_F(DensityFittingOptionsTest, InternalsToKvt)
{
    // stores the default internal options
    DensityFittingOptions densityFittingOptions;
    KeyValueTreeBuilder   builder;
    densityFittingOptions.writeInternalParametersToKvt(builder.rootObject());
    const auto kvtTree = builder.build();
    EXPECT_TRUE(kvtTree.keyExists("density-guided-simulation-group"));
    EXPECT_TRUE(kvtTree["density-guided-simulation-group"].isArray());
    auto storedIndex = kvtTree["density-guided-simulation-group"].asArray().values();

    EXPECT_EQ(0, storedIndex.size());
}

TEST_F(DensityFittingOptionsTest, KvtToInternal)
{
    setFromMdpValues(densityFittingSetActiveAsMdpValues());

    KeyValueTreeBuilder builder;
    auto                addedArray =
            builder.rootObject().addUniformArray<std::int64_t>("density-guided-simulation-group");
    addedArray.addValue(1);
    addedArray.addValue(15);
    const auto tree = builder.build();

    densityFittingOptions_.readInternalParametersFromKvt(tree);

    EXPECT_EQ(2, densityFittingOptions_.buildParameters().indices_.size());
    EXPECT_EQ(1, densityFittingOptions_.buildParameters().indices_[0]);
    EXPECT_EQ(15, densityFittingOptions_.buildParameters().indices_[1]);
}

TEST_F(DensityFittingOptionsTest, RoundTripForInternalsIsIdempotent)
{
    setFromMdpValues(densityFittingSetActiveAsMdpValues());
    {
        const IndexGroupsAndNames indexGroupAndNames = genericIndexGroupsAndNames();
        densityFittingOptions_.setFitGroupIndices(indexGroupAndNames);
    }

    DensityFittingParameters parametersBefore = densityFittingOptions_.buildParameters();

    KeyValueTreeBuilder builder;
    densityFittingOptions_.writeInternalParametersToKvt(builder.rootObject());
    const auto inputTree = builder.build();

    mangleInternalParameters();

    DensityFittingParameters parametersAfter = densityFittingOptions_.buildParameters();
    EXPECT_NE(parametersBefore, parametersAfter);

    densityFittingOptions_.readInternalParametersFromKvt(inputTree);

    parametersAfter = densityFittingOptions_.buildParameters();
    EXPECT_EQ(parametersBefore, parametersAfter);
}

} // namespace

} // namespace gmx
