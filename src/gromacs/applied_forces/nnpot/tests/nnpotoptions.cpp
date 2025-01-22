/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for functionality of the NNPotOptions
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpotoptions.h"

#include <gtest/gtest.h>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

class NNPotOptionsTest : public ::testing::Test
{
public:
    void setFromMdpValues(const KeyValueTreeObject& qmmmMdpValues)
    {
        // Setup options
        Options qmmmModuleOptions;
        nnpotOptions_.initMdpOptions(&qmmmModuleOptions);

        // Add rules to transform mdp inputs to densityFittingModule data
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);

        nnpotOptions_.initMdpTransform(transform.rules());

        // Execute the transform on the mdpValues
        auto transformedMdpValues = transform.transform(qmmmMdpValues, nullptr);
        assignOptionsFromKeyValueTree(&qmmmModuleOptions, transformedMdpValues.object(), nullptr);
    }

    static KeyValueTreeObject nnpotBuildDefaultMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-active", std::string("true"));
        return mdpValueBuilder.build();
    }

    static KeyValueTreeObject nnpotBuildInputMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-active", std::string("true"));
        mdpValueBuilder.rootObject().addValue(
                c_nnpotModuleName + "-modelfile",
                gmx::test::TestFileManager::getInputFilePath("model.pt").string());
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-model-input1",
                                              std::string("atom-positions"));
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-model-input2",
                                              std::string("atom-numbers"));
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-model-input3", std::string("box"));
        mdpValueBuilder.rootObject().addValue(c_nnpotModuleName + "-model-input4", std::string("pbc"));
        return mdpValueBuilder.build();
    }

    static IndexGroupsAndNames indexGroupsAndNamesGeneric()
    {
        // System group is default for QM atoms
        std::vector<IndexGroup> indexGroups;
        indexGroups.push_back({ "A", { 1 } });
        indexGroups.push_back({ "System", { 1, 2, 3 } });
        indexGroups.push_back({ "C", { 2, 3 } });

        return IndexGroupsAndNames(indexGroups);
    }

protected:
    NNPotOptions nnpotOptions_;
};

TEST_F(NNPotOptionsTest, DefaultParameters)
{
    const NNPotParameters&          defaultParams = nnpotOptions_.parameters();
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkBoolean(defaultParams.active_, "active");
    checker.checkString(defaultParams.modelFileName_, "modelFileName");
    checker.checkString(defaultParams.inputGroup_, "inputGroup");
    checker.checkString(defaultParams.modelInput_[0], "modelInput1");
    checker.checkString(defaultParams.modelInput_[1], "modelInput2");
    checker.checkString(defaultParams.modelInput_[2], "modelInput3");
    checker.checkString(defaultParams.modelInput_[3], "modelInput4");
}

TEST_F(NNPotOptionsTest, OptionSetsActive)
{
    EXPECT_FALSE(nnpotOptions_.parameters().active_);
    setFromMdpValues(nnpotBuildDefaultMdpValues());
    EXPECT_TRUE(nnpotOptions_.parameters().active_);
}

TEST_F(NNPotOptionsTest, OutputNoDefaultValuesWhenInactive)
{
    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    nnpotOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(NNPotOptionsTest, OutputDefaultValuesWhenActive)
{

    // Set nnpot-active = true
    setFromMdpValues(nnpotBuildDefaultMdpValues());

    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    nnpotOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(NNPotOptionsTest, InternalsToKvtAndBack)
{
    // Set nnpot-active = true
    setFromMdpValues(nnpotBuildInputMdpValues());

    // Set indices
    const IndexGroupsAndNames indexGroupAndNames = indexGroupsAndNamesGeneric();
    nnpotOptions_.setInputGroupIndices(indexGroupAndNames);

    // Set dummy logger and warning handler
    MDLogger logger;
    nnpotOptions_.setLogger(logger);
    WarningHandler warninp(true, 0);
    nnpotOptions_.setWarninp(&warninp);

    // Copy internal parameters
    const NNPotParameters& params           = nnpotOptions_.parameters();
    auto                   inpIndicesBefore = params.inpIndices_;
    auto                   mmIndicesBefore  = params.mmIndices_;

    KeyValueTreeBuilder builder;
    if (GMX_TORCH)
    {
        EXPECT_NO_THROW(nnpotOptions_.writeParamsToKvt(builder.rootObject()));
        const auto inputTree = builder.build();

        EXPECT_NO_THROW(nnpotOptions_.readParamsFromKvt(inputTree));

        // Check Internal parameters taken back from KVT
        const NNPotParameters& params2 = nnpotOptions_.parameters();
        EXPECT_EQ(inpIndicesBefore, params2.inpIndices_);
        EXPECT_EQ(mmIndicesBefore, params2.mmIndices_);
    }
    else
    {
        EXPECT_ANY_THROW(nnpotOptions_.writeParamsToKvt(builder.rootObject()));
    }
}

} // namespace test

} // namespace gmx
