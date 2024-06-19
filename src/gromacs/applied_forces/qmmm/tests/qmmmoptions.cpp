/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for QMMMOptions class of QMMM MDModule.
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/qmmm/qmmmoptions.h"

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/qmmm/qmmmtypes.h"
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
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

class QMMMOptionsTest : public ::testing::Test
{
public:
    void setFromMdpValues(const KeyValueTreeObject& qmmmMdpValues)
    {
        // Setup options
        Options qmmmModuleOptions;
        qmmmOptions_.initMdpOptions(&qmmmModuleOptions);

        // Add rules to transform mdp inputs to densityFittingModule data
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);

        qmmmOptions_.initMdpTransform(transform.rules());

        // Execute the transform on the mdpValues
        auto transformedMdpValues = transform.transform(qmmmMdpValues, nullptr);
        assignOptionsFromKeyValueTree(&qmmmModuleOptions, transformedMdpValues.object(), nullptr);
    }

    static KeyValueTreeObject qmmmBuildDefaulMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(c_qmmmCP2KModuleName + "-active", std::string("true"));
        return mdpValueBuilder.build();
    }

    static KeyValueTreeObject qmmmBuildMethodInputMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(c_qmmmCP2KModuleName + "-active", std::string("true"));
        mdpValueBuilder.rootObject().addValue(c_qmmmCP2KModuleName + "-qmmethod", std::string("INPUT"));
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


    static IndexGroupsAndNames indexGroupsAndNamesNoQM()
    {
        // System group is default for QM atoms (not present here)
        std::vector<IndexGroup> indexGroups;
        indexGroups.push_back({ "A", { 1 } });
        indexGroups.push_back({ "B", { 2 } });
        indexGroups.push_back({ "C", { 3 } });

        return IndexGroupsAndNames(indexGroups);
    }

    static IndexGroupsAndNames indexGroupsAndNamesEmptyQM()
    {
        // System group is default for QM atoms and empty in this case
        std::vector<IndexGroup> indexGroups;
        indexGroups.push_back({ "A", { 1 } });
        indexGroups.push_back({ "System", {} });
        indexGroups.push_back({ "C", { 2, 3 } });

        return IndexGroupsAndNames(indexGroups);
    }

protected:
    QMMMOptions qmmmOptions_;
};

TEST_F(QMMMOptionsTest, DefaultParameters)
{
    const QMMMParameters&           defaultParameters = qmmmOptions_.parameters();
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // Tolerance of all coordinates and vectors should be 1E-3 (as in gro or pdb files)
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    checker.checkBoolean(defaultParameters.active_, "Active");
    checker.checkInteger(defaultParameters.atomNumbers_.size(), "Size of atomNumbers");
    checker.checkInteger(defaultParameters.link_.size(), "Size of link");
    checker.checkInteger(defaultParameters.mmIndices_.size(), "Size of mmIndices");
    checker.checkInteger(defaultParameters.qmIndices_.size(), "Size of qmIndices");
    checker.checkInteger(defaultParameters.qmCharge_, "QM charge");
    checker.checkInteger(defaultParameters.qmMultiplicity_, "QM multiplicity");
    checker.checkInteger(static_cast<int>(defaultParameters.qmMethod_), "QM method");
    checker.checkVector(defaultParameters.qmTrans_, "QM Translation");
    checker.checkVector(defaultParameters.qmBox_[0], "QM Box Vector 1");
    checker.checkVector(defaultParameters.qmBox_[1], "QM Box Vector 2");
    checker.checkVector(defaultParameters.qmBox_[2], "QM Box Vector 3");
    checker.checkString(defaultParameters.qmFileNameBase_, "QM filename base");
    checker.checkString(defaultParameters.qmInput_, "Input");
    checker.checkString(defaultParameters.qmPdb_, "PDB");
}

TEST_F(QMMMOptionsTest, OptionSetsActive)
{
    EXPECT_FALSE(qmmmOptions_.parameters().active_);
    setFromMdpValues(qmmmBuildDefaulMdpValues());
    EXPECT_TRUE(qmmmOptions_.parameters().active_);
}

TEST_F(QMMMOptionsTest, OutputNoDefaultValuesWhenInactive)
{
    // Transform module data into a flat key-value tree for output.

    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    qmmmOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(QMMMOptionsTest, OutputDefaultValuesWhenActive)
{

    // Set qmmm-active = true
    setFromMdpValues(qmmmBuildDefaulMdpValues());

    // Transform module data into a flat key-value tree for output.

    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    qmmmOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(QMMMOptionsTest, CanConvertGroupStringToIndexGroup)
{
    // Set qmmm-active = true
    setFromMdpValues(qmmmBuildDefaulMdpValues());

    // Generic index data
    const auto indexGroupAndNames = indexGroupsAndNamesGeneric();
    qmmmOptions_.setQMMMGroupIndices(indexGroupAndNames);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkInteger(qmmmOptions_.parameters().qmIndices_.size(), "Size of qmIndices");
    checker.checkInteger(qmmmOptions_.parameters().mmIndices_.size(), "Size of mmIndices");
}

TEST_F(QMMMOptionsTest, NoQMGroupConvertGroupStringToIndexGroup)
{
    // Set qmmm-active = true
    setFromMdpValues(qmmmBuildDefaulMdpValues());

    const auto indexGroupAndNames = indexGroupsAndNamesNoQM();
    EXPECT_ANY_THROW(qmmmOptions_.setQMMMGroupIndices(indexGroupAndNames));
}

TEST_F(QMMMOptionsTest, EmptyQMGroupConvertGroupStringToIndexGroup)
{
    // Set qmmm-active = true
    setFromMdpValues(qmmmBuildDefaulMdpValues());

    const auto indexGroupAndNames = indexGroupsAndNamesEmptyQM();
    EXPECT_ANY_THROW(qmmmOptions_.setQMMMGroupIndices(indexGroupAndNames));
}

TEST_F(QMMMOptionsTest, InternalsToKvtAndBack)
{
    // Set qmmm-active = true
    setFromMdpValues(qmmmBuildDefaulMdpValues());

    // Set indices
    const IndexGroupsAndNames indexGroupAndNames = indexGroupsAndNamesGeneric();
    qmmmOptions_.setQMMMGroupIndices(indexGroupAndNames);

    // Copy internal parameters
    const QMMMParameters& params            = qmmmOptions_.parameters();
    auto                  qmIndicesBefore   = params.qmIndices_;
    auto                  mmIndicesBefore   = params.mmIndices_;
    auto                  atomNumbersBefore = params.atomNumbers_;
    auto                  linkBefore        = params.link_;
    auto                  qmInputBefore     = params.qmInput_;
    auto                  qmPdbBefore       = params.qmPdb_;


    KeyValueTreeBuilder builder;
    qmmmOptions_.writeInternalParametersToKvt(builder.rootObject());
    const auto inputTree = builder.build();

    qmmmOptions_.readInternalParametersFromKvt(inputTree);

    // Check Internal parameters taken back from KVT
    const QMMMParameters& params2 = qmmmOptions_.parameters();
    EXPECT_EQ(qmIndicesBefore, params2.qmIndices_);
    EXPECT_EQ(mmIndicesBefore, params2.mmIndices_);
    EXPECT_EQ(atomNumbersBefore, params2.atomNumbers_);
    EXPECT_EQ(linkBefore.size(), params2.link_.size());
    EXPECT_EQ(qmInputBefore, params2.qmInput_);
    EXPECT_EQ(qmPdbBefore, params2.qmPdb_);
}

TEST_F(QMMMOptionsTest, CP2KInputProcessing)
{
    // Set qmmm-active = true and qmmm-qmmethod = INPUT
    setFromMdpValues(qmmmBuildMethodInputMdpValues());

    // Path to the sample CP2K input file
    std::string cp2kInput =
            gmx::test::TestFileManager::getInputFilePath("sample_cp2k_input.inp").string();

    // Process input file
    qmmmOptions_.setQMExternalInputFile({ true, cp2kInput });

    const QMMMParameters& params = qmmmOptions_.parameters();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(params.qmInput_, "CP2K input");
    checker.checkInteger(params.qmCharge_, "QM charge");
    checker.checkInteger(params.qmMultiplicity_, "QM multiplicity");
}

} // namespace gmx
