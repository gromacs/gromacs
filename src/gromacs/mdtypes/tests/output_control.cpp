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
 * Tests for OutputControl module's MDP output functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/output_control.h"

#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/inputrecstrings.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

class OutputControlModuleBuildMdpOutputTest : public ::testing::Test
{
public:
    OutputControlModuleBuildMdpOutputTest() : module_(OutputControlModuleInfo::create())
    {
        // Set up the OutputControl target
        outputControl_ = std::make_unique<OutputControl>();
        setOutputControlTarget(module_.get(), outputControl_.get());
    }

    //! Convert buildMdpOutput() result to text string
    std::string buildMdpOutputAsString()
    {
        KeyValueTreeBuilder builder;
        auto                builderObject = builder.rootObject();

        auto* mdpProvider = module_->mdpOptionProvider();
        mdpProvider->buildMdpOutput(&builderObject);

        StringOutputStream stream;
        {
            TextWriter writer(&stream);
            writeKeyValueTreeAsMdp(&writer, builder.build());
        }
        stream.close();

        return stream.toString();
    }

    std::unique_ptr<IMDModule>     module_;
    std::unique_ptr<OutputControl> outputControl_;
};

TEST_F(OutputControlModuleBuildMdpOutputTest, BuildsNoOutputByDefault)
{
    // By default, buildMdpOutput() should not write anything
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();

    auto* mdpProvider = module_->mdpOptionProvider();
    ASSERT_NE(mdpProvider, nullptr);

    // Build output without enabling writeToMdpOutput
    mdpProvider->buildMdpOutput(&builderObject);

    auto tree = builder.build();

    // Tree should be empty - no fields written
    EXPECT_FALSE(tree.keyExists("nstlog"));
    EXPECT_FALSE(tree.keyExists("nstxout"));
}

TEST_F(OutputControlModuleBuildMdpOutputTest, BuildsOutputWithDefaultValues)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing for this test
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Create KeyValueTree builder
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();

    // Get the MDP option provider
    auto* mdpProvider = module_->mdpOptionProvider();
    ASSERT_NE(mdpProvider, nullptr);

    // Build output
    mdpProvider->buildMdpOutput(&builderObject);

    // Get the resulting tree
    auto tree = builder.build();

    // Verify all 8 fields are present with default values
    EXPECT_EQ(tree["nstlog"].cast<int>(), 1000);
    EXPECT_EQ(tree["nstxout"].cast<int>(), 0);
    EXPECT_EQ(tree["nstvout"].cast<int>(), 0);
    EXPECT_EQ(tree["nstfout"].cast<int>(), 0);
    EXPECT_EQ(tree["nstenergy"].cast<int>(), 1000);
    EXPECT_EQ(tree["nstxout-compressed"].cast<int>(), 0);
    EXPECT_EQ(tree["compressed-x-precision"].cast<real>(), 1000.0);
    EXPECT_EQ(tree["nstcalcenergy"].cast<int>(), 100);
}

TEST_F(OutputControlModuleBuildMdpOutputTest, BuildsOutputWithCustomValues)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing for this test
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Set custom values
    outputControl_->nstlog                  = 500;
    outputControl_->nstxout                 = 1000;
    outputControl_->nstvout                 = 2000;
    outputControl_->nstfout                 = 3000;
    outputControl_->nstenergy               = 100;
    outputControl_->nstxout_compressed      = 5000;
    outputControl_->x_compression_precision = 2500.0;
    outputControl_->nstcalcenergy           = 50;

    // Create KeyValueTree builder
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();

    // Get the MDP option provider and build output
    auto* mdpProvider = module_->mdpOptionProvider();
    mdpProvider->buildMdpOutput(&builderObject);

    // Get the resulting tree
    auto tree = builder.build();

    // Verify all fields have the custom values
    EXPECT_EQ(tree["nstlog"].cast<int>(), 500);
    EXPECT_EQ(tree["nstxout"].cast<int>(), 1000);
    EXPECT_EQ(tree["nstvout"].cast<int>(), 2000);
    EXPECT_EQ(tree["nstfout"].cast<int>(), 3000);
    EXPECT_EQ(tree["nstenergy"].cast<int>(), 100);
    EXPECT_EQ(tree["nstxout-compressed"].cast<int>(), 5000);
    EXPECT_EQ(tree["compressed-x-precision"].cast<real>(), 2500.0);
    EXPECT_EQ(tree["nstcalcenergy"].cast<int>(), 50);
}

TEST_F(OutputControlModuleBuildMdpOutputTest, FieldNamesMatchMdpFormat)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing for this test
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Create KeyValueTree builder
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();

    // Build output
    auto* mdpProvider = module_->mdpOptionProvider();
    mdpProvider->buildMdpOutput(&builderObject);

    // Get the resulting tree
    auto tree = builder.build();

    // Verify field names use hyphens (MDP format) not underscores (C++ format)
    EXPECT_TRUE(tree.keyExists("nstxout-compressed"));
    EXPECT_TRUE(tree.keyExists("compressed-x-precision"));
    EXPECT_TRUE(tree.keyExists("nstcalcenergy"));

    // These keys should NOT exist (would be programming errors)
    EXPECT_FALSE(tree.keyExists("nstxout_compressed"));
    EXPECT_FALSE(tree.keyExists("x_compression_precision"));
}

TEST_F(OutputControlModuleBuildMdpOutputTest, IncludesCommentsInOutput)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing for this test
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Create KeyValueTree builder
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();

    // Build output
    auto* mdpProvider = module_->mdpOptionProvider();
    mdpProvider->buildMdpOutput(&builderObject);

    // Get the resulting tree
    auto tree = builder.build();

    // Verify all comment fields are present
    EXPECT_TRUE(tree.keyExists("comment-output-control-section"));
    EXPECT_TRUE(tree.keyExists("comment-output-control-coords"));
    EXPECT_TRUE(tree.keyExists("comment-output-control-energy"));
    EXPECT_TRUE(tree.keyExists("comment-output-control-compressed"));

    // Verify comment contents match legacy format
    EXPECT_EQ(tree["comment-output-control-section"].cast<std::string>(),
              "\n; OUTPUT CONTROL OPTIONS");
    EXPECT_EQ(tree["comment-output-control-coords"].cast<std::string>(),
              "; Output frequency for coords (x), velocities (v) and forces (f)");
    EXPECT_EQ(tree["comment-output-control-energy"].cast<std::string>(),
              "; Output frequency for energies to log file and energy file");
    EXPECT_EQ(tree["comment-output-control-compressed"].cast<std::string>(),
              "; Output frequency and precision for .xtc file");
}

TEST_F(OutputControlModuleBuildMdpOutputTest, NoTextOutputByDefault)
{
    // Do NOT enable output writing - test default behavior

    // Get the text output (should be empty)
    std::string output = buildMdpOutputAsString();

    // Check against reference data
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(output, "MdpOutput");
}

TEST_F(OutputControlModuleBuildMdpOutputTest, OutputOrderWithDefaultValues)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Get the text output
    std::string output = buildMdpOutputAsString();

    // Check against reference data
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(output, "MdpOutput");
}

TEST_F(OutputControlModuleBuildMdpOutputTest, OutputOrderWithCustomValues)
{
    // Setup preprocessing strings (required for MDP output in tests)
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);

    // Enable output writing
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Set custom values
    outputControl_->nstlog                  = 500;
    outputControl_->nstxout                 = 1000;
    outputControl_->nstvout                 = 2000;
    outputControl_->nstfout                 = 3000;
    outputControl_->nstenergy               = 100;
    outputControl_->nstxout_compressed      = 5000;
    outputControl_->x_compression_precision = 2500.0;
    outputControl_->nstcalcenergy           = 50;

    // Get the text output
    std::string output = buildMdpOutputAsString();

    // Check against reference data
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(output, "MdpOutput");
}

TEST_F(OutputControlModuleBuildMdpOutputTest, WritesGroupParameters)
{
    // Setup preprocessing strings
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "Protein";
    preprocessingStrings.energyGroups      = "Protein Non-Protein";

    // Configure module
    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);
    setOutputControlWriteToMdpOutput(module_.get(), true);

    // Build MDP output (accesses preprocessingStrings directly)
    std::string output = buildMdpOutputAsString();

    // Verify group parameters are present
    EXPECT_THAT(output, testing::HasSubstr("compressed-x-grps"));
    EXPECT_THAT(output, testing::HasSubstr("= Protein"));
    EXPECT_THAT(output, testing::HasSubstr("energygrps"));
    EXPECT_THAT(output, testing::HasSubstr("= Protein Non-Protein"));

    // Verify comments are present
    EXPECT_THAT(output, testing::HasSubstr("This selects the subset of atoms for the compressed"));
    EXPECT_THAT(output, testing::HasSubstr("Selection of energy groups"));

    // Check against reference data
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(output, "MdpOutput");
}

TEST_F(OutputControlModuleBuildMdpOutputTest, WritesEmptyGroupParameters)
{
    // Setup empty preprocessing strings
    gmx_inputrec_strings preprocessingStrings;
    preprocessingStrings.compressedXGroups = "";
    preprocessingStrings.energyGroups      = "";

    setOutputControlPreprocessingStrings(module_.get(), &preprocessingStrings);
    setOutputControlWriteToMdpOutput(module_.get(), true);

    std::string output = buildMdpOutputAsString();

    // Empty strings should still be written (means "all atoms")
    EXPECT_THAT(output, testing::HasSubstr("compressed-x-grps"));
    EXPECT_THAT(output, testing::HasSubstr("energygrps"));

    // Check against reference data
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(output, "MdpOutput");
}

} // namespace

// Death tests to verify assertions
TEST(OutputControlModuleDeathTest, SetPreprocessingStringsWithNull)
{
    auto module = OutputControlModuleInfo::create();

    // Should fail: passing null pointer
    GMX_EXPECT_DEATH_IF_SUPPORTED(setOutputControlPreprocessingStrings(module.get(), nullptr),
                                  "null pointer");
}

TEST(OutputControlModuleDeathTest, BuildMdpOutputWithoutPreprocessingStrings)
{
    auto module = OutputControlModuleInfo::create();

    // Set up OutputControl but NOT preprocessingStrings
    OutputControl outputControl;
    setOutputControlTarget(module.get(), &outputControl);

    // Enable MDP output (test-only mode)
    setOutputControlWriteToMdpOutput(module.get(), true);

    // Try to build MDP output - should fail assertion
    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();
    auto*               mdpProvider   = module->mdpOptionProvider();

    GMX_EXPECT_DEATH_IF_SUPPORTED(mdpProvider->buildMdpOutput(&builderObject),
                                  "preprocessingStrings not set");
}

} // namespace test
} // namespace gmx
