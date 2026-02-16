/*! \internal \file
 * \brief
 * Tests for RAMD module options.
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/ramd/ramdoptions.h"

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/ramd/ramd.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/logger.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace
{

class RAMDOptionsTest : public ::testing::Test
{
public:
    static KeyValueTreeObject ramdBuildDefaultMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        return mdpValueBuilder.build();
    }

    static KeyValueTreeObject ramdBuildMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-seed",
                                              std::string("42"));
        return mdpValueBuilder.build();
    }
};

TEST_F(RAMDOptionsTest, DefaultParameters)
{
    RAMDOptions ramdOptions;
    const auto  defaultParameters = ramdOptions.parameters();
    EXPECT_FALSE(defaultParameters.active_);
    EXPECT_EQ(1234, defaultParameters.seed_);
    EXPECT_EQ(0, defaultParameters.ngroups_);
}

TEST_F(RAMDOptionsTest, OptionSetsActive)
{
    RAMDOptions ramdOptions;
    test::fillOptionsFromMdpValues(ramdBuildMdpValues(), &ramdOptions);

    EXPECT_TRUE(ramdOptions.active());
    EXPECT_TRUE(ramdOptions.parameters().active_);
    EXPECT_EQ(42, ramdOptions.parameters().seed_);

    // Write parameters to the KVT
    KeyValueTreeBuilder builder;
    MDLogger            logger;
    ramdOptions.setLogger(logger);
    ramdOptions.writeInternalParametersToKvt(builder.rootObject());
    const auto inputTree = builder.build();

    // Retrieve parameters from the KVT
    ramdOptions.readInternalParametersFromKvt(inputTree);

    EXPECT_TRUE(ramdOptions.active());
    EXPECT_TRUE(ramdOptions.parameters().active_);
    EXPECT_EQ(42, ramdOptions.parameters().seed_);
}

} // namespace
} // namespace gmx
