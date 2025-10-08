/*! \internal \file
 * \brief
 * Tests for RAMD module options.
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/ramd/ramd.h"
#include "gromacs/applied_forces/ramd/ramdoptions.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

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
};

TEST_F(RAMDOptionsTest, OptionSetsActive)
{
    RAMDOptions ramdOptions;
    EXPECT_FALSE(ramdOptions.parameters().active_);
    test::fillOptionsFromMdpValues(ramdBuildDefaultMdpValues(), &ramdOptions);
    EXPECT_TRUE(ramdOptions.parameters().active_);
}

TEST_F(RAMDOptionsTest, DefaultParameters)
{
    RAMDOptions ramdOptions;
    const auto defaultParameters = ramdOptions.parameters();
    EXPECT_FALSE(defaultParameters.active_);
    EXPECT_EQ(1234, defaultParameters.seed_);
    EXPECT_EQ(0, defaultParameters.ngroups_);

    // EXPECT_REAL_EQ(0.0025, defaultParameters.groups_[0].r_min_dist_);
}


} // namespace anonymous
} // namespace gmx
