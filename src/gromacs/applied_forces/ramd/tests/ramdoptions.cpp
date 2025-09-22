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

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

TEST(RAMDOptionsTest, DefaultParameters)
{
    RAMDOptions ramdOptions;
    const auto defaultParameters = ramdOptions.parameters();
    EXPECT_FALSE(defaultParameters.active_);
    EXPECT_EQ(1234, defaultParameters.seed_);
    EXPECT_EQ(0, defaultParameters.ngroups_);

    // EXPECT_REAL_EQ(0.0025, defaultParameters.groups_[0].r_min_dist_);
}


} // namespace

} // namespace gmx
