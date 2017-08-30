#include "gmxapi/version.h"
#include <gtest/gtest.h>


using namespace gmxapi;

const int current_major = GMXAPI_MAJOR;
const int current_minor = GMXAPI_MINOR;
const int current_patch = GMXAPI_PATCH;

namespace
{

// The is_at_least functin uses major(), minor(), and patch()
// so testing them is maybe superfluous

TEST(VersionTest, SaneComparisons)
{
    ASSERT_TRUE(Version::is_at_least(0, 0, 0));
    ASSERT_FALSE(Version::is_at_least(-1, -1, -1));
    ASSERT_TRUE(Version::is_at_least(current_major, current_minor, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major + 1, current_minor, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major, current_minor + 1, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major, current_minor, current_patch + 1));
}

} // end anonymous namespace
