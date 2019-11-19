#include <vector>

#include "gromacs/nblib/util.h"
#include "gromacs/math/vectypes.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(NBlibTest, checkNumericValues)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, true);
}

TEST(NBlibTest, checkNumericValuesHasNan)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    vec.push_back({ NAN, NAN, NAN });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}

TEST(NBlibTest, checkNumericValuesHasInf)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    vec.push_back({ INFINITY, INFINITY, INFINITY });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}


}  // namespace
}  // namespace test
}  // namespace gmx
