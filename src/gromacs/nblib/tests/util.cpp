#include <vector>

#include "gromacs/nblib/util.h"
#include "gromacs/math/vectypes.h"

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

    vec.push_back({ (real)NAN, (real)NAN, (real)NAN });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}

TEST(NBlibTest, checkNumericValuesHasInf)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    vec.push_back({ (real)INFINITY, (real)INFINITY, (real)INFINITY });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}


TEST(NBlibTest, generateVelocity)
{
    constexpr size_t N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);

    std::vector<gmx::RVec> expected = {
        { 2.698104, 0.553971 ,-0.669996 },
        { 1.208262, -0.947503 ,-0.393945 },
        { -0.256397, 0.150944 ,-1.902301 },
        { -2.665339, 1.028487, 1.863356 },
        { -0.519059, -1.580847, 0.596605 },
        { -1.535892, -4.004550, 2.329542 },
        { 2.046137, -0.657188 ,-0.847896 },
        { 0.524716, 2.047179, 1.075778 },
        { -0.530676, 1.008563, 1.509182 },
        { 0.710458, -1.426227, 2.217572 }
    };

    for(size_t i = 0; i < N; i++) {
        for (int m = 0; (m < DIM); m++)
        {
            EXPECT_REAL_EQ_TOL(out[XX][XX], expected[XX][XX], defaultRealTolerance());
        }
    }
}
TEST(NBlibTest, generateVelocitySize)
{
    constexpr int N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);
    EXPECT_EQ(out.size(), N);
}

TEST(NBlibTest, generateVelocityCheckNumbers)
{
    constexpr int N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);
    bool ret = checkNumericValues(out);
    EXPECT_EQ(ret, true);
}

}  // namespace
}  // namespace test
}  // namespace gmx
