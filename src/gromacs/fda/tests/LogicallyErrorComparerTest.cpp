/*
 * LogicallyErrorComparerTest.cpp
 *
 *  Created on: Jun 26, 2017
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <gtest/gtest.h>
#include "testutils/LogicallyErrorComparer.h"
#include <tuple>
#include <vector>

namespace fda
{

TEST(LogicallyErrorComparerTest, DISABLED_Base)
{
    LogicallyEqualComparer<true, true> comparer(1.0);

    std::vector<std::tuple<float, float, bool>> input{
        std::make_tuple(1000000, 1000001, true),
        std::make_tuple(1000, 1001, false),
        std::make_tuple(1.0000001, 1.0000002, true),
        std::make_tuple(1.0002, 1.0001, false),
        std::make_tuple(0.000000001000001, 0.000000001000002, true),
        std::make_tuple(0.000000000001002, 0.000000000001001, false)
    };

    for (auto const& e : input) {
        EXPECT_TRUE(comparer(std::get<0>(e), std::get<1>(e)) == std::get<2>(e)) << std::get<0>(e) << " " << std::get<1>(e) << " " << std::get<2>(e);
        EXPECT_TRUE(comparer(std::get<1>(e), std::get<0>(e)) == std::get<2>(e)) << std::get<0>(e) << " " << std::get<1>(e) << " " << std::get<2>(e);
    }
}

TEST(LogicallyErrorComparerTest, SpecialTest1)
{
    LogicallyEqualComparer<true, true> comparer(1e2);

    EXPECT_TRUE(comparer(1.516308e-06, -1.140742e-06));
}

TEST(LogicallyErrorComparerTest, SpecialTest2)
{
    LogicallyEqualComparer<true, true> comparer(1e3);

    EXPECT_TRUE(comparer(-4.506144e+00, -4.505985e+00));
}

TEST(LogicallyErrorComparerTest, SpecialTest3)
{
    LogicallyEqualComparer<true, true> comparer(1e4);

    EXPECT_TRUE(comparer(5.264869e-01, 5.262685e-01));
}

} // namespace fda
