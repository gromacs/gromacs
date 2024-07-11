/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/ramd/randomsphericaldirectiongenerator.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

struct RandomSphericalDirectionGeneratorTest : public ::testing::Test
{
    int number_of_directions = 1000000;
    int number_of_buckets    = 32;

    std::vector<int> create_histogram(int number_of_buckets, std::vector<double> const& data)
    {
        const double     bucket_size = 2 * M_PI / number_of_buckets;
        std::vector<int> histogram(number_of_buckets, 0);

        for (auto const& d : data)
        {
            ++histogram[floor((d + M_PI) / bucket_size)];
        }

        return histogram;
    }

    void check_distribution(int x, int y, std::vector<DVec> const& random_directions)
    {
        std::vector<double> thetas(random_directions.size());
        for (size_t i = 0; i < random_directions.size(); ++i)
        {
            thetas[i] = std::atan2(random_directions[i][y], random_directions[i][x]);
        }

        auto histogram   = create_histogram(32, thetas);
        int  mean        = number_of_directions / number_of_buckets;
        int  lower_bound = mean - mean * 0.1;
        int  upper_bound = mean + mean * 0.1;

        for (auto const& h : histogram)
        {
            EXPECT_GT(h, lower_bound);
            EXPECT_LT(h, upper_bound);
        }
    }
};

TEST_F(RandomSphericalDirectionGeneratorTest, angle_distribution)
{
    RandomSphericalDirectionGenerator random_spherical_direction_generator(1234);

    std::vector<DVec> random_directions(number_of_directions);
    for (auto& d : random_directions)
    {
        d = random_spherical_direction_generator();
    }

    for (auto const& d : random_directions)
    {
        EXPECT_NEAR(1.0, d.norm2(), 1e-6);
    }

    check_distribution(0, 1, random_directions); // xy
    check_distribution(0, 2, random_directions); // xz
    check_distribution(2, 1, random_directions); // zy
}

} // namespace
} // namespace test
} // namespace gmx
