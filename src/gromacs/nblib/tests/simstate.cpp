/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * This implements SimState tests
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 */

#include <vector>

#include "gmxpre.h"

#include "gromacs/nblib/box.h"
#include "gromacs/nblib/setup.h"
#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

class SimStateTester
{
public:
    std::vector<gmx::RVec> coords_;
    std::vector<gmx::RVec> vel_;

    ::nblib::Box box_;
    ::nblib::TopologyBuilder topologyBuilder_;

    SimStateTester()
            :box_(6.05449)
    {
        coords_ =
        {

            { 0.794,  1.439,  0.610 },
            { 1.397,  0.673,  1.916 },
            { 0.659,  1.080,  0.573 },
            { 1.105,  0.090,  3.431 },
            { 1.741,  1.291,  3.432 },
            { 1.936,  1.441,  5.873 },
            { 0.960,  2.246,  1.659 },
            { 0.382,  3.023,  2.793 },
            { 0.053,  4.857,  4.242 },
            { 2.655,  5.057,  2.211 },
            { 4.114,  0.737,  0.614 },
            { 5.977,  5.104,  5.217 }
        };
        vel_    =
        {
            { 0.326183,  -0.119182,  0.0135174 },
            { 0.355938,  0.0374831,  0.252879 },
            { 0.383617,  0.285418,  0.426147 },
            { -0.0766866,  0.551282,  -0.37828 },
            { 0.0502158,  0.263403,  0.0817824 },
            { 0.286086,  -0.0559638,  0.0653588 },
            { 0.310655,  -0.282913,  -0.08545 },
            { 0.154439,  -0.0347022,  0.299215 },
            { 0.14032,  0.232671,  0.278278 },
            { 0.391716,  -0.0599792,  0.0041444 },
            { -0.0910481,  -0.0146716,  -0.25681 },
            { -0.140364,  -0.371435,  -0.256161 }
        };


        topologyBuilder_.setAtomTypes({0});
        topologyBuilder_.setCharges({1.});
        topologyBuilder_.setMasses({1.});
        topologyBuilder_.setExclusions({},{}, 1);
        topologyBuilder_.setNonbondedParameters({0,1});
    }

    void setupSimState()
    {
        auto topology = topologyBuilder_.buildTopology(coords_.size());
        ::nblib::SimState(coords_, box_, topology, vel_);
    }
};

TEST(NBlibTest, SimStateArgonBox)
{
    SimStateTester simStateTester;
    EXPECT_NO_THROW(simStateTester.setupSimState());
}

TEST(NBlibTest, SimStateArgonBoxCoordThrowNAN)
{
    SimStateTester simStateTester;
    simStateTester.coords_[2][0] = NAN;
    EXPECT_THROW(simStateTester.setupSimState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimStateArgonBoxCoordThrowINF)
{

    SimStateTester simStateTester;
    simStateTester.coords_[2][0] = INFINITY;
    EXPECT_THROW(simStateTester.setupSimState(), gmx::InvalidInputError);
}

}  // namespace
}  // namespace test
}  // namespace gmx