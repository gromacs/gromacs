/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * This implements position restraints tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "gromacs/listed_forces/position_restraints.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "nblib/box.h"
#include "nblib/listed_forces/positionrestraints.hpp"
#include "nblib/vector.h"

namespace nblib
{
namespace test
{
namespace
{

class PosreDxTest : public ::testing::TestWithParam<std::tuple<RefCoordScaling, PbcType>>
{
protected:
    t_pbc   pbc_;
    PbcType pbcType_;
    Box     box_;
    Vec3    posres_com_ = { 0, 0.5, 0 };
    Vec3    zero_       = { 0, 0, 0 };
    int     npbcdim_;

    RefCoordScaling refCoordScaling_;

    PosreDxTest() : box_(0.9, 1.0, 1.1)
    {
        refCoordScaling_ = std::get<0>(GetParam());
        pbcType_         = std::get<1>(GetParam());

        Vec3 x    = { 0.4, 0.5, 0.6 };
        Vec3 pos0 = { 0.5, 0.6, 0.0 };
        set_pbc(&pbc_, pbcType_, box_.legacyMatrix());
        PbcHolder pbcHolder(pbcType_, box_);
        npbcdim_ = numPbcDimensions(pbcType_);

        Vec3 gmxDx, gmxRDist, gmxDpdl;
        gmx::posres_dx(
                x, pos0, zero_, posres_com_, zero_, 0, &pbc_, refCoordScaling_, npbcdim_, gmxDx, gmxRDist, gmxDpdl);
        {
            auto [nblibDx, nblibRDist, nblibDpdl] = nblib::posres_dx(
                    x, pos0, zero_, posres_com_, zero_, 0, pbcHolder, box_, refCoordScaling_, npbcdim_);
            EXPECT_TRUE(gmxDx == nblibDx);
            EXPECT_TRUE(gmxRDist == nblibRDist);
            EXPECT_TRUE(gmxDpdl == nblibDpdl);
        }
        {
            auto [nblibDx, nblibRDist] =
                    nblib::posres_dx(x, pos0, posres_com_, pbcHolder, box_, refCoordScaling_, npbcdim_);
            EXPECT_TRUE(gmxDx == nblibDx);
            EXPECT_TRUE(gmxRDist == nblibRDist);
        }
    }
};

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };
//! Reference Coordinate Scaling values for testing
std::vector<RefCoordScaling> c_refCoordScalingForTests = { RefCoordScaling::No,
                                                           RefCoordScaling::Com,
                                                           RefCoordScaling::All };

TEST_P(PosreDxTest, PosreDx)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
}

INSTANTIATE_TEST_SUITE_P(PosReDxBasicTest,
                         PosreDxTest,
                         ::testing::Combine(::testing::ValuesIn(c_refCoordScalingForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));


} // namespace
} // namespace test
} // namespace nblib
