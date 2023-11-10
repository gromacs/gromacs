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

#include "nblib/listed_forces/positionrestraints.hpp"

#include "gromacs/utility/arrayref.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "nblib/box.h"
#include "nblib/listed_forces/dataflowrestraints.hpp"
#include "nblib/tests/testhelpers.h"

#include "listedtesthelpers.h"

namespace nblib
{
namespace test
{
namespace
{

class PositionRestraintsTest : public ::testing::TestWithParam<std::tuple<RefCoordScaling, PbcType>>
{
protected:
    Box             box_;
    PbcType         pbcType_;
    RefCoordScaling refCoordScaling_;

    PositionRestraintsTest() : box_(0.9, 1.0, 1.1)
    {
        refCoordScaling_ = std::get<0>(GetParam());
        pbcType_         = std::get<1>(GetParam());
        PbcHolder                         pbcHolder(pbcType_, box_);
        auto                              indices = indexVector<PositionRestraints>();
        std::vector<util::array<real, 3>> coords  = { { 0.4, 0.5, 0.6 } };
        std::vector<util::array<real, 3>> forces(coords.size(), util::array<real, 3>{ 0, 0, 0 });
        std::vector<real>                 virials(9, 0);
        util::array<real, 3>              com = { 0, 0.5, 0 };
        test::RefDataChecker              refDataChecker(1e-4);
        gmx::ArrayRef<const util::array<real, 3>> coordinates             = coords;
        std::vector<PositionRestraints>           InputPositionRestraints = { PositionRestraints(
                0.5, 0.6, 0.0, 0, 200, 400) };

        auto energy = computeForcesRestraints(
                gmx::ArrayRef<const InteractionIndex<PositionRestraints>>(indices),
                gmx::ArrayRef<const PositionRestraints>(InputPositionRestraints),
                gmx::ArrayRef<const PositionRestraints>(InputPositionRestraints),
                coordinates,
                &forces,
                gmx::ArrayRef<real>(virials),
                pbcHolder,
                box_,
                pbcType_,
                refCoordScaling_,
                com);

        refDataChecker.testReal(energy.potentialEnergy(), "Epot");

        refDataChecker.test3DVectors<util::array<real, 3>>(forces, "forces");

        refDataChecker.testArrays<real>(virials, "virials");
    }
};

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };
//! Reference Coordinate Scaling values for testing
std::vector<RefCoordScaling> c_refCoordScalingForTests = { RefCoordScaling::No,
                                                           RefCoordScaling::Com,
                                                           RefCoordScaling::All };

TEST_P(PositionRestraintsTest, BasicPosResNoFreeEnergy)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
}

INSTANTIATE_TEST_SUITE_P(PosResBasicTest,
                         PositionRestraintsTest,
                         ::testing::Combine(::testing::ValuesIn(c_refCoordScalingForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));


} // namespace
} // namespace test
} // namespace nblib
