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
/*!\internal
 * \file
 * \brief
 * Tests for gmx::SetStartTime
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include <memory>

#include "gromacs/coordinateio/outputadapters/setstarttime.h"
#include "gromacs/trajectory/trajectoryframe.h"

#include "gromacs/coordinateio/tests/coordinate_test.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Test fixture to prepare a setatoms object to pass data through.
 */
class SetStartTimeTest : public gmx::test::CommandLineTestBase
{
public:
    SetStartTimeTest() { clear_trxframe(frame(), true); }
    /*! \brief
     * Get access to the method for changing frame time information.
     *
     * \param[in] startTime User supplied start time to test.
     */
    SetStartTime* setStartTime(real startTime)
    {
        if (!setStartTime_)
        {
            setStartTime_ = std::make_unique<SetStartTime>(startTime);
        }
        return setStartTime_.get();
    }
    //! Get access to trajectoryframe to mess with.
    t_trxframe* frame() { return &frame_; }

private:
    //! Object to use for tests
    SetStartTimePointer setStartTime_;
    //! Storage of trajectoryframe.
    t_trxframe frame_;
};

TEST_F(SetStartTimeTest, WorksWithNonZeroStart)
{
    frame()->bTime = false;
    frame()->time  = 5;
    // Set step to nonsense value to check that it is ignored.
    SetStartTime* method = setStartTime(42);
    EXPECT_NO_THROW(method->processFrame(0, frame()));
    EXPECT_TRUE(frame()->bTime);
    EXPECT_EQ(frame()->time, 42);
    // Set frame time again to simulate advancing to new time.
    frame()->time = 15;
    EXPECT_NO_THROW(method->processFrame(1, frame()));
    EXPECT_EQ(frame()->time, 52);
    // Expect to use next frame time to get correct time step again.
    frame()->time = 20;
    EXPECT_NO_THROW(method->processFrame(2, frame()));
    EXPECT_EQ(frame()->time, 57);
}

TEST_F(SetStartTimeTest, WorksWithZeroStart)
{
    frame()->time        = 42;
    SetStartTime* method = setStartTime(0);
    EXPECT_NO_THROW(method->processFrame(0, frame()));
    EXPECT_EQ(frame()->time, 0);
    // No matter what the next time in the frame is, ignore it.
    frame()->time = 65;
    EXPECT_NO_THROW(method->processFrame(1, frame()));
    EXPECT_EQ(frame()->time, 23);
    // And so on for more frames.
    frame()->time = 72;
    EXPECT_NO_THROW(method->processFrame(2, frame()));
    EXPECT_EQ(frame()->time, 30);
}

} // namespace test

} // namespace gmx
