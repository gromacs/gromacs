/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Tests for the mdrun signalling functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/simulationsignal.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace test
{

//! Check that a null signaller can be called without problems
TEST(NullSignalTest, NullSignallerWorks)
{
    SimulationSignaller signaller(nullptr, nullptr, nullptr, false, false);
    EXPECT_EQ(0, signaller.getCommunicationBuffer().size());
    signaller.finalizeSignals();
}

//! Test fixture for mdrun signalling
class SignalTest : public ::testing::Test
{
public:
    SignalTest() : signals_{}
    {
        signals_[0].sig = 1;
        signals_[1].sig = -1;
    }
    //! Default object to hold signals
    SimulationSignals signals_;
};

TEST_F(SignalTest, NoSignalPropagatesIfNoSignallingTakesPlace)
{
    SimulationSignaller signaller(&signals_, nullptr, nullptr, false, false);
    EXPECT_EQ(0, signaller.getCommunicationBuffer().size());
    signaller.finalizeSignals();
    EXPECT_EQ(1, signals_[0].sig);
    EXPECT_EQ(-1, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(0, signals_[0].set);
    EXPECT_EQ(0, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, LocalIntraSimSignalPropagatesWhenIntraSimSignalTakesPlace)
{
    SimulationSignaller signaller(&signals_, nullptr, nullptr, false, true);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    signaller.finalizeSignals();
    EXPECT_EQ(0, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(1, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, LocalIntraSimSignalPropagatesWhenInterSimTakesPlace)
{
    SimulationSignaller signaller(&signals_, nullptr, nullptr, true, false);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    // Can't call finalizeSignals without a full commrec
    signaller.setSignals();
    EXPECT_EQ(0, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(1, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, LocalIntraSimSignalPropagatesWhenBothTakePlace)
{
    SimulationSignaller signaller(&signals_, nullptr, nullptr, true, true);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    // Can't call finalizeSignals without a full commrec
    signaller.setSignals();
    EXPECT_EQ(0, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(1, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, NonLocalSignalDoesntPropagateWhenIntraSimSignalTakesPlace)
{
    signals_[0].isLocal = false;
    SimulationSignaller signaller(&signals_, nullptr, nullptr, false, true);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    signaller.finalizeSignals();
    EXPECT_EQ(1, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(0, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, NonLocalSignalPropagatesWhenInterSimSignalTakesPlace)
{
    signals_[0].isLocal = false;
    SimulationSignaller signaller(&signals_, nullptr, nullptr, true, false);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    // Can't call finalizeSignals without a full commrec
    signaller.setSignals();
    EXPECT_EQ(0, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(1, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

TEST_F(SignalTest, NonLocalSignalPropagatesWhenBothTakePlace)
{
    signals_[0].isLocal = false;
    SimulationSignaller signaller(&signals_, nullptr, nullptr, true, true);
    EXPECT_NE(0, signaller.getCommunicationBuffer().size());
    // Can't call finalizeSignals without a full commrec
    signaller.setSignals();
    EXPECT_EQ(0, signals_[0].sig);
    EXPECT_EQ(0, signals_[1].sig);
    EXPECT_EQ(0, signals_[2].sig);
    EXPECT_EQ(1, signals_[0].set);
    EXPECT_EQ(-1, signals_[1].set);
    EXPECT_EQ(0, signals_[2].set);
}

} // namespace test
} // namespace gmx
