/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Tests for the EnergyDriftTacker class
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/energydrifttracker.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(EnergyDriftTracker, emptyWorks)
{
    EnergyDriftTracker tracker(1);

    EXPECT_EQ(tracker.timeInterval(), 0);
    EXPECT_EQ(tracker.energyDrift(), 0);
}

TEST(EnergyDriftTracker, onePointWorks)
{
    EnergyDriftTracker tracker(1);

    tracker.addPoint(1.5, -3.5_real);
    EXPECT_EQ(tracker.timeInterval(), 0);
    EXPECT_EQ(tracker.energyDrift(), 0);
}

TEST(EnergyDriftTracker, manyPointsWorks)
{
    EnergyDriftTracker tracker(10);

    tracker.addPoint(1.5, 2.5_real);
    tracker.addPoint(3.5, 4.0_real);
    tracker.addPoint(5.5, -5.5_real);
    EXPECT_FLOAT_EQ(tracker.timeInterval(), 4.0_real);
    EXPECT_FLOAT_EQ(tracker.energyDrift(), -0.2_real);
}

} // namespace
} // namespace test
} // namespace gmx
