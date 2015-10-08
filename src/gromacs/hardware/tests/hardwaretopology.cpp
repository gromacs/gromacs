/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests for gmx::HardwareTopology
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/hardware/hardwaretopology.h"

#include <gtest/gtest.h>

namespace
{

TEST(HardwareTopologyTest, Execute)
{
    // There is no way we can compare to any reference data since that
    // depends on the architecture, but we can at least make sure that it
    // works to execute the tests

    gmx::HardwareTopology hwTop;

    // If we cannot even find the number of logical processors we want to flag it
    EXPECT_GT(hwTop.supportLevel(), gmx::HardwareTopology::SupportLevel::None);

    if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
    {
        int socketsInSystem  = hwTop.machine().sockets.size();
        int coresPerSocket   = hwTop.machine().sockets[0].cores.size();
        int hwThreadsPerCore = hwTop.machine().sockets[0].cores[0].hwThreads.size();

        // Check that logical processor information is reasonable
        for (auto &l : hwTop.machine().logicalProcessors)
        {
            EXPECT_TRUE(l.socket >= 0 && l.socket < socketsInSystem);
            EXPECT_TRUE(l.core >= 0 && l.core < coresPerSocket);
            EXPECT_TRUE(l.hwThread >= 0 && l.hwThread < hwThreadsPerCore);
        }

        // Double-check that the tree is self-consistent with logical processor info
        for (int s = 0; s < socketsInSystem; s++)
        {
            for (int c = 0; c < coresPerSocket; c++)
            {
                for (int t = 0; t < hwThreadsPerCore; t++)
                {
                    int idx = hwTop.machine().sockets[s].cores[c].hwThreads[t].logicalProcessorId;
                    EXPECT_LT(idx, hwTop.machine().logicalProcessorCount);
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].socket, s);
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].core, c);
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].hwThread, t);
                }
            }
        }
    }
}

} // namespace
