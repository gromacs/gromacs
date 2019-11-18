/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2019, by the GROMACS development team, led by
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
 * Tests for infrastructure for running tests under MPI.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/mpitest.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"

namespace
{

class MpiSelfTest : public ::testing::Test
{
public:
    MpiSelfTest() : reached{ 0, 0 } {}

    int reached[2];
};

TEST_F(MpiSelfTest, Runs)
{
    GMX_MPI_TEST(2);
#if GMX_THREAD_MPI
    reached[gmx_node_rank()] = 1;
    MPI_Barrier(MPI_COMM_WORLD);
#else
    int value = 1;
    MPI_Gather(&value, 1, MPI_INT, reached, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (gmx_node_rank() == 0)
    {
        EXPECT_EQ(1, reached[0]);
        EXPECT_EQ(1, reached[1]);
    }
}

} // namespace
