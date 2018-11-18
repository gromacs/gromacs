/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests for hardwareOptions
 *
 * \author Kevin Boyd <kevin.boyd@uconn.edu>
 * \ingroup module_hardware
 */

#include "gmxpre.h"

#include "gromacs/hardware/hw_info.h"

#include <gtest/gtest.h>

namespace test
{
namespace
{

TEST(hardwareOptionsManagerTest, assignmentAndAccessingWorks)
{
    // construct default hw_opt_t, change some but not all of the options
    gmx_hw_opt_t userOptions;
    userOptions.nthreads_tmpi   = 4;
    userOptions.gpuIdsAvailable = "01";

    hardwareOptionsManager hwOpts(userOptions);
    // user set options
    EXPECT_EQ(4, hwOpts.nthreads_tmpi());
    EXPECT_EQ("01", hwOpts.gpuIdsAvailable());

    // default options
    EXPECT_EQ(0, hwOpts.nthreads_tot());
    EXPECT_EQ("", hwOpts.gpuTaskAssignment());

    // check boolean queries
    EXPECT_TRUE(hwOpts.nthreads_tmpi.isSetByUser());
    EXPECT_FALSE(hwOpts.nthreads_tot.isSetByUser());

}

TEST(hardwareOptionsManagerTest, setOptionsWork)
{
    gmx_hw_opt_t           userOptions;
    userOptions.nthreads_tmpi = 4;
    hardwareOptionsManager hwOpts(userOptions);

    // amend a value not defined by a user
    hwOpts.nthreads_tot.set(2);
    EXPECT_EQ(2, hwOpts.nthreads_tot());

    // Unset value
    EXPECT_FALSE(hwOpts.nthreads_omp.isSet());
    EXPECT_FALSE(hwOpts.nthreads_omp.isSetByUser());
    // Value set after instantiation
    EXPECT_TRUE(hwOpts.nthreads_tot.isSet());
    EXPECT_FALSE(hwOpts.nthreads_tot.isSetByUser());
    // Value set at instantiation
    EXPECT_TRUE(hwOpts.nthreads_tmpi.isSet());
    EXPECT_TRUE(hwOpts.nthreads_tmpi.isSetByUser());

    // try to amend a user-set value, should trip exception
    EXPECT_THROW(hwOpts.nthreads_tmpi.set(2), gmx::APIError);
}


} // namespace
} // namespace test
