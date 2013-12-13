/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
 * Tests utilities for "Computational Electrophysiology" setups.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun
 */
#include <gtest/gtest.h>
#include "moduletest.h"
#include "gromacs/options/filenameoption.h"
#include "programs/mdrun/tests/moduletest.h"

namespace
{

//! Test fixture for mdrun with "Computational Electrophysiology" settings,
// i.e. double membrane sandwich with ion/water exchange protocol
typedef gmx::test::SwapTestFixture CompelTest;

/* This test ensures that the compel protocol can be run, and that all of
 * the swapcoords parameters from the .mdp file are understood */
TEST_F(CompelTest, SwapCanRun)
{
    std::string name = "OctaneSandwich";
    topFileName      = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    groInputFileName = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    mdpInputFileName = fileManager_.getInputFilePath((std::string(name) + ".mdp").c_str());
    ndxFileName      = fileManager_.getInputFilePath((std::string(name) + ".ndx").c_str());
    
    swapFileName     = fileManager_.getTemporaryFilePath("swap.xvg");
            
    EXPECT_EQ(0, callGrompp());

    ASSERT_EQ(0, callMdrun());

    ASSERT_EQ(0, callMdrunAppend()); // with checkpoint
}



/*! \todo Add other tests for the compel module, e.g.
 *
 *  - a test that checks that actually ion/water swaps have been done, by
 *    calling gmxcheck on the swap output file and a reference file
 */

} // namespace
