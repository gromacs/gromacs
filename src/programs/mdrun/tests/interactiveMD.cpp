/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Tests utilities for interactive molecular dynamics (IMD) setups.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class ImdTestFixture : public MdrunTestFixture
{
    protected:
        ImdTestFixture();
        ~ImdTestFixture();
};


ImdTestFixture::ImdTestFixture()
{
}

ImdTestFixture::~ImdTestFixture()
{
}


//! Test fixture for mdrun with IMD settings
typedef gmx::test::ImdTestFixture ImdTest;

/* If GROMACS was compiled with IMD support, this test checks
 * - whether the IMD-group parameter from the .mdp file is understood,
 * - whether mdrun understands the IMD-related command line parameters -imdpull, -imdwait, -imdterm,
 * - whether mdrun finishes without error when IMD is enabled.
 */
TEST_F(ImdTest, ImdCanRun)
{
    std::string name = "spc2";
    runner_.useTopGroAndNdxFromDatabase(name.c_str());
    runner_.mdpInputFileName_ = fileManager_.getInputFilePath((name + "-IMD.mdp").c_str());

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine imdCaller;
    imdCaller.append("mdrun");

    imdCaller.addOption("-imdport", 0); // automatically assign a free port
    imdCaller.append("-imdpull");
    imdCaller.append("-noimdwait");     // cannot use -imdwait: then mdrun would not return control ...
    imdCaller.append("-noimdterm");

    // Do an mdrun with IMD enabled
    ASSERT_EQ(0, runner_.callMdrun(imdCaller));
}



} // namespace test
} // namespace gmx
