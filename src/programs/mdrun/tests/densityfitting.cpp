/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests molecular dynamics with a density fitting potential.
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/grompp.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class DensfitTestFixture : public MdrunTestFixture
{
    protected:
        DensfitTestFixture();
        ~DensfitTestFixture();
};


DensfitTestFixture::DensfitTestFixture()
{
}

DensfitTestFixture::~DensfitTestFixture()
{
}



//! Test fixture for mdrun with density fitting switched on
typedef gmx::test::DensfitTestFixture DensfitTest;

/* This basic test ensures that the mdrun with an added density fitting potential
 * run, and that the related .mdp input parameters are understood.
 */
TEST_F(DensfitTest, DensfitCanRun)
{
    // Prepare .tpr file:
    std::string name = "1AKE";
    runner_.topFileName_      = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    runner_.groFileName_      = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    runner_.mdpInputFileName_ = fileManager_.getInputFilePath((std::string(name) + ".mdp").c_str());
    auto        mapFileName          = fileManager_.getInputFilePath("4AKE.ccp4");
    auto        densfitOutFileName   = fileManager_.getTemporaryFilePath((name + "_densfit.xvg").c_str());

    CommandLine gromppCaller;
    gromppCaller.append("grompp");
    gromppCaller.addOption("-f", runner_.mdpInputFileName_);
    gromppCaller.addOption("-p", runner_.topFileName_);
    gromppCaller.addOption("-c", runner_.groFileName_);
    gromppCaller.addOption("-mi", mapFileName);

    gromppCaller.addOption("-po", runner_.mdpOutputFileName_);
    gromppCaller.addOption("-o", runner_.tprFileName_);

    EXPECT_EQ(0, gmx_grompp(gromppCaller.argc(), gromppCaller.argv()));

    // Now run the .tpr file:
    CommandLine mdrunCaller;
    mdrunCaller.append("mdrun");
    mdrunCaller.append("-v");
    mdrunCaller.addOption("-d", densfitOutFileName);

    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
}


} // namespace test
} // namespace gmx
