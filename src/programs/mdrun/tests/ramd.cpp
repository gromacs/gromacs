/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright None- The GROMACS Authors
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
/**
 * \internal \file
 *
 * \brief
 * Implements test of RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 *
 */

#include "gmxpre.h"

#include "gromacs/mdlib/sighandler.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "programs/mdrun/tests/moduletest.h"

namespace gmx
{
namespace test
{

class RAMDTestFixture : public MdrunTestFixture
{
protected:
    RAMDTestFixture();
    ~RAMDTestFixture() override;
};


RAMDTestFixture::RAMDTestFixture() {}

RAMDTestFixture::~RAMDTestFixture() {}


//! Test fixture for mdrun with RAMD settings
typedef gmx::test::RAMDTestFixture RAMDTest;

const std::string water4_mdp_base = R"(
    integrator               = md
    dt                       = 0.001
    nsteps                   = 100000
    nstlog                   = 10
    rlist                    = 1.0
    coulombtype              = Cut-off
    rcoulomb-switch          = 0
    rcoulomb                 = 1.0
    epsilon-r                = 1
    epsilon-rf               = 1
    vdw-type                 = Cut-off
    rvdw-switch              = 0
    rvdw                     = 1.0
    DispCorr                 = no
    Tcoupl                   = no
    Pcoupl                   = no

    ramd                     = yes
    ramd-seed                = 1234
    ramd-eval-freq           = 10
    ramd-force-out-freq      = 10
    ramd-old-angle-dist      = no
    ramd-ngroups             = 2
    ramd-group1-receptor     = 1SOL
    ramd-group1-ligand       = 2SOL
    ramd-group1-force        = 100
    ramd-group1-max-dist     = 1.0
    ramd-group1-r-min-dist   = 0.0025
    ramd-group2-receptor     = 1SOL
    ramd-group2-ligand       = 3SOL
    ramd-group2-force        = 100
    ramd-group2-max-dist     = 1.0
    ramd-group2-r-min-dist   = 0.0025
)";

TEST_F(RAMDTest, RAMD_connected_ligands)
{
    runner_.useTopGroAndNdxFromDatabase("4water");
    auto mdpContents = water4_mdp_base + R"(
        ramd-connected-ligands = no
    )";
    runner_.useStringAsMdpFile(mdpContents);

    CommandLine caller;
    caller.addOption("-ramd");
    caller.addOption("-reprod");

    EXPECT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(2, runner_.callMdrun(caller));
    gmx_reset_stop_condition();

    TextReader  reader_log(runner_.logFileName_);
    std::string line;
    int         number_of_steps = -1;
    while (reader_log.readLine(&line))
    {
        if (line.find("==== RAMD ==== GROMACS will be stopped after") != std::string::npos)
        {
            number_of_steps = stoi(gmx::splitString(line)[8]);
        }
    }
    EXPECT_EQ(number_of_steps, 630);

    TextReader reader_pullx(fileManager_.getTemporaryFilePath("state_pullx.xvg"));
    // std::cout << reader_pullx.readAll();
    while (reader_pullx.readLine(&line))
    {
        if (line.rfind("0.000", 0) != std::string::npos)
        {
            EXPECT_EQ(std::string("0.0593702"), gmx::splitString(line)[1]);
        }
    }

    TextReader reader_ramd(fileManager_.getTemporaryFilePath("state.xvg"));
    // std::cout << reader_ramd.readAll();
    while (reader_ramd.readLine(&line))
    {
        if (line.rfind("0.000", 0) != std::string::npos)
        {
            EXPECT_EQ(std::string("0.423961"), gmx::splitString(line)[1]);
        }
    }
}

} // namespace test
} // namespace gmx
