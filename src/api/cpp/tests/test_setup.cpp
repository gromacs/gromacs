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
//
// Created by Eric Irrgang on 10/18/17.
//

#include "gmxpre.h"

#include <gtest/gtest.h>
#include <fstream>

#include "moduletest.h"

namespace
{

static const std::string atomicConfiguration {
    R"delimiter(spc-and-methanol
    6
    1MeOH   Me1    1   1.965   1.457   1.206 -0.8621 -0.4370 -0.5161
    1MeOH    O2    2   1.978   1.412   1.078 -0.0083 -0.4996 -0.6302
    1MeOH    H3    3   1.903   1.451   1.025 -0.2944 -1.5927 -1.1407
    2SOL     OW    4   1.559   1.516   0.709  0.7434  0.8949  1.0702
    2SOL    HW1    5   1.502   1.495   0.789  0.5852 -0.0293  0.7108
    2SOL    HW2    6   1.501   1.532   0.630  0.9009  1.8670  1.1441
   30.0000   30.0000   30.0000
)delimiter"
};

//! Test fixture for grompp
class GromppTest :
    public gmx::test::MdrunTestFixture
{
    public:
        void runGrompp()
        {
            runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
            {
                std::fstream groFile(runner_.groFileName_, groFile.trunc | groFile.out);
                groFile << atomicConfiguration;
            }
            EXPECT_EQ(0, runner_.callGrompp());
        }

        void runMD()
        {
            EXPECT_EQ(0, runner_.callMdrun());
        }
};

/* Initial target distance is `pull-coord1-init`. Rate is in nm / ps.
 *
 * Apply stochastic dynamics (Langevin?) at very low temperature. Pull from 0.6 to 10.6nm
 * distance over 10ps. This is a small system, but use a big box.
 */
std::string theMdpFile {
    R"""(cutoff-scheme = Verlet
integrator = sd
nsteps = 10000
tc-grps = System
tau-t = 0.1
ref-t = 10
pull = yes
pull-ngroups = 2
pull-ncoords = 1
pull-group1-name = SOL
pull-group2-name = Methanol
pull-coord1-groups = 1 2
pull-coord1-type = umbrella
pull-coord1-geometry = distance
pull-coord1-k = 1000
pull-coord1-init = 0.6
pull-coord1-rate = 1
)"""
};

TEST_F(GromppTest, MdpFileWorks)
{
    runner_.useStringAsMdpFile(theMdpFile);
    runGrompp();
}

TEST_F(GromppTest, RunPull)
{
    runner_.useStringAsMdpFile(theMdpFile);
    runGrompp();
    runMD();
}



} // end anonymous namespace
