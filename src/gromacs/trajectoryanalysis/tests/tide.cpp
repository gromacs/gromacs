/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
 * Tests for functionality of the "tide" time-averaged density trajectory analysis module.
 *
 * \author Camilo Aponte <ca.aponte.uniandes.edu.co>
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/tide.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::ExactTextMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Tide.
 */

//! Test fixture for the angle analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::TimeAveragedDensityInfo>
    DensityModuleTest;

TEST_F(DensityModuleTest, SingleLipidYieldsCorrectDensity)
{
    const char *const cmdline[] = {
        "tide",
        "-size", "40",  "40",  "80",
        "-min",  "79",  "80",  "61",
        "-max", "119", "120", "141",
        "-box",   "2",   "2",   "4",
        "-cutoff", "2.5",
        "-Ngaussians", "2",
        "-select", "resname DPPC"
    };

    setTopology  ("DppcLipid.gro");
    setTrajectory("DppcLipid.gro");
    setInputFile ("-Gauss_coef", "DppcGaussians.dat");
    setOutputFile("-map_av", "av.xplor", ExactTextMatch());

    runTest(CommandLine(cmdline));
}

TEST_F(DensityModuleTest, TrajectoryYieldsCorrectAverages)
{
    const char *const cmdline[] = {
        "tide",
        "-size", "40", "40",  "80",
        "-min",  "95", "47", "-22",
        "-max", "135", "87",  "58",
        "-box",   "2",  "2",   "4",
        "-cutoff", "0.4",
        "-Ngaussians", "5",
        "-select", "resname DMP",
        "-nonorm",
        "-time_window", "10",
        "-corr"
    };

    setTopology  ("DmpcLipid.gro");
    setTrajectory("DmpcLipid.xtc");
    setInputFile ("-Gauss_coef", "DmpcGaussians.dat");
    setInputFile ("-mi", "av_0-10ps.xplor");
    setOutputFile("-map_av", "av.xplor", ExactTextMatch());
    setOutputFile("-map_sd", "sigma.xplor", ExactTextMatch());
    setOutputFile("-rho", "rho.dat", ExactTextMatch());

    runTest(CommandLine(cmdline));
}
} // namespace
