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
 * Tests for the mdrun -rerun functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

//! Test fixture for mdrun -rerun
class MdrunRerun : public gmx::test::MdrunTestFixture,
                   public ::testing::WithParamInterface<const char *>
{
};

/* Among other things, this test ensures mdrun can read a trajectory. */
TEST_P(MdrunRerun, WithDifferentInputFormats)
{
    runner_.useEmptyMdpFile();
    runner_.useTopGroAndNdxFromDatabase("spc2");
    EXPECT_EQ(0, runner_.callGrompp());

    std::string rerunFileName = fileManager_.getInputFilePath(GetParam());

    ::gmx::test::CommandLine rerunCaller;
    rerunCaller.append("mdrun");
    rerunCaller.addOption("-rerun", rerunFileName);
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
}

/*! \brief Helper array of input files present in the source repo
 * database. These all have two identical frames of two SPC water
 * molecules, which were generated via trjconv from the .gro
 * version. */
const char *trajectoryFileNames[] = {
    "../../../gromacs/gmxana/legacytests/spc2-traj.trr",
#if defined GMX_USE_TNG && defined HAVE_ZLIB
    "../../../gromacs/gmxana/legacytests/spc2-traj.tng",
#endif
    "../../../gromacs/gmxana/legacytests/spc2-traj.xtc",
    "../../../gromacs/gmxana/legacytests/spc2-traj.gro",
    "../../../gromacs/gmxana/legacytests/spc2-traj.pdb",
    "../../../gromacs/gmxana/legacytests/spc2-traj.g96"
};
// TODO later. Find a better way to manage this file database and
// these string arrays that index it

#ifdef __INTEL_COMPILER
#pragma warning( disable : 177 )
#endif

INSTANTIATE_TEST_CASE_P(NoFatalErrorFrom,
                        MdrunRerun,
                            ::testing::ValuesIn(gmx::ArrayRef<const char*>(trajectoryFileNames)));

/*! \todo Add other tests for mdrun -rerun, e.g.
 *
 * - RerunReproducesRunWhenRunOnlyWroteEnergiesOnNeighborSearchSteps
 *   (e.g. do such a run, do a rerun, call gmxcheck)
 * - TpiExitsNormally (since it uses the -rerun machinery)
 */

} // namespace
