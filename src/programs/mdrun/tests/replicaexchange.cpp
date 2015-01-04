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
 * Tests for the mdrun replica exchange functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "config.h"

#include <math.h>

#include <algorithm>
#include <iostream>

#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "programs/mdrun/mdrun_main.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

/*! \brief
 * Test fixture for replica exchange
 */
class ReplicaExchangeTest : public gmx::test::ParameterizedMdrunTestFixture
{
    public:
        //! Constructor
        ReplicaExchangeTest() : size_(gmx_node_num()),
                                rank_(gmx_node_rank())
        {
            mdrunCaller_.append("mdrun_mpi");
        }

        /*! \brief Organize the .mdp file for this rank
         *
         * \param controlVariable Allows parameterization to work with
         * T, P or (later) lamda as the control variable, by passing a
         * string with "mdp-param = value" such that different paths
         * in init_replica_exchange() are followed.
         */
        void organizeMdpFile(const char *controlVariable)
        {
            const real  baseTemperature = 298;
            const real  basePressure    = 1;
            std::string mdpFileContents =
                gmx::formatString
                    ("nsteps = 2\n"
                    "nstlog = 1\n"
                    "nstcalcenergy = 1\n"
                    "tcoupl = v-rescale\n"
                    "tc-grps = System\n"
                    "tau-t = 1\n"
                    "ref-t = %f\n"
                    // pressure coupling (if active)
                    "tau-p = 1\n"
                    "ref-p = %f\n"
                    "compressibility = 4.5e-5\n"
                    // velocity generation
                    "gen-vel = yes\n"
                    "gen-temp = %f\n"
                    // control variable specification
                    "%s\n",
                    baseTemperature + 0.0001*rank_,
                    basePressure * pow(1.01, rank_),
                    /* Set things up so that the initial KE decreases
                       with increasing replica number, so that the
                       (identical) starting PE decreases on the first
                       step more for the replicas with higher number,
                       which will tend to force replica exchange to
                       occur. */
                    std::max(baseTemperature - 10 * rank_, real(0)),
                    controlVariable);
            runner_.useStringAsMdpFile(mdpFileContents);
        }

        //! Number of MPI ranks
        int                    size_;
        //! MPI rank of this process
        int                    rank_;
        //! Object for building the mdrun command line
        gmx::test::CommandLine mdrunCaller_;
};

/* This test ensures mdrun can run NVT REMD under the supported
 * conditions. It runs one replica per MPI rank.
 *
 * TODO Preferably, we could test that mdrun correctly refuses to run
 * replica exchange unless under real MPI with more than one rank
 * available. However, if we just call mdrun blindly, those cases
 * trigger an error that is currently fatal to mdrun and also to the
 * test binary. So, in the meantime we must not test those cases. If
 * there is no MPI, we disable the test, so that there is a reminder
 * that it is disabled. There's no elegant way to conditionally
 * disable a test at run time, so currently there is no feedback if
 * only one rank is available. However, the test harness knows
 * to run this test with more than one rank. */
TEST_P(ReplicaExchangeTest, ExitsNormally)
{
    if (size_ <= 1)
    {
        /* Can't test replica exchange without multiple ranks. */
        return;
    }

    runner_.mdpInputFileName_  = fileManager_.getTemporaryFilePath(gmx::formatString("input%d.mdp", rank_));
    runner_.mdpOutputFileName_ = fileManager_.getTemporaryFilePath(gmx::formatString("output%d.mdp", rank_));

    /* grompp needs to name the .tpr file so that when mdrun appends
       the MPI rank, it will find the right file. If we just used
       "%d.tpr" then \c TestFileManager prefixes that with an
       underscore. Then, there is no way for mdrun to be told the
       right name, because if you add the underscore manually, you get
       a second one from \c TestFileManager. However, it's easy to
       just start the suffix with "topol" in both cases. */
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath(gmx::formatString("topol%d.tpr", rank_));

    const char *pcoupl = GetParam();
    organizeMdpFile(pcoupl);
    runner_.useTopGroAndNdxFromDatabase("spc2");
    /* Call grompp on every rank - the standard callGrompp() only runs
       grompp on rank 0. */
    EXPECT_EQ(0, runner_.callGromppOnThisRank());

    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("topol.tpr");
    mdrunCaller_.addOption("-multi", size_);
    mdrunCaller_.addOption("-replex", 1);
    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller_));
}

#ifdef GMX_LIB_MPI
INSTANTIATE_TEST_CASE_P(WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#endif

} // namespace
