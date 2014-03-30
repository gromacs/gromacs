/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "moduletest.h"

#include <math.h>

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/fileio/path.h"
#include "gromacs/utility/stringutil.h"

#include "../mdrun_main.h"

#include "testutils/cmdlinetest.h"

namespace
{

/*! \brief
 * Test fixture for replica exchange
 */
class ReplicaExchangeTest : public gmx::test::ParameterizedMdrunTestFixture
{
    public:
        //! Constructor
        ReplicaExchangeTest() : size(gmx_node_num()),
                                rank(gmx_node_rank())
        {
            mdrunCaller.append("mdrun_mpi");
        }

        /*! \brief Organize the -multidir directories for the test.
         *
         * These are created inside the temporary directory for the
         * test case, and added to the eventual mdrun command
         * line. The temporary directory for the call to grompp by
         * this rank is set to the appropriate -multidir directory, so
         * the grompp output files go to the right place. */
        void organizeMultidir()
        {
            mdrunCaller.append("-multidir");

            std::string futureTestTempDirectory;
            for (int i = 0; i != size; ++i)
            {
                std::string replicaTempDirectory =
                    gmx::formatString("%s/multidir_%d",
                                      fileManager_.getOutputTempDirectory(), i);
                mdrunCaller.append(replicaTempDirectory);

                if (rank == i)
                {
                    gmx::Directory::create(replicaTempDirectory);
                    futureTestTempDirectory = std::string(replicaTempDirectory);
                }
            }
            fileManager_.setOutputTempDirectory(futureTestTempDirectory);

            /* Prepare grompp output filenames inside the new
               temporary directory */
            mdpInputFileName  = fileManager_.getTemporaryFilePath("input.mdp");
            mdpOutputFileName = fileManager_.getTemporaryFilePath("output.mdp");
            tprFileName       = fileManager_.getTemporaryFilePath(".tpr");

            mdrunCaller.addOption("-deffnm", fileManager_.getTestSpecificFileNameRoot());
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
                    baseTemperature + 0.0001*rank,
                    basePressure * pow(1.01, rank),
                    /* Set things up so that the initial KE decreases
                       with increasing replica number, so that the
                       (identical) starting PE decreases on the first
                       step more for the replicas with higher number,
                       which will tend to force replica exchange to
                       occur. */
                    std::max(baseTemperature - 10 * rank, real(0)),
                    controlVariable);
            useStringAsMdpFile(mdpFileContents);
        }

        //! MPI process set size
        int                    size;
        //! MPI rank of this process
        int                    rank;
        //! Object for building the mdrun command line
        gmx::test::CommandLine mdrunCaller;
};

/* This test ensures mdrun can run NVT REMD under the supported
 * conditions.
 *
 * TODO Preferably, we could test that mdrun correctly refuses to run
 * replica exchange unless under real MPI with more than one rank
 * available. However, those cases trigger an error that is currently
 * fatal to mdrun and also to the test binary. So, in the meantime we
 * must not test those cases. This is done via disabling the test, so
 * that there is a reminder that it is disabled. There's no elegant
 * way to conditionally disable a test. */
TEST_P(ReplicaExchangeTest, ExitsNormally)
{
    if (size <= 1)
    {
        /* Can't test replica exchange without multiple ranks. */
        return;
    }

    organizeMultidir();
    const char *pcoupl = GetParam();
    organizeMdpFile(pcoupl);
    useTopGroAndNdxFromDatabase("spc2");
    /* Call grompp on every rank - the standard callGrompp() only runs
       grompp on rank 0. */
    EXPECT_EQ(0, callGromppOnThisRank());

    mdrunCaller.addOption("-replex", 1);
    ASSERT_EQ(0, gmx_mdrun(mdrunCaller.argc(), mdrunCaller.argv()));
}

#ifdef GMX_LIB_MPI
INSTANTIATE_TEST_CASE_P(WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#endif

} // namespace
