/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Tests for a simple run using Plumed MDModule.
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include <limits>
#include <numeric>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/mdrun_main.h"
#include "programs/mdrun/tests/moduletest.h"
#include "programs/mdrun/tests/trajectorycomparison.h"

#include "plumedTestUtils.h"

namespace gmx
{
namespace test
{

namespace
{
class PlumedRun :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
protected:
    void SetUp() override
    {
        if (!plumedKernelAvailable())
        {
            const char* requirePlumedKernel = getenv("GMX_TEST_PLUMED_KERNEL_IS_AVAILABLE");
            if (requirePlumedKernel == nullptr)
            {
                GTEST_SKIP() << "The GMX_TEST_PLUMED_KERNEL_IS_AVAILABLE environment variable is "
                                "not set, and the PLUMED kernel was not found.";
            }
            GTEST_FAIL()
                    << "GMX_TEST_PLUMED_KERNEL_IS_AVAILABLE was set but PLUMED kernel not found";
        }
        gmxSetenv("PLUMED_MAXBACKUP", "0", 1);
    }

public:
    static constexpr std::array<const char*, 3> component_names_ = { "x", "y", "z" };
};

TEST_P(PlumedRun, PlumedSees)
{
    auto [simulationName, integrator] = GetParam();

    SCOPED_TRACE(
            formatString("Checking if we can run with plumed using "
                         "simulation '%s' with integrator '%s'",
                         simulationName.c_str(),
                         integrator.c_str()));

    auto mdpFieldValues =
            prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no");
    auto tmpdir = fileManager_.getOutputTempDirectory();
    // we'll need these later
    auto plumedxyz  = tmpdir / "plumed.xyz";
    auto cell_data  = tmpdir / "cell_data";
    auto time       = tmpdir / "time";
    auto plumed_dat = tmpdir / "plumed.dat";

    {
        std::string plumed_string;

        plumed_string += "DUMPATOMS PRECISION=10  ATOMS=@mdatoms FILE=" + plumedxyz.string() + "\n";
        plumed_string += "c: CELL\n";
        plumed_string += "PRINT ARG=c.* FILE=" + cell_data.string() + "\n";
        plumed_string += "t1: TIME\n";
        plumed_string += "PRINT ARG=t1 FILE=" + time.string() + "\n";
        gmx::TextWriter::writeFileFromString(plumed_dat, plumed_string);
    }
    // settings are stolen from simple_mdrun.cpp
    mdpFieldValues["nsteps"]        = "10";
    mdpFieldValues["nstxout"]       = "1";
    mdpFieldValues["nstfout"]       = "1";
    mdpFieldValues["constraints"]   = "none";
    mdpFieldValues["nstcalcenergy"] = "1";
    mdpFieldValues["coulombtype"]   = "Cut-off";
    mdpFieldValues["vdwtype"]       = "Cut-off";
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        ASSERT_EQ(0, runner_.callGrompp(caller));
    }

    // Do mdrun
    {
        CommandLine mdrunCaller;
        mdrunCaller.addOption("-reprod");
        mdrunCaller.addOption("-plumed", plumed_dat.string());
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }
    // Check if expected files are present
    for (const auto& file : { time, plumedxyz, cell_data })
    {
        EXPECT_TRUE(File::exists(file, File::returnFalseOnError)) << "File " << file << " was not found.";
    }
    {
        TextInputFile cell(cell_data);

        TextInputFile timestep(time);
        TextInputFile plumedTraj(plumedxyz);

        std::string plumedBoxStr;
        // values from trajectorycomparison.cpp
        FloatingPointTolerance tol_box = defaultRealTolerance();
        FloatingPointTolerance tol_pos = defaultRealTolerance();
        // removing comments
        cell.readLine(&plumedBoxStr);
        timestep.readLine(&plumedBoxStr);
        TrajectoryFrameReader reader(runner_.fullPrecisionTrajectoryFileName_);
        do
        {
            TrajectoryFrame frame = reader.frame();
            std::string     timeStr;
            { // check box
                cell.readLine(&plumedBoxStr);
                std::stringstream ss(plumedBoxStr);
                ss >> timeStr;

                for (int d = 0; d < DIM; ++d)
                {
                    for (int dd = 0; dd < DIM; ++dd)
                    {
                        real val;
                        ss >> val;
                        EXPECT_REAL_EQ_TOL(frame.box()[d][dd], val, tol_box);
                    }
                }
            }
            { // check time
                timestep.readLine(&timeStr);

                std::stringstream ss(timeStr);
                ss >> timeStr;
                double mytime;
                ss >> mytime;
                EXPECT_REAL_EQ_TOL(frame.time(), mytime, tol_box);
            }
            { // check positions (this actually is checking data tha has did the path gmx_prec->double->string->gmx_prec)
                std::string       line;
                std::stringstream atoms;
                auto              pos = frame.x();

                plumedTraj.readLine(&line);
                auto nat = std::stoi(line);
                // comments
                plumedTraj.readLine(&line);

                for (int i = 0; i < nat; ++i)
                {
                    plumedTraj.readLine(&line);
                    atoms.str(line);

                    std::string tmp;
                    real        x, y, z;
                    atoms >> tmp >> x >> y >> z;

                    EXPECT_REAL_EQ_TOL(pos[i][0], x, tol_pos);
                    EXPECT_REAL_EQ_TOL(pos[i][1], y, tol_pos);
                    EXPECT_REAL_EQ_TOL(pos[i][2], z, tol_pos);
                }
            }


        } while (reader.readNextFrame());
    }
}

TEST_P(PlumedRun, PlumedDoes)
{
    // some settings fast to change (or eventually to set up as test inputs)
    real constexpr force           = 10.0;
    int constexpr target_atom      = 2;
    int constexpr target_component = 0;

    auto [simulationName, integrator] = GetParam();

    auto mdpFieldValues =
            prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no");
    auto tmpdir = fileManager_.getOutputTempDirectory();
    // we'll need these later
    auto dist       = tmpdir / "dist";
    auto plumed_dat = tmpdir / "plumed_dd.dat";
    {
        std::string plumed_string;
        // using target atom+1 since in input plumed uses 1-based indexing
        plumed_string += "p: POSITION ATOM=" + std::to_string(target_atom + 1) + "\n";
        // I am applying a constant force of 10 on thex component of atom 3
        plumed_string += "res: RESTRAINT ARG=p." + std::string(component_names_[target_component])
                         + " KAPPA=0 SLOPE=" + std::to_string(force) + " AT=0\n";
        plumed_string += "PRINT ARG=p.*,res.* FILE=" + dist.string() + "\n";
        /*
        `plumed_string` should look like:
        ```plumed
        p: POSITION ATOM=3
        res: RESTRAINT ARG=p.x KAPPA=0 SLOPE=10.0 AT=0
        PRINT ARG=p.*,res.* FILE=/path/to/test/plumed_dd.dat
        ```
        */
        gmx::TextWriter::writeFileFromString(plumed_dat, plumed_string);
    }
    // settings are stolen from simple_mdrun.cpp
    mdpFieldValues["nsteps"]        = "10";
    mdpFieldValues["nstxout"]       = "1";
    mdpFieldValues["nstfout"]       = "1";
    mdpFieldValues["constraints"]   = "none";
    mdpFieldValues["nstcalcenergy"] = "1";
    mdpFieldValues["coulombtype"]   = "Cut-off";
    mdpFieldValues["vdwtype"]       = "Cut-off";

    auto baseTrajectoryFileName   = fileManager_.getTemporaryFilePath("plainSim.trr");
    auto baseEdrFileName          = fileManager_.getTemporaryFilePath("plainSim.edr");
    auto plumedTrajectoryFileName = fileManager_.getTemporaryFilePath("plumedSim.trr");
    auto plumedEdrFileName        = fileManager_.getTemporaryFilePath("plumedSim.edr");

    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        ASSERT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;

        runner_.fullPrecisionTrajectoryFileName_ = baseTrajectoryFileName.string();
        runner_.edrFileName_                     = baseEdrFileName.string();
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

        runner_.fullPrecisionTrajectoryFileName_ = plumedTrajectoryFileName.string();
        runner_.edrFileName_                     = plumedEdrFileName.string();
        mdrunCaller.addOption("-rerun", baseTrajectoryFileName.string());
        mdrunCaller.addOption("-plumed", plumed_dat.string());

        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    TrajectoryFrameReader  reader(baseTrajectoryFileName.string());
    TrajectoryFrameReader  readerPlumed(plumedTrajectoryFileName.string());
    FloatingPointTolerance tol = defaultRealTolerance();

    do
    {
        TrajectoryFrame frame       = reader.frame();
        TrajectoryFrame plumedframe = readerPlumed.frame();
        auto            nat         = frame.x().size();
        // only  a components of the force should have changed
        for (auto at = 0u; at < nat; at++)
        {
            for (auto i = 0; i < 3; i++)
            {
                if (i == target_component && at == target_atom)
                {
                    EXPECT_REAL_EQ_TOL(force, frame.f()[at][i] - plumedframe.f()[at][i], tol);
                }
                else
                {
                    EXPECT_REAL_EQ_TOL(frame.f()[at][i], plumedframe.f()[at][i], tol);
                }
            }
        }
        readerPlumed.readNextFrame();
    } while (reader.readNextFrame());
}

INSTANTIATE_TEST_SUITE_P(SimplePlumedMD,
                         PlumedRun,
                         ::testing::Combine(::testing::Values("argon12"), ::testing::Values("md")));


} // namespace
} // namespace test
} // namespace gmx
