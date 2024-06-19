/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Test for Ewald 3DC and epsilon-surface terms.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstdio>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/trajectoryreader.h"

#include "energycomparison.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! Helper type
using MdpField = MdpFieldValues::value_type;

/*! \brief Test fixture base for simple mdrun systems
 *
 * This test ensures mdrun can run a simulation, reaching reproducible
 * energies and forces.
 * The starting coordinates are set up such that both molecules are
 * broken over PBC and the PBC treatment is tested.
 */
class EwaldSurfaceTermTest : public MdrunTestFixture, public ::testing::WithParamInterface<std::string>
{
};

TEST_P(EwaldSurfaceTermTest, WithinTolerances)
{
    auto simulationName = GetParam();
    SCOPED_TRACE(formatString("Comparing simple mdrun for '%s'", simulationName.c_str()));

    int numRanksAvailable = getNumberOfTestMpiRanks();
    /* For epsilon-surface we need whole molecules.
     * Without constraints we can make molecules whole on a sinlge rank.
     * With constraints molecules are whole with update groups with DD.
     *
     * TODO: Remove the rank=1 check when DD also runs single rank.
     */
    if ((simulationName == "epsilon-surface" && numRanksAvailable > 1)
        || (simulationName == "epsilon-surface-constraint" && numRanksAvailable == 1))
    {
        fprintf(stdout,
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are %s 1.\n",
                simulationName.c_str(),
                numRanksAvailable,
                numRanksAvailable == 1 ? ">" : "=");
        return;
    }

    std::string theMdpFile =
            "coulombtype     = PME\n"
            "nstcalcenergy   = 1\n"
            "nstenergy       = 4\n"
            "nstfout         = 4\n"
            "rcoulomb        = 1.5\n"
            "rvdw            = 1.5\n"
            "pme-order       = 4\n"
            "fourier-spacing = 0.2\n"
            "dt              = 0.0025\n"
            "nsteps          = 20\n";

    if (simulationName == "3DC")
    {
        theMdpFile +=
                "ewald-geometry = 3DC\n"
                "pbc             = xy\n"
                "nwall           = 2\n"
                "wall-type       = 12-6\n"
                "wall-atomtype   = C C\n"
                "wall-ewald-zfac = 2\n";
    }
    else
    {
        theMdpFile += "epsilon-surface = 1\n";
        if (simulationName == "epsilon-surface-constraint")
        {
            theMdpFile += "constraints = all-bonds";
        }
    }

    // Prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase("dipoles");
        runner_.useStringAsMdpFile(theMdpFile);
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
        EnergyTermsToCompare energyTermsToCompare{
            { { interaction_function[F_EPOT].longname, absoluteTolerance(1e-3) },
              { interaction_function[F_ETOT].longname, absoluteTolerance(1e-3) } }
        };
        TestReferenceData refData;
        auto checker = refData.rootChecker().checkCompound("Simulation", simulationName);
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
        // Now check the forces
        TrajectoryFrameReader reader(runner_.fullPrecisionTrajectoryFileName_);
        checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1, 1e-3));
        do
        {
            auto frame = reader.frame();
            auto force = frame.f();
            int  atom  = 0;
            for (const auto& f : force)
            {
                std::string forceName = frame.frameName() + " F[" + toString(atom) + "]";

                checker.checkVector(f, forceName.c_str());
                atom++;
            }
        } while (reader.readNextFrame());
    }
}

//! Containers of systems to test.
//! \{
std::vector<std::string> surfaceTerm = { "3DC", "epsilon-surface-constraint", "epsilon-surface" };
//! \}

INSTANTIATE_TEST_SUITE_P(EwaldSurfaceTerm, EwaldSurfaceTermTest, ::testing::ValuesIn(surfaceTerm));

} // namespace
} // namespace test
} // namespace gmx
