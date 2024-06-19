/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for mdrun enforced rotation (rotational pulling) functionality.
 *
 * These tests check the mdrun results for cases where one of the rotation
 * potentials (see the corresponding section in the manual) is turned on.
 * For a small rotation group of 4 atoms with the following atomic positions
 * (.gro format):
 *
 *     4
 *     1AR      AR    1   1.500   1.500   3.200
 *     2AR      AR    2   3.900   5.800   3.500
 *     3AR      AR    3   4.700   6.200   2.100
 *     4AR      AR    4   5.500   4.000   6.600
 *     8.000 8.000 8.000
 *
 * and these positions of the reference group:
 *
 *     4
 *     1AR      AR    1   1.000   2.000   3.000
 *     2AR      AR    2   4.000   5.500   5.000
 *     3AR      AR    3   4.500   6.000   2.500
 *     4AR      AR    4   5.000   4.500   7.000
 *     8.000 8.000 8.000
 *
 * the tests for the iso, iso-pf, pm, pm-pf, rm, rm-pf, rm2, rm2-pf, flex, flex-t,
 * flex2 and flex2-t potential types check the energies and forces from each
 * rotational potential for a 5 degree rotation against the results obtained
 * from a Mathematica notebook. As the cut-offs are small, the forces calculated
 * at step zero are from the rotational potential only. The forces in the first
 * frame of the .trr output files should therefore be identical to the Mathematica
 * results.
 *
 * In addition to the test described above, which checks energies and forces
 * at time step 0, a second test compares a short 25-step trajectory to reference
 * data with slightly looser tolerances.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstddef>

#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/tests/comparison_helpers.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"
#include "trajectorycomparison.h"

namespace gmx
{
namespace test
{
namespace
{

//! Convenience type to provide test parameters as a single argument
// 1st argument: mdp file entry for rotation potential type, e.g. 'iso', 'iso-pf', ...
// 2nd argument: expected value for potential energy E_rot at first step
// 3rd argument: expected value for torque tau at first step
typedef std::tuple<std::string, double, double> Parameters;

//! Reference forces per potential resulting from a 5 degree rotation
std::map<std::string, std::vector<std::vector<double>>> referenceForces = {
    { { "iso",
        { { -131.51692692, 127.41828586, -54.439881601 },
          { 90.187124112, -258.41723581, 1208.8824492 },
          { -410.64765335, -162.01496986, 484.89253102 },
          { -522.65423367, 810.87229117, 540.30321711 } } },
      { "iso-pf",
        { { -34.053757936, 75.632448726, -272.40371317 },
          { 285.11346208, -361.98891008, 772.95478603 },
          { -118.2581464, -317.37248127, -168.99896369 },
          { -132.80155774, 603.72894263, -331.55210917 } } },
      { "pm",
        { { -128.65978406, 133.13257158, -45.868453029 },
          { -138.38430446, -715.56009296, 523.16816346 },
          { -462.07622478, -264.87211272, 330.60681674 },
          { -716.93994796, 422.3008626, -42.553925743 } } },
      { "pm-pf",
        { { 15.946242064, 175.63244873, -122.40371317 },
          { 150.82774779, -630.56033866, 370.09764317 },
          { -28.2581464, -137.37248127, 101.00103631 },
          { -138.51584346, 592.3003712, -348.69496631 } } },
      { "rm",
        { { -132.43053309, 130.32663443, -42.740911919 },
          { 335.27526601, -551.67578239, 256.02543292 },
          { -315.35595461, 130.96336469, 17.809741744 },
          { -577.37525264, -45.713857761, 222.93432272 } } },
      { "rm-pf",
        { { -60.682586583, 131.62304711, -67.521169211 },
          { -214.87938414, -191.7299335, 199.44641705 },
          { 97.203055624, 23.46213548, -48.042442195 },
          { 178.3589151, 36.644750908, -83.88280564 } } },
      { "rm2",
        { { -143.73628978, 108.23867091, -24.247017347 },
          { 31.674798399, -2.5514862071, -8.8572753284 },
          { -189.23935665, 102.24209625, -5.0816119495 },
          { -387.59084854, -89.705802628, 189.00081793 } } },
      { "rm2-pf",
        { { -72.570497152, 79.439407569, -28.769439329 },
          { -30.440888716, -27.798490561, 28.679289946 },
          { 31.555426196, -23.378441432, 5.0671522227 },
          { 71.455959672, -28.262475577, -4.9770028397 } } },
      { "flex",
        { { -42.305472698, 25.397494447, -22.79279844 },
          { 243.49534318, -260.46845668, 99.087010409 },
          { -71.573582075, 200.8909235, 79.762542508 },
          { 12.712812823, 318.83824119, 270.93054921 } } },
      { "flex-t",
        { { -36.90697396, 27.986236644, -30.97153705 },
          { 127.85407917, -146.63332095, 140.35716757 },
          { -92.458053794, 120.04725225, -7.5156601118 },
          { 73.958305837, 143.49454656, 115.47210135 } } },
      { "flex2",
        { { -41.50124682, 8.0976554981, -28.190729908 },
          { 182.26944646, -173.94075441, 112.45880245 },
          { -125.42736594, 164.47230806, -68.461771568 },
          { 7.8675958351, 47.787649935, 53.818987649 } } },
      { "flex2-t",
        { { -50.140956713, 6.6349690252, -40.465489987 },
          { 156.89909499, -170.72723015, 140.66803971 },
          { -127.41036756, 173.38786648, -52.41251982 },
          { 48.021017959, 45.441971984, 34.316336103 } } } }
};


//! Compare the forces at the first step to the reference
void checkRotForcesAtStepZero(const std::string& fn, const std::vector<std::vector<double>>& reference)
{
    auto reader = TrajectoryFrameReader(fn);
    auto frame  = reader.frame();
    auto f      = frame.f();

    // Loop over all 4 atoms of the system
    for (size_t i = 0; i < 4; i++)
    {
        // Loop over x, y, and z entry of the forces
        for (size_t j = 0; j < 3; j++)
        {
            EXPECT_REAL_EQ_TOL(f[i][j],
                               reference[i][j],
                               absoluteTolerance(std::is_same_v<real, double> ? 1e-7 : 1e-3));
        }
    }
}

//! Return the first entry for the rotational energy stored in the .edr file
real getFirstRotEnergyValue(const std::string& fn)
{
    auto E   = 0.0;
    auto efr = openEnergyFileToReadTerms(fn, { interaction_function[F_COM_PULL].longname });
    if (efr->readNextFrame())
    {
        auto fr = efr->frame();
        E       = fr.at(interaction_function[F_COM_PULL].longname);
    }
    else
    {
        gmx_fatal(FARGS, "Cannot read first energy frame of file %s\n", fn.c_str());
    }
    return E;
}

//! Parameterized test fixture for enforced rotation
class RotationTest : public MdrunTestFixture, public ::testing::WithParamInterface<Parameters>
{
};

//! Tests to see if the rotation potentials are giving the correct energies and forces
//
// All tests are for a four atom argon system.
INSTANTIATE_TEST_SUITE_P(RotationWorks,
                         RotationTest,
                         ::testing::Values(Parameters("iso", 1567.0459109, -1033.7307213),
                                           Parameters("iso-pf", 820.93248453, -123.87809352),
                                           Parameters("pm", 929.18876803, -1033.7307213),
                                           Parameters("pm-pf", 572.0039131, -123.87809352),
                                           Parameters("rm", 515.95322462, -1595.3684936),
                                           Parameters("rm-pf", 143.15040559, 165.78560092),
                                           Parameters("rm2", 196.69333887, -1189.9052573),
                                           Parameters("rm2-pf", 42.345322053, -232.66325156),
                                           Parameters("flex", 255.58915604, -111.76617807),
                                           Parameters("flex-t", 114.13605288, 186.27588869),
                                           Parameters("flex2", 170.6496809, 142.13808477),
                                           Parameters("flex2-t", 101.11490074, 270.30162398)));


TEST_P(RotationTest, CheckEnergiesForcesAndTraj)
{
    Parameters testPar = GetParam();

    // Unpack:
    auto rotTypeString  = std::get<0>(testPar);
    auto expectedEnergy = std::get<1>(testPar);
    auto expectedTorque = std::get<2>(testPar);

    SCOPED_TRACE(formatString("Checking enforced rotation for potential type '%s'", rotTypeString.c_str()));

    auto mdpStub = formatString(
            "integrator     = md          \n"
            "tinit          = 0.002       \n"
            "dt             = 0.002       \n"
            "nsteps         = 25          \n"
            "coulombtype    = Cut-off     \n"
            "nstxout        = 4           \n"
            "nstfout        = 4           \n"
            "nstcalcenergy  = 5           \n"
            "nstenergy      = 5           \n"
            "nstlist        = 5           \n"
            "rlist          = 0.4         \n"
            "rvdw           = 0.4         \n"
            "rcoulomb       = 0.4         \n"
            "Tcoupl         = no          \n"
            "Pcoupl         = no          \n"
            "gen-vel        = no          \n"
            "rotation       = yes         \n"
            "rot_nstrout    = 1           \n"
            "rot-nstsout    = 1           \n"
            "rot-group0     = system      \n"
            "rot-massw0     = yes         \n"
            "rot-vec0       = 1.0 2.0 3.0 \n"
            "rot-pivot0     = 4 5 4       \n"
            "rot-rate0      = 2500        \n"
            "rot-k0         = 1000        \n"
            "rot-min-gauss0 = 1.0e-12     \n"
            "rot-eps0       = 0.01        \n");

    // Prepare the simulation input .tpr file
    {
        runner_.useStringAsMdpFile(mdpStub + "rot_type0 = " + rotTypeString);
        runner_.useTopGroAndNdxFromDatabase("argon4");
        auto gromppCaller = CommandLine();
        auto rotrefFileName = std::filesystem::path(TestFileManager::getTestSimulationDatabaseDirectory())
                                      .append("rotref.trr");
        gromppCaller.addOption("-ref", rotrefFileName.string());
        ASSERT_EQ(0, runner_.callGrompp(gromppCaller));
    }

    auto fn_xvg = fileManager_.getTemporaryFilePath("rotation.xvg");

    // Do a short MD simulation
    {
        auto mdrunCaller = CommandLine();
        mdrunCaller.addOption("-ro", fn_xvg.string());
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }


    /******************************************************/
    /* Check computed values at step 0 against references */
    /* obtained from a Mathematica notebook               */
    /******************************************************/

    // Check energies
    {
        auto E_rot = getFirstRotEnergyValue(runner_.edrFileName_);

        EXPECT_REAL_EQ_TOL(
                expectedEnergy, E_rot, absoluteTolerance(std::is_same_v<real, double> ? 1e-8 : 0.001));
    }

    // Check forces
    checkRotForcesAtStepZero(runner_.fullPrecisionTrajectoryFileName_, referenceForces[rotTypeString]);

    // Check the torques and other values written to .xvg file
    {
        auto        timeSeriesData = readXvgTimeSeries(fn_xvg, 0.0, 1.0);
        const auto& step0Data      = timeSeriesData.asConstView()[0];
        EXPECT_EQ(step0Data[0], 0.002); // time (fs)
        EXPECT_EQ(step0Data[1], 5.000); // theta_ref (degrees)
        // As we only write out a few digits to the .xvg file, we cannot use the same tight
        // tolerances we get from .edr
        EXPECT_REAL_EQ_TOL(
                expectedTorque, step0Data[3], relativeToleranceAsFloatingPoint(expectedEnergy, 1e-3));
        // Next test is somewhat redundant, but it does check whether the correct values actually end up in the .xvg file
        EXPECT_REAL_EQ_TOL(
                expectedEnergy, step0Data[4], relativeToleranceAsFloatingPoint(expectedEnergy, 1e-3));
    }


    /******************************************************/
    /* Now compare computed values to a short reference   */
    /* trajectory                                         */
    /******************************************************/
    TestReferenceData refData;
    auto              checker = refData.rootChecker();

    // Compare energies for .edr time series
    {
        auto energyTolerance = absoluteTolerance(std::is_same_v<real, double> ? 1e-8 : 0.01);

        EnergyTermsToCompare energyTermsToCompare{
            { { interaction_function[F_COM_PULL].longname, energyTolerance },
              { interaction_function[F_EPOT].longname, energyTolerance } }
        };
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }

    // Compare positions and forces for .trr time series
    {
        // Specify how trajectory frame matching must work
        const TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,  // box
                                                                    false, // no need to handle PBC
                                                                    false, // no need to handle PBC
                                                                    ComparisonConditions::MustCompare, // x
                                                                    ComparisonConditions::NoComparison, // v
                                                                    ComparisonConditions::MustCompare, // f
                                                                    MaxNumFrames::compareAllFrames() };
        TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
        TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

        checkTrajectoryAgainstReferenceData(
                runner_.fullPrecisionTrajectoryFileName_, trajectoryComparison, &checker);
    }
}
} // namespace
} // namespace test
} // namespace gmx
