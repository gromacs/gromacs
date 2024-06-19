/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Tests that mdrun restarts are exact, that is that two
 * successive runs without appending reproduce a single-part run.
 *
 * \todo Extend the coverage to the appending case.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "programs/mdrun/tests/comparison_helpers.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{


/*! \brief Manages comparing a pair of matching energy frames from a
 * normal run and the same run in two parts.
 *
 * The passed frame readers must contain logically equivalent
 * contents, with no extra or missing frames from either the
 * one- or two-part run.
 *
 * \todo This class has a similar interface with one in
 * mdcomparison.h, but not enough to warrant extracting an interface
 * class. Perhaps parts of this should be cast as a class that returns
 * iterators to successive frames, so that the fancy footwork for
 * pretending a two-part run is one sequence is separate from the
 * comparison of that sequene with that of the one-part run?
 *
 * \tparam FrameReader  Has readNextFrame() and frame() methods
 *                      useful for returning successive Frame objects */
template<class FrameReader>
class ContinuationFramePairManager
{
public:
    //! Convenience typedef
    typedef std::unique_ptr<FrameReader> FrameReaderPtr;
    //! Constructor
    ContinuationFramePairManager(FrameReaderPtr full, FrameReaderPtr firstPart, FrameReaderPtr secondPart) :
        full_(std::move(full)),
        firstPart_(std::move(firstPart)),
        secondPart_(std::move(secondPart)),
        isFirstPart_(true)
    {
    }
    /*! \brief Probe for a pair of valid frames, and return true if both are found.
     *
     * Gives a test failure if exactly one frame is found, because
     * it is an error for either run to have missing or extra
     * frames.  Note that the frame where the two-part run ends
     * and begins is duplicated between the two output files by
     * mdrun, and the test accommodates this.
     *
     * \todo This would be straightforward if velocity Verlet
     * behaved like other integrators. */
    bool shouldContinueComparing()
    {
        if (full_->readNextFrame())
        {
            if (isFirstPart_)
            {
                if (firstPart_->readNextFrame())
                {
                    // Two valid next frames exist, so we should continue comparing.
                    return true;
                }
                else
                {
                    // First part ran out of frames, move on to the second part
                    isFirstPart_ = false;
                    if (secondPart_->readNextFrame())
                    {
                        // Skip a second-part frame so the one we will
                        // read can compare with the next full-run
                        // frames.
                        secondPart_->frame();
                        if (secondPart_->readNextFrame())
                        {
                            // Two valid next frames exist, so we should continue comparing.
                            return true;
                        }
                        else
                        {
                            ADD_FAILURE() << "Second-part energy file had no (new) frames";
                        }
                    }
                    else
                    {
                        ADD_FAILURE() << "Second-part energy file had no frames";
                    }
                }
            }
            else
            {
                if (secondPart_->readNextFrame())
                {
                    // Two valid next frames exist, so we should continue comparing.
                    return true;
                }
                else
                {
                    ADD_FAILURE() << "Full run energy file had at least one more frame than "
                                     "two-part run energy file";
                }
            }
        }
        else
        {
            if (isFirstPart_)
            {
                ADD_FAILURE() << "Full-run energy file ran out of frames before the first part of "
                                 "the two-part run completed";
            }
            else
            {
                if (secondPart_->readNextFrame())
                {
                    ADD_FAILURE() << "Two-part run energy file had at least one more frame than "
                                     "full-run energy file";
                }
                else
                {
                    // Both files ran out of frames at the same time, which is the expected behaviour.
                }
            }
        }
        // At least one file is out of frames, so should not continue comparing.
        return false;
    }
    /*! \brief Compare all possible pairs of frames using \c compareTwoFrames.
     *
     * \tparam Frame  The type of frame used in the comparison (returned
     *                by FrameReader and used by compareTwoFrames). */
    template<class Frame>
    void compareAllFramePairs(std::function<void(const Frame&, const Frame&)> compareTwoFrames)
    {
        while (shouldContinueComparing())
        {
            EnergyFrame firstFrame  = full_->frame();
            EnergyFrame secondFrame = isFirstPart_ ? firstPart_->frame() : secondPart_->frame();
            SCOPED_TRACE("Comparing frames from two runs '" + firstFrame.frameName() + "' and '"
                         + secondFrame.frameName() + "'");
            compareTwoFrames(firstFrame, secondFrame);
        }
    }

private:
    EnergyFrameReaderPtr full_;
    EnergyFrameReaderPtr firstPart_;
    EnergyFrameReaderPtr secondPart_;
    bool                 isFirstPart_;
};

/*! \brief Run grompp for a normal mdrun, the same mdrun stopping part
 * way, doing a continuation, and compare the results. */
void runTest(TestFileManager*            fileManager,
             SimulationRunner*           runner,
             const std::string&          simulationName,
             const MdpFieldValues&       mdpFieldValues,
             const bool                  requestReproducibility,
             const EnergyTermsToCompare& energyTermsToCompare)
{
    int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout,
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(),
                numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    // prepare some names for files to use with the two mdrun calls
    std::string fullRunTprFileName      = fileManager->getTemporaryFilePath("full.tpr").string();
    std::string firstPartRunTprFileName = fileManager->getTemporaryFilePath("firstpart.tpr").string();
    std::string fullRunEdrFileName      = fileManager->getTemporaryFilePath("full.edr").string();
    std::string firstPartRunEdrFileName = fileManager->getTemporaryFilePath("firstpart.edr").string();
    std::string firstPartRunCheckpointFileName =
            fileManager->getTemporaryFilePath("firstpart.cpt").string();
    std::string secondPartRunEdrFileName = fileManager->getTemporaryFilePath("secondpart").string();

    // prepare the full run .tpr file, which will be used for the full
    // run, and for the second part of the two-part run.
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        runner->useTopGroAndNdxFromDatabase(simulationName);
        runner->useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        runner->tprFileName_ = fullRunTprFileName;
        EXPECT_EQ(0, runner->callGrompp(caller));
    }

    const std::string splitPoint = std::to_string(std::stoi(mdpFieldValues.at("nsteps")) / 2);

    // prepare the .tpr file for the first part of the two-part run
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        runner->useTopGroAndNdxFromDatabase(simulationName);
        auto firstPartMdpFieldValues      = mdpFieldValues;
        firstPartMdpFieldValues["nsteps"] = splitPoint;
        runner->useStringAsMdpFile(prepareMdpFileContents(firstPartMdpFieldValues));
        runner->tprFileName_ = firstPartRunTprFileName;
        EXPECT_EQ(0, runner->callGrompp(caller));
    }

    // do a normal mdrun
    {
        runner->tprFileName_ = fullRunTprFileName;
        runner->edrFileName_ = fullRunEdrFileName;
        CommandLine fullRunCaller;
        fullRunCaller.append("mdrun");
        if (requestReproducibility)
        {
            fullRunCaller.append("-reprod");
            /* Force neighborlist update at the beginning of the second half of the trajectory.
             * Doing so through CLI options prevents pairlist tuning from changing it. */
            fullRunCaller.addOption("-nstlist", splitPoint);
        }
        ASSERT_EQ(0, runner->callMdrun(fullRunCaller));
    }

    // do a repeat of the first part of the same mdrun
    {
        runner->tprFileName_ = firstPartRunTprFileName;
        runner->edrFileName_ = firstPartRunEdrFileName;
        CommandLine firstPartRunCaller;
        firstPartRunCaller.append("mdrun");
        if (requestReproducibility)
        {
            firstPartRunCaller.append("-reprod");
            firstPartRunCaller.addOption("-nstlist", splitPoint);
        }
        runner->cptOutputFileName_ = firstPartRunCheckpointFileName;
        ASSERT_EQ(0, runner->callMdrun(firstPartRunCaller));
    }

    // do a continuation (without appending) from the first part of
    // that same mdrun
    {
        runner->tprFileName_ = fullRunTprFileName;
        runner->edrFileName_ = secondPartRunEdrFileName;
        CommandLine secondPartRunCaller;
        secondPartRunCaller.append("mdrun");
        if (requestReproducibility)
        {
            secondPartRunCaller.append("-reprod");
            secondPartRunCaller.addOption("-nstlist", splitPoint);
        }
        // TODO We could test with appending but it would need a
        // different implementation.
        secondPartRunCaller.append("-noappend");
        secondPartRunCaller.addOption("-cpi", firstPartRunCheckpointFileName);
        ASSERT_EQ(0, runner->callMdrun(secondPartRunCaller));
        // Cope with how -noappend works
        secondPartRunEdrFileName += ".part0002.edr";
    }

    // Build the functor that will compare energy frames on the chosen
    // energy terms.
    EnergyComparison energyComparison(energyTermsToCompare, MaxNumFrames::compareAllFrames());

    // Build the manager that will present matching pairs of frames to compare.
    //
    // TODO Here is an unnecessary copy of keys (ie. the energy term
    // names), for convenience. In the future, use a range.
    auto namesOfEnergiesToMatch = energyComparison.getEnergyNames();
    ContinuationFramePairManager<EnergyFrameReader> energyManager(
            openEnergyFileToReadTerms(fullRunEdrFileName, namesOfEnergiesToMatch),
            openEnergyFileToReadTerms(firstPartRunEdrFileName, namesOfEnergiesToMatch),
            openEnergyFileToReadTerms(secondPartRunEdrFileName, namesOfEnergiesToMatch));
    // Compare the energy frames.
    energyManager.compareAllFramePairs<EnergyFrame>(energyComparison);
}

/*! \brief Test fixture for mdrun exact continuations
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce to within a tolerance the same
 * energies from runs that stopped part of the way, and restarted from
 * the checkpoint.
 *
 * \todo Is there value in testing with mdrun -reprod? As well as
 * without it?
 *
 * \todo Add FEP case. */
class MdrunNoAppendContinuationIsExact :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string, std::string, std::string, MdpParameterDatabase, bool>>
{
public:
    //! Constructor
    MdrunNoAppendContinuationIsExact() {}
};

/* This test setup runs a simulation and then runs the same
 * simulation as two half simulations with continuation from checkpoint.
 * The test can be ran in binary-reproducible mode which should result in
 * binary exact continuation. It can also be run in non binary-reproducible
 * mode, then some tolerances are needed, but we can then check continuation
 * with GPU and at non pair-search steps.
 */
TEST_P(MdrunNoAppendContinuationIsExact, WithinTolerances)
{
    auto params                  = GetParam();
    auto simulationName          = std::get<0>(params);
    auto integrator              = std::get<1>(params);
    auto temperatureCoupling     = std::get<2>(params);
    auto pressureCoupling        = std::get<3>(params);
    auto additionalMdpParameters = std::get<4>(params);
    auto binaryReproducible      = std::get<5>(params);

    // Check for unimplemented functionality
    // TODO: Update this as modular simulator gains functionality
    const bool isModularSimulatorExplicitlyDisabled = (getenv("GMX_DISABLE_MODULAR_SIMULATOR") != nullptr);
    const bool isTCouplingCompatibleWithModularSimulator =
            (temperatureCoupling == "no" || temperatureCoupling == "v-rescale");
    if (integrator == "md-vv" && pressureCoupling == "parrinello-rahman"
        && (isModularSimulatorExplicitlyDisabled || !isTCouplingCompatibleWithModularSimulator))
    {
        // Under md-vv, Parrinello-Rahman is only implemented for the modular simulator
        return;
    }
    if (additionalMdpParameters == MdpParameterDatabase::ExpandedEnsemble && isModularSimulatorExplicitlyDisabled)
    {
        // Checkpointing is disabled in the legacy simulator (#4629),
        // so exact continuation is impossible,
        // so we skip the test if it would run in the legacy simulator.
        // With the current test system, this only happens if modular simulator
        // is explicitly disabled, or if GPU update was requested (see #4711).
        return;
    }

    SCOPED_TRACE(
            formatString("Comparing normal and two-part run of simulation '%s' "
                         "with integrator '%s'",
                         simulationName.c_str(),
                         integrator.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(),
                                                integrator.c_str(),
                                                temperatureCoupling.c_str(),
                                                pressureCoupling.c_str(),
                                                additionalMdpParameters);
    // The exact lambda state choice is unimportant, so long as there
    // is one when using an FEP input.
    mdpFieldValues["init-lambda-state"] = "3";
    mdpFieldValues["nsteps"]            = "16";
    if (!binaryReproducible)
    {
        // With reproducible==false we don't generate the pair list
        // at the same step, so we should avoid missing non-bonded
        // interactions causing significant differences
        mdpFieldValues["verlet-buffer-tolerance"] = "1e-5";
    }

    int ulpToleranceInMixed  = 0;
    int ulpToleranceInDouble = 0;
    if (!binaryReproducible)
    {
        // Pair list generation at different steps results in different
        // sumation order.
        // Forces and update on GPUs are generally result in different
        // sumation order
        ulpToleranceInMixed  = 64;
        ulpToleranceInDouble = 128;
    }

    EnergyTermsToCompare energyTermsToCompare{
        { { interaction_function[F_EPOT].longname,
            relativeToleranceAsPrecisionDependentUlp(10.0, ulpToleranceInMixed, ulpToleranceInDouble) },
          { interaction_function[F_EKIN].longname,
            relativeToleranceAsPrecisionDependentUlp(10.0, ulpToleranceInMixed, ulpToleranceInDouble) } }
    };

    if (temperatureCoupling != "no" || pressureCoupling != "no")
    {
        if (simulationName == "alanine_vacuo")
        {
            // This is slightly less reproducible
            // Wider tolerance is needed because of sum of negative PE
            // with positive KE producing small conserved energy.
            const int factor = 4;
            energyTermsToCompare.insert(
                    { interaction_function[F_ECONSERVED].longname,
                      relativeToleranceAsPrecisionDependentUlp(
                              10.0, factor * ulpToleranceInMixed, factor * ulpToleranceInDouble) });
        }
        else
        {
            energyTermsToCompare.insert({ interaction_function[F_ECONSERVED].longname,
                                          relativeToleranceAsPrecisionDependentUlp(
                                                  10.0, ulpToleranceInMixed, ulpToleranceInDouble) });
        }
    }

    if (mdpFieldValues.count("free-energy") > 0 && mdpFieldValues.at("free-energy") != "no")
    {
        energyTermsToCompare.insert({ interaction_function[F_DVDL_COUL].longname,
                                      relativeToleranceAsPrecisionDependentUlp(
                                              10.0, ulpToleranceInMixed, ulpToleranceInDouble) });
        energyTermsToCompare.insert({ interaction_function[F_DVDL_VDW].longname,
                                      relativeToleranceAsPrecisionDependentUlp(
                                              10.0, ulpToleranceInMixed, ulpToleranceInDouble) });
        energyTermsToCompare.insert({ interaction_function[F_DVDL_BONDED].longname,
                                      relativeToleranceAsPrecisionDependentUlp(
                                              10.0, ulpToleranceInMixed, ulpToleranceInDouble) });
        energyTermsToCompare.insert({ interaction_function[F_DVDL_RESTRAINT].longname,
                                      relativeToleranceAsPrecisionDependentUlp(
                                              10.0, ulpToleranceInMixed, ulpToleranceInDouble) });
    }

    if (pressureCoupling == "parrinello-rahman")
    {
        energyTermsToCompare.insert({ "Box-Vel-XX",
                                      relativeToleranceAsPrecisionDependentUlp(
                                              1e-12, ulpToleranceInMixed, ulpToleranceInDouble) });
        energyTermsToCompare.insert({ "Box-Vel-YY",
                                      relativeToleranceAsPrecisionDependentUlp(
                                              1e-12, ulpToleranceInMixed, ulpToleranceInDouble) });
        energyTermsToCompare.insert({ "Box-Vel-ZZ",
                                      relativeToleranceAsPrecisionDependentUlp(
                                              1e-12, ulpToleranceInMixed, ulpToleranceInDouble) });
    }

    runTest(&fileManager_, &runner_, simulationName, mdpFieldValues, binaryReproducible, energyTermsToCompare);
}

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if !GMX_GPU_OPENCL

INSTANTIATE_TEST_SUITE_P(
        NormalIntegrators,
        MdrunNoAppendContinuationIsExact,
        ::testing::Combine(::testing::Values("argon12", "spc2", "alanine_vsite_vacuo"),
                           ::testing::Values("md", "md-vv", "bd", "sd"),
                           ::testing::Values("no"),
                           ::testing::Values("no"),
                           ::testing::Values(MdpParameterDatabase::Default),
                           ::testing::Values(true, false)));

INSTANTIATE_TEST_SUITE_P(NormalIntegratorsWithFEP,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                            ::testing::Values("md", "md-vv", "bd", "sd"),
                                            ::testing::Values("no"),
                                            ::testing::Values("no"),
                                            ::testing::Values(MdpParameterDatabase::Default),
                                            ::testing::Values(true, false)));

INSTANTIATE_TEST_SUITE_P(NVT,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("argon12"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("v-rescale", "nose-hoover"),
                                            ::testing::Values("no"),
                                            ::testing::Values(MdpParameterDatabase::Default),
                                            ::testing::Values(true, false)));

INSTANTIATE_TEST_SUITE_P(NPH,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("argon12"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("no"),
                                            // C-rescale temporarily removed because a reference temperature is needed
                                            //                           ::testing::Values("parrinello-rahman", "C-rescale"),
                                            ::testing::Values("parrinello-rahman"),
                                            ::testing::Values(MdpParameterDatabase::Default),
                                            ::testing::Values(true)));

INSTANTIATE_TEST_SUITE_P(NPT,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("argon12"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("v-rescale", "nose-hoover"),
                                            ::testing::Values("parrinello-rahman", "C-rescale"),
                                            ::testing::Values(MdpParameterDatabase::Default),
                                            ::testing::Values(true)));

INSTANTIATE_TEST_SUITE_P(MTTK,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("argon12"),
                                            ::testing::Values("md-vv"),
                                            ::testing::Values("nose-hoover"),
                                            ::testing::Values("mttk"),
                                            ::testing::Values(MdpParameterDatabase::Default),
                                            ::testing::Values(true)));

INSTANTIATE_TEST_SUITE_P(Pull,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("spc2"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("no"),
                                            ::testing::Values("no"),
                                            ::testing::Values(MdpParameterDatabase::Pull),
                                            ::testing::Values(true)));

INSTANTIATE_TEST_SUITE_P(Awh,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("alanine_vacuo"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("v-rescale"),
                                            ::testing::Values("no"),
                                            ::testing::Values(MdpParameterDatabase::Awh),
                                            ::testing::Values(true)));

INSTANTIATE_TEST_SUITE_P(ExpandedEnsemble,
                         MdrunNoAppendContinuationIsExact,
                         ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                            ::testing::Values("md-vv"),
                                            ::testing::Values("v-rescale"),
                                            ::testing::Values("no"),
                                            ::testing::Values(MdpParameterDatabase::ExpandedEnsemble),
                                            ::testing::Values(true)));

#else
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(MdrunNoAppendContinuationIsExact);
#endif

} // namespace
} // namespace test
} // namespace gmx
