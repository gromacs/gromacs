/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Tests to verify that a simulator that only does some actions
 * periodically produces the expected results.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <tuple>

#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Mdp parameters that determine the manner of simulation
 * propagation. */
using PropagationParameters = MdpFieldValues;

/*! \brief Mdp parameters that should only affect the observations,
 *  not the simulation propagation. */
using PeriodicOutputParameters = MdpFieldValues;

//! Helper type of output file names for the reference mdrun call
struct ReferenceFileNames
{
    //! Name of energy file
    std::string edrFileName_;
};

//! Function type to produce sets of .mdp parameters for testing periodic output
using OutputParameterGeneratorFunction = std::vector<PeriodicOutputParameters> (*)();

/*! \brief Test fixture base for comparing a simulator with one that
 * does everything every step
 *
 * This test ensures that two simulator code paths called via
 * different mdp options yield identical energy trajectories,
 * up to some (arbitrary) precision.
 *
 * These tests are useful to check that periodic actions implemented
 * in simulators are correct, and that different code paths expected
 * to yield identical results are equivalent.
 */
class PeriodicActionsTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<PropagationParameters, OutputParameterGeneratorFunction>>
{
public:
    // PeriodicActionsTest();
    //! Run mdrun with given output parameters
    void doMdrun(const PeriodicOutputParameters& output);
    //! Generate reference data from mdrun writing everything every step.
    void prepareReferenceData();
    //! Names for the output files from the reference mdrun call
    ReferenceFileNames referenceFileNames_ = { fileManager_.getTemporaryFilePath("reference.edr") };
    //! Functor for energy comparison
    EnergyComparison energyComparison_{ EnergyComparison::defaultEnergyTermsToCompare() };
    //! Names of energies compared by energyComparison_
    std::vector<std::string> namesOfEnergiesToMatch_ = energyComparison_.getEnergyNames();
};

void PeriodicActionsTest::doMdrun(const PeriodicOutputParameters& output)
{
    auto propagation = std::get<0>(GetParam());
    SCOPED_TRACE(
            formatString("Doing %s simulation with %s integrator, %s tcoupling and %s pcoupling\n",
                         propagation["simulationName"].c_str(), propagation["integrator"].c_str(),
                         propagation["tcoupl"].c_str(), propagation["pcoupl"].c_str()));
    auto mdpFieldValues = prepareMdpFieldValues(propagation["simulationName"], propagation["integrator"],
                                                propagation["tcoupl"], propagation["pcoupl"]);
    mdpFieldValues.insert(propagation.begin(), propagation.end());
    mdpFieldValues.insert(output.begin(), output.end());

    // prepare the tpr file
    {
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", propagation["maxGromppWarningsTolerated"]);
        runner_.useTopGroAndNdxFromDatabase(propagation["simulationName"]);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // do a normal mdrun
    {
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }
}

void PeriodicActionsTest::prepareReferenceData()
{
    SCOPED_TRACE("Preparing reference data");

    // Configure the SimulationRunner to write output to files we can
    // use in comparing each test run.
    std::swap(runner_.edrFileName_, referenceFileNames_.edrFileName_);

    // Run the reference simulation with everything output every step
    PeriodicOutputParameters outputEveryStep = {
        { "nstenergy", "1" },
        { "nstlog", "1" },
        { "nstdhdl", "1" },
        { "description", "output everything every step" },
    };
    doMdrun(outputEveryStep);

    // Restore the standard filenames we'll use for the runs whose behavior we are testing.
    std::swap(runner_.edrFileName_, referenceFileNames_.edrFileName_);
}

/*! \brief Compare the next frame returned by \c testFrameReader with the
 * matching frame from \c referenceTestReader using \c comparator.
 *
 * \returns  True when a successful comparison was made
 */
template<typename Reader, typename Comparator>
bool compareFrames(Reader referenceFrameReader, Reader testFrameReader, const Comparator& comparator)
{
    if (!testFrameReader->readNextFrame())
    {
        return false;
    }

    auto        testFrame = testFrameReader->frame();
    std::string frameName = testFrame.frameName();
    SCOPED_TRACE("Found frame from test file named " + frameName);
    bool foundMatch = false;
    while (!foundMatch)
    {
        if (referenceFrameReader->readNextFrame())
        {
            auto referenceFrame = referenceFrameReader->frame();
            // Can the reference and test frames be compared, because they have the same name?
            if (frameName == referenceFrame.frameName())
            {
                SCOPED_TRACE("Found frame from reference file named " + frameName);
                comparator(referenceFrame, testFrame);
                foundMatch = true;
            }
        }
        else
        {
            ADD_FAILURE() << "Ran out of reference frames to compare";
            return false;
        }
    }
    return foundMatch;
}

TEST_P(PeriodicActionsTest, PeriodicActionsAgreeWithReference)
{
    auto propagation = std::get<0>(GetParam());
    SCOPED_TRACE(formatString("Comparing two simulations of '%s' with integrator '%s'",
                              propagation["simulationName"].c_str(), propagation["integrator"].c_str()));

    prepareReferenceData();

    auto outputParametersGenerator = std::get<1>(GetParam());
    for (const PeriodicOutputParameters& output : outputParametersGenerator())
    {
        SCOPED_TRACE("Comparing to observe " + output.at("description"));

        // Run the test simulation
        doMdrun(output);

        // Prepare readers for the reference and test output files
        auto referenceEnergyFrameReader =
                openEnergyFileToReadTerms(referenceFileNames_.edrFileName_, namesOfEnergiesToMatch_);
        auto testEnergyFrameReader =
                openEnergyFileToReadTerms(runner_.edrFileName_, namesOfEnergiesToMatch_);

        const bool shouldCompareEnergies   = fromString<int>(output.at("nstenergy")) > 0;
        bool       shouldContinueComparing = shouldCompareEnergies;
        while (shouldContinueComparing)
        {
            if (shouldCompareEnergies)
            {
                SCOPED_TRACE("Comparing energy frames from reference '" + referenceFileNames_.edrFileName_
                             + "' and test '" + runner_.edrFileName_ + "'");
                shouldContinueComparing = shouldContinueComparing
                                          && compareFrames(referenceEnergyFrameReader.get(),
                                                           testEnergyFrameReader.get(), energyComparison_);
            }
        }
    }
}

/*! \brief Some common choices of periodic output mdp parameters to
 * simplify defining values for the combinations under test */
PeriodicOutputParameters g_basicPeriodicOutputParameters = {
    { "nstenergy", "0" }, { "nstlog", "0" },           { "nstxout", "0" },
    { "nstvout", "0" },   { "nstfout", "0" },          { "nstxout-compressed", "0" },
    { "nstdhdl", "0" },   { "description", "unknown" }
};

/*! \brief Return vector of mdp parameter sets to test
 *
 * These are constructed to observe the mdp parameter choices that
 * only affect when output is written agree with those that were
 * written from a reference run where output was done every step. The
 * numbers are chosen in the context of the defaults in
 * prepareDefaultMdpFieldValues().
 *
 * \todo Test nstlog, nstdhdl, nstxout-compressed */
std::vector<PeriodicOutputParameters> outputParameters()
{
    std::vector<PeriodicOutputParameters> parameterSets;

    parameterSets.push_back(g_basicPeriodicOutputParameters);
    parameterSets.back()["nstenergy"]   = "1";
    parameterSets.back()["description"] = "energies every step works";

    parameterSets.push_back(g_basicPeriodicOutputParameters);
    parameterSets.back()["nstenergy"]   = "4";
    parameterSets.back()["description"] = "energies every fourth step works";

    parameterSets.push_back(g_basicPeriodicOutputParameters);
    parameterSets.back()["nstxout"]     = "4";
    parameterSets.back()["description"] = "coordinates every fourth step works";

    parameterSets.push_back(g_basicPeriodicOutputParameters);
    parameterSets.back()["nstvout"]     = "4";
    parameterSets.back()["description"] = "velocities every fourth step works";

    parameterSets.push_back(g_basicPeriodicOutputParameters);
    parameterSets.back()["nstfout"]     = "4";
    parameterSets.back()["description"] = "forces every fourth step works";

    return parameterSets;
}

//! Returns sets of simple simulation propagators
std::vector<PropagationParameters> simplePropagationParameters()
{
    return {
        { { "simulationName", "argon12" },
          { "integrator", "md" },
          { "comm-mode", "linear" },
          { "nstcomm", "1" },
          { "tcoupl", "no" },
          { "nsttcouple", "0" },
          { "pcoupl", "no" },
          { "nstpcouple", "0" },
          { "maxGromppWarningsTolerated", "0" } },
    };
}

/*! \brief Returns sets of simulation propagators including coupling
 *
 * These are chosen to cover the commonly used space of propagation
 * algorithms togther with the perdiods between their actions. The
 * periods tested are chosen to be mutually co-prime and distinct from
 * the pair search and user output period (4), so that over the
 * duration of a short simulation many kinds of simulation step
 * behavior are tested. */
std::vector<PropagationParameters> propagationParametersWithCoupling()
{
    std::string nsttcouple = "2";
    std::string nstpcouple = "3";
    std::string comm_mode  = "linear";
    std::string nstcomm    = "5";

    std::vector<PropagationParameters> parameterSets;
    for (std::string simulationName : { "argon12" })
    {
        for (std::string integrator : { "md", "sd", "md-vv" })
        {
            for (std::string tcoupl : { "no", "v-rescale", "Nose-Hoover" })
            {
                // SD doesn't support temperature-coupling algorithms,
                if (integrator == "sd" && tcoupl != "no")
                {
                    continue;
                }
                for (std::string pcoupl : { "no", "Berendsen", "Parrinello-Rahman" })
                {
                    // VV supports few algorithm combinations
                    if (integrator == "md-vv")
                    {
                        // P-R with VV is called MTTK
                        if (pcoupl == "Parrinello-Rahman")
                        {
                            pcoupl = "MTTK";
                        }
                        if ((tcoupl == "Nose-Hoover" && pcoupl == "Berendsen")
                            || (tcoupl != "Nose-Hoover" && pcoupl == "MTTK"))
                        {
                            continue;
                        }
                    }
                    if (pcoupl == "MTTK" && simulationName != "argon12")
                    {
                        // MTTK does not support constraints
                        continue;
                    }

                    int maxGromppWarningsTolerated = 0;
                    if (pcoupl == "Berendsen")
                    {
                        ++maxGromppWarningsTolerated;
                    }
                    parameterSets.emplace_back(PropagationParameters{
                            { "simulationName", simulationName },
                            { "integrator", integrator },
                            { "comm-mode", comm_mode },
                            { "nstcomm", nstcomm },
                            { "tcoupl", tcoupl },
                            { "nsttcouple", nsttcouple },
                            { "pcoupl", pcoupl },
                            { "nstpcouple", nstpcouple },
                            { "maxGromppWarningsTolerated", toString(maxGromppWarningsTolerated) } });
                }
            }
        }
    }
    return parameterSets;
}

/*! \brief Returns sets of simulation propagators including coupling
 *
 * These are chosen to cover the commonly used space of propagation
 * algorithms on systems with constraints. */
std::vector<PropagationParameters> propagationParametersWithConstraints()
{
    std::string nsttcouple = "2";
    std::string nstpcouple = "3";
    std::string comm_mode  = "linear";
    std::string nstcomm    = "5";

    std::vector<PropagationParameters> parameterSets;
    for (std::string simulationName : { "tip3p5" })
    {
        for (std::string integrator : { "md", "sd", "md-vv" })
        {
            for (std::string tcoupl : { "no", "v-rescale" })
            {
                // SD doesn't support temperature-coupling algorithms,
                if (integrator == "sd" && tcoupl != "no")
                {
                    continue;
                }
                for (std::string pcoupl : { "no", "Parrinello-Rahman" })
                {
                    // VV supports few algorithm combinations
                    if (integrator == "md-vv")
                    {
                        // P-R with VV is called MTTK
                        if (pcoupl == "Parrinello-Rahman")
                        {
                            pcoupl = "MTTK";
                        }
                        if ((tcoupl == "Nose-Hoover" && pcoupl == "Berendsen")
                            || (tcoupl != "Nose-Hoover" && pcoupl == "MTTK"))
                        {
                            continue;
                        }
                    }
                    if (pcoupl == "MTTK" && simulationName != "argon12")
                    {
                        // MTTK does not support constraints
                        continue;
                    }

                    int maxGromppWarningsTolerated = 0;
                    parameterSets.emplace_back(PropagationParameters{
                            { "simulationName", simulationName },
                            { "integrator", integrator },
                            { "comm-mode", comm_mode },
                            { "nstcomm", nstcomm },
                            { "tcoupl", tcoupl },
                            { "nsttcouple", nsttcouple },
                            { "pcoupl", pcoupl },
                            { "nstpcouple", nstpcouple },
                            { "maxGromppWarningsTolerated", toString(maxGromppWarningsTolerated) } });
                }
            }
        }
    }
    return parameterSets;
}

using ::testing::Combine;
using ::testing::Values;
using ::testing::ValuesIn;

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if GMX_GPU != GMX_GPU_OPENCL
INSTANTIATE_TEST_CASE_P(BasicPropagators,
                        PeriodicActionsTest,
                        Combine(ValuesIn(simplePropagationParameters()), Values(outputParameters)));
INSTANTIATE_TEST_CASE_P(PropagatorsWithCoupling,
                        PeriodicActionsTest,
                        Combine(ValuesIn(propagationParametersWithCoupling()), Values(outputParameters)));
INSTANTIATE_TEST_CASE_P(PropagatorsWithConstraints,
                        PeriodicActionsTest,
                        Combine(ValuesIn(propagationParametersWithConstraints()), Values(outputParameters)));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_BasicPropagators,
                        PeriodicActionsTest,
                        Combine(ValuesIn(simplePropagationParameters()), Values(outputParameters)));
INSTANTIATE_TEST_CASE_P(DISABLED_PropagatorsWithCoupling,
                        PeriodicActionsTest,
                        Combine(ValuesIn(propagationParametersWithCoupling()), Values(outputParameters)));
INSTANTIATE_TEST_CASE_P(DISABLED_PropagatorsWithConstraints,
                        PeriodicActionsTest,
                        Combine(ValuesIn(propagationParametersWithConstraints()), Values(outputParameters)));
#endif

} // namespace
} // namespace test
} // namespace gmx
