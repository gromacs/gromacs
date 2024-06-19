/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Utility functions for tests to verify that a simulator that only does some actions
 * periodically produces the expected results.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "periodicactions.h"

#include "config.h"

#include <initializer_list>
#include <map>
#include <utility>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"

#include "programs/mdrun/tests/energycomparison.h"
#include "programs/mdrun/tests/energyreader.h"
#include "programs/mdrun/tests/moduletest.h"

namespace gmx
{
namespace test
{

/*! \brief Mdp parameters that determine the manner of simulation
 * propagation. */
using PropagationParameters = MdpFieldValues;

/*! \brief Mdp parameters that should only affect the observations,
 *  not the simulation propagation. */
using PeriodicOutputParameters = MdpFieldValues;

//! Function type to produce sets of .mdp parameters for testing periodic output
using OutputParameterGeneratorFunction = std::vector<PeriodicOutputParameters> (*)();

void PeriodicActionsTest::doMdrun(const PeriodicOutputParameters& output)
{
    auto propagation = std::get<0>(GetParam());
    SCOPED_TRACE(
            formatString("Doing %s simulation with %s integrator, %s tcoupling and %s pcoupling\n",
                         propagation["simulationName"].c_str(),
                         propagation["integrator"].c_str(),
                         propagation["tcoupl"].c_str(),
                         propagation["pcoupl"].c_str()));
    auto mdpFieldValues = prepareMdpFieldValues(propagation["simulationName"],
                                                propagation["integrator"],
                                                propagation["tcoupl"],
                                                propagation["pcoupl"]);

    // This lambda writes all mdp options in `source` into `target`, overwriting options already
    // present in `target`. It also filters out non-mdp option entries in the source maps
    auto overWriteMdpMapValues = [](const MdpFieldValues& source, MdpFieldValues& target) {
        for (auto const& [key, value] : source)
        {
            if (key == "simulationName" || key == "maxGromppWarningsTolerated" || key == "description")
            {
                // Remove non-mdp entries used in propagation and output
                continue;
            }
            target[key] = value;
        }
    };
    // Add options in propagation and output to the mdp options
    overWriteMdpMapValues(propagation, mdpFieldValues);
    overWriteMdpMapValues(output, mdpFieldValues);

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
                              propagation["simulationName"].c_str(),
                              propagation["integrator"].c_str()));

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
                                                           testEnergyFrameReader.get(),
                                                           energyComparison_);
            }
        }
    }
}

/*! \brief Some common choices of periodic output mdp parameters to
 * simplify defining values for the combinations under test */
static PeriodicOutputParameters g_basicPeriodicOutputParameters = {
    { "nstenergy", "0" }, { "nstlog", "0" },           { "nstxout", "0" },
    { "nstvout", "0" },   { "nstfout", "0" },          { "nstxout-compressed", "0" },
    { "nstdhdl", "0" },   { "description", "unknown" }
};

// \todo Test nstlog, nstdhdl, nstxout-compressed
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

std::vector<PropagationParameters> propagationParametersWithCoupling()
{
    std::string nsttcouple = "2";
    std::string nstpcouple = "3";
    std::string comm_mode  = "linear";
    std::string nstcomm    = "5";

    std::vector<PropagationParameters> parameterSets;
    std::vector<std::string>           simulations = { "argon12" };
    for (const std::string& simulationName : simulations)
    {
        std::vector<std::string> integrators{ "md", "sd", "md-vv" };
        for (const std::string& integrator : integrators)
        {
            std::vector<std::string> tcouplValues{ "no", "v-rescale", "Nose-Hoover" };
            for (const std::string& tcoupl : tcouplValues)
            {
                // SD doesn't support temperature-coupling algorithms,
                if (integrator == "sd" && tcoupl != "no")
                {
                    continue;
                }
                for (std::string pcoupl : { "no", "Berendsen", "Parrinello-Rahman", "C-rescale" })
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

std::vector<PropagationParameters> propagationParametersWithConstraints()
{
    std::string nsttcouple = "2";
    std::string nstpcouple = "3";
    std::string comm_mode  = "linear";
    std::string nstcomm    = "5";

    std::vector<PropagationParameters> parameterSets;
    std::vector<std::string>           simulations = { "tip3p5" };
    for (const std::string& simulationName : simulations)
    {
        std::vector<std::string> integrators{ "md", "sd", "md-vv" };
        for (const std::string& integrator : integrators)
        {
            std::vector<std::string> tcouplValues{ "no", "v-rescale" };
            for (const std::string& tcoupl : tcouplValues)
            {
                // SD doesn't support temperature-coupling algorithms,
                if (integrator == "sd" && tcoupl != "no")
                {
                    continue;
                }
                for (std::string pcoupl : { "no", "Parrinello-Rahman", "C-rescale" })
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

} // namespace test
} // namespace gmx
