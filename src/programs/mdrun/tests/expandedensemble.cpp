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

#include <cmath>

#include <filesystem>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <regex>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! The verifiable quantities reported per lambda
struct PerLambdaInformation
{
    long   index;
    double count, G, dG;
};

/*! \brief Helper struct for parsing the log-file MC-lambda output
 *
 * This contains both verifiable quantities, and a marker in the log
 * file at the end of the MC-lambda table to start the search for the
 * next table. */
struct LambdaInformation
{
    std::optional<std::size_t>        currentState;
    std::vector<PerLambdaInformation> perLambdaInformation;
    std::size_t                       finalPos;
};

/*! \brief Test fixture for expanded ensemble continuations
 *
 * This test ensures mdrun can run an expanded-ensemble simulation,
 * obtaining the histogram count from both .mdp fields or the
 * checkpoint. */
class ExpandedEnsembleTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string, std::string, std::string, MdpParameterDatabase>>
{
public:
    //! Expects that the Wang-Landau incrementor will not increase
    static double incrementorDoesNotIncrease(const std::string& logFileContents,
                                             double             lastValue,
                                             const int          expectedNumMatches)
    {
        std::regex incrementorRegex("Wang-Landau incrementor is: *(.*)");
        std::sregex_iterator logBegin(logFileContents.begin(), logFileContents.end(), incrementorRegex);
        std::sregex_iterator logEnd;
        int                  numMatches = 0;
        for (std::sregex_iterator i = logBegin; i != logEnd; ++i)
        {
            const std::smatch& match = *i;
            EXPECT_EQ(match.size(), 2);
            double value = std::stod(match[1].str());
            EXPECT_LE(value, lastValue);
            lastValue = value;
            numMatches++;
        }
        EXPECT_EQ(numMatches, expectedNumMatches);
        return lastValue;
    }
    //! Extract information from an MC-lambda table in the log file next following \c pos
    static LambdaInformation extractLambdaInformation(const std::size_t  numLambdas,
                                                      const std::string& logFileContents,
                                                      std::size_t        pos)
    {
        // Find the next MC-lambda table header, starting from pos
        pos = logFileContents.find("MC-lambda information", pos);
        if (pos == std::string::npos)
        {
            GMX_THROW(InvalidInputError("log file contents too short"));
        }
        // Now inspect the first table of lambda output to check
        // count and G. It can look a bit like
        //
        //             MC-lambda information
        //  Wang-Landau incrementor is:          10
        //  Nmass-lambdascoul-lambdasvdw-lambdasbonded-lambdasrestraint-lambdas    Count   G(in kT) dG(in kT)
        //  1  0.000  0.000  0.000  0.000  0.000   29.000    2.00000    1.00000
        //  2  0.500  0.000  0.000  0.000  0.000   30.000    3.00000    1.00000
        //  3  1.000  0.000  0.000  0.000  0.000   31.000    4.00000    1.00000
        //  4  1.000  0.000  0.000  0.500  0.000   32.000    5.00000    1.00000 <<
        //  5  1.000  0.000  0.000  1.000  0.000   33.000    6.00000    1.00000
        //  6  1.000  0.000  0.000  1.000  0.500   34.000    7.00000    1.00000
        //  7  1.000  0.000  0.000  1.000  1.000   35.000    8.00000    1.00000
        //  8  1.000  0.000  0.500  1.000  1.000   36.000    9.00000    1.00000
        //  9  1.000  0.000  1.000  1.000  1.000   37.000   10.00000    1.00000
        // 10  1.000  0.500  1.000  1.000  1.000   38.000   11.00000    1.00000
        // 11  1.000  1.000  1.000  1.000  1.000   39.000   12.00000

        // Skip some header lines
        std::size_t endOfLine = logFileContents.find('\n', pos);
        if (endOfLine == std::string::npos)
        {
            GMX_THROW(InvalidInputError("log file contents too short"));
        }
        endOfLine = logFileContents.find('\n', endOfLine + 1);
        if (endOfLine == std::string::npos)
        {
            GMX_THROW(InvalidInputError("log file contents too short"));
        }
        endOfLine = logFileContents.find('\n', endOfLine + 1);
        if (endOfLine == std::string::npos)
        {
            GMX_THROW(InvalidInputError("log file contents too short"));
        }
        LambdaInformation lambdaInformation;
        // Parse the rows of the table to fill lambdaInformation
        for (std::size_t i = 0; i < numLambdas; ++i)
        {
            const std::size_t beginningOfLine = endOfLine + 1;
            endOfLine                         = logFileContents.find('\n', beginningOfLine);
            if (endOfLine == std::string::npos)
            {
                GMX_THROW(InvalidInputError("log file contents too short"));
            }
            const std::string line = logFileContents.substr(beginningOfLine, endOfLine - beginningOfLine);
            std::vector<std::string> columns = splitString(line);
            if (columns.back() == "<<")
            {
                lambdaInformation.currentState = i;
                // Get rid of the indicator of the current lambda
                // state so we can index uniformly below.
                columns.pop_back();
            }
            if (i == numLambdas - 1)
            {
                // There is no dG entry in the final row of the table.
                // Append a zero to make the logic uniform across rows.
                columns.push_back("0.0");
            }
            lambdaInformation.perLambdaInformation.emplace_back(
                    PerLambdaInformation{ std::stol(columns[0]),
                                          std::stod(columns[columns.size() - 3]),
                                          std::stod(columns[columns.size() - 2]),
                                          std::stod(columns[columns.size() - 1]) });
            // Keep the final position to know where to start the
            // search for the next MC-lambda table.
            lambdaInformation.finalPos = endOfLine;
        }
        return lambdaInformation;
    }
};

TEST_F(ExpandedEnsembleTest, ContinuationPreservesExpandedEnsembleState)
{
    const char* simulationName          = "nonanol_vacuo";
    const char* integrator              = "md-vv";
    const char* temperatureCoupling     = "v-rescale";
    const char* pressureCoupling        = "";
    const auto  additionalMdpParameters = MdpParameterDatabase::ExpandedEnsemble;

    auto mdpFieldValues = prepareMdpFieldValues(
            simulationName, integrator, temperatureCoupling, pressureCoupling, additionalMdpParameters);
    // The exact choice of initial lambda state choice is unimportant,
    // so long as there is one.
    const int initLambdaState           = 3;
    mdpFieldValues["init-lambda-state"] = std::to_string(initLambdaState);
    const int numStepsInFirstPart       = 16;
    mdpFieldValues["nsteps"]            = std::to_string(numStepsInFirstPart);

    // The number of lambdas is set by the contents of the mdp files
    // for nonanol_vacuo in the simulation database.
    const int numLambdas = splitString(mdpFieldValues["mass-lambdas"]).size();

    // Set some arbitary unique values into the expanded-ensemble mdp fields
    std::vector<float> initLambdaWeightValues(numLambdas);
    std::iota(initLambdaWeightValues.begin(), initLambdaWeightValues.end(), 2);
    mdpFieldValues["init-lambda-weights"] = formatAndJoin(
            initLambdaWeightValues.begin(), initLambdaWeightValues.end(), " ", StringFormatter("%g"));

    std::vector<int> initLambdaCountValues(numLambdas);
    std::iota(initLambdaCountValues.begin(), initLambdaCountValues.end(), 14);
    mdpFieldValues["init-lambda-counts"] = formatAndJoin(
            initLambdaCountValues.begin(), initLambdaCountValues.end(), " ", StringFormatter("%d"));

    std::vector<float> initWlHistogramCountsValues(numLambdas);
    std::iota(initWlHistogramCountsValues.begin(), initWlHistogramCountsValues.end(), 29);
    mdpFieldValues["init-wl-histogram-counts"] = formatAndJoin(initWlHistogramCountsValues.begin(),
                                                               initWlHistogramCountsValues.end(),
                                                               " ",
                                                               StringFormatter("%g"));

    // Don't use domain decompositions that won't work
    const int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        GTEST_SKIP() << formatString(
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName,
                numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
    }

    // Make .tpr file
    {
        CommandLine caller;
        caller.append("grompp");
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Run first simulation part
    {
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    double            lastValueOfIncrementor = std::stod(mdpFieldValues["init-wl-delta"]);
    LambdaInformation finalLambdaInformation;
    {
        auto logFileContents = TextReader::readFileToString(runner_.logFileName_);

        {
            SCOPED_TRACE("Verify that Wang-Landau incrementor does not increase");
            lastValueOfIncrementor = incrementorDoesNotIncrease(
                    logFileContents, lastValueOfIncrementor, numStepsInFirstPart + 1);
        }

        {
            SCOPED_TRACE("Check that contents of mdp vectors are transferred to mdrun");

            // The format string matches that in pr_rvec() used when
            // dumping the inputrec to the mdrun log file.
            for (int i = 0; i < numLambdas; ++i)
            {
                const std::string expectedValueFromMdpFile = formatString(
                        "%s[%d]=%12.5e", "init-lambda-weights", i, initLambdaWeightValues[i]);
                EXPECT_NE(logFileContents.find(expectedValueFromMdpFile), std::string::npos);
            }
            for (int i = 0; i < numLambdas; ++i)
            {
                const std::string expectedValueFromMdpFile =
                        formatString("%s[%d]=%12.5e",
                                     "init-lambda-counts",
                                     i,
                                     static_cast<float>(initLambdaCountValues[i]));
                EXPECT_NE(logFileContents.find(expectedValueFromMdpFile), std::string::npos);
            }
            for (int i = 0; i < numLambdas; ++i)
            {
                const std::string expectedValueFromMdpFile = formatString(
                        "%s[%d]=%12.5e", "init-wl-histogram-counts", i, initWlHistogramCountsValues[i]);
                EXPECT_NE(logFileContents.find(expectedValueFromMdpFile), std::string::npos);
            }
        }

        std::size_t pos = 0;
        {
            SCOPED_TRACE("Check that mdp input values end up in the first histogram output");

            // The first report of MC-lambda information contains the
            // Wang-Landau histogram counts that were passed from the
            // .mdp file, so we can check that. Values at subsequent
            // steps are not checked in this test.
            LambdaInformation lambdaInformation =
                    extractLambdaInformation(numLambdas, logFileContents, pos);
            ASSERT_TRUE(lambdaInformation.currentState.has_value());
            EXPECT_EQ(lambdaInformation.currentState, initLambdaState);
            EXPECT_EQ(lambdaInformation.perLambdaInformation.size(), numLambdas);
            for (std::size_t i = 0; i != lambdaInformation.perLambdaInformation.size(); ++i)
            {
                EXPECT_EQ(i + 1, lambdaInformation.perLambdaInformation[i].index);
                EXPECT_EQ(initWlHistogramCountsValues[i],
                          std::lround(lambdaInformation.perLambdaInformation[i].count));
                EXPECT_EQ(initLambdaWeightValues[i], lambdaInformation.perLambdaInformation[i].G);
            }
            pos = lambdaInformation.finalPos;
        }

        {
            SCOPED_TRACE(
                    "Capture the final MC-lambda information to verify it is faithfully restored "
                    "from the checkpoint");
            std::size_t lastPos = std::string::npos;
            // Search for the final MC-lambda information section in the log file
            do
            {
                lastPos = pos;
                pos     = logFileContents.find("MC-lambda information", lastPos + 1);
            } while (pos != std::string::npos);

            finalLambdaInformation = extractLambdaInformation(numLambdas, logFileContents, lastPos);
            ASSERT_NE(finalLambdaInformation.finalPos, std::string::npos)
                    << "Expected more than one MC-information output";
        }
    }

    // Check for unimplemented functionality
    // TODO: Update this as modular simulator gains functionality
    const bool isModularSimulatorExplicitlyDisabled = (getenv("GMX_DISABLE_MODULAR_SIMULATOR") != nullptr);
    if (additionalMdpParameters == MdpParameterDatabase::ExpandedEnsemble && isModularSimulatorExplicitlyDisabled)
    {
        // Checkpointing is disabled in the legacy simulator (#4629),
        // so continuation is impossible, so we skip the remainder of
        // the test if it would run in the legacy simulator. With the
        // current test system, this only happens if modular simulator
        // is explicitly disabled, or if GPU update was requested (see
        // #4711).
        return;
    }

    const int numStepsInSecondPart = 8;
    {
        SCOPED_TRACE("Expanded ensemble can do checkpoint restart");
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        mdrunCaller.addOption("-cpi", runner_.cptOutputFileName_);
        mdrunCaller.addOption("-noappend");
        mdrunCaller.addOption("-nsteps", numStepsInSecondPart);
        runner_.logFileName_ = fileManager_.getTemporaryFilePath(".part0002.log").string();
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    {
        auto logFileContents = TextReader::readFileToString(runner_.logFileName_);

        {
            SCOPED_TRACE(
                    "Verify that Wang-Landau incrementor continues to not increase starting from "
                    "the value from the first simulation part");
            incrementorDoesNotIncrease(logFileContents, lastValueOfIncrementor, numStepsInSecondPart + 1);
        }

        {
            SCOPED_TRACE(
                    "Check that checkpointed values end up in the first histogram output after the "
                    "restart");

            // The first report of MC-lambda information contains the
            // Wang-Landau histogram counts that were passed from the
            // checkpoint file.
            LambdaInformation lambdaInformation =
                    extractLambdaInformation(numLambdas, logFileContents, 0);
            EXPECT_EQ(lambdaInformation.currentState, finalLambdaInformation.currentState);
            EXPECT_EQ(lambdaInformation.perLambdaInformation.size(), numLambdas);
            ASSERT_EQ(lambdaInformation.perLambdaInformation.size(),
                      finalLambdaInformation.perLambdaInformation.size());
            for (std::size_t i = 0; i != lambdaInformation.perLambdaInformation.size(); ++i)
            {
                EXPECT_EQ(lambdaInformation.perLambdaInformation[i].index,
                          finalLambdaInformation.perLambdaInformation[i].index);
                EXPECT_EQ(lambdaInformation.perLambdaInformation[i].count,
                          finalLambdaInformation.perLambdaInformation[i].count);
                EXPECT_EQ(lambdaInformation.perLambdaInformation[i].G,
                          finalLambdaInformation.perLambdaInformation[i].G);
                EXPECT_EQ(lambdaInformation.perLambdaInformation[i].dG,
                          finalLambdaInformation.perLambdaInformation[i].dG);
            }
        }
    }
}

TEST_F(ExpandedEnsembleTest, WeightEquilibrationReported)
{
    const char* simulationName          = "nonanol_vacuo";
    const char* integrator              = "md-vv";
    const char* temperatureCoupling     = "v-rescale";
    const char* pressureCoupling        = "";
    const auto  additionalMdpParameters = MdpParameterDatabase::ExpandedEnsemble;

    auto mdpFieldValues = prepareMdpFieldValues(
            simulationName, integrator, temperatureCoupling, pressureCoupling, additionalMdpParameters);
    // The exact choice of initial lambda state choice is unimportant,
    // so long as there is one.
    mdpFieldValues["init-lambda-state"] = "3";
    mdpFieldValues["nsteps"]            = "1";
    // Values chosen to provoke logfile output at the first nstexpanded step
    mdpFieldValues["lmc-weights-equil"]              = "number-all-lambda";
    mdpFieldValues["weight-equil-number-all-lambda"] = "4";
    mdpFieldValues["init-lambda-counts"]             = "5 5 5 5 5 5 5 5 5 5 5";

    // Don't use domain decompositions that won't work
    const int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        GTEST_SKIP() << formatString(
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName,
                numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
    }

    // Make .tpr file
    {
        CommandLine caller;
        caller.append("grompp");
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Run simulation
    {
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    auto logFileContents = TextReader::readFileToString(runner_.logFileName_);
    EXPECT_NE(logFileContents.find(formatString(
                      "Step %s: Weights have equilibrated, using criteria: number-all-lambda",
                      mdpFieldValues["nstexpanded"].c_str())),
              std::string::npos);
}

} // namespace
} // namespace test
} // namespace gmx
