/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Regression style tests for mdrun functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energyreader.h"
#include "mdruncomparison.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

// TODO does this belong in a more central place?
//! Convenience function to get std::string keys from some kind of \c Map.
template <typename Map>
std::vector<std::string> getKeys(const Map &m)
{
    std::vector<std::string> keys;
    for (const auto &it : m)
    {
        keys.push_back(it.first);
    }
    return keys;
}

//! Run grompp and then mdrun
void runTest(TestFileManager        *fileManager,
             SimulationRunner       *runner,
             const std::string      &simulationName,
             int                     maxWarningsTolerated,
             const MdpFieldValues   &mdpFieldValues,
             const EnergyTolerances &energiesToMatch)
{
    auto trajectoryFileName = fileManager->getTemporaryFilePath(".trr");
    auto edrFileName        = fileManager->getTemporaryFilePath(".edr");

    // prepare the .tpr file
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", maxWarningsTolerated);
        runner->useTopGroAndNdxFromDatabase(simulationName);
        runner->useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner->callGrompp(caller));
    }

    // run mdrun
    //
    // TODO Implement a way to loop over things like -nb, -pme, -npme
    // in a stable way, comparing with the same reference data each
    // time.
    {
        runner->fullPrecisionTrajectoryFileName_ = trajectoryFileName;
        runner->edrFileName_                     = edrFileName;
        CommandLine normalRunCaller;
        normalRunCaller.append("mdrun");
        ASSERT_EQ(0, runner->callMdrun(normalRunCaller));
    }

    // Check the energies agree with the refdata within \c tolerance.

    auto namesOfEnergiesToMatch = getKeys(energiesToMatch);
    auto energyReader           = openEnergyFileToReadFields(edrFileName,
                                                             namesOfEnergiesToMatch);

    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());

    int64_t              frameNumber = 0;
    while (energyReader->readNextFrame())
    {
        auto frameChecker = rootChecker.checkCompound("Frame", int64ToString(frameNumber));
        auto frame        = energyReader->frame();
        frameChecker.checkString(frame.frameName(), "Name");
        for (const auto &energyToMatch : energiesToMatch)
        {
            auto &energyName = energyToMatch.first;
            SCOPED_TRACE("Comparing " +  energyName + " with reference data");
            auto  energyIterator = frame.find(energyName);
            if (energyIterator != frame.end())
            {
                auto &energyValue = energyIterator->second;
                frameChecker.checkReal(energyValue, energyName.c_str());
                // energyToMatch->second);
            }
            else
            {
                ADD_FAILURE() << "Could not find required energy component";
            }
        }
        ++frameNumber;
    }
}

/*! \brief Test fixture base for nbnxn module in mdrun
 *
 * \todo Use a struct for the parameterization so we know what the
 * different strings are.
 */
class NbnxnRegressionTest : public MdrunTestFixture,
                            public ::testing::WithParamInterface <
                            std::tuple < std::string>>
{
};

TEST_P(NbnxnRegressionTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    SCOPED_TRACE(formatString("Checking run of simulation '%s' ",
                              simulationName.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(),
                                                "md", "no", "no");

    using MdpField = MdpFieldValues::value_type;
    // TODO Work out how to test the range of things the nbnxn module
    // supports, e.g. PME (including different orders and LJ-PME), RF
    // (and flavors), VDW switch behaviors, combination rules
    //
    // TODO This doesn't work, but something like it probably could
    // with some more work.
    mdpFieldValues.insert(MdpField("other",
                                   R"(coulombtype = pme)"));

    // TODO Decide some useful tolerances, and how to decide them.
    EnergyTolerances energiesToMatch
    {{
         {
             interaction_function[F_COUL_SR].longname,
             relativeToleranceAsPrecisionDependentUlp(-20000, 5, 10)
         },
         // TODO work out how to avoid inserting this when not testing electrostatic PME
         {
             interaction_function[F_COUL_RECIP].longname,
             relativeToleranceAsPrecisionDependentUlp(-19980, 5, 10)
         },
         {
             interaction_function[F_EPOT].longname,
             relativeToleranceAsPrecisionDependentUlp(-22000, 5, 10)
         },
         {
             interaction_function[F_PRES].longname,
             relativeToleranceAsPrecisionDependentUlp(1350.0, 5, 10)
         },
     }};

    int numWarningsToTolerate = 0;
    runTest(&fileManager_, &runner_,
            simulationName, numWarningsToTolerate, mdpFieldValues,
            energiesToMatch);
}

INSTANTIATE_TEST_CASE_P(IsReproduced, NbnxnRegressionTest,
                            ::testing::Values("hexane-and-water"));

} // namespace
} // namespace
} // namespace
