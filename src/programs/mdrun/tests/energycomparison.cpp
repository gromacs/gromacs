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
 * \brief Implementions of related classes for tests that want to
 * inspect energies produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "energycomparison.h"

#include <cstdint>

#include <map>
#include <memory>
#include <utility>

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "programs/mdrun/tests/comparison_helpers.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{

EnergyTermsToCompare EnergyComparison::defaultEnergyTermsToCompare()
{
    return {
        { interaction_function[F_EPOT].longname, relativeToleranceAsPrecisionDependentUlp(60.0, 200, 160) },
        { interaction_function[F_EKIN].longname, relativeToleranceAsPrecisionDependentUlp(60.0, 200, 160) },
        // The pressure is very strongly affected by summation errors,
        // so we need a large tolerance.
        // The value of 18000 is calibrated for running a small water box for 16 steps.
        // For a single frame for a water box a value of 150 could work.
        { interaction_function[F_PRES].longname, relativeToleranceAsUlp(10.0, 18000) },
    };
};

EnergyComparison::EnergyComparison(const EnergyTermsToCompare& energyTermsToCompare,
                                   MaxNumFrames                maxNumFrames) :
    energyTermsToCompare_(energyTermsToCompare), maxNumFrames_(maxNumFrames)
{
}

std::vector<std::string> EnergyComparison::getEnergyNames() const
{
    std::vector<std::string> keys;
    keys.reserve(energyTermsToCompare_.size());
    for (const auto& it : energyTermsToCompare_)
    {
        keys.push_back(it.first);
    }
    return keys;
}

void EnergyComparison::operator()(const EnergyFrame& reference, const EnergyFrame& test) const
{
    if (numComparedFrames_ >= maxNumFrames_)
    {
        // Nothing should be compared
        return;
    }

    SCOPED_TRACE("Comparing energy reference frame " + reference.frameName() + " and test frame "
                 + test.frameName());
    for (auto referenceIt = reference.begin(); referenceIt != reference.end(); ++referenceIt)
    {
        const auto& energyName = referenceIt->first;
        SCOPED_TRACE("Comparing " + energyName + " between frames");
        auto testIt = test.find(energyName);
        if (testIt != test.end())
        {
            const auto& energyValueInReference = referenceIt->second;
            const auto& energyValueInTest      = testIt->second;
            EXPECT_REAL_EQ_TOL(
                    energyValueInReference, energyValueInTest, energyTermsToCompare_.at(energyName));
        }
        else
        {
            ADD_FAILURE() << "Could not find energy component from reference frame in test frame";
        }
    }
    numComparedFrames_++;
}

void checkEnergiesAgainstReferenceData(const std::string&          energyFilename,
                                       const EnergyTermsToCompare& energyTermsToCompare,
                                       TestReferenceChecker*       checker,
                                       MaxNumFrames                maxNumEnergyFrames)
{
    const bool thisRankChecks = (gmx_node_rank() == 0);

    if (thisRankChecks)
    {
        EnergyComparison energyComparison(energyTermsToCompare, maxNumEnergyFrames);
        auto energyReader = openEnergyFileToReadTerms(energyFilename, energyComparison.getEnergyNames());

        std::unordered_map<std::string, TestReferenceChecker> checkers;
        for (const auto& energyTermToCompare : energyTermsToCompare)
        {
            const auto& energyName      = energyTermToCompare.first;
            checkers[energyName]        = checker->checkCompound("Energy", energyName.c_str());
            const auto& energyTolerance = energyTermToCompare.second;
            checkers[energyName].setDefaultTolerance(energyTolerance);
        }

        // We can't assume that frame names based purely on frame
        // contents are unique. For example, CG can write multiple
        // frames with the same step number. But we need a unique
        // identifier so we match the intended reference data, so we
        // keep track of the number of the frame read from the file.
        unsigned int frameNumber = 0;
        while (frameNumber < maxNumEnergyFrames && energyReader->readNextFrame())
        {
            const EnergyFrame& frame = energyReader->frame();
            const std::string  frameName =
                    frame.frameName() + " in frame " + toString(static_cast<int64_t>(frameNumber));

            SCOPED_TRACE("Comparing frame " + frameName);
            for (const auto& energyTermToCompare : energyTermsToCompare)
            {
                const std::string& energyName  = energyTermToCompare.first;
                const real         energyValue = frame.at(energyName);

                SCOPED_TRACE("Comparing energy " + energyName);
                checkers[energyName].checkReal(energyValue, frameName.c_str());
            }
            ++frameNumber;
        }
        if (frameNumber == maxNumEnergyFrames && energyReader->readNextFrame())
        {
            // There would have been at least one more frame!
            checker->disableUnusedEntriesCheck();
        }
    }
    else
    {
        EXPECT_NONFATAL_FAILURE(checker->checkUnusedEntries(), ""); // skip checks on other ranks
    }
}

void checkEnergiesAgainstReferenceData(const std::string&          energyFilename,
                                       const EnergyTermsToCompare& energyTermsToCompare,
                                       TestReferenceChecker*       checker)
{
    checkEnergiesAgainstReferenceData(
            energyFilename, energyTermsToCompare, checker, MaxNumFrames::compareAllFrames());
}

} // namespace test
} // namespace gmx
