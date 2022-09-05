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
 * This implements TPR reading tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/tpr.h"

#include <gtest/gtest.h>

#include "gromacs/tools/convert_tpr.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"

#include "programs/mdrun/tests/moduletest.h"

#include "nblib/gmxcalculatorcpu.h"

#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{

//! Test fixture for running the utils needed for reading TPR files from simulation database
class TprReaderTest : public gmx::test::MdrunTestFixture
{
public:
    //! Returns a TprReader object populated with the given files in the simulation database
    TprReader makeTPRfromSimulationDatabase(const std::string& filename)
    {
        gmx::test::CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(filename);
        auto mdpFieldValues = gmx::test::prepareMdpFieldValues(filename.c_str(), "md", "no", "no");
        mdpFieldValues["init-lambda-state"] = "0";
        mdpFieldValues["constraints"]       = "none"; // constraints not supported in NBLIB
        runner_.useStringAsMdpFile(gmx::test::prepareMdpFileContents(mdpFieldValues));
        runner_.callGrompp(caller);

        return TprReader(runner_.tprFileName_);
    }
};

TEST_F(TprReaderTest, SimDBTprIsCreated)
{
    EXPECT_NO_THROW(makeTPRfromSimulationDatabase("argon12"));
}

TEST_F(TprReaderTest, Spc2Reads)
{
    TprReader tprReader = makeTPRfromSimulationDatabase("spc2");

    EXPECT_EQ(tprReader.coordinates_.size(), 6);
}

TEST_F(TprReaderTest, ArgonImportedDataIsCorrect)
{
    TprReader      tprReader = makeTPRfromSimulationDatabase("argon12");
    RefDataChecker dataImportTest(5e-5);

    // check number of coordinates in the system
    EXPECT_EQ(tprReader.coordinates_.size(), 12);

    // check coordinates
    dataImportTest.testArrays<Vec3>(tprReader.coordinates_, "coordinates");

    // check number of velocities in the system
    EXPECT_EQ(tprReader.velocities_.size(), 12);

    // check velocities
    dataImportTest.testArrays<Vec3>(tprReader.velocities_, "velocities");

    // check box
    EXPECT_TRUE(tprReader.getBox() == Box(6.05449));

    // check non-bonded params
    dataImportTest.testArrays<real>(tprReader.nonbondedParameters_, "nbparams");

    // check exclusion elements
    dataImportTest.testArrays<int>(tprReader.exclusionListElements_, "exclusion elements");

    // check exclusion ranges
    dataImportTest.testArrays<int>(tprReader.exclusionListRanges_, "exclusion ranges");
}

TEST_F(TprReaderTest, FCfromTprDataWorks)
{
    auto options        = NBKernelOptions();
    options.nbnxmSimd   = SimdKernels::SimdNo;
    options.coulombType = CoulombType::Cutoff;

    TprReader tprReader = makeTPRfromSimulationDatabase("argon12");

    auto forceCalculatorTpr = setupGmxForceCalculatorCpu(tprReader, options);

    std::vector<Vec3> forceOutputTpr(tprReader.coordinates_.size(), Vec3{ 0, 0, 0 });

    forceCalculatorTpr->updatePairlist(tprReader.coordinates_, tprReader.getBox());
    forceCalculatorTpr->compute(tprReader.coordinates_, tprReader.getBox(), forceOutputTpr);

    RefDataChecker forcesOutputTest;
    forcesOutputTest.testArrays<Vec3>(forceOutputTpr, "Argon forces");
}

} // namespace test

} // namespace nblib
