/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * This implements basic PME sanity test (using single-rank mdrun).
 * It runs the input system with PME for several steps (on CPU and GPU, if available),
 * and checks the reciprocal and conserved energies.
 * TODO: implement multi-rank tests as well.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <algorithm>
#include <list>

#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! A basic PME runner
using PmeTest = MdrunTestFixture;

TEST_F(PmeTest, ReproducesEnergies)
{
    int         nsteps     = 40;
    std::string theMdpFile = formatString("coulombtype     = PME\n"
                                          "nstcalcenergy   = 1\n"
                                          "nstenergy       = 1\n"
                                          "pme-order       = 4\n"
                                          "tcoupl          = v-rescale\n"
                                          "tau_t           = 0.1\n"
                                          "ref_t           = 300\n"
                                          "tc_grps         = system\n"
                                          "nsteps          = %d\n",
                                          nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "OctaneSandwich";
    runner_.useTopGroAndNdxFromDatabase(inputFile.c_str());

    EXPECT_EQ(0, runner_.callGrompp());

    std::list<std::string> runModes;
    runModes.push_back("-pme cpu");
#if GMX_GPU == GMX_GPU_CUDA
    runModes.push_back("-pme gpu");
#endif
    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());

    for (const auto &mode : runModes)
    {
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + "_" + mode + ".edr");

        // The following 4 lines just turn std::string into CommandLine - is there a better way?
        std::vector<std::string>  argList = splitString(mode);
        std::vector<const char *> argPtrList;
        std::transform(argList.begin(), argList.end(), std::back_inserter(argPtrList), [](const std::string &s){return s.c_str(); });
        CommandLine               commandLine(argPtrList);

        ASSERT_EQ(0, runner_.callMdrun(commandLine));

        EnergyFrameReaderPtr energyReader      = openEnergyFileToReadFields(runner_.edrFileName_, {{"Coul. recip."}, {"Conserved En."}});
        TestReferenceChecker conservedChecker  = rootChecker.checkCompound("Energy", "Conserved");
        TestReferenceChecker reciprocalChecker = rootChecker.checkCompound("Energy", "Reciprocal");
        for (int i = 0; i < nsteps; i++)
        {
            EnergyFrame frame            = energyReader->frame();
            std::string stepNum          = gmx::formatString("%d", i);
            const real  conservedEnergy  = frame.at("Conserved En.");
            const real  reciprocalEnergy = frame.at("Coul. recip.");
            if (i == 0)
            {
                // use first step values as references for tolerance
                const auto conservedTolerance  = relativeToleranceAsFloatingPoint(conservedEnergy, 4e-3);
                const auto reciprocalTolerance = relativeToleranceAsFloatingPoint(reciprocalEnergy, 4e-3);
                reciprocalChecker.setDefaultTolerance(reciprocalTolerance);
                conservedChecker.setDefaultTolerance(conservedTolerance);
            }
            conservedChecker.checkReal(conservedEnergy, stepNum.c_str());
            reciprocalChecker.checkReal(reciprocalEnergy, stepNum.c_str());
        }
    }
}

}
}
}
