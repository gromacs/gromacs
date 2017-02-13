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
 * A basic PME GPU sanity test.
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

//! A basic PME energy conservation test
class PmeTest :
    public MdrunTestFixture //,
    //public ::testing::WithParamInterface<MdrunInputParameters>
{
//SimulationRunner runner_;
//   PmeTest(): runner_(this){};
};

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
                                          //"verlet-buffer-tolerance=-1\n"
                                          //"rlist           = 0.6\n"
                                          //"rcoulomb        = 0.6\n" //FIXME
                                          //"rvdw             = 0.6\n"
                                          "nsteps          = %d\n",
                                          nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "OctaneSandwich"; //"w48k"; // //spc2
    runner_.useTopGroAndNdxFromDatabase(inputFile.c_str());

    EXPECT_EQ(0, runner_.callGrompp());

    std::list<std::string> runParameters;
    runParameters.push_back("-pme cpu");
#if GMX_GPU == GMX_GPU_CUDA
    runParameters.push_back("-pme gpu");
    //FIXME runParameters.push_back("-pme gpu -npme 1 -ntmpi 4");
#endif
    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());

    for (const auto &runIt : runParameters)
    {
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + "_" + runIt + ".edr"); //FIXME!
        // TODO jsut use teh test name generation? also add the options


        // PMECommandLine.addOption("-ntmpi", 1);
        //commandLine.addOption("-pme", codePathIt);

        std::vector<std::string>  temp = splitString(runIt);
        std::vector<const char *> temp2;
        std::transform(temp.begin(), temp.end(), std::back_inserter(temp2), [](const std::string &s){return s.c_str(); });
        CommandLine               commandLine;
        commandLine.initFromArray(temp2);             //TODO betetr way?

        ASSERT_EQ(0, runner_.callMdrun(commandLine)); //FIXME - some stuff already appended inside?


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
                const auto conservedTolerance  = relativeToleranceAsFloatingPoint(conservedEnergy, 4e-3);
                const auto reciprocalTolerance = relativeToleranceAsFloatingPoint(reciprocalEnergy, 1e-2);
                reciprocalChecker.setDefaultTolerance(reciprocalTolerance);
                conservedChecker.setDefaultTolerance(conservedTolerance);
            }
            conservedChecker.checkReal(conservedEnergy, stepNum.c_str());
            reciprocalChecker.checkReal(reciprocalEnergy, stepNum.c_str());
        }
        //compareFrames(std::make_pair(energyReadersByMode[pmeModes[0]]->frame(),
        //                           energyReadersByMode[pmeModes[1]]->frame()),
        //          tolerance);
    }
}

//INSTANTIATE_TEST_CASE_P(What, PmeTest, ::testing::Values("-pme cpu", " -pme gpu", "-pme gpu -npme 1"));
//INSTANTIATE_TEST_F(What, PmeTest, ::testing::Values("-pme cpu", " -pme gpu", "-pme gpu -npme 1"));

}
}
}
