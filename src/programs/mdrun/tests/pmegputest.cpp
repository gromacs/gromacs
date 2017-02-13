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

#include "gromacs/utility/stringutil.h"

#include "energyreader.h"
#include "moduletest.h"

#include <list>
#include <map>

namespace gmx
{
namespace test
{
namespace
{

//! A basic PME GPU test
class PmeGpuTest :
    public test::MdrunTestFixture,
    public ::testing::WithParamInterface<const char *>
{

};

/* Ensure 2 mdruns with CPU and GPU PME produce similar reciprocal and conserved energies. */
TEST_F(PmeGpuTest, ReproducesEnergies)
{
    int         nsteps     = 200;
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

    const std::string inputFile = "spc2";
    runner_.useTopGroAndNdxFromDatabase(inputFile.c_str());

    /* For spc2 Coulomb reciprocal energy is ~7.5, conserved is ~1.3 => abs. tolerance of 3e-4 is OK for both */
    //const FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, 1e-2);

    EXPECT_EQ(0, runner_.callGrompp());


    enum class CodePath //fixme - borrow
    {
        CPU,
        CUDA
    };

    //std::map<CodePath, std::list<std::string> > modesByCodePath;

    //std::vector<std::string> pmeModes;
    //pmeModes.push_back(CodePath::CPU, {""});
//    /pmeModes.push_back(CodePath::CUDA, );
    /*

       std::map<std::string, EnergyFrameReaderPtr> energyReadersByMode;

       for (auto &it : pmeModes)
       {
        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + "_" + it + ".edr");

        CommandLine commandLine;
        // PMECommandLine.addOption("-ntmpi", 1);
        commandLine.addOption("-pme", it);

        ASSERT_EQ(0, runner_.callMdrun(commandLine));

        energyReadersByMode[it] = openEnergyFileToReadFields(runner_.edrFileName_, {{"Coul. recip."}, {"Conserved En."}});
       }

       for (int i = 0; i <= nsteps; i++)
       {
        for (auto &it : pmeModes)
        {
            energyReadersByMode[it]->readNextFrame();
        }

        compareFrames(std::make_pair(energyReadersByMode[pmeModes[0]]->frame(),
                                     energyReadersByMode[pmeModes[1]]->frame()),
                      tolerance);
       }
     */
}

//TODO remove #ifdef __INTEL_COMPILER
//#pragma warning( disable : 177 )
//#endif

}
}
}
