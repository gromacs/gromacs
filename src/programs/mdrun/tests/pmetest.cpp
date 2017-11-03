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

#include <map>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/cstringutil.h"
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
class PmeTest : public MdrunTestFixture
{
    public:
        //! Before any test is run, work out whether any compatible GPUs exist.
        static void SetUpTestCase();
        //! Store whether any compatible GPUs exist.
        static bool s_hasCompatibleCudaGpus;
};

bool PmeTest::s_hasCompatibleCudaGpus = false;

void PmeTest::SetUpTestCase()
{
    gmx_gpu_info_t gpuInfo {};
    char           detection_error[STRLEN];
    GMX_UNUSED_VALUE(detection_error); //TODO
    // It would be nicer to do this detection once and have mdrun
    // re-use it, but this is OK. Note that this also caters for when
    // there is no GPU support in the build.
    if (GMX_GPU == GMX_GPU_CUDA &&
        (detect_gpus(&gpuInfo, detection_error) >= 0) &&
        gpuInfo.n_dev_compatible > 0)
    {
        s_hasCompatibleCudaGpus = true;
    }
    free_gpu_info(&gpuInfo);
}

TEST_F(PmeTest, ReproducesEnergies)
{
    const int   nsteps     = 20;
    std::string theMdpFile = formatString("coulombtype     = PME\n"
                                          "nstcalcenergy   = 1\n"
                                          "nstenergy       = 1\n"
                                          "pme-order       = 4\n"
                                          "nsteps          = %d\n",
                                          nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "spc-and-methanol";
    runner_.useTopGroAndNdxFromDatabase(inputFile.c_str());

    EXPECT_EQ(0, runner_.callGrompp());

    //TODO test all proper/improper combinations in more thorough way?
    std::map < std::string, std::vector < const char *>> runModes;
    runModes["PmeOnCpu"] = {"-pme", "cpu"};
    runModes["PmeAuto"]  = {"-pme", "auto"};
    // TODO uncomment this when functionality gets activated.
    //runModes["PmeOnGpuFftOnCpu"] = {"-pme", "gpu", "-pmefft", "cpu"};
    runModes["PmeOnGpuFftOnGpu"] = {"-pme", "gpu", "-pmefft", "gpu"};
    runModes["PmeOnGpuFftAuto"] = {"-pme", "gpu", "-pmefft", "auto"};
    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());

    for (const auto &mode : runModes)
    {
        auto modeTargetsGpus = (mode.first.find("Gpu") != std::string::npos);
        if (modeTargetsGpus && !s_hasCompatibleCudaGpus)
        {
            // This run mode will cause a fatal error from mdrun when
            // it can't find GPUs, which is not something we're trying
            // to test here.
            continue;
        }

        runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + "_" + mode.first + ".edr");

        CommandLine commandLine(mode.second);
        commandLine.append("-notunepme"); // for reciprocal energy reproducibility
        ASSERT_EQ(0, runner_.callMdrun(commandLine));

        auto energyReader      = openEnergyFileToReadFields(runner_.edrFileName_, {"Coul. recip.", "Total Energy", "Kinetic En."});
        auto conservedChecker  = rootChecker.checkCompound("Energy", "Conserved");
        auto reciprocalChecker = rootChecker.checkCompound("Energy", "Reciprocal");
        for (int i = 0; i <= nsteps; i++)
        {
            EnergyFrame frame            = energyReader->frame();
            std::string stepNum          = gmx::formatString("%d", i);
            const real  conservedEnergy  = frame.at("Total Energy");
            const real  reciprocalEnergy = frame.at("Coul. recip.");
            if (i == 0)
            {
                // use first step values as references for tolerance
                const real startingKineticEnergy = frame.at("Kinetic En.");
                const auto conservedTolerance    = relativeToleranceAsFloatingPoint(startingKineticEnergy, 2e-5);
                const auto reciprocalTolerance   = relativeToleranceAsFloatingPoint(reciprocalEnergy, 2e-5);
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
