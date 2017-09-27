/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests molecular dynamics with a density fitting potential.
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <math.h>

#include "gromacs/gmxpreprocess/grompp.h"

#include "energyreader.h"
#include "moduletest.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{

//! Field in .edr file that contains the correlation coefficient c.c.
const std::vector<std::string> &nameOfcorrelationCoefficientInEdrFile = { "densfit c.c." };


class DensfitIntegrationTest : public MdrunTestFixture
{
    protected:
        DensfitIntegrationTest();
        ~DensfitIntegrationTest();

    public:
        /*! \brief Convenience method to get the first energy entry from the .edr file
         *
         * Any following frames of the .edr file are ignored!
         */
        real getFirstEnergyValue(const std::string fn, const std::vector<std::string> &namesOfRequiredEnergyFields) const;

        //! \brief Compares all frames of two trajectories
        gmx_bool compareTrajectories(const std::string fnA, const std::string fnB, double absTolerance) const;

        //! \brief Compares energy on a frame-by-frame basis for two energy files
        gmx_bool compareEnergies(const std::string fnA, const std::string fnB, double absTolerance) const;

};


DensfitIntegrationTest::DensfitIntegrationTest()
{
}

DensfitIntegrationTest::~DensfitIntegrationTest()
{
}

real DensfitIntegrationTest::getFirstEnergyValue(const std::string fn, const std::vector<std::string> &namesOfRequiredEnergyFields) const
{
    auto E   = 0.0;
    auto efr = openEnergyFileToReadFields(fn, namesOfRequiredEnergyFields );
    if (efr->readNextFrame())
    {
        auto fr = efr->frame();
        for (auto &element : namesOfRequiredEnergyFields)
        {
            E += fr.at(element);
        }
    }
    else
    {
        // What TODO if we cannot read the .edr file?
    }
    return E;
}


gmx_bool DensfitIntegrationTest::compareEnergies(const std::string fnA, const std::string fnB, double absTolerance) const
{
    auto efrA = openEnergyFileToReadFields(fnA, nameOfcorrelationCoefficientInEdrFile);
    auto efrB = openEnergyFileToReadFields(fnB, nameOfcorrelationCoefficientInEdrFile);

    // Loop over both trajectory files and compare the contained data frame-by-frame
    while (efrA->readNextFrame() && efrB->readNextFrame() )
    {
        auto ccA   = 0.0; // correlation coefficient A
        auto ccB   = 0.0;
        auto frPme = efrA->frame();
        auto frFmm = efrB->frame();
        for (auto &element : nameOfcorrelationCoefficientInEdrFile)
        {
            ccA += frPme.at(element);
            ccB += frFmm.at(element);
        }
        fprintf(stderr, "\nRead .edr frame A '%s' with c.c. = %g and frame B '%s' with c.c. = %g (difference %g)\n",
                frPme.getFrameName().c_str(), ccA,
                frFmm.getFrameName().c_str(), ccB,
                fabs(ccA-ccB));

        EXPECT_NEAR(ccA, ccB, absTolerance);
    }

    return 0;
}


gmx_bool DensfitIntegrationTest::compareTrajectories(const std::string fnA, const std::string fnB, double absTolerance) const
{
    fprintf(stderr, "\nComparing frames of %s and %s\n", fnA.c_str(), fnB.c_str());
    auto tfr1(new TrajectoryFrameReader(fnA));
    auto tfr2(new TrajectoryFrameReader(fnB));

    // Loop over both trajectory files and compare the contained data frame-by-frame
    while (tfr1->readNextFrame() && tfr2->readNextFrame() )
    {
        auto fr1 = tfr1->frame();
        auto fr2 = tfr2->frame();
        fprintf(stderr, "\nRead trajectory frames '%s' and '%s'\n", fr1.getFrameName().c_str(), fr2.getFrameName().c_str());

        auto PairOfFrames = std::pair<TrajectoryFrame, TrajectoryFrame> (fr1, fr2);
        compareFrames(PairOfFrames, absoluteTolerance(absTolerance));
    }

    return 0;
}


/* This basic test ensures that the mdrun with an added density fitting potential
 * run, and that the related .mdp input parameters are understood.
 */
TEST_F(DensfitIntegrationTest, DensfitForcesAreCorrect)
{
    // Prepare .tpr file:
    std::string name = "1AKE";
    runner_.topFileName_                     = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    runner_.groFileName_                     = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    runner_.mdpInputFileName_                = fileManager_.getInputFilePath((std::string(name) + ".mdp").c_str());
    runner_.fullPrecisionTrajectoryFileName_ = fileManager_.getTemporaryFilePath((name + ".trr").c_str());
    auto        mapFileName                  = fileManager_.getInputFilePath("4AKE.ccp4");
    auto        densfitOutFileName           = fileManager_.getTemporaryFilePath((name + "_densfit.xvg"));
    auto        edrReferenceFileName         = fileManager_.getInputFilePath("1AKE_with_4AKE_map_forces.trr");

    CommandLine gromppCaller;
    gromppCaller.append("grompp");
    gromppCaller.addOption("-f", runner_.mdpInputFileName_);
    gromppCaller.addOption("-p", runner_.topFileName_);
    gromppCaller.addOption("-c", runner_.groFileName_);
    gromppCaller.addOption("-mi", mapFileName);

    gromppCaller.addOption("-po", runner_.mdpOutputFileName_);
    gromppCaller.addOption("-o", runner_.tprFileName_);

    EXPECT_EQ(0, gmx_grompp(gromppCaller.argc(), gromppCaller.argv()));

    // Now run the .tpr file:
    CommandLine mdrunCaller;
    mdrunCaller.append("mdrun");
    mdrunCaller.append("-v");
    mdrunCaller.addOption("-d", densfitOutFileName);

    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

    // ... and compare the output
    auto cc = getFirstEnergyValue(runner_.edrFileName_, nameOfcorrelationCoefficientInEdrFile );
    EXPECT_REAL_EQ_TOL(0.6914835, cc, absoluteTolerance(1e-6));

    ASSERT_EQ(0, compareTrajectories(edrReferenceFileName, runner_.fullPrecisionTrajectoryFileName_, 0.01));
}


} // namespace test
} // namespace gmx
