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
 * Tests molecular dynamics with a density fitting potential.
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <math.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxpreprocess/grompp.h"

#include "energyreader.h"
#include "moduletest.h"
#include "trajectoryreader.h"


namespace gmx
{
namespace test
{


class DensfitIntegrationTest : public MdrunTestFixture
{
    protected:
        DensfitIntegrationTest();
        ~DensfitIntegrationTest();

    public:
        //! \brief Compares all frames of two trajectories
        bool trajectoriesAreEqualWithinTolerance(const std::string &fnA, const std::string &fnB, double absTolerance) const;

        //! \brief Read the correlation coefficient from  the densfit output file
        real correlationFromXvgFile(const std::string &fn) const;
};


DensfitIntegrationTest::DensfitIntegrationTest()
{
}

DensfitIntegrationTest::~DensfitIntegrationTest()
{
}


real DensfitIntegrationTest::correlationFromXvgFile(const std::string &fn) const
{
    real   firstTimeToRead = -1.0;
    real   lastTimeToRead  = 1.0;
    int    nset            = 0, n = 0;
    real  *t               = nullptr, dt = -99;

    real **val = read_xvg_time(fn.c_str(), true,
                               true, firstTimeToRead,
                               true, lastTimeToRead,
                               1, &nset, &n, &dt, &t);
    printf("Read %d sets of %d points, dt = %g\n\n", nset, n, dt);

    real cc = val[3][0];

    return cc;
}


bool DensfitIntegrationTest::trajectoriesAreEqualWithinTolerance(const std::string &fnA, const std::string &fnB, double absTolerance) const
{
    fprintf(stderr, "\nComparing frames of %s and %s\n", fnA.c_str(), fnB.c_str());
    TrajectoryFrameReader tfr1(fnA);
    TrajectoryFrameReader tfr2(fnB);

    // Loop over both trajectory files and compare the contained data frame-by-frame
    while (tfr1.readNextFrame() && tfr2.readNextFrame() )
    {
        const auto &fr1 = tfr1.frame();
        const auto &fr2 = tfr2.frame();
        fprintf(stderr, "\nRead trajectory frames '%s' and '%s'\n", fr1.getFrameName().c_str(), fr2.getFrameName().c_str());

        auto PairOfFrames = std::pair<TrajectoryFrame, TrajectoryFrame> (fr1, fr2);
        compareFrames(PairOfFrames, absoluteTolerance(absTolerance));
    }

    return true;
}


/* This basic test ensures that the mdrun with an added density fitting potential
 * run, and that the related .mdp input parameters are understood.
 */
TEST_F(DensfitIntegrationTest, DensfitForcesAreCorrect)
{
    std::string name = "1AKE";
    runner_.useTopGroAndNdxFromDatabase(name.c_str());
    auto        inputMapFileName             = fileManager_.getInputFilePath("4AKE-aligned-Calpha.ccp4");
    auto        forcesReferenceFileName      = fileManager_.getInputFilePath("1AKE_with_4AKE_map_forces.trr");
    auto        outputMapFileName            = fileManager_.getTemporaryFilePath(".ccp4");
    auto        densfitOutFileName           = fileManager_.getTemporaryFilePath(".xvg");
    runner_.fullPrecisionTrajectoryFileName_ = fileManager_.getTemporaryFilePath(".trr");

    // Prepare .tpr file:
    std::string mdpString =
        "dt                       = 0.004e-10\n"
        "nsteps                   = 0        \n"
        "nstfout                  = 1        \n"
        "cutoff-scheme            = Verlet   \n"
        "coulombtype              = pme      \n"
        "fourierspacing           = 0.135    \n"
        "ewald_rtol               = 1e-06    \n"
        "constraints              = all-bonds\n"
        "continuation             = yes      \n"
        "DensityFitting           = yes      \n"
        "densfit-potential        = yes      \n"
        "reference-density        = " + inputMapFileName + "\n"
        "spread-group             = C-alpha  \n"
        "densfit-time             = 0        \n"
        "densfit-sigma            = 0.3      \n"
        "densfit-k                = 10000e10 \n";

    CommandLine gromppCaller;
    runner_.useStringAsMdpFile(mdpString);

    ASSERT_EQ(0, runner_.callGrompp(gromppCaller));

    // Now run the .tpr file:
    CommandLine mdrunCaller;
    mdrunCaller.addOption("-d", densfitOutFileName);
    mdrunCaller.addOption("-mo", outputMapFileName);
    mdrunCaller.addOption("-noconfout");

    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

    EXPECT_REAL_EQ_TOL(0.58290, correlationFromXvgFile(densfitOutFileName), absoluteTolerance(1e-6));

    // densfit forces were multiplied by 1e10 before (see mdp input file), therefore the factor 1e10 below
#if GMX_DOUBLE
    ASSERT_EQ(true, trajectoriesAreEqualWithinTolerance(forcesReferenceFileName, runner_.fullPrecisionTrajectoryFileName_, 0.00001e10));
#else
    ASSERT_EQ(true, trajectoriesAreEqualWithinTolerance(forcesReferenceFileName, runner_.fullPrecisionTrajectoryFileName_, 0.00050e10));
#endif
}


} // namespace test
} // namespace gmx
