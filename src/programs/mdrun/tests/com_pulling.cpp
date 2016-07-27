/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Test for center of mass pulling
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"

#include "energyreader.h"
#include "moduletest.h"

namespace
{

//! Test center of mass pulling by checking potential energy in
//! steepest-descent energy minimization of two water molecules.
class CenterOfMassPullingTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<const char *>
{
    public:
        //! The file name of the MDP file
        std::string                     theMdpFile;
        //! To read the potential energy during energy minimization
        gmx::test::EnergyFrameReaderPtr energy;
        //! run the test
        void runTest()
        {
            runner_.useStringAsMdpFile(theMdpFile);
            runner_.useTopGroAndNdxFromDatabase("spc2");

            // grompp should run without error
            EXPECT_EQ(0, runner_.callGrompp());

            runner_.edrFileName_ = fileManager_.getTemporaryFilePath("spc2-ener.edr");

            // assert that mdrun is finished without error
            ASSERT_EQ(0, runner_.callMdrun());

            energy = gmx::test::openEnergyFileToReadFields(runner_.edrFileName_, {{"Potential"}});

            gmx::test::TestReferenceData    data;
            gmx::test::TestReferenceChecker checker(data.rootChecker());

            // Tolerance is very lax, because steepest descent energy minimisation uses random numbers.
            // Due to the virtual lack of local minima here, energy minimisation is expected to
            // converge always in the same manner.
            checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-2));

            // Check that the potential energy behaves consistently in all frames.
            int i_frame = 0;
            while (energy->readNextFrame())
            {
                checker.checkReal(energy->frame().at("Potential"),
                                  std::string("PotentialEnergyFrame"+std::to_string(i_frame++)).c_str());
            }
        }
};

//! Helper typedef for naming test cases like sentences
typedef CenterOfMassPullingTest CenterOfMassPulling;

/* Ensure grompp and mdrun run and potential energy converges in an expected way.*/
TEST_F(CenterOfMassPulling, ReproducesEnergiesInEnergyMinimization)
{
    theMdpFile = gmx::formatString("integrator = steep\n"
                                   "nsteps              = 5 \n"
                                   "pull                = yes                \n"
                                   "pull-ngroups        = 2                  \n"
                                   "pull-ncoords        = 1                  \n"
                                   "pull-group1-name    = FirstWaterMolecule \n"
                                   "pull-group2-name    = SecondWaterMolecule\n"
                                   "pull-coord1-k       = 100                \n"
                                   "pull-coord1-groups  = 1 2                \n");
    runTest();
}

/* Wether FristWaterMolecule is pulled towards SecondWatermolecule or vice versa
 * may not make a difference in potential energies.*/
TEST_F(CenterOfMassPulling, ReproducesEnergiesInEnergyMinimizationSwappingPulledGroups)
{
    theMdpFile = gmx::formatString("integrator = steep\n"
                                   "nsteps              = 5 \n"
                                   "pull                = yes                \n"
                                   "pull-ngroups        = 2                  \n"
                                   "pull-ncoords        = 1                  \n"
                                   "pull-group1-name    = SecondWaterMolecule\n"
                                   "pull-group2-name    = FirstWaterMolecule \n"
                                   "pull-coord1-k       = 100                \n"
                                   "pull-coord1-groups  = 1 2                \n");
    runTest();
}

#ifdef __INTEL_COMPILER
#pragma warning( disable : 177 )
#endif

} // namespace
