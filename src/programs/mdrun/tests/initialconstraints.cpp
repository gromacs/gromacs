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
 * This implements basic initial constrains test (using single-rank mdrun).
 * It runs the input system for 1 step (no continuation), and compares the total energy.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "gromacs/utility/stringutil.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! A convenience typedef
using InitialConstraintsTest = MdrunTestFixture;

TEST_F(InitialConstraintsTest, Works)
{
    const int nsteps = 1;
    const float timestep = 0.001;
    const std::string theMdpFile = formatString("nstcalcenergy   = 1\n"
                                                "nstenergy       = 1\n"
                                                "continuation    = no\n"
                                                "nsteps          = %d\n"
                                                "dt              = %f\n",
                                                nsteps, timestep);

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "spc-and-methanol";
    runner_.useTopGroAndNdxFromDatabase(inputFile.c_str());
    EXPECT_EQ(0, runner_.callGrompp());

    runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + ".edr");
    ASSERT_EQ(0, runner_.callMdrun());

    EnergyFrameReaderPtr energyReader = openEnergyFileToReadFields(runner_.edrFileName_, {{"Total Energy"}});
    std::vector<real> totalEnergies;
    for (int i = 0; i <= nsteps; i++)
    {
        EnergyFrame frame = energyReader->frame();
        totalEnergies.push_back(frame.at("Total Energy"));
    }
    EXPECT_REAL_EQ_TOL(totalEnergies[0], totalEnergies[1], relativeToleranceAsFloatingPoint(totalEnergies[0], 1e-5));
}

}
}
}
