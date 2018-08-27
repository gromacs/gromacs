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
 * Tests whether the class MetropolisStepMehlig is working as it should.
 *
 * \author Sebastian Wingbermuehle
 */
#include "gmxpre.h"

#include "gromacs/hybridMCMD/metropolis.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/topology/idef.h"

namespace gmx
{
namespace test
{

class MetropolisStepTest : public ::testing::Test
{
    public:
        MetropolisStepTest()
        {
        }
        ~MetropolisStepTest()
        {
        }
        /*! \brief Test if the Metropolis criterion is evaluated correctly and the correct boolean is returned */
        void runTest(IMetropolisStep &metropolisStep, IMetropolisStepVelocities &metropolisStepVelocities)
        {
            // set initial potential energy
            enerd_.term[F_EPOT] = -1;
            accepted_           = metropolisStep.accept(0, &enerd_);
            EXPECT_TRUE(accepted_);

            // 1) total energy decreases => accept
            metropolisStepVelocities.setInitialKineticEnergy(5);
            enerd_.term[F_EPOT] = -3;
            enerd_.term[F_EKIN] = 4;
            accepted_           = metropolisStep.accept(1, &enerd_);
            EXPECT_TRUE(accepted_);

            // 2) increase in total energy that is rejected with this random seed (test for proper management of energy values)
            metropolisStepVelocities.setInitialKineticEnergy(1);
            enerd_.term[F_EPOT] = -1;
            enerd_.term[F_EKIN] = 3;
            accepted_           = metropolisStep.accept(1, &enerd_);
            EXPECT_FALSE(accepted_);

            // 3) slight increase in total energy => accepted with this random seed
            metropolisStepVelocities.setInitialKineticEnergy(2);
            enerd_.term[F_EPOT] = -1;
            enerd_.term[F_EKIN] = 3;
            accepted_           = metropolisStep.accept(1, &enerd_);
            EXPECT_TRUE(accepted_);

            // 4) total energy increases (beyond PROBABILITYCUTOFF = -100) => reject regardless of random seed
            metropolisStepVelocities.setInitialKineticEnergy(2);
            enerd_.term[F_EPOT] = 200;
            enerd_.term[F_EKIN] = 200;
            accepted_           = metropolisStep.accept(1, &enerd_);
            EXPECT_FALSE(accepted_);
        }

        gmx_enerdata_t  enerd_;
        bool            accepted_;
};

TEST_F(MetropolisStepTest, MetropolisStepMehligWorks)
{
    MetropolisStepMehlig metropolisStepMehlig(300, 300, 1, 0, false);
    runTest(metropolisStepMehlig, metropolisStepMehlig);
}

} // namespace
} // namespace
