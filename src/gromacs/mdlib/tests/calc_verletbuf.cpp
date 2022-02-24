/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * \brief Tests for the Verlet buffer calculation algorithm.
 *
 * \author Berk Hess <hess@kth.se>
 */
#include "gmxpre.h"

#include "gromacs/mdlib/calc_verletbuf.h"

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace
{

class VerletBufferConstraintTest : public ::testing::Test
{
};

/* This test covers the displacement correction for constrained atoms.
 * This test does not check exact values, but rather checks that the MSD
 * estimate for a constrained atom is smaller than that of a free atom
 * and checks that the MSD is not smaller and also not much larger
 * than the maximum of the exact value for rotational MSD beyond
 * the location of the maximum. Furthermore, we check that the MSD estimate
 * never decreases, as this is a requirement for the Verlet buffer size
 * estimation. Together these criteria provide tight margins on
 * the shape and values of the estimate.
 *
 * Additionally we check the 3D MSD for the COM of the two atoms.
 */
TEST_F(VerletBufferConstraintTest, EqualMasses)
{
    // The location and value of the MSD maximum for the exact displacement
    // is described in the source file. We need to divide the maximum given
    // there by 2, since sigma2 is per DOF for the 2 DOF constraint rotation.
    const real sigma2RelMaxLocation = 4.5119;
    const real sigma2RelMaxValue    = 2.5695 / 2;

    // Our max of our current estimate is 3% above the exact value.
    const real sigma2RelMaxMargin = 1.04;

    // The exact parameter values here don't actually matter.
    real mass = 10;
    real arm  = 0.1;

    atom_nonbonded_kinetic_prop_t prop;
    prop.mass     = mass;
    prop.type     = -1;
    prop.q        = 0;
    prop.bConstr  = TRUE;
    prop.con_mass = mass;
    prop.con_len  = 2 * arm;

    // We scan a range of rotation distributions by scanning over T.
    int  numPointsBeforeMax = 0;
    int  numPointsAfterMax  = 0;
    real sigma2_2d_prev     = 0;
    for (int i = 0; i <= 200; i++)
    {
        real ktFac = i * 0.01;
        // The rotational displacement is Gaussian with a sigma^2 of:
        real sigma2_rot = ktFac / (2 * mass);

        // Get the estimate for the Cartesian displacement.
        real sigma2_2d, sigma2_3d;
        constrained_atom_sigma2(ktFac, &prop, &sigma2_2d, &sigma2_3d);

        // Check we are not decreasing sigma2_2d
        EXPECT_EQ(std::max(sigma2_2d_prev, sigma2_2d), sigma2_2d);
        // Check that sigma2_2d is not larger than sigma2 for free motion.
        EXPECT_EQ(std::min(sigma2_rot, sigma2_2d), sigma2_2d);

        // Check that we don't underestimate sigma2_rot beyond the real maximum
        // and that our overestimate is tight.
        real sigma2Rel = sigma2_rot / gmx::square(arm);
        if (sigma2Rel >= sigma2RelMaxLocation)
        {
            EXPECT_EQ(std::max(sigma2_2d, sigma2RelMaxValue * gmx::square(arm)), sigma2_2d);
            EXPECT_EQ(std::min(sigma2_2d, sigma2RelMaxMargin * sigma2RelMaxValue * gmx::square(arm)),
                      sigma2_2d);

            numPointsAfterMax++;
        }
        else
        {
            numPointsBeforeMax++;
        }

        // Also check sigma2 for the COM of the two atoms
        EXPECT_REAL_EQ_TOL(sigma2_rot, sigma2_3d, test::defaultRealTolerance());
    }

    GMX_RELEASE_ASSERT(
            numPointsBeforeMax >= 20 && numPointsAfterMax >= 20,
            "This test only provides full coverage when we test a sufficient number of points "
            "before and after the location of the maximum value for the exact formula.");
}

} // namespace

} // namespace gmx
