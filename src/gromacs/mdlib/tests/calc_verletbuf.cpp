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

#include <cstdlib>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

// Effective density test with 8 cells with 2 atoms in both a tight and large unit cell
TEST(EffectiveAtomDensity, VolumeIndependence)
{
    const std::vector<RVec> coordinates      = { { 1, 1, 1 }, { 1, 1, 1 }, { 1, 3, 1 }, { 1, 3, 1 },
                                            { 3, 1, 1 }, { 3, 1, 1 }, { 3, 3, 1 }, { 3, 3, 1 },
                                            { 1, 1, 3 }, { 1, 1, 3 }, { 1, 3, 3 }, { 1, 3, 3 },
                                            { 3, 1, 3 }, { 3, 1, 3 }, { 3, 3, 3 }, { 3, 3, 3 } };
    const matrix            tightBox         = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    const matrix            largeBox         = { { 40, 0, 0 }, { 0, 40, 0 }, { 0, 0, 4 } };
    const real              cutoff           = 2;
    const real              referenceDensity = coordinates.size() / real(4 * 4 * 4);

    const real tightBoxDensity = computeEffectiveAtomDensity(coordinates, tightBox, cutoff, MPI_COMM_NULL);
    EXPECT_FLOAT_EQ(tightBoxDensity, referenceDensity);

    const real largeBoxDensity = computeEffectiveAtomDensity(coordinates, largeBox, cutoff, MPI_COMM_NULL);
    EXPECT_FLOAT_EQ(largeBoxDensity, referenceDensity);
}

// Effective density test with 1 cell with 3 atoms, 7 cells with 1 atom
TEST(EffectiveAtomDensity, WeightingWorks)
{
    const std::vector<RVec> coordinates = { { 1, 1, 1 }, { 1, 1, 1 }, { 1, 1, 1 }, { 1, 3, 1 },
                                            { 3, 1, 1 }, { 3, 3, 1 }, { 1, 1, 3 }, { 1, 3, 3 },
                                            { 3, 1, 3 }, { 3, 3, 3 } };
    const matrix            box         = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    const real              cutoff      = 2;
    const real referenceDensity         = (3 * 3 + 7 * 1) / (coordinates.size() * real(2 * 2 * 2));

    const real density = computeEffectiveAtomDensity(coordinates, box, cutoff, MPI_COMM_NULL);
    EXPECT_FLOAT_EQ(density, referenceDensity);
}

TEST(AtomNonbondedAndKineticProperties, IsAccurate)
{
    const real c_resolution = 0.1;

    std::vector<real> invMasses = { 8.0, 10.149, 20.051 };
    std::vector<real> charges   = { -2.0, 0.149, -0.951 };

    for (Index i = 0; i < gmx::ssize(invMasses); i++)
    {
        AtomNonbondedAndKineticProperties props({ c_resolution, c_resolution, c_resolution });
        props.setMassTypeCharge(1 / invMasses[i], 0, charges[i]);

        EXPECT_LT(std::abs(props.invMass() - invMasses[i]), 0.5 * c_resolution);
        EXPECT_LT(std::abs(props.charge() - charges[i]), 0.5 * c_resolution);
    }
}

TEST(AtomNonbondedAndKineticProperties, ConstraintsWork)
{
    const real c_resolution = 0.1;

    std::vector<real> invMasses = { 8.0, 10.149, 20.051 };
    std::vector<real> lengths   = { 0.0, 3.2499, 4.9501 };

    for (Index i = 0; i < gmx::ssize(invMasses); i++)
    {
        AtomNonbondedAndKineticProperties props({ c_resolution, c_resolution, c_resolution });
        props.setMassTypeCharge(1 / invMasses[i], 0, 0);
        if (i > 0)
        {
            // Add a first constraint with low mass
            props.addConstraint(1 / invMasses[i], 1.0);
            // Add a second one with higher mass that is the one that should be used
            props.addConstraint(2 / invMasses[i], lengths[i]);
        }

        EXPECT_EQ(props.hasConstraint(), i > 0);
        if (props.hasConstraint())
        {
            EXPECT_LT(std::abs(props.constraintInvMass() - 0.5 * invMasses[i]), 0.5 * c_resolution);
            EXPECT_LT(std::abs(props.constraintLength() - lengths[i]), 0.5 * c_resolution);
        }
    }
}

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

    AtomNonbondedAndKineticPropertiesResolutions resolutions;

    resolutions.invMassResolution          = 1 / mass * 0.01;
    resolutions.chargeResolution           = 1;
    resolutions.constraintLengthResolution = 2 * arm * 0.01;

    AtomNonbondedAndKineticProperties prop(resolutions);
    prop.setMassTypeCharge(mass, -1, 0);
    prop.addConstraint(mass, 2 * arm);

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
        constrained_atom_sigma2(ktFac, prop, &sigma2_2d, &sigma2_3d);

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

// Issue #5002
TEST(EffectiveAtomDensity, LargeValuesHandledWell)
{
    const std::vector<RVec> coordinates = { { 13.132, -8.229, -2.700 } };
    const matrix            box         = { { 6.2, 0, 0 }, { 0, 6.2, 0 }, { 0, 0, 6.2 } };
    const real              cutoff      = 1;
    const real referenceDensity         = (1) / (coordinates.size() * gmx::power3<real>(6.2 / 6));

    const real density = computeEffectiveAtomDensity(coordinates, box, cutoff, MPI_COMM_NULL);
    EXPECT_FLOAT_EQ(density, referenceDensity);
}

} // namespace
} // namespace test
} // namespace gmx
