/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Integrator tests
 *
 * Test for CUDA implementation of the Leap-Frog integrator.
 *
 * \todo Connect to the CPU-based version.
 * \todo Prepare for temperature and pressure controlled integrators.
 * \todo Add PBC handling test.
 * \todo Take over coordinates/velocities/forces handling - this infrastructure
 *       will be removed from integrator.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/leapfrog_cuda.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

#if GMX_GPU == GMX_GPU_CUDA

/*! \brief The parameter space for the test.
 *
 * The test will run for all possible combinations of accessible
 * values of the:
 * 1. Number of atoms
 * 2. Timestep
 * 3-5. Velocity components
 * 6-8. Force components
 * 9. Number of steps
 */
typedef std::tuple<int, real, real, real, real, real, real, real, int> IntegratorTestParameters;

/*! \brief Test fixture for integrator.
 *
 *  Creates a system of independent particles exerting constant external forces,
 *  makes several numerical integration timesteps and compares the result
 *  with analytical solution.
 *
 */
class IntegratorTest : public ::testing::TestWithParam<IntegratorTestParameters>
{
    public:
        //! Number of atoms in the system
        int  numAtoms_;
        //! Integration timestep
        real timestep_;

        //! Initial coordinates
        std::vector<RVec> x0_;
        //! Current coordinates
        std::vector<RVec> x_;
        //! Coordinates after integrator update
        std::vector<RVec> xPrime_;
        //! Initial velocities
        std::vector<RVec> v0_;
        //! Current velocities
        std::vector<RVec> v_;
        //! External forces
        std::vector<RVec> f_;
        //! Inverse masses of the particles
        std::vector<real> inverseMasses_;
        //! MD atoms structure in which inverse masses will be passed to the integrator
        t_mdatoms         mdAtoms_;

        /*! \brief Test initialization function.
         *
         * \param[in]  numAtoms  Number of atoms in the system
         * \param[in]  timestep  Integration timestep
         * \param[in]  v0        Initial velocity (same for all particles)
         * \param[in]  f0        External constant force, acting on all particles
         */
        void init(int numAtoms, real timestep, rvec v0, rvec f0)
        {
            numAtoms_    = numAtoms;
            timestep_    = timestep;

            x0_.resize(numAtoms_);
            x_.resize(numAtoms_);
            xPrime_.resize(numAtoms_);
            v0_.resize(numAtoms_);
            v_.resize(numAtoms_);
            f_.resize(numAtoms_);
            inverseMasses_.resize(numAtoms_);
            for (unsigned i = 0; i < x_.size(); i++)
            {
                // Typical PBC box size is tens of nanometers
                x_[i][XX] = (i%21)*1.0;
                x_[i][YY] = 6.5 + (i%13)*(-1.0);
                x_[i][ZZ] = (i%32)*(0.0);

                for (int d = 0; d < DIM; d++)
                {
                    xPrime_[i][d] = 0.0;
                    // Thermal velocity is ~1 nm/ps (|v0| = 1-2 nm/ps)
                    v_[i][d] = v0[d];
                    // TODO Check what value typical MD forces have (now ~ 1 kJ/mol/nm)
                    f_[i][d] = f0[d];

                    x0_[i][d] = x_[i][d];
                    v0_[i][d] = v_[i][d];
                }
                // Atom masses are ~1-100 g/mol
                inverseMasses_[i] = 1.0/(1.0 + i%100);
            }

            mdAtoms_.invmass = inverseMasses_.data();
        }

};

TEST_P(IntegratorTest, SimpleIntegration){
    // Do nothing if it is a CUDA build but there are no CUDA-capable GPUs
    std::string errorMessage;
    if (!canDetectGpus(&errorMessage))
    {
        return;
    }
    int  numAtoms;    // 1. Number of atoms
    real timestep;    // 2. Timestep
    rvec v0;          // 3. Velocity
    rvec f0;          // 4. Force
    int  nStep;       // 5. Number of steps
    std::tie(numAtoms, timestep, v0[XX], v0[YY], v0[ZZ], f0[XX], f0[YY], f0[ZZ], nStep) = GetParam();

    std::string testDescription = formatString("while testing %d atoms for %d timestep (dt = %f, v0=(%f, %f, %f), f0=(%f, %f, %f))",
                                               numAtoms, nStep, timestep,
                                               v0[XX], v0[YY], v0[ZZ],
                                               f0[XX], f0[YY], f0[ZZ]);


    init(numAtoms, timestep, v0, f0);

    std::unique_ptr<LeapFrogCuda> integrator = std::make_unique<LeapFrogCuda>(numAtoms);

    integrator->set(mdAtoms_);
    integrator->copyCoordinatesToGpu((rvec*)(x_.data()));
    integrator->copyVelocitiesToGpu((rvec*)(v_.data()));
    integrator->copyForcesToGpu((rvec*)(f_.data()));

    for (int step = 0; step < nStep; step++)
    {
        integrator->integrate(timestep_);
        integrator->copyCoordinatesFromGpu((rvec*)(xPrime_.data()));
        integrator->copyCoordinatesToGpu((rvec*)(xPrime_.data()));
    }

    integrator->copyCoordinatesFromGpu((rvec*)(xPrime_.data()));
    integrator->copyVelocitiesFromGpu((rvec*)(v_.data()));

    real totalTime = nStep*timestep;
    // TODO For the case of constant force, the numerical scheme is exact and
    //      the only source of errors is floating point arithmetic. Hence,
    //      the tolerance can be calculated.
    FloatingPointTolerance tolerance = absoluteTolerance(nStep*0.000005);

    for (unsigned i = 0; i < x_.size(); i++)
    {
        rvec xAnalytical;
        rvec vAnalytical;
        for (int d = 0; d < DIM; d++)
        {
            // Analytical solution for constant-force particle movement
            xAnalytical[d] = x0_[i][d] + v0_[i][d]*totalTime + 0.5*f_[i][d]*totalTime*totalTime*inverseMasses_[i];
            vAnalytical[d] = v0_[i][d] + f_[i][d]*totalTime*inverseMasses_[i];

            EXPECT_REAL_EQ_TOL(xPrime_[i][d], xAnalytical[d], tolerance)
            << gmx::formatString("Coordinate %d of atom %d is different from analytical solution", d, i)
            << testDescription;

            EXPECT_REAL_EQ_TOL(v_[i][d], vAnalytical[d], tolerance)
            << gmx::formatString("Velocity component %d of atom %d is different from analytical solution", d, i)
            << testDescription;
        }
    }
}

INSTANTIATE_TEST_CASE_P(WithParameters, IntegratorTest,
                            ::testing::Combine(
                                    ::testing::Values(1, 10, 300),    // Number of atoms
                                    ::testing::Values(0.001, 0.0005), // Timestep
                                    ::testing::Values(-2.0, 0.0),     // vx
                                    ::testing::Values( 0.0, 2.0),     // vy
                                    ::testing::Values( 0.0),          // vz
                                    ::testing::Values(-1.0, 0.0),     // fx
                                    ::testing::Values( 0.0, 1.0),     // fy
                                    ::testing::Values( 2.0),          // fz
                                    ::testing::Values(1, 10)          // Number of steps
                                ));
#endif                                                                // GMX_GPU == GMX_GPU_CUDA

}                                                                     // namespace test
}                                                                     // namespace gmx
