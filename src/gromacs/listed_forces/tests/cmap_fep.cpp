/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Tests for CMAP free energy (lambda) interpolation in cmap_dihs().
 *
 * Verifies that energy and forces are linearly interpolated between A and B
 * state CMAP grids as a function of lambda, and that dH/dlambda equals
 * the energy difference E_B - E_A.
 *
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include <cmath>

#include <array>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Build a non-flat CMAP grid with random values and numerically derived derivatives
static CmapGrid buildRandomGrid(int gridExtent, uint32_t seed)
{
    std::mt19937                         rng(seed);               // seeded RNG for reproducibility
    std::uniform_real_distribution<real> dist(-10.0, 10.0);       // random energy range
    const real                           dx = 360.0 / gridExtent; // grid spacing in degrees

    CmapGrid grid(gridExtent, gridExtent); // allocate grid

    for (int i = 0; i < gridExtent; i++) // fill grid values first
    {
        for (int j = 0; j < gridExtent; j++)
        {
            grid(i, j)[0] = dist(rng); // random energy value
        }
    }

    for (int i = 0; i < gridExtent; i++) // compute derivatives via central differences
    {
        const int im1 = (i - 1 + gridExtent) % gridExtent; // periodic index i-1
        const int ip1 = (i + 1) % gridExtent;              // periodic index i+1
        for (int j = 0; j < gridExtent; j++)
        {
            const int jm1 = (j - 1 + gridExtent) % gridExtent;                // periodic index j-1
            const int jp1 = (j + 1) % gridExtent;                             // periodic index j+1
            grid(i, j)[1] = (grid(ip1, j)[0] - grid(im1, j)[0]) / (2.0 * dx); // d/dphi1
            grid(i, j)[2] = (grid(i, jp1)[0] - grid(i, jm1)[0]) / (2.0 * dx); // d/dphi2
            grid(i, j)[3] = (grid(ip1, jp1)[0] - grid(ip1, jm1)[0]            // cross
                             - grid(im1, jp1)[0] + grid(im1, jm1)[0])
                            / (4.0 * dx * dx);
        }
    }
    return grid;
}

//! Build a flat (constant-value) CMAP grid of given extent for testing
static CmapGrid buildFlatGrid(int gridExtent, real value)
{
    CmapGrid grid(gridExtent, gridExtent); // allocate grid
    for (int i = 0; i < gridExtent; i++)
    {
        for (int j = 0; j < gridExtent; j++)
        {
            grid(i, j)[0] = value; // grid value
            grid(i, j)[1] = 0.0;   // first derivative in phi1
            grid(i, j)[2] = 0.0;   // first derivative in phi2
            grid(i, j)[3] = 0.0;   // cross derivative
        }
    }
    return grid;
}

//! Atom coordinates for five atoms forming two consecutive dihedrals
static std::array<rvec, 5> buildAtomPositions()
{
    std::array<rvec, 5> x; // five-atom coordinates for phi/psi pair
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 0.0;
    x[1][0] = 0.1;
    x[1][1] = 0.0;
    x[1][2] = 0.0;
    x[2][0] = 0.2;
    x[2][1] = 0.1;
    x[2][2] = 0.0;
    x[3][0] = 0.3;
    x[3][1] = 0.1;
    x[3][2] = 0.1;
    x[4][0] = 0.4;
    x[4][1] = 0.2;
    x[4][2] = 0.0;
    return x;
}

TEST(CmapFepTest, EnergyInterpolatesLinearly)
{
    // Two flat grids with different constant values; energy must interpolate linearly
    constexpr int  gridExtent = 24;   // grid size matching CHARMM default
    constexpr real valueA     = 10.0; // A-state energy value
    constexpr real valueB     = 20.0; // B-state energy value

    CmapGrids grids;                                    // grid container
    grids.push_back(buildFlatGrid(gridExtent, valueA)); // grid 0: A-state
    grids.push_back(buildFlatGrid(gridExtent, valueB)); // grid 1: B-state

    t_iparams forceparams[1];      // interaction parameters for one CMAP type
    forceparams[0].cmap.cmapA = 0; // A-state grid index
    forceparams[0].cmap.cmapB = 1; // B-state grid index

    std::array<rvec, 5> xArr = buildAtomPositions(); // five-atom coordinate array with ownership
    const rvec*         x    = xArr.data();          // atom positions

    t_pbc pbc; // PBC object (no PBC)
    set_pbc(&pbc, PbcType::No, nullptr);

    // forceatoms: type index followed by five atom indices
    t_iatom   forceatoms[6] = { 0, 0, 1, 2, 3, 4 }; // one CMAP interaction
    const int nbonds        = 6;                    // length of forceatoms array

    for (int iLambda = 0; iLambda <= 4; iLambda++) // test lambda = 0, 0.25, 0.5, 0.75, 1.0
    {
        const real lambda    = iLambda * 0.25; // current lambda value
        real       dvdlambda = 0.0;            // dH/dlambda accumulator

        rvec4 f[5]                      = {}; // forces zeroed each iteration
        rvec  fshift[c_numShiftVectors] = {}; // shift forces zeroed each iteration

        const real energy = cmap_dihs(
                nbonds, forceatoms, forceparams, grids, x, f, fshift, &pbc, lambda, &dvdlambda, {}, nullptr, nullptr, nullptr, nullptr);

        const real expectedEnergy = (1.0 - lambda) * valueA + lambda * valueB; // linear interpolation
        EXPECT_REAL_EQ_TOL(energy, expectedEnergy, defaultRealTolerance())
                << "energy not linearly interpolated at lambda=" << lambda;

        const real expectedDvdl = valueB - valueA; // dH/dlambda = E_B - E_A for flat grids
        EXPECT_REAL_EQ_TOL(dvdlambda, expectedDvdl, defaultRealTolerance())
                << "dvdlambda incorrect at lambda=" << lambda;
    }
}

TEST(CmapFepTest, AtLambdaZeroMatchesAState)
{
    // At lambda=0 result must be identical to pure A-state evaluation
    constexpr int  gridExtent = 24;   // grid size
    constexpr real valueA     = 5.0;  // A-state energy value
    constexpr real valueB     = 15.0; // B-state energy value (must not affect result)

    CmapGrids grids;
    grids.push_back(buildFlatGrid(gridExtent, valueA)); // grid 0: A-state
    grids.push_back(buildFlatGrid(gridExtent, valueB)); // grid 1: B-state

    t_iparams forceparams[1];
    forceparams[0].cmap.cmapA = 0; // A-state grid index
    forceparams[0].cmap.cmapB = 1; // B-state grid index

    std::array<rvec, 5> xArr = buildAtomPositions(); // five-atom coordinate array with ownership
    const rvec*         x    = xArr.data();          // atom positions
    rvec4               f[5] = {};
    rvec                fshift[c_numShiftVectors] = {};
    t_pbc               pbc;
    set_pbc(&pbc, PbcType::No, nullptr);
    t_iatom forceatoms[6] = { 0, 0, 1, 2, 3, 4 };
    real    dvdlambda     = 0.0;

    const real energy = cmap_dihs(
            6, forceatoms, forceparams, grids, x, f, fshift, &pbc, 0.0, &dvdlambda, {}, nullptr, nullptr, nullptr, nullptr);

    // at lambda=0 energy equals A-state value
    EXPECT_REAL_EQ_TOL(energy, valueA, defaultRealTolerance());
}

TEST(CmapFepTest, AtLambdaOneMatchesBState)
{
    // At lambda=1 result must be identical to pure B-state evaluation
    constexpr int  gridExtent = 24;   // grid size
    constexpr real valueA     = 5.0;  // A-state energy value (must not affect result)
    constexpr real valueB     = 15.0; // B-state energy value

    CmapGrids grids;
    grids.push_back(buildFlatGrid(gridExtent, valueA)); // grid 0: A-state
    grids.push_back(buildFlatGrid(gridExtent, valueB)); // grid 1: B-state

    t_iparams forceparams[1];
    forceparams[0].cmap.cmapA = 0; // A-state grid index
    forceparams[0].cmap.cmapB = 1; // B-state grid index

    std::array<rvec, 5> xArr = buildAtomPositions(); // five-atom coordinate array with ownership
    const rvec*         x    = xArr.data();          // atom positions
    rvec4               f[5] = {};
    rvec                fshift[c_numShiftVectors] = {};
    t_pbc               pbc;
    set_pbc(&pbc, PbcType::No, nullptr);
    t_iatom forceatoms[6] = { 0, 0, 1, 2, 3, 4 };
    real    dvdlambda     = 0.0;

    const real energy = cmap_dihs(
            6, forceatoms, forceparams, grids, x, f, fshift, &pbc, 1.0, &dvdlambda, {}, nullptr, nullptr, nullptr, nullptr);

    // at lambda=1 energy equals B-state value
    EXPECT_REAL_EQ_TOL(energy, valueB, defaultRealTolerance());
}

TEST(CmapFepTest, UnperturbedCmapHasZeroDvdl)
{
    // When cmapA == cmapB (unperturbed), dvdlambda must be zero
    constexpr int  gridExtent = 24;  // grid size
    constexpr real value      = 7.0; // grid value (same for A and B)

    CmapGrids grids;
    grids.push_back(buildFlatGrid(gridExtent, value)); // single grid used for both states

    t_iparams forceparams[1];
    forceparams[0].cmap.cmapA = 0; // A-state grid index
    forceparams[0].cmap.cmapB = 0; // B-state same as A: unperturbed

    std::array<rvec, 5> xArr = buildAtomPositions(); // five-atom coordinate array with ownership
    const rvec*         x    = xArr.data();          // atom positions
    rvec4               f[5] = {};
    rvec                fshift[c_numShiftVectors] = {};
    t_pbc               pbc;
    set_pbc(&pbc, PbcType::No, nullptr);
    t_iatom forceatoms[6] = { 0, 0, 1, 2, 3, 4 };
    real    dvdlambda     = 0.0;

    cmap_dihs(6, forceatoms, forceparams, grids, x, f, fshift, &pbc, 0.5, &dvdlambda, {}, nullptr, nullptr, nullptr, nullptr);

    EXPECT_EQ(dvdlambda, 0.0); // unperturbed CMAP must produce zero dH/dlambda
}

TEST(CmapFepTest, EnergyInterpolatesLinearlyNonFlatGrids)
{
    // Non-flat grids with numerically derived derivatives; interpolation must still be linear
    constexpr int gridExtent = 24; // grid size

    CmapGrids grids;
    grids.push_back(buildRandomGrid(gridExtent, 42));  // grid 0: A-state
    grids.push_back(buildRandomGrid(gridExtent, 137)); // grid 1: B-state

    t_iparams forceparams[1];
    forceparams[0].cmap.cmapA = 0; // A-state grid index
    forceparams[0].cmap.cmapB = 1; // B-state grid index

    std::array<rvec, 5> xArr = buildAtomPositions(); // five-atom coordinate array with ownership
    const rvec*         x    = xArr.data();          // atom positions
    t_pbc               pbc;
    set_pbc(&pbc, PbcType::No, nullptr);
    t_iatom forceatoms[6] = { 0, 0, 1, 2, 3, 4 };

    // Obtain eA and eB by evaluating at the endpoints
    rvec4 fA[5]                      = {};
    rvec  fshiftA[c_numShiftVectors] = {};
    real  dvdlA                      = 0.0;
    // energy at lambda=0
    const real eA = cmap_dihs(
            6, forceatoms, forceparams, grids, x, fA, fshiftA, &pbc, 0.0, &dvdlA, {}, nullptr, nullptr, nullptr, nullptr);

    rvec4 fB[5]                      = {};
    rvec  fshiftB[c_numShiftVectors] = {};
    real  dvdlB                      = 0.0;
    // energy at lambda=1
    const real eB = cmap_dihs(
            6, forceatoms, forceparams, grids, x, fB, fshiftB, &pbc, 1.0, &dvdlB, {}, nullptr, nullptr, nullptr, nullptr);

    for (int iLambda = 1; iLambda <= 3; iLambda++) // test lambda = 0.25, 0.5, 0.75
    {
        const real lambda                    = iLambda * 0.25; // current lambda value
        real       dvdl                      = 0.0;
        rvec4      f[5]                      = {};
        rvec       fshift[c_numShiftVectors] = {};

        const real energy = cmap_dihs(
                6, forceatoms, forceparams, grids, x, f, fshift, &pbc, lambda, &dvdl, {}, nullptr, nullptr, nullptr, nullptr);

        const real expectedEnergy = (1.0 - lambda) * eA + lambda * eB; // linear interpolation
        EXPECT_REAL_EQ_TOL(energy, expectedEnergy, defaultRealTolerance())
                << "energy not linearly interpolated at lambda=" << lambda;

        const real expectedDvdl = eB - eA; // dH/dlambda = E_B - E_A
        EXPECT_REAL_EQ_TOL(dvdl, expectedDvdl, defaultRealTolerance())
                << "dvdlambda incorrect at lambda=" << lambda;
    }
}

} // namespace
} // namespace test
} // namespace gmx
