/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018, by the GROMACS development team, led by
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
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "gromacs/applied-forces/densityfitting/densfitforces.h"
#include "gromacs/applied-forces/densityfitting/measures.h"
#include "gromacs/math/griddata/grid.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/operations/densityspreader.h"

#include "gromacs/fileio/griddataio.h"


#include <gtest/gtest.h>

#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

/*!\brief Genereate densities and calculate forces on an atom.
 *
 * Genereate a reference map from an atom. Then calculate the forces on an
 * atom at another position from a simulated map.
 */
class DensfitForcesTest : public ::testing::Test
{
    public:
        DensfitForcesTest()
        {
        }

        std::vector<RVec> calculateForces(
            ArrayRef<const RVec>   referenceAtomCoordinates,
            ArrayRef<const RVec>   fitAtomCoordinates,
            const DensityPotential fitMethod)
        {
            // Setup a grid
            const CanonicalVectorBasis<DIM> cell(cellLength_);
            const Grid<DIM>                 grid(cell, numberGridPoints_);

            // generate a reference map on that grid
            DensitySpreader densitySpreader(grid, nThreads_, nSigmaSpreading_, sigma_, integrateDensityWhenSpreading_);

            densitySpreader.spreadLocalAtoms(referenceAtomCoordinates, weights_);
            const auto referenceMap = densitySpreader.getSpreadGrid();
            densitySpreader.zero();

            densitySpreader.spreadLocalAtoms(fitAtomCoordinates, weights_);
            const auto simulatedMap = densitySpreader.getSpreadGrid();

            const auto measure = createDensityDensityMeasure(fitMethod, referenceMap);
            measure->evaluateDensityDensityDerivative(simulatedMap, forceConstant_);

            const auto        derivativeMap = measure->densityDensityDerivative();

            DensfitForces     densfitForces(grid, sigma_, nSigmaSpreading_);

            std::vector<RVec> forces;
            for (const auto &x : fitAtomCoordinates)
            {
                forces.push_back(densfitForces.force(x, derivativeMap));
            }
            return forces;
        };

        TestReferenceData data_;

    private:

        CanonicalVectorBasis<DIM>::NdVector cellLength_ = {{1, 1, 1}};

        ColumnMajorLattice<DIM>::MultiIndex numberGridPoints_ = {{20, 20, 20}};

        const std::vector<real>             weights_ = {1};

        const int  nThreads_                      = 1;
        const real sigma_                         = 0.2;
        const real nSigmaSpreading_               = 5;
        const bool integrateDensityWhenSpreading_ = true;

        const real forceConstant_ = 1.;
};

TEST_F(DensfitForcesTest, PullAtomToGeneratedDensity)
{
    TestReferenceChecker   checker(data_.rootChecker());
    FloatingPointTolerance tolerance( gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-2));
    checker.setDefaultTolerance(tolerance);

    const std::vector<RVec> referenceAtomCoordinates = {{0.25, 0.25, 0.25}};
    const std::vector<RVec> fitAtomCoordinates       = {{0.75, 0.75, 0.75}};

    auto                    forces = calculateForces(referenceAtomCoordinates, fitAtomCoordinates, DensityPotential::CrossCorrelation);
    checker.checkSequence(forces.begin(), forces.end(), "CrossCorrelationForces");

    forces = calculateForces(referenceAtomCoordinates, fitAtomCoordinates, DensityPotential::CrossEntropy);
    checker.checkSequence(forces.begin(), forces.end(), "CrossEntropyForces");

    forces = calculateForces(referenceAtomCoordinates, fitAtomCoordinates, DensityPotential::MeanSquareDeviation);
    checker.checkSequence(forces.begin(), forces.end(), "MsdForces");

}

} // namespace test

} // namespace gmx
