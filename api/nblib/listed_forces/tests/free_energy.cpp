/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief
 * This generates free energy
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include <vector>

#include "testutils/testasserts.h"

#include "nblib/listed_forces/dataflow.hpp"

namespace nblib
{

TEST(FEP, minimumViableProduct)
{
    real kA = 2.1;
    real kB = 2.2;

    real xA = 1.0;
    real xB = 1.1;

    std::vector<HarmonicBondType>                   parametersA{ { kA, xA } };
    std::vector<HarmonicBondType>                   parametersB{ { kB, xB } };
    std::vector<InteractionIndex<HarmonicBondType>> indices{ { 0, 1, 0 } };

    std::vector<gmx::RVec> x{ { 0.0, 0.0, 0.0 }, { 0.1, 0.1, 0.1 } };
    std::vector<gmx::RVec> f(2, gmx::RVec{ 0, 0, 0 });

    real lambda = 0.1;

    auto energies = computeForces(gmx::ArrayRef<const InteractionIndex<HarmonicBondType>>(indices),
                                  gmx::ArrayRef<const HarmonicBondType>(parametersA),
                                  gmx::ArrayRef<const HarmonicBondType>(parametersB),
                                  gmx::ArrayRef<const gmx::RVec>(x),
                                  lambda,
                                  &f,
                                  gmx::ArrayRef<std::nullptr_t>{},
                                  NoPbc{});

    real kk = (real(1.0) - lambda) * kA + lambda * kB;
    real x0 = (real(1.0) - lambda) * xA + lambda * xB;

    real dx = std::sqrt(real(0.03)) - x0;

    real dvdlRef = real(0.5) * (kB - kA) * dx * dx + (xA - xB) * kk * dx;

    EXPECT_REAL_EQ_TOL(dvdlRef, energies.freeEnergyDerivative(), gmx::test::absoluteTolerance(1e-6));
}

} // namespace nblib
