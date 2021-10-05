/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_AWH_TEST_SETUP_H
#define GMX_AWH_TEST_SETUP_H

#include "gmxpre.h"

#include <memory>
#include <vector>

#include "gromacs/applied_forces/awh/bias.h"
#include "gromacs/mdtypes/awh_params.h"

namespace gmx
{

template<typename>
class ArrayRef;

namespace test
{

/*! \internal \brief
 * Prepare a memory buffer with serialized AwhDimParams.
 */
std::vector<char> awhDimParamSerialized(
        AwhCoordinateProviderType inputCoordinateProvider = AwhCoordinateProviderType::Pull,
        int                       inputCoordIndex         = 0,
        double                    inputOrigin             = 0.5,
        double                    inputEnd                = 1.5,
        double                    inputPeriod             = 0,
        // Correction for removal of GaussianGeometryFactor/2 in histogram size
        double inputDiffusion = 0.1 / (0.144129616073222 * 2));

/*! \internal \brief
 * Struct that gathers all input for setting up and using a Bias
 */
struct AwhTestParameters
{
    explicit AwhTestParameters(ISerializer* serializer);
    //! Move constructor
    AwhTestParameters(AwhTestParameters&& o) noexcept :
        beta(o.beta), awhParams(std::move(o.awhParams)), dimParams(std::move(o.dimParams))
    {
    }
    //! 1/(kB*T).
    double beta;

    //! AWH parameters, this is the struct to actually use.
    AwhParams awhParams;
    //! Dimension parameters for setting up Bias.
    std::vector<DimParams> dimParams;
};

/*! \brief
 * Helper function to set up the C-style AWH parameters for the test.
 *
 * Builds the test input data from serialized data.
 */
AwhTestParameters getAwhTestParameters(AwhHistogramGrowthType            eawhgrowth,
                                       AwhPotentialType                  eawhpotential,
                                       ArrayRef<const std::vector<char>> dimensionParameterBuffers,
                                       bool                              inputUserData,
                                       double                            beta,
                                       bool                              useAwhFep,
                                       double                            inputErrorScaling,
                                       int                               numFepLambdaStates,
                                       int                               biasShareGroup = 0);

} // namespace test
} // namespace gmx

#endif
