/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Implements density fitting forces.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/densityfittingforce.h"

#include <array>
#include <memory>

#include "gromacs/math/functions.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/********************************************************************
 * DensityFittingForce::Impl
 */

/*! \internal \brief
 * Private implementation class for DensityFittingForce.
 */
class DensityFittingForce::Impl
{
public:
    /*! \brief Construct densityfitting force implementation from
     * spread of function used to generate the density and spread range.
     * \param[in] kernelShapeParameters determine the shape of the spreading kernel
     */
    explicit Impl(const GaussianSpreadKernelParameters::Shape& kernelShapeParameters);
    //! \copydoc DensityFittingForce::evaluateForce
    RVec evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters,
                       basic_mdspan<const float, dynamicExtents3D> densityDerivative);
    //! The width of the Gaussian in lattice spacing units
    DVec sigma_;
    //! The spread range in lattice points
    IVec latticeSpreadRange_;
    //! The three one-dimensional Gaussians that are used in the force calculation
    std::array<GaussianOn1DLattice, DIM> gauss1d_;
    //! The outer product of a Gaussian along the z and y dimension
    OuterProductEvaluator outerProductZY_;
};

DensityFittingForce::Impl::Impl(const GaussianSpreadKernelParameters::Shape& kernelShapeParameters) :
    sigma_{ kernelShapeParameters.sigma_ },
    latticeSpreadRange_{ kernelShapeParameters.latticeSpreadRange()[XX],
                         kernelShapeParameters.latticeSpreadRange()[YY],
                         kernelShapeParameters.latticeSpreadRange()[ZZ] },
    gauss1d_({ GaussianOn1DLattice(latticeSpreadRange_[XX], sigma_[XX]),
               GaussianOn1DLattice(latticeSpreadRange_[YY], sigma_[YY]),
               GaussianOn1DLattice(latticeSpreadRange_[ZZ], sigma_[ZZ]) })
{
}

RVec DensityFittingForce::Impl::evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters,
                                              basic_mdspan<const float, dynamicExtents3D> densityDerivative)
{
    const IVec closestLatticePoint(roundToInt(localParameters.coordinate_[XX]),
                                   roundToInt(localParameters.coordinate_[YY]),
                                   roundToInt(localParameters.coordinate_[ZZ]));
    const auto spreadRange = spreadRangeWithinLattice(
            closestLatticePoint, densityDerivative.extents(), latticeSpreadRange_);

    // do nothing if the added Gaussian will never reach the lattice
    if (spreadRange.empty())
    {
        return { 0., 0., 0. };
    }

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        // multiply with amplitude so that Gauss3D = (amplitude * Gauss_x) * Gauss_y * Gauss_z
        const float gauss1DAmplitude = dimension > XX ? 1.0 : localParameters.amplitude_;
        gauss1d_[dimension].spread(
                gauss1DAmplitude, localParameters.coordinate_[dimension] - closestLatticePoint[dimension]);
    }

    const auto spreadZY = outerProductZY_(gauss1d_[ZZ].view(), gauss1d_[YY].view());
    const auto spreadX  = gauss1d_[XX].view();
    const IVec spreadGridOffset(latticeSpreadRange_[XX] - closestLatticePoint[XX],
                                latticeSpreadRange_[YY] - closestLatticePoint[YY],
                                latticeSpreadRange_[ZZ] - closestLatticePoint[ZZ]);

    const DVec differenceVectorScale  = { 1. / (square(sigma_[XX])),
                                         1. / (square(sigma_[YY])),
                                         1. / (square(sigma_[ZZ])) };
    const DVec differenceVectorOffset = scaleByVector(
            spreadRange.begin().toDVec() - localParameters.coordinate_.toDVec(), differenceVectorScale);

    DVec differenceVector = differenceVectorOffset;

    DVec force = { 0., 0., 0. };

    for (int zLatticeIndex = spreadRange.begin()[ZZ]; zLatticeIndex < spreadRange.end()[ZZ];
         ++zLatticeIndex, differenceVector[ZZ] += differenceVectorScale[ZZ])
    {
        auto zSliceOfDerivative = densityDerivative[zLatticeIndex];

        differenceVector[YY] = differenceVectorOffset[YY];
        for (int yLatticeIndex = spreadRange.begin()[YY]; yLatticeIndex < spreadRange.end()[YY];
             ++yLatticeIndex, differenceVector[YY] += differenceVectorScale[YY])
        {
            auto       ySliceOfDerivative = zSliceOfDerivative[yLatticeIndex];
            const auto zyPrefactor        = spreadZY(zLatticeIndex + spreadGridOffset[ZZ],
                                              yLatticeIndex + spreadGridOffset[YY]);

            differenceVector[XX] = differenceVectorOffset[XX];
            for (int xLatticeIndex = spreadRange.begin()[XX]; xLatticeIndex < spreadRange.end()[XX];
                 ++xLatticeIndex, differenceVector[XX] += differenceVectorScale[XX])
            {
                const double preFactor = zyPrefactor * spreadX[xLatticeIndex + spreadGridOffset[XX]]
                                         * ySliceOfDerivative[xLatticeIndex];
                force += preFactor * differenceVector;
            }
        }
    }
    return localParameters.amplitude_ * force.toRVec();
}

/********************************************************************
 * DensityFittingForce
 */

DensityFittingForce::DensityFittingForce(const GaussianSpreadKernelParameters::Shape& kernelShapeParameters) :
    impl_(new Impl(kernelShapeParameters))
{
}

RVec DensityFittingForce::evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters,
                                        basic_mdspan<const float, dynamicExtents3D> densityDerivative)
{
    return impl_->evaluateForce(localParameters, densityDerivative);
}

DensityFittingForce::~DensityFittingForce() {}

DensityFittingForce::DensityFittingForce(const DensityFittingForce& other) :
    impl_(new Impl(*other.impl_))
{
}

DensityFittingForce& DensityFittingForce::operator=(const DensityFittingForce& other)
{
    *impl_ = *other.impl_;
    return *this;
}

DensityFittingForce::DensityFittingForce(DensityFittingForce&&) noexcept = default;

DensityFittingForce& DensityFittingForce::operator=(DensityFittingForce&&) noexcept = default;

} // namespace gmx
