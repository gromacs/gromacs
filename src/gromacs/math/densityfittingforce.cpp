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
 * \brief
 * Implements density fitting forces.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "densityfittingforce.h"

#include "gromacs/math/functions.h"

#include "multidimarray.h"

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
        Impl(const GaussianSpreadKernelParameters::Shape &kernelShapeParameters);
        ~Impl() {}
        //! \copydoc DensityFittingForce::evaluateForce(const DensitySpreadKernelParameters::PositionAndAmplitude & localParameters, basic_mdspan<const float, dynamicExtents3D> densityDerivative)
        RVec evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude &localParameters, basic_mdspan<const float, dynamicExtents3D> densityDerivative);
        //! The width of the Gaussian in lattice spacing units
        double sigma_;
        //! The spread range in lattice points
        int    latticeSpreadRange_;
        //! The three one-dimensional Gaussians that are used in the force calculation
        std::array<GaussianOn1DLattice, DIM> gauss1d_;
        //! The outer product of a Gaussian along the z and y dimension
        OuterProductEvaluator                outerProductZY_;
};

DensityFittingForce::Impl::Impl(const GaussianSpreadKernelParameters::Shape &kernelShapeParameters)
    : sigma_ {kernelShapeParameters.sigma_},
latticeSpreadRange_ {
    kernelShapeParameters.latticeSpreadRange()
},
gauss1d_({GaussianOn1DLattice(latticeSpreadRange_, sigma_),
          GaussianOn1DLattice(latticeSpreadRange_, sigma_),
          GaussianOn1DLattice(latticeSpreadRange_, sigma_)})
{
}

RVec DensityFittingForce::Impl::evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude &localParameters,
                                              basic_mdspan<const float, dynamicExtents3D> densityDerivative)
{
    const IVec closestLatticePoint(roundToInt(localParameters.coordinate_[XX]), roundToInt(localParameters.coordinate_[YY]), roundToInt(localParameters.coordinate_[ZZ]));
    const auto spreadRange = spreadRangeWithinLattice(closestLatticePoint, densityDerivative.extents(), latticeSpreadRange_);

    // do nothing if the added Gaussian will never reach the lattice
    if (spreadRange.empty())
    {
        return {};
    }

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        // multiply with amplitude so that Gauss3D = (amplitude * Gauss_x) * Gauss_y * Gauss_z
        const float gauss1DAmplitude = dimension > XX ? 1.0 : localParameters.amplitude_;
        gauss1d_[dimension].spread(gauss1DAmplitude, localParameters.coordinate_[dimension] - closestLatticePoint[dimension]);
    }

    const auto  spreadZY         = outerProductZY_(gauss1d_[ZZ].view(), gauss1d_[YY].view());
    const auto  spreadX          = gauss1d_[XX].view();
    const IVec  spreadGridOffset = {
        latticeSpreadRange_ - closestLatticePoint[XX],
        latticeSpreadRange_ - closestLatticePoint[YY],
        latticeSpreadRange_ - closestLatticePoint[ZZ]
    };

    const RVec  differenceVectorOffset =
        (RVec(spreadRange.begin()[XX], spreadRange.begin()[YY], spreadRange.begin()[ZZ]) - localParameters.coordinate_) / square(sigma_);
    RVec        differenceVector       = differenceVectorOffset;
    const auto  differenceVectorScale  = 1./(square(sigma_));

    RVec        force = {0., 0., 0.};

    for (int zLatticeIndex = spreadRange.begin()[ZZ]; zLatticeIndex < spreadRange.end()[ZZ]; ++zLatticeIndex)
    {
        auto zSliceOfDerivative = densityDerivative[zLatticeIndex];
        for (int yLatticeIndex = spreadRange.begin()[YY]; yLatticeIndex < spreadRange.end()[YY]; ++yLatticeIndex)
        {
            auto       ySliceOfDerivative = zSliceOfDerivative[yLatticeIndex];
            const auto zyPrefactor        = spreadZY(zLatticeIndex + spreadGridOffset[ZZ],
                                                     yLatticeIndex + spreadGridOffset[YY]);

            for (int xLatticeIndex = spreadRange.begin()[XX]; xLatticeIndex < spreadRange.end()[XX]; ++xLatticeIndex)
            {
                const real preFactor = zyPrefactor *
                    spreadX[xLatticeIndex + spreadGridOffset[XX]] * ySliceOfDerivative[xLatticeIndex];
                force += preFactor * differenceVector;

                differenceVector[XX] += differenceVectorScale;
            }
            differenceVector[XX]  = differenceVectorOffset[XX];
            differenceVector[YY] += differenceVectorScale;
        }
        differenceVector[YY]  = differenceVectorOffset[YY];
        differenceVector[ZZ] += differenceVectorScale;
    }
    return localParameters.amplitude_*force;
}

/********************************************************************
 * DensityFittingForce
 */

DensityFittingForce::DensityFittingForce(const GaussianSpreadKernelParameters::Shape &kernelShapeParameters) :
    impl_(new Impl(kernelShapeParameters))
{
}

RVec DensityFittingForce::evaluateForce(const GaussianSpreadKernelParameters::PositionAndAmplitude &localParameters,
                                        basic_mdspan<const float, dynamicExtents3D> densityDerivative)
{
    return impl_->evaluateForce(localParameters, densityDerivative);
}

DensityFittingForce::~DensityFittingForce()
{
}

} // namespace gmx
