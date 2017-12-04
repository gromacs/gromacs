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

#include "gmxpre.h"

#include <algorithm>
#include <numeric>
#include "measures.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/container/containeroperation.h"


namespace gmx
{

namespace
{

#if GMX_SIMD_HAVE_REAL

SimdReal times(SimdReal a, SimdReal b) { return a * b; }
typedef SimdReal (*binary_real_op)(SimdReal, SimdReal);
template <binary_real_op F>
real innerProduct(PaddedArrayRef<const real> reference,
                  PaddedArrayRef<const real> comparand)
{
    SimdReal   simdSum = 0.;

    const auto simdAlignedSize =
        (reference.size() / GMX_SIMD_REAL_WIDTH) * GMX_SIMD_REAL_WIDTH;
    for (long i = 0; i < simdAlignedSize; i += GMX_SIMD_REAL_WIDTH)
    {
        const SimdReal sx = load<SimdReal>(reference.data() + i);
        const SimdReal sy = load<SimdReal>(comparand.data() + i);
        simdSum = simdSum + F(sx, sy);
    }
    std::vector < real, AlignedAllocator < real>> summationArray(GMX_SIMD_REAL_WIDTH);
    store(summationArray.data(), simdSum);
    real result =
        std::accumulate(std::begin(summationArray), std::end(summationArray), 0.);

    std::vector < real, AlignedAllocator < real>> referenceExtra(GMX_SIMD_REAL_WIDTH,
                                                                 0.);
    std::vector < real, AlignedAllocator < real>> comparandExtra(GMX_SIMD_REAL_WIDTH,
                                                                 0.);
    for (long i = simdAlignedSize; i < reference.size(); ++i)
    {
        referenceExtra[i - simdAlignedSize] = reference[i];
        comparandExtra[i - simdAlignedSize] = comparand[i];
    }
    std::vector < real, AlignedAllocator < real>> lastvalues(GMX_SIMD_REAL_WIDTH, 0.);
    store(lastvalues.data(), F(load<SimdReal>(referenceExtra.data()),
                               load<SimdReal>(comparandExtra.data())));
    for (long i = simdAlignedSize; i < reference.size(); ++i)
    {
        result += lastvalues[i - simdAlignedSize];
    }
    return result;
}

struct gaussianLogLikelihoodFunctor {
    static SimdReal value(SimdReal reference, SimdReal comparand);
    static real     derivative(real reference, real comparand);
};


#else   // GMX_SIMD_HAVE_REAL

real times(real a, real b) { return a * b; }
typedef real (*binary_real_op)(real, real);
template <binary_real_op F>
real innerProduct(PaddedArrayRef<const real> reference,
                  PaddedArrayRef<const real> comparand)
{
    return std::inner_product(std::begin(reference), std::end(reference),
                              std::begin(comparand), 0., std::plus<real>(), F);
}
struct gaussianLogLikelihoodFunctor {
    static real value(real reference, real comparand);
    static real derivative(real reference, real comparand);
};

#endif

struct poissonLogLikelihoodFunctor {
#if GMX_SIMD_HAVE_REAL
    static SimdReal value(SimdReal reference, SimdReal comparand);
#else
    static real value(real reference, real comparand);
#endif
    static real derivative(real reference, real comparand);
};

#if GMX_SIMD_HAVE_REAL
SimdReal gaussianLogLikelihoodFunctor::value(SimdReal reference,
                                             SimdReal comparand)
#else
real gaussianLogLikelihoodFunctor::value(real reference, real comparand)
#endif  // GMX_SIMD_HAVE_REAL
{
    return 0.5 * (reference - comparand) * (reference - comparand);
}

real gaussianLogLikelihoodFunctor::derivative(real reference, real comparand)
{
    return comparand - reference;
}

#if GMX_SIMD_HAVE_REAL
SimdReal poissonLogLikelihoodFunctor::value(SimdReal reference,
                                            SimdReal comparand)
#else
real poissonLogLikelihoodFunctor::value(real reference, real comparand)
#endif  // GMX_SIMD_HAVE_REAL
{
    return -reference * log(comparand + 1e-6);
}

real poissonLogLikelihoodFunctor::derivative(real reference, real comparand)
{
    return -reference / (comparand + 1e-6);
}

real normalizeMap(GridDataReal3D *map)
{
    const auto mean = containermeasure::mean(*map);
    containeroperation::addScalar(map, -mean);
    const auto squaredNorm = containermeasure::normSquared(*map);
    containeroperation::divide(map, sqrt(squaredNorm));
    return squaredNorm;
}

}   // namespace

/*************************DensityCrossCorrelationMeasure***********************/

DensityCrossCorrelationMeasure::DensityCrossCorrelationMeasure(const GridDataReal3D &referenceMap)
    :  hasNormalizedSimulatedMap(false),
      normalizedReferenceMap_(referenceMap),
      normlizationFactorReferenceMap_(normalizeMap(&normalizedReferenceMap_)),
      derivativeMap_(referenceMap.getGrid().duplicate())
{

}

const GridDataReal3D &DensityCrossCorrelationMeasure::densityDensityDerivative()
{
    return derivativeMap_;
}

void DensityCrossCorrelationMeasure::evaluateDensityDensityDerivative( const GridDataReal3D &simulatedMap, real forceConstant)
{
    if (!hasNormalizedSimulatedMap)
    {
        normalizedSimulatedMap_         = simulatedMap;
        normlizationFactorSimulatedMap_ =  normalizeMap(&normalizedSimulatedMap_);
    }
    hasNormalizedSimulatedMap = false;
    auto scale = -1.f/normlizationFactorSimulatedMap_;
    std::transform(
            std::begin(normalizedReferenceMap_), std::end(normalizedReferenceMap_),
            std::begin(simulatedMap),
            std::begin(derivativeMap_),
            [scale](real simulated, real reference){return scale * simulated * reference; }
            );
    containeroperation::multiply(&derivativeMap_, forceConstant);

}

real DensityCrossCorrelationMeasure::goodnessOfFit(
        const GridDataReal3D &simulatedMap)
{
    normalizedSimulatedMap_         = simulatedMap;
    normlizationFactorSimulatedMap_ = normalizeMap(&normalizedSimulatedMap_);
    hasNormalizedSimulatedMap       = true;
    return innerProduct<times>(normalizedReferenceMap_, normalizedSimulatedMap_);
}

/*************************DensityCrossEntropyMeasure***************************/

DensityCrossEntropyMeasure::DensityCrossEntropyMeasure(const GridDataReal3D
                                                       &referenceMap) : referenceMap_(referenceMap),
                                                                        derivativeMap_(referenceMap.getGrid().duplicate())
{

}

const GridDataReal3D &DensityCrossEntropyMeasure::densityDensityDerivative()
{
    return derivativeMap_;
}

void DensityCrossEntropyMeasure::evaluateDensityDensityDerivative(const GridDataReal3D &simulatedMap, real forceConstant)
{
    std::transform(
            std::begin(referenceMap_), std::end(referenceMap_),
            std::begin(simulatedMap),
            std::begin(derivativeMap_),
            poissonLogLikelihoodFunctor::derivative
            );
    containeroperation::multiply(&derivativeMap_, forceConstant);
}

real DensityCrossEntropyMeasure::goodnessOfFit(const GridDataReal3D &simulatedMap)
{
    return innerProduct<poissonLogLikelihoodFunctor::value>(referenceMap_, simulatedMap);
}

/*************************DensityMeanSquareDeviationMeasure********************/

DensityMeanSquareDeviationMeasure::DensityMeanSquareDeviationMeasure(const
                                                                     GridDataReal3D &referenceMap) : referenceMap_(referenceMap),
                                                                                                     derivativeMap_(referenceMap.getGrid().duplicate()) {}

const GridDataReal3D &DensityMeanSquareDeviationMeasure::densityDensityDerivative()
{
    return derivativeMap_;
}

void DensityMeanSquareDeviationMeasure::evaluateDensityDensityDerivative(
        const GridDataReal3D &simulatedMap, real forceConstant)
{
    std::transform(
            std::begin(referenceMap_), std::end(referenceMap_),
            std::begin(simulatedMap),
            std::begin(derivativeMap_),
            gaussianLogLikelihoodFunctor::derivative
            );
    containeroperation::multiply(&derivativeMap_, forceConstant);

}

real DensityMeanSquareDeviationMeasure::goodnessOfFit(
        const GridDataReal3D &simulatedMap)
{
    return innerProduct<gaussianLogLikelihoodFunctor::value>(referenceMap_, simulatedMap);
}

std::unique_ptr<IDensityDensityMeasure> createDensityDensityMeasure(DensityPotential potential, const GridDataReal3D &referenceMap)
{
    switch (potential)
    {
        case DensityPotential::CrossCorrelation:
            return std::unique_ptr<IDensityDensityMeasure>(new DensityCrossCorrelationMeasure(referenceMap));
        case DensityPotential::CrossEntropy:
            return std::unique_ptr<IDensityDensityMeasure>(new DensityCrossEntropyMeasure(referenceMap));
        case DensityPotential::MeanSquareDeviation:
            return std::unique_ptr<IDensityDensityMeasure>(new DensityMeanSquareDeviationMeasure(referenceMap));
        default:
            return nullptr;
    }
}


}  // namespace gmx
