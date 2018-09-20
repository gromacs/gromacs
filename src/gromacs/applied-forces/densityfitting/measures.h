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
#ifndef GMX_APPLIEDFORCES_DENSITYFITTING_MEASURES_H_
#define GMX_APPLIEDFORCES_DENSITYFITTING_MEASURES_H_

#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/griddata/griddata.h"


namespace gmx
{

//! Defines how simulated and reference density are compared.
enum class DensityPotential
{
    CrossCorrelation,
    CrossEntropy,
    MeanSquareDeviation,
    NumberOfPotentials
};

class IDensityDensityMeasure;

/*!\brief Make a new density-density measure, depending on the type of comparison. */
std::unique_ptr<IDensityDensityMeasure> createDensityDensityMeasure(DensityPotential potential, const GridDataFloat3D &referenceMap);

/**! \brief Abstract base class for measures that compare densities */
class IDensityDensityMeasure
{
    public:
        virtual void
        evaluateDensityDensityDerivative(const GridDataFloat3D &simulatedMap, real forceConstant) = 0;
        virtual const GridDataFloat3D &densityDensityDerivative()       = 0;
        virtual real goodnessOfFit(const GridDataFloat3D &simulatedMap) = 0;
};

class DensityCrossCorrelationMeasure : public IDensityDensityMeasure
{
    public:
        DensityCrossCorrelationMeasure(const GridDataFloat3D &referenceMap);
        void
        evaluateDensityDensityDerivative(const GridDataFloat3D &simulatedMap, real forceConstant) override;
        const GridDataFloat3D &densityDensityDerivative() override;
        real goodnessOfFit(const GridDataFloat3D &simulatedMap) override;
    private:
        bool            hasNormalizedSimulatedMap;

        GridDataFloat3D normalizedSimulatedMap_;
        GridDataFloat3D normalizedReferenceMap_; //< Reference map with zero mean and unit variance for speed-up of cross correlation calculations

        real            normlizationFactorReferenceMap_;
        real            normlizationFactorSimulatedMap_;

        GridDataFloat3D derivativeMap_;
};

class DensityMeanSquareDeviationMeasure : public IDensityDensityMeasure
{
    public:
        DensityMeanSquareDeviationMeasure(const GridDataFloat3D &referenceMap);

        void
        evaluateDensityDensityDerivative(const GridDataFloat3D &simulatedMap, real forceConstant) override;

        const GridDataFloat3D &densityDensityDerivative() override;

        real goodnessOfFit( const GridDataFloat3D &simulatedMap) override;
    private:
        const GridDataFloat3D &referenceMap_;
        GridDataFloat3D        derivativeMap_;
};

class DensityCrossEntropyMeasure : public IDensityDensityMeasure
{
    public:
        DensityCrossEntropyMeasure(const GridDataFloat3D &referenceMap);
        void
        evaluateDensityDensityDerivative(const GridDataFloat3D &simulatedMap, real forceConstant) override;
        const GridDataFloat3D &densityDensityDerivative() override;
        real goodnessOfFit(const GridDataFloat3D &simulatedMap) override;
    private:
        const GridDataFloat3D &referenceMap_;
        GridDataFloat3D        derivativeMap_;
};


}

#endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTING_MEASURES_H_ */
