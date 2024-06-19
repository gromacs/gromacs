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
 * Implements density similarity measures and their derivatives.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/densityfit.h"

#include <cmath>

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class DensitySimilarityMeasureImpl
{
public:
    virtual ~DensitySimilarityMeasureImpl();
    //! convenience typedef
    using density = DensitySimilarityMeasure::density;
    //! \copydoc DensitySimilarityMeasure::gradient(DensitySimilarityMeasure::density comparedDensity)
    virtual density gradient(density comparedDensity) = 0;
    //! \copydoc DensitySimilarityMeasure::similarity(density comparedDensity)
    virtual real similarity(density comparedDensity) = 0;
    //! clone to allow copy operations
    virtual std::unique_ptr<DensitySimilarityMeasureImpl> clone() = 0;
};
DensitySimilarityMeasureImpl::~DensitySimilarityMeasureImpl() = default;

namespace
{

/****************** Inner Product *********************************************/

/*! \internal
 * \brief Implementation for DensitySimilarityInnerProduct.
 *
 * The similarity measure itself is documented in DensitySimilarityMeasureMethod::innerProduct.
 */
class DensitySimilarityInnerProduct final : public DensitySimilarityMeasureImpl
{
public:
    //! Construct similarity measure by setting the reference density
    DensitySimilarityInnerProduct(density referenceDensity);
    //! The gradient for the inner product similarity measure is the reference density divided by the number of voxels
    density gradient(density comparedDensity) override;
    //! Clone this
    std::unique_ptr<DensitySimilarityMeasureImpl> clone() override;
    //! The similarity between reference density and compared density
    real similarity(density comparedDensity) override;

private:
    //! A view on the reference density
    const density referenceDensity_;
    //! Stores the gradient of the similarity measure in memory
    MultiDimArray<std::vector<float>, dynamicExtents3D> gradient_;
};

DensitySimilarityInnerProduct::DensitySimilarityInnerProduct(density referenceDensity) :
    referenceDensity_{ referenceDensity }, gradient_{ referenceDensity.extents() }
{
    const auto numVoxels = gradient_.asConstView().mapping().required_span_size();
    /* the gradient for the inner product measure of fit is constant and does not
     * depend on the compared density, so it is pre-computed here */
    std::transform(begin(referenceDensity_), end(referenceDensity_), begin(gradient_), [numVoxels](float x) {
        return x / numVoxels;
    });
}

real DensitySimilarityInnerProduct::similarity(density comparedDensity)
{
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }
    /* the similarity measure uses the gradient instead of the reference,
     * here, because it is the reference density divided by the number of voxels */
    return std::inner_product(begin(gradient_), end(gradient_), begin(comparedDensity), 0.);
}

DensitySimilarityMeasure::density DensitySimilarityInnerProduct::gradient(density comparedDensity)
{
    /* even though the gradient density does not depend on the compad density,
     * still checking the extents to make sure we're consistent */
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }

    return gradient_.asConstView();
}

std::unique_ptr<DensitySimilarityMeasureImpl> DensitySimilarityInnerProduct::clone()
{
    return std::make_unique<DensitySimilarityInnerProduct>(referenceDensity_);
}

/****************** Relative Entropy *****************************************/

//! Calculate a single summand in the relative entropy sum.
real relativeEntropyAtVoxel(real reference, real comparison)
{
    if ((reference > 0) && (comparison > 0))
    {
        return reference * (std::log(comparison / reference));
    }
    return 0.;
}

//! Calculate a single relative entropy gradient entry at a voxel.
real relativeEntropyGradientAtVoxel(real reference, real comparison)
{
    if ((reference > 0) && (comparison > 0))
    {
        return reference / comparison;
    }
    return 0.;
}

/*! \internal
 * \brief Implementation for DensitySimilarityRelativeEntropy.
 *
 * The similarity measure itself is documented in DensitySimilarityMeasureMethod::RelativeEntropy.
 */
class DensitySimilarityRelativeEntropy final : public DensitySimilarityMeasureImpl
{
public:
    //! Construct similarity measure by setting the reference density
    DensitySimilarityRelativeEntropy(density referenceDensity);
    //! The gradient for the relative entropy similarity measure
    density gradient(density comparedDensity) override;
    //! Clone this
    std::unique_ptr<DensitySimilarityMeasureImpl> clone() override;
    //! The similarity between reference density and compared density
    real similarity(density comparedDensity) override;

private:
    //! A view on the reference density
    const density referenceDensity_;
    //! Stores the gradient of the similarity measure in memory
    MultiDimArray<std::vector<float>, dynamicExtents3D> gradient_;
};

DensitySimilarityRelativeEntropy::DensitySimilarityRelativeEntropy(density referenceDensity) :
    referenceDensity_{ referenceDensity }, gradient_(referenceDensity.extents())
{
}

real DensitySimilarityRelativeEntropy::similarity(density comparedDensity)
{
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }
    return std::inner_product(begin(referenceDensity_),
                              end(referenceDensity_),
                              begin(comparedDensity),
                              0.,
                              std::plus<>(),
                              relativeEntropyAtVoxel);
}

DensitySimilarityMeasure::density DensitySimilarityRelativeEntropy::gradient(density comparedDensity)
{
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }
    std::transform(begin(referenceDensity_),
                   end(referenceDensity_),
                   begin(comparedDensity),
                   begin(gradient_),
                   relativeEntropyGradientAtVoxel);
    return gradient_.asConstView();
}

std::unique_ptr<DensitySimilarityMeasureImpl> DensitySimilarityRelativeEntropy::clone()
{
    return std::make_unique<DensitySimilarityRelativeEntropy>(referenceDensity_);
}

/****************** Cross Correlation *****************************************/

//! Helper values for evaluating the cross correlation
struct CrossCorrelationEvaluationHelperValues
{
    //! The mean of the reference density
    real meanReference = 0;
    //! The mean of the compared density
    real meanComparison = 0;
    //! The sum of the squared reference density voxel values
    real referenceSquaredSum = 0;
    //! The sum of the squared compared density voxel values
    real comparisonSquaredSum = 0;
    //! The covariance of the reference and the compared density
    real covariance = 0;
};

/*! \brief Calculate helper values for the cross-correlation.

 * Enables numerically stable single-pass cross-correlation evaluation algorithm
 * as described in Bennett, J., Grout, R. , Pebay, P., Roe D., Thompson D.
 * "Numerically Stable, Single-Pass, Parallel Statistics Algorithms"
 * and implemented in boost's correlation coefficient
 */
CrossCorrelationEvaluationHelperValues evaluateHelperValues(DensitySimilarityMeasure::density reference,
                                                            DensitySimilarityMeasure::density compared)
{
    CrossCorrelationEvaluationHelperValues helperValues;

    Index i = 0;

    const auto* referenceIterator = begin(reference);
    for (const real comp : compared)
    {
        const real refHelper        = *referenceIterator - helperValues.meanReference;
        const real comparisonHelper = comp - helperValues.meanComparison;
        helperValues.referenceSquaredSum += (i * square(refHelper)) / (i + 1);
        helperValues.comparisonSquaredSum += (i * square(comparisonHelper)) / (i + 1);
        helperValues.covariance += i * refHelper * comparisonHelper / (i + 1);
        helperValues.meanReference += refHelper / (i + 1);
        helperValues.meanComparison += comparisonHelper / (i + 1);

        ++referenceIterator;
        ++i;
    }

    return helperValues;
}

//! Calculate a single cross correlation gradient entry at a voxel.
class CrossCorrelationGradientAtVoxel
{
public:
    //! Set up the gradient calculation with pre-computed values
    CrossCorrelationGradientAtVoxel(const CrossCorrelationEvaluationHelperValues& preComputed) :
        prefactor_(evaluatePrefactor(preComputed.comparisonSquaredSum, preComputed.referenceSquaredSum)),
        comparisonPrefactor_(preComputed.covariance / preComputed.comparisonSquaredSum),
        meanReference_(preComputed.meanReference),
        meanComparison_(preComputed.meanComparison)
    {
    }
    //! Evaluate the cross correlation gradient at a voxel
    real operator()(real reference, real comparison) const
    {
        return prefactor_
               * (reference - meanReference_ - comparisonPrefactor_ * (comparison - meanComparison_));
    }

private:
    static real evaluatePrefactor(real comparisonSquaredSum, real referenceSquaredSum)
    {
        GMX_ASSERT(comparisonSquaredSum > 0,
                   "Squared sum of comparison values needs to be larger than zero.");
        GMX_ASSERT(referenceSquaredSum > 0,
                   "Squared sum of reference values needs to be larger than zero.");
        return 1.0 / (std::sqrt(comparisonSquaredSum) * std::sqrt(referenceSquaredSum));
    }
    const real prefactor_;
    const real comparisonPrefactor_;
    const real meanReference_;
    const real meanComparison_;
};

/*! \internal
 * \brief Implementation for DensitySimilarityCrossCorrelation.
 *
 * The similarity measure itself is documented in DensitySimilarityMeasureMethod::crossCorrelation.
 */
class DensitySimilarityCrossCorrelation final : public DensitySimilarityMeasureImpl
{
public:
    //! Construct similarity measure by setting the reference density
    DensitySimilarityCrossCorrelation(density referenceDensity);
    //! The gradient for the cross correlation similarity measure
    density gradient(density comparedDensity) override;
    //! Clone this
    std::unique_ptr<DensitySimilarityMeasureImpl> clone() override;
    //! The similarity between reference density and compared density
    real similarity(density comparedDensity) override;

private:
    //! A view on the reference density
    const density referenceDensity_;
    //! Stores the gradient of the similarity measure in memory
    MultiDimArray<std::vector<float>, dynamicExtents3D> gradient_;
};

DensitySimilarityCrossCorrelation::DensitySimilarityCrossCorrelation(density referenceDensity) :
    referenceDensity_{ referenceDensity }, gradient_(referenceDensity.extents())
{
}

real DensitySimilarityCrossCorrelation::similarity(density comparedDensity)
{
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }

    CrossCorrelationEvaluationHelperValues helperValues =
            evaluateHelperValues(referenceDensity_, comparedDensity);

    if ((helperValues.referenceSquaredSum == 0) || (helperValues.comparisonSquaredSum == 0))
    {
        return 0;
    }

    // To avoid numerical instability due to large squared density value sums
    // division is re-written to avoid multiplying two large numbers
    // as product of two separate divisions of smaller numbers
    const real covarianceSqrt = std::sqrt(std::fabs(helperValues.covariance));
    const int  sign           = helperValues.covariance > 0 ? 1 : -1;
    return sign * (covarianceSqrt / std::sqrt(helperValues.referenceSquaredSum))
           * (covarianceSqrt / std::sqrt(helperValues.comparisonSquaredSum));
}

DensitySimilarityMeasure::density DensitySimilarityCrossCorrelation::gradient(density comparedDensity)
{
    if (comparedDensity.extents() != referenceDensity_.extents())
    {
        GMX_THROW(RangeError("Reference density and compared density need to have same extents."));
    }

    CrossCorrelationEvaluationHelperValues helperValues =
            evaluateHelperValues(referenceDensity_, comparedDensity);

    std::transform(begin(referenceDensity_),
                   end(referenceDensity_),
                   begin(comparedDensity),
                   begin(gradient_),
                   CrossCorrelationGradientAtVoxel(helperValues));

    return gradient_.asConstView();
}

std::unique_ptr<DensitySimilarityMeasureImpl> DensitySimilarityCrossCorrelation::clone()
{
    return std::make_unique<DensitySimilarityCrossCorrelation>(referenceDensity_);
}


} // namespace


DensitySimilarityMeasure::DensitySimilarityMeasure(DensitySimilarityMeasureMethod method,
                                                   density                        referenceDensity)
{
    // chose the implementation depending on the method of density comparison
    // throw an error if the method is not known
    switch (method)
    {
        case DensitySimilarityMeasureMethod::innerProduct:
            impl_ = std::make_unique<DensitySimilarityInnerProduct>(referenceDensity);
            break;
        case DensitySimilarityMeasureMethod::relativeEntropy:
            impl_ = std::make_unique<DensitySimilarityRelativeEntropy>(referenceDensity);
            break;
        case DensitySimilarityMeasureMethod::crossCorrelation:
            impl_ = std::make_unique<DensitySimilarityCrossCorrelation>(referenceDensity);
            break;
        default: GMX_THROW(NotImplementedError("Similarity measure not implemented."));
    }
}

DensitySimilarityMeasure::density DensitySimilarityMeasure::gradient(density comparedDensity)
{
    return impl_->gradient(comparedDensity);
}

real DensitySimilarityMeasure::similarity(density comparedDensity)
{
    return impl_->similarity(comparedDensity);
}

DensitySimilarityMeasure::~DensitySimilarityMeasure() = default;

DensitySimilarityMeasure::DensitySimilarityMeasure(const DensitySimilarityMeasure& other) :
    impl_(other.impl_->clone())
{
}

DensitySimilarityMeasure& DensitySimilarityMeasure::operator=(const DensitySimilarityMeasure& other)
{
    impl_ = other.impl_->clone();
    return *this;
}

DensitySimilarityMeasure::DensitySimilarityMeasure(DensitySimilarityMeasure&&) noexcept = default;

DensitySimilarityMeasure& DensitySimilarityMeasure::operator=(DensitySimilarityMeasure&&) noexcept = default;

void normalizeSumPositiveValuesToUnity(ArrayRef<float> data)
{
    const double sumDataLargerZero =
            std::accumulate(std::begin(data), std::end(data), 0., [](double sum, float value) {
                return value > 0 ? sum + value : sum;
            });

    // leave the data untouched if there are no values larger than zero
    if (sumDataLargerZero == 0.)
    {
        return;
    }

    std::transform(std::begin(data), std::end(data), std::begin(data), [sumDataLargerZero](float& datum) {
        return datum / sumDataLargerZero;
    });
}

} // namespace gmx
