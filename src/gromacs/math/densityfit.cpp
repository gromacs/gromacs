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
 * Implements density similarity measures and their derivatives.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "densityfit.h"

#include <algorithm>
#include <numeric>

#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"

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
        //! \copydoc DensitySimilarityMeasure::gradient(density comparedDensity)
        virtual float similarity(density comparedDensity) = 0;
        //! clone to allow copy operations
        virtual std::unique_ptr<DensitySimilarityMeasureImpl> clone() = 0;
};
DensitySimilarityMeasureImpl::~DensitySimilarityMeasureImpl() = default;

namespace
{

/*! \internal \brief
 * Implementation for DensitySimilarityInnerProduct.
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
        float similarity(density comparedDensity) override;
    private:
        //! A view on the reference density
        const density referenceDensity_;
        //! Stores the gradient of the similarity measure in memory
        MultiDimArray<std::vector<float>, dynamicExtents3D> gradient_;
};

DensitySimilarityInnerProduct::DensitySimilarityInnerProduct(density referenceDensity) :
    referenceDensity_ {referenceDensity },
gradient_ {
    referenceDensity.extents()
}
{
    const auto numVoxels = gradient_.asConstView().mapping().required_span_size();
    /* the gradient for the inner product measure of fit is constant and does not
     * depend on the compared density, so it is pre-computed here */
    std::transform(begin(referenceDensity_), end(referenceDensity), begin(gradient_),
                   [numVoxels](float x){return x/numVoxels; });
}

float DensitySimilarityInnerProduct::similarity(density comparedDensity)
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

}   // namespace


DensitySimilarityMeasure::DensitySimilarityMeasure(DensitySimilarityMeasureMethod method, density referenceDensity)
{
    // chose the implementation depending on the method of density comparison
    switch (method)
    {
        case DensitySimilarityMeasureMethod::innerProduct:
            impl_ = std::make_unique<DensitySimilarityInnerProduct>(referenceDensity);
            break;

        default:
            break;
    }
}

DensitySimilarityMeasure::density DensitySimilarityMeasure::gradient(density comparedDensity)
{
    return impl_->gradient(comparedDensity);
}

float DensitySimilarityMeasure::similarity(density comparedDensity)
{
    return impl_->similarity(comparedDensity);
}

DensitySimilarityMeasure::~DensitySimilarityMeasure()  = default;

DensitySimilarityMeasure::DensitySimilarityMeasure(const DensitySimilarityMeasure &other)
    : impl_(other.impl_->clone())
{
}

DensitySimilarityMeasure &DensitySimilarityMeasure::operator=(const DensitySimilarityMeasure &other)
{
    impl_ = other.impl_->clone();
    return *this;
}

DensitySimilarityMeasure::DensitySimilarityMeasure(DensitySimilarityMeasure &&) noexcept = default;

DensitySimilarityMeasure &DensitySimilarityMeasure::operator=(DensitySimilarityMeasure &&) noexcept = default;

} // namespace gmx
