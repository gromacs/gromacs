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
 * \brief Implements coordinate transformation routines.
 *
 * \author Christian Blau <blau@kth.se>
 */
#include "gmxpre.h"

#include "gromacs/math/coordinatetransformation.h"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/matrix.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/********************************************************************
 * ScaleCoordinates::Impl
 */

class ScaleCoordinates::Impl
{
public:
    Impl(const RVec& scale);
    void inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const;
    void scale(ArrayRef<RVec> coordinates) const;

private:
    RVec scale_;
};
ScaleCoordinates::Impl::Impl(const RVec& scale) : scale_{ scale } {}
void ScaleCoordinates::Impl::scale(ArrayRef<RVec> coordinates) const
{
    for (auto& coordinate : coordinates)
    {
        coordinate[XX] *= scale_[XX];
        coordinate[YY] *= scale_[YY];
        coordinate[ZZ] *= scale_[ZZ];
    }
}

void ScaleCoordinates::Impl::inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const
{
    RVec inverseScale;
    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        inverseScale[dimension] = scale_[dimension] != 0 ? 1. / scale_[dimension] : 1.;
    }

    for (auto& coordinate : coordinates)
    {
        coordinate[XX] *= inverseScale[XX];
        coordinate[YY] *= inverseScale[YY];
        coordinate[ZZ] *= inverseScale[ZZ];
    }
}

/********************************************************************
 * ScaleCoordinates
 */

ScaleCoordinates::ScaleCoordinates(const RVec& scale) : impl_{ new Impl(scale) } {}

void ScaleCoordinates::operator()(ArrayRef<RVec> coordinates) const
{
    impl_->scale(coordinates);
}

void ScaleCoordinates::operator()(RVec* coordinate) const
{
    impl_->scale({ coordinate, coordinate + 1 });
}

void ScaleCoordinates::inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const
{
    impl_->inverseIgnoringZeroScale(coordinates);
}

void ScaleCoordinates::inverseIgnoringZeroScale(RVec* coordinate) const
{
    impl_->inverseIgnoringZeroScale({ coordinate, coordinate + 1 });
}

ScaleCoordinates::~ScaleCoordinates() = default;

ScaleCoordinates::ScaleCoordinates(const ScaleCoordinates& other) : impl_(new Impl(*other.impl_)) {}

ScaleCoordinates& ScaleCoordinates::operator=(const ScaleCoordinates& other)
{
    *impl_ = *other.impl_;
    return *this;
}

ScaleCoordinates::ScaleCoordinates(ScaleCoordinates&&) noexcept = default;

ScaleCoordinates& ScaleCoordinates::operator=(ScaleCoordinates&&) noexcept = default;

/********************************************************************
 * TranslateAndScale::Impl
 */

class TranslateAndScale::Impl
{
public:
    Impl(const RVec& scale, const RVec& translation);
    void transform(ArrayRef<RVec> coordinates) const;
    RVec scale_;
    RVec translation_;
};

TranslateAndScale::Impl::Impl(const RVec& scale, const RVec& translation) :
    scale_{ scale }, translation_{ translation }
{
}

void TranslateAndScale::Impl::transform(ArrayRef<RVec> coordinates) const
{
    for (auto& coordinate : coordinates)
    {
        coordinate += translation_;
        coordinate[XX] *= scale_[XX];
        coordinate[YY] *= scale_[YY];
        coordinate[ZZ] *= scale_[ZZ];
    }
}


/********************************************************************
 * TranslateAndScale
 */

TranslateAndScale::TranslateAndScale(const RVec& scale, const RVec& translation) :
    impl_(new Impl(scale, translation))
{
}

void TranslateAndScale::operator()(ArrayRef<RVec> coordinates) const
{
    impl_->transform(coordinates);
}

void TranslateAndScale::operator()(RVec* coordinate) const
{
    impl_->transform({ coordinate, coordinate + 1 });
}

ScaleCoordinates TranslateAndScale::scaleOperationOnly() const
{
    return ScaleCoordinates{ impl_->scale_ };
}

TranslateAndScale::~TranslateAndScale() = default;

TranslateAndScale::TranslateAndScale(const TranslateAndScale& other) : impl_(new Impl(*other.impl_))
{
}

TranslateAndScale& TranslateAndScale::operator=(const TranslateAndScale& other)
{
    *impl_ = *other.impl_;
    return *this;
}

TranslateAndScale::TranslateAndScale(TranslateAndScale&&) noexcept = default;

TranslateAndScale& TranslateAndScale::operator=(TranslateAndScale&&) noexcept = default;

/********************************************************************
 * AffineTransformation
 */

AffineTransformation::AffineTransformation(Matrix3x3ConstSpan mat, const RVec& translation) :
    translation_{ translation }
{
    std::copy(begin(mat), end(mat), begin(matrix_));
}

void AffineTransformation::operator()(ArrayRef<RVec> vectors) const
{
    for (RVec& vector : vectors)
    {
        matrixVectorMultiply(matrix_.asConstView(), &vector);
        vector += translation_;
    }
}

void AffineTransformation::operator()(RVec* vector) const
{
    (*this)({ vector, vector + 1 });
}

Matrix3x3 AffineTransformation::gradient() const
{
    return transpose(matrix_);
}

} // namespace gmx
