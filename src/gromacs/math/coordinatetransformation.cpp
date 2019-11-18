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
 * \brief Implements coordinate transformation routines.
 *
 * \author Christian Blau <blau@kth.se>
 */
#include "gmxpre.h"

#include "coordinatetransformation.h"

#include <vector>

#include "gromacs/math/vec.h"

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

void ScaleCoordinates::inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const
{
    impl_->inverseIgnoringZeroScale(coordinates);
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
    scale_{ scale },
    translation_{ translation }
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


} // namespace gmx
