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
/*! \libinternal \file
 * \brief Declares coordinate transformation routines.
 *
 * \author Christian Blau <blau@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_COORDINATETRANSFORMATION_H
#define GMX_MATH_COORDINATETRANSFORMATION_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class ScaleCoordinates
{
public:
    //! Set up coordinate scaling with the scaling factor in each dimension
    explicit ScaleCoordinates(const RVec& scale);
    ~ScaleCoordinates();

    //! Copy constructor
    ScaleCoordinates(const ScaleCoordinates& other);
    //! Copy assignment
    ScaleCoordinates& operator=(const ScaleCoordinates& other);
    //! Move constructor
    ScaleCoordinates(ScaleCoordinates&& other) noexcept;
    //! Move assignment
    ScaleCoordinates& operator=(ScaleCoordinates&& other) noexcept;

    /*! \brief Perform a coordinate transformation on input coordinates.
     * \param[in] coordinates to be transformed
     */
    void operator()(ArrayRef<RVec> coordinates) const;
    /*! \brief Apply the inverse scale to coordinates, ignoring dimensions for which scale is zero.
     * \param[in] coordinates to be transformed
     */
    void inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const;

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief Transform coordinates in three dimensions by first
 * translating, then scaling them.
 */
class TranslateAndScale
{
public:
    /*! \brief Construct a three-dimensional coordinate transformation.
     * Coordinates are first translated, then scaled.
     * \param[in] translation to be performed on the coordinates
     * \param[in] scale to be applied to the coordinates
     */
    TranslateAndScale(const RVec& scale, const RVec& translation);

    ~TranslateAndScale();

    //! Copy constructor
    TranslateAndScale(const TranslateAndScale& other);
    //! Copy assignment
    TranslateAndScale& operator=(const TranslateAndScale& other);
    //! Move constructor
    TranslateAndScale(TranslateAndScale&& other) noexcept;
    //! Move assignment
    TranslateAndScale& operator=(TranslateAndScale&& other) noexcept;

    /*! \brief Perform a coordinate transformation on input coordinates.
     * \param[in] coordinates to be transformed
     */
    void operator()(ArrayRef<RVec> coordinates) const;

    /*! \brief Returns the scaling operation, discarding the translation.
     */
    ScaleCoordinates scaleOperationOnly() const;

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};
} // namespace gmx
#endif // CoordinateTransformation
