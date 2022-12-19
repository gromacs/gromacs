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

#include <memory>

#include "gromacs/math/vectypes.h"

#include "matrix.h"

namespace gmx
{

template<typename>
class ArrayRef;
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

    /*! \brief Perform a coordinate transformation on an input coordinate.
     * \param[in] coordinate to be transformed
     */
    void operator()(RVec* coordinate) const;

    /*! \brief Apply the inverse scale to coordinates, ignoring dimensions for which scale is zero.
     * \param[in] coordinates to be transformed
     */
    void inverseIgnoringZeroScale(ArrayRef<RVec> coordinates) const;

    /*! \brief Apply the inverse scale to a coordinate, ignoring dimensions for which scale is zero.
     * \param[in] coordinate to be transformed
     */
    void inverseIgnoringZeroScale(RVec* coordinate) const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
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

    /*! \brief Perform a coordinate transformation on a coordinate.
     * \param[in] coordinate to be transformed
     */
    void operator()(RVec* coordinate) const;

    /*! \brief Returns the scaling operation, discarding the translation.
     */
    ScaleCoordinates scaleOperationOnly() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/*! \libinternal
 * \brief Affine transformation of three-dimensional coordinates.
 *
 * Perfoms in-place coordinate transformations.
 *
 * Coordinates are first multiplied by a matrix, then translated.
 */
class AffineTransformation
{
public:
    /*! \brief Construct a three-dimensional affine transformation.
     * \param[in] matrix to be applied to the vectors
     * \param[in] translation to be performed on the vectors
     */
    AffineTransformation(Matrix3x3ConstSpan matrix, const RVec& translation);

    /*! \brief Perform an affine transformation on input vectors.
     * \param[in,out] vectors to be transformed in-place
     */
    void operator()(ArrayRef<RVec> vectors) const;

    /*! \brief Perform an affine transformation on a vector.
     * \param[in,out] vector to be transformed in-place
     */
    void operator()(RVec* vector) const;

    /*! \brief Return the gradient of this affine transformation.
     *
     * The gradient of an affine transformation is the transpose of the
     * linear transformation matrix. See, eg, here
     * https://inst.eecs.berkeley.edu/~ee127/sp21/livebook/def_gradient.html
     *
     * \returns matrix that acts as gradient
     */
    Matrix3x3 gradient() const;

private:
    //! The matrix describing the affine transformation A(x) = matrix_ * x + translation_
    Matrix3x3 matrix_;
    //! The translation vector describing the affine transformation A(x) = matrix * x + translation
    RVec translation_;
};

} // namespace gmx
#endif // GMX_MATH_COORDINATETRANSFORMATION_H
