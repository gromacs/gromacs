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
/*! \internal
 * \file
 * \brief
 * Define coordinate transformations.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_COORDINATETRANSFORMATIONS_H
#define GMX_MATH_COORDINATETRANSFORMATIONS_H

#include <algorithm>
#include <functional>
#include <numeric>

namespace gmx
{

/*! \internal \brief Coordinate transformation base class.
 * Provides cooridnate transformation, it's inverse and the transformations volume element at the origin.
 * \tparam Transformation The coordinate transformation
 * \tparam CoordinateType The coordinate type used in the transformation, used explicitely
 *                        here to avoid introducing traits since this is the only type that
 *                        needs to be known to the base class apart from the actual Transformation
 */
template<class Transformation, class CoordinateType>
class CoordinateTransformationBase
{
    public:
        //! The base constructor, asserts that general conditions for all coordinate transformations are met.
        CoordinateTransformationBase()
        {
            static_assert(std::tuple_size<CoordinateType>::value,
                          "Coordinate transformations requires a coordinate type with fixed compile-time size.");
        }
        /*! \brief Apply coordinate transformation to input coordinate.
         * \param[in] x The input coordinate
         * \returns the transformed coordinate
         */
        CoordinateType operator()(const CoordinateType &x) const
        {
            return Transformation(x);
        }
        /*! \brief Apply inverse coordinate transformation to input coordinate.
         * inverse(operator()(x)) shall be equivalent to x.
         * \param[in] x The input coordinate
         * \returns the inverse transformed coordinate
         */
        CoordinateType inverse(const CoordinateType &x) const
        {
            return Transformation::inverse(x);
        }
        /*!\brief The volume element of the transformation at zero.
         * Usually the Jacobian determinant of the transformation at the origin.
         * For linear scaling, this is the unit cell volume.
         * \note a more general volumeElement at an arbitrary coordinate is not implemented, because it is not needed for now.
         * \retuns the volume element at the origin
         */
        double volumeElementAtZero() const
        {
            return Transformation::volumeElementAtZero();
        }
};
/*!\internal \brief Compose two coordiante transformations.
 * \tparam CoordinateType The type of coordinate to be transformed
 * \tparam FirstTransformation The first transformation to be carried out
 * \tparam SecondTransformation The second transformation to be carried out.
 */
template <class FirstTransformation, class SecondTransformation>
class CompositeTransform :
    public CoordinateTransformationBase<
        CompositeTransform<FirstTransformation, SecondTransformation>,
        typename FirstTransformation::coordinate_type>
{
    public:
        using coordinate_type = typename FirstTransformation::coordinate_type;
        //! Construct the composite transformation by copying the two transformations to be combined.
        CompositeTransform(const FirstTransformation &first, const SecondTransformation &second) :
            first_ {first},
        second_ {second}
        { }

        //! \copydoc CoordinateTransformationBase::operator()
        constexpr coordinate_type operator()(const coordinate_type &x) const
        {
            return second_(first_(x));
        }
        //! \copydoc CoordinateTransformationBase::inverse()
        constexpr coordinate_type inverse(const coordinate_type &x) const
        {
            return first_.inverse(second_.inverse(x));
        }
        //! \copydoc CoordinateTransformationBase::volumeElementAtZero()
        constexpr double volumeElementAtZero() const
        {
            return first_.volumeElementAtZero() * second_.volumeElementAtZero();
        }

    private:
        //! The first transformation to be applied
        FirstTransformation  first_;
        //! The second transformation to be applied
        SecondTransformation second_;
};

/*!\internal \brief Builder function for composite coordinate transforms.
 * Enables automated template parameter deduction for composite transforms.
 */
template<class FirstTransformation, class SecondTransformation>
CompositeTransform<FirstTransformation, SecondTransformation>
composeTransforms(const FirstTransformation &first, const SecondTransformation &second)
{
    return {first, second};
}

/*! \internal \brief Scale input coordiantes.
 * Equivalent to multiplication with diagonal matrix.
 * \tparam CoordinateType the type of the input coordinate and scale.
 */
template <typename CoordinateType>
class ScaleTransform :
    public CoordinateTransformationBase<ScaleTransform<CoordinateType>, CoordinateType>
{
    public:
        using coordinate_type = CoordinateType;

        /*!\brief Construct by setting a scale.
         * NOTE Undefined behaviour if any scale parameter is zero.
         */
        ScaleTransform(const coordinate_type &scale) : scale_ {scale},
        scaleInverse_ {calculateInverse(scale_)}
        {}
        //! \copydoc CoordinateTransformationBase::operator()
        constexpr coordinate_type operator()(const coordinate_type &x) const
        {
            return elementwiseProduct(x, scale_);
        }
        //! \copydoc CoordinateTransformationBase::inverse()
        constexpr coordinate_type inverse(const coordinate_type &x) const
        {
            return elementwiseProduct(x, scaleInverse_);
        }
        //! \copydoc CoordinateTransformationBase::volumeElementAtZero()
        double volumeElementAtZero() const
        {
            return std::accumulate(std::begin(scale_), std::end(scale_),
                                   1.0, std::multiplies<double>());
        }

    private:
        //! Multiply two coordiantes element by element.
        coordinate_type elementwiseProduct(const coordinate_type &x, const coordinate_type &y) const
        {
            coordinate_type result;
            std::transform(std::begin(x), std::end(x), std::begin(y), std::begin(result),
                           std::multiplies<typename coordinate_type::value_type>() );
            return result;
        }
        /*! \brief Invert coordinate values.
         * Does not check on the validity of inverting the values in x.
         */
        coordinate_type calculateInverse(const coordinate_type &x)
        {
            coordinate_type result;
            std::transform(std::begin(x), std::end(x), std::begin(result),
                           [](const typename  coordinate_type::value_type &xValue)
                           {return 1./xValue; } );
            return result;
        }
        //! The scale for each input coordinate
        coordinate_type scale_;
        //! Storing the inverse scale for each input coordinate avoids repeat 1/x operation for inverse scaling
        coordinate_type scaleInverse_;
};

/*!\internal \brief Builder function for scaling coordinate transforms.
 * Enables automated template parameter deduction for scaling transforms.
 */
template< typename CoordinateType>
ScaleTransform<CoordinateType>
makeScaleTransform(const CoordinateType &scale)
{
    return scale;
}

/*!\internal \brief Translational coordinate system transformation.
 * \tparam CoordinateType The type of coordiante to be transformed, as well as
 * for defining the transformation.
 */
template <typename CoordinateType>
class TranslationTransform :
    public CoordinateTransformationBase<TranslationTransform<CoordinateType>, CoordinateType>
{
    public:
        using coordinate_type = CoordinateType;

        //! The default constructor initializes the translation to zero.
        constexpr TranslationTransform() : translation_ {}
        { }
        //! Construct with translation vector.
        TranslationTransform(const coordinate_type &translation) : translation_ {translation}
        {
        }
        //! \copydoc CoordinateTransformationBase::operator()
        coordinate_type operator()(const coordinate_type &x) const
        {
            coordinate_type result;
            std::transform(std::begin(x), std::end(x), std::begin(translation_),
                           std::begin(result), std::plus<typename coordinate_type::value_type>());
            return result;
        }
        //! \copydoc CoordinateTransformationBase::inverse()
        coordinate_type inverse(const coordinate_type &x) const
        {
            coordinate_type result;
            std::transform(std::begin(x), std::end(x), std::begin(translation_),
                           std::begin(result), std::minus<typename coordinate_type::value_type>());
            return result;
        }
        //! \copydoc CoordinateTransformationBase::volumeElementAtZero()
        constexpr double volumeElementAtZero() const
        {
            // Translation transformations do not change the volume element.
            return 1;
        }
    private:
        //! The external coordinates of (0,...,0)
        coordinate_type translation_;
};

/*!\internal \brief Builder function for translation transforms.
 * Enables automated template parameter deduction for translation transforms.
 */
template< typename CoordinateType>
TranslationTransform<CoordinateType>
makeTranslationTransform(const CoordinateType &translation)
{
    return translation;
}

} // namespace gmx

#endif
