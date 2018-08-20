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

#ifndef GROMACS_RESTRAINT_VECTORTYPE_H
#define GROMACS_RESTRAINT_VECTORTYPE_H

/*! \file
 * \brief Template header for 3D vector types and operations.
 *
 * Insulate client code from changes to vector data type conventions in the
 * library.
 *
 * Motivation:
 *
 * 1. Make data types and precision explicit and unambiguous.
 * 2. Provide an abstraction from storage method.
 * 3. Avoid pointer dereferencing.
 * 4. Avoid indexing and bounds-checking.
 * 5. Allow more compiler optimizations.
 *
 * These types should map easily to float3 (or float4) as in CUDA and other libraries,
 * as well as to arrays or even non-contiguous structures, at least insofar as
 * the compiler should be able to optimize away copies. When it can't,
 * a reinterpret_cast<>() in an inline helper function should do the trick.
 * But the template nature of the type and inline operators should generally
 * optimize away.
 *
 * Along these lines, the structures are intended to be short-lived handles
 * for convenience and strong typing of operations. Arrays of vec3 should not
 * be necessary and are probably not desirable, at least across C++ translation units.
 *
 * \ingroup module_restraint
 */

// TODO: It's worth taking a look at how compilers handle iterative extraction of vec3 from Nx3 data in practice.

#include <cassert>
#include <cmath>

#include <ostream>

#include "gromacs/math/vectypes.h"

namespace gmx
{
namespace detail
{

/*!
 * \brief 3-dimensional vector types.
 *
 * Provide strongly-typed vector type for unambiguous operations. Maps easily to CUDA float3.
 *
 * Conversions to other types with different semantics should be explicit and avoid accidental loss of precision.
 * \tparam Scalar
 *
 * \inlibraryapi
 */
template<typename Scalar>
class vec3
{
    public:
        /// magnitude in first dimension
        Scalar x;
        /// magnitude in second dimension
        Scalar y;
        /// magnitude in third dimension
        Scalar z;

        /*!
         * \brief Require type matching for direct construction.
         *
         * \param X
         * \param Y
         * \param Z
         */
        constexpr explicit vec3(const Scalar &X, const Scalar &Y, const Scalar &Z) : x {X}, y {Y}, z {Z} {};

        /// \brief Default construct a zero-vector.
        vec3() : vec3{Scalar(0), Scalar(0), Scalar(0)}
        {};

        /// \brief Default copy construct from another vector.
        vec3(const vec3 &) = default;
        /*!
         * \brief Copy assign a vec3.
         *
         * \return reference
         */
        vec3 &operator=(const vec3 &) = default;
        /*!
         * \brief Move assign a vec3.
         *
         * \return reference
         */
        vec3 &operator=(vec3 &&) noexcept = default;

        /*!
         * \brief Implicit non-narrowing conversion between vec3<>
         *
         * \tparam T
         * \return converted vec3 constructed from the original x, y, z members.
         */
        template<typename T>
        explicit operator vec3<T>() { return vec3<T>(x, y, z); }

        /*!
         * \brief Templated constructor.
         *
         * If vec3<Scalar> can be constructed from a list of three numbers of type T, construct from vec3<T>.
         * \tparam T source vector type.
         * \param a source vector.
         */
        template<typename T>
        explicit vec3(const T &a) : vec3(a.x, a.y, a.z) {};
};

//
// Arithmetic
//
// A common idiom in vector math libraries is to overload operator*(), operator/(), and operator%(),
// but in the context of matrices and tensor algebra, it is not unambiguous whether multiplication
// should imply dot product. Let's stick to explicit free functions.
//

/*!
 * \brief Unary negation operator
 *
 * \tparam Scalar underlying vector type
 * \param a input vector
 * \return (-a.x, -a.y, -a.z)
 * \inlibraryapi
 */
template<typename Scalar>
inline vec3<Scalar> operator-(const vec3<Scalar> &a)
{
    return vec3<Scalar>(-a.x, -a.y, -a.z);
}

/*!
 * \brief Binary addition operator.
 *
 * \tparam Scalar
 * \param a first vector
 * \param b second vector
 * \return resulting vector
 * \inlibraryapi
 */
template<typename Scalar>
inline vec3<Scalar> operator+(const vec3<Scalar> &a, const vec3<Scalar> &b)
{
    return vec3<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z);
};

/*!
 * \brief Binary subtraction operator.
 *
 * \tparam Scalar
 * \param a first vector
 * \param b second vector
 * \return resulting vector
 * \inlibraryapi
 */
template<typename Scalar>
inline vec3<Scalar> operator-(const vec3<Scalar> &a, const vec3<Scalar> &b)
{
    return vec3<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z);
}


/*!
 * \brief Multiply vector by scalar
 *
 * \tparam Scalar vector type
 * \param a input vector
 * \param s input scalar
 * \return (a.x, a.y, a.z) * s
 *
 * Note that input scalar may be implicitly narrowed if higher precision than input vector.
 * \inlibraryapi
 *
 * \{
 */
template<typename T1, typename T2>
inline vec3<T1> operator*(const vec3<T1> &a, const T2 &s)
{
    return vec3<T1>(T1(a.x * s), T1(a.y * s), T1(a.z * s));
}
template<typename T1, typename T2>
inline vec3<T1> operator*(const T2 &s, const vec3<T1> &a)
{
    return a * s;
};
/*! \} */

/*!
 * \brief Vector division by scalar
 * \tparam Scalar underlying vector type
 * \param a input vector
 * \param s input scalar
 * \return
 *
 * Note that input scalar may be implicitly narrowed if higher precision than input vector.
 * \inlibraryapi
 */
template<typename T1, typename T2>
inline vec3<T1> operator/(const vec3<T1> &a, const T2 &s)
{
    assert(s != T2(0.0));
    const T2 inv_s {
        T2(1.0)/s
    };
    return vec3<T1>(T1(a.x * inv_s), T1(a.y * inv_s), T1(a.z * inv_s));
}

/*!
 * \brief Scalar vector product
 * \tparam Scalar underlying vector type
 * \param a first vector
 * \param b second vector
 * \return (a.x * b.x) + (a.y * b.y) + (a.z * b.z)
 * \inlibraryapi
 */
template<typename Scalar>
inline Scalar dot(const vec3<Scalar> &a, const vec3<Scalar> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/*!
 * \brief Norm or magnitude of vector
 *
 * \param v input vector
 * \return magnitude of v
 * To specify the precision of the calculation and result, choose the output type template parameter.
 *
 * Example
 *
 *     constexpr vec3<float> v{1,0,0};
 *     float magnitude = norm(v);
 *     double magnitude = norm<double>(v);
 *
 * \inlibraryapi
 */
template<typename Scalar, typename Tout = Scalar>
inline Tout norm(const vec3<Scalar> &v)
{
    return sqrt(dot(vec3<Tout>(v), vec3<Tout>(v)));
}

//
// Comparisons
//

/*!
 * \brief Equality comparison operator
 *
 * \tparam Scalar
 * \param a
 * \param b
 * \return true if all elements of vectors are arithmetically equal.
 * \inlibraryapi
 */
template<typename S1, typename S2>
bool operator==(const vec3<S1> &a, const vec3<S2> &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

//
// Conversions
//

//
// Helpers
//

/*!
 * \brief Flexibly produce vector of a given type.
 *
 * \tparam Scalar underlying output vector type
 * \tparam T1 x input type
 * \tparam T2 y input type
 * \tparam T3 z input type
 * \param x
 * \param y
 * \param z
 * \return vector (x, y, z) of type vec3<Scalar>
 *
 * Helper function allows narrowing and mismatched types. Constructs a vec3<Scalar> from any x, y, and z that can be
 * implicitly converted to type Scalar.
 * \inlibraryapi
 */
template<typename Scalar, typename T1, typename T2, typename T3>
inline constexpr vec3<Scalar> make_vec3(T1 x, T2 y, T3 z)
{
    return vec3<Scalar>(Scalar(x), Scalar(y), Scalar(z));
};

/*!
 * \brief Helper operator overload to format string output.
 *
 * \tparam T a valid type with which to template vec3<>
 * \param stream reference to output stream
 * \param vec vector to format
 * \return string representation of a vector as "(x,y,z)"
 */
template<class T>
std::ostream &operator<<(std::ostream &stream, const vec3<T> &vec)
{
    stream << "(" << vec.x << "," << vec.y << "," << vec.z << ")";
    return stream;
}

}      // end namespace gmx::detail
}      // end namespace gmx

#endif //GROMACS_VECTORTYPE_H
