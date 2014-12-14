/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifndef GMX_MATH_VECTYPES_H
#define GMX_MATH_VECTYPES_H

#include "gromacs/utility/real.h"

#define XX      0 /* Defines for indexing in */
#define YY      1 /* vectors                 */
#define ZZ      2
#define DIM     3 /* Dimension of vectors    */

typedef real    rvec[DIM];

typedef double  dvec[DIM];

typedef real    matrix[DIM][DIM];

typedef real    tensor[DIM][DIM];

typedef int     ivec[DIM];

typedef int     imatrix[DIM][DIM];

#ifdef __cplusplus

namespace gmx
{

/*! \brief
 * C++ class for 3D vectors.
 *
 * \tparam ValueType  Type
 *
 * This class provides a C++ version of rvec/dvec/ivec that can be put into STL
 * containers etc.  It is more or less a drop-in replacement for `rvec` and
 * friends: it can be used in most contexts that accept the equivalent C type.
 * However, there are two cases where explicit conversion is necessary:
 *  - An array of these objects needs to be converted with as_vec_array() (or
 *    convenience methods like as_rvec_array()).
 *  - Passing an RVec as a `const rvec &` parameter to a function needs an
 *    explicit call to as_vec().  The implicit conversion should work for this
 *    as well, but cppcheck parses the necessary implicit conversion operator
 *    incorrectly and MSVC fails to compile code that relies on the implicit
 *    conversion, so the explicit method is necessary.
 *
 * For the array conversion to work, the compiler should not add any extra
 * alignment/padding in the layout of this class;  that this actually works as
 * intended is tested in the unit tests.
 *
 * \inpublicapi
 */
template <typename ValueType>
class BasicVector
{
    public:
        //! Underlying raw C array type (rvec/dvec/ivec).
        typedef ValueType RawArray[DIM];

        //! Constructs default (uninitialized) vector.
        BasicVector() {}
        //! Constructs a vector from given values.
        BasicVector(ValueType x, ValueType y, ValueType z)
        {
            x_[XX] = x;
            x_[YY] = y;
            x_[ZZ] = z;
        }
        /*! \brief
         * Constructs a vector from given values.
         *
         * This constructor is not explicit to support implicit conversions
         * that allow, e.g., calling `std::vector<RVec>:``:push_back()` directly
         * with an `rvec` parameter.
         */
        BasicVector(const RawArray x)
        {
            x_[XX] = x[XX];
            x_[YY] = x[YY];
            x_[ZZ] = x[ZZ];
        }
        //! Indexing operator to make the class work as the raw array.
        ValueType &operator[](int i) { return x_[i]; }
        //! Indexing operator to make the class work as the raw array.
        ValueType operator[](int i) const { return x_[i]; }
        // The conversion functions below could more accurately return
        // RawArray &, but this fails with cppcheck and does not solve the
        // issue with MSVC, so as_vec() should be used instead.
        //! Makes BasicVector usable in contexts where a raw C array is expected.
        operator ValueType *() { return x_; }
        //! Makes BasicVector usable in contexts where a raw C array is expected.
        operator const ValueType *() const { return x_; }

        //! Converts to a raw C array where implicit conversion does not work.
        RawArray &as_vec() { return x_; }
        //! Converts to a raw C array where implicit conversion does not work.
        const RawArray &as_vec() const { return x_; }

    private:
        RawArray x_;
};

/*! \brief
 * Casts a gmx::BasicVector array into an equivalent raw C array.
 */
template <typename ValueType> static inline
typename BasicVector<ValueType>::RawArray *
as_vec_array(BasicVector<ValueType> *x)
{
    return reinterpret_cast<typename BasicVector<ValueType>::RawArray *>(x);
}

/*! \brief
 * Casts a gmx::BasicVector array into an equivalent raw C array.
 */
template <typename ValueType> static inline
const typename BasicVector<ValueType>::RawArray *
as_vec_array(const BasicVector<ValueType> *x)
{
    return reinterpret_cast<const typename BasicVector<ValueType>::RawArray *>(x);
}

//! Shorthand for C++ `rvec`-equivalent type.
typedef BasicVector<real> RVec;
//! Casts a gmx::RVec array into an `rvec` array.
static inline rvec *as_rvec_array(RVec *x)
{
    return as_vec_array(x);
}
//! Casts a gmx::RVec array into an `rvec` array.
static inline const rvec *as_rvec_array(const RVec *x)
{
    return as_vec_array(x);
}

} // namespace gmx

#endif

#endif
