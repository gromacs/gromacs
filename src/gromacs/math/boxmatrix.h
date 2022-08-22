/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * \brief
 * Declares routines to operate on 3x3 matrices that are boxes,
 * ie. have zero values in their upper right triangle.
 *
 * This header is distinct from matrix.h so that we can still inline
 * these operations while not requiring headers that merely need the
 * name Matrix3x3 to pull in the dependencies needed here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 * \inlibraryapi
 */
#ifndef GMX_MATH_BOXMATRIX_H
#define GMX_MATH_BOXMATRIX_H

#include "gromacs/math/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \brief Assert that the matrix \c m describes a simulation box
 *
 * The GROMACS convention is that all simulation box descriptions are
 * normalized to have zero entries in the upper triangle. This function
 * asserts if that is not true. */
static inline void assertMatrixIsBoxMatrix(const Matrix3x3 gmx_used_in_debug& m)
{
    GMX_ASSERT((m(XX, YY) == 0.0) && (m(XX, ZZ) == 0.0) && (m(YY, ZZ) == 0.0),
               gmx::formatString(
                       "Box matrix should contain zero in the upper triangle, but "
                       "was\n%10.6g %10.6g %10.6g\n%10.6g %10.6g %10.6g\n%10.6g %10.6g %10.6g\n",
                       m(0, 0),
                       m(0, 1),
                       m(0, 2),
                       m(1, 0),
                       m(1, 1),
                       m(1, 2),
                       m(2, 0),
                       m(2, 1),
                       m(2, 2))
                       .c_str());
}

/*! \brief Assert that the matrix \c m describes a simulation box
 *
 * The GROMACS convention is that all simulation box descriptions are
 * normalized to have zero entries in the upper triangle. This function
 * asserts if that is not true. */
static inline void assertMatrixIsBoxMatrix(const matrix gmx_used_in_debug m)
{
    GMX_ASSERT((m[XX][YY] == 0.0) && (m[XX][ZZ] == 0.0) && (m[YY][ZZ] == 0.0),
               gmx::formatString(
                       "Box matrix should contain zero in the upper triangle, but was\n%10.6g "
                       "%10.6g %10.6g\n%10.6g %10.6g %10.6g\n%10.6g %10.6g %10.6g\n",
                       m[0][0],
                       m[0][1],
                       m[0][2],
                       m[1][0],
                       m[1][1],
                       m[1][2],
                       m[2][0],
                       m[2][1],
                       m[2][2])
                       .c_str());
}

/*! \brief Multiply vector \c v by the transpose of the box matrix \c
 * m, relying on the fact that the upper triangle of the matrix is
 * zero. */
static inline void multiplyVectorByTransposeOfBoxMatrix(const Matrix3x3& m, const rvec v, rvec result)
{
    assertMatrixIsBoxMatrix(m);
    result[XX] = m(XX, XX) * v[XX] + m(YY, XX) * v[YY] + m(ZZ, XX) * v[ZZ];
    result[YY] = m(YY, YY) * v[YY] + m(ZZ, YY) * v[ZZ];
    result[ZZ] = m(ZZ, ZZ) * v[ZZ];
}

/*! \brief Multiply vector \c v by the transpose of the box matrix \c
 * m, relying on the fact that the upper triangle of the matrix is
 * zero. */
static inline RVec multiplyVectorByTransposeOfBoxMatrix(const Matrix3x3& m, const RVec& v)
{
    assertMatrixIsBoxMatrix(m);

    RVec dest;
    dest[XX] = m(XX, XX) * v[XX] + m(YY, XX) * v[YY] + m(ZZ, XX) * v[ZZ];
    dest[YY] = m(YY, YY) * v[YY] + m(ZZ, YY) * v[ZZ];
    dest[ZZ] = m(ZZ, ZZ) * v[ZZ];
    return dest;
}

//! Multiply box matrices, relying on the fact that their upper triangle is zero
static inline Matrix3x3 multiplyBoxMatrices(const Matrix3x3& a, const Matrix3x3& b)
{
    assertMatrixIsBoxMatrix(a);
    assertMatrixIsBoxMatrix(b);
    Matrix3x3 result;
    result(XX, XX) = a(XX, XX) * b(XX, XX);
    result(XX, YY) = 0.0;
    result(XX, ZZ) = 0.0;
    result(YY, XX) = a(YY, XX) * b(XX, XX) + a(YY, YY) * b(YY, XX);
    result(YY, YY) = a(YY, YY) * b(YY, YY);
    result(YY, ZZ) = 0.0;
    result(ZZ, XX) = a(ZZ, XX) * b(XX, XX) + a(ZZ, YY) * b(YY, XX) + a(ZZ, ZZ) * b(ZZ, XX);
    result(ZZ, YY) = a(ZZ, YY) * b(YY, YY) + a(ZZ, ZZ) * b(ZZ, YY);
    result(ZZ, ZZ) = a(ZZ, ZZ) * b(ZZ, ZZ);
    return result;
}

//! Multiply box matrices, relying on the fact that their upper triangle is zero
static inline Matrix3x3 multiplyBoxMatrices(const Matrix3x3& a, const matrix b)
{
    assertMatrixIsBoxMatrix(a);
    assertMatrixIsBoxMatrix(b);
    Matrix3x3 result;
    result(XX, XX) = a(XX, XX) * b[XX][XX];
    result(XX, YY) = 0.0;
    result(XX, ZZ) = 0.0;
    result(YY, XX) = a(YY, XX) * b[XX][XX] + a(YY, YY) * b[YY][XX];
    result(YY, YY) = a(YY, YY) * b[YY][YY];
    result(YY, ZZ) = 0.0;
    result(ZZ, XX) = a(ZZ, XX) * b[XX][XX] + a(ZZ, YY) * b[YY][XX] + a(ZZ, ZZ) * b[ZZ][XX];
    result(ZZ, YY) = a(ZZ, YY) * b[YY][YY] + a(ZZ, ZZ) * b[ZZ][YY];
    result(ZZ, ZZ) = a(ZZ, ZZ) * b[ZZ][ZZ];
    return result;
}

/*! \brief Invert a simulation-box matrix
 *
 * This routine assumes that src is a simulation-box matrix, i.e. has
 * zeroes in the upper-right triangle.
 *
 * \throws RangeError if the product of the leading diagonal is too small.
 */
static inline Matrix3x3 invertBoxMatrix(const Matrix3x3& src)
{
    assertMatrixIsBoxMatrix(src);

    double tmp = src(XX, XX) * src(YY, YY) * src(ZZ, ZZ);
    if (std::fabs(tmp) <= 100 * GMX_REAL_MIN)
    {
        GMX_THROW(RangeError("Cannot invert matrix, determinant is too close to zero"));
    }

    Matrix3x3 dest;
    dest(XX, XX) = 1 / src(XX, XX);
    dest(YY, YY) = 1 / src(YY, YY);
    dest(ZZ, ZZ) = 1 / src(ZZ, ZZ);
    dest(ZZ, XX) = (src(YY, XX) * src(ZZ, YY) * dest(YY, YY) - src(ZZ, XX)) * dest(XX, XX) * dest(ZZ, ZZ);
    dest(YY, XX) = -src(YY, XX) * dest(XX, XX) * dest(YY, YY);
    dest(ZZ, YY) = -src(ZZ, YY) * dest(YY, YY) * dest(ZZ, ZZ);
    dest(XX, YY) = 0.0;
    dest(XX, ZZ) = 0.0;
    dest(YY, ZZ) = 0.0;
    return dest;
}

/*! \brief Invert a simulation-box matrix in \c src, return in \c dest
 *
 * This routine assumes that src is a simulation-box matrix, i.e. has
 * zeroes in the upper-right triangle. A fatal error occurs if the
 * product of the leading diagonal is too small. The inversion can be
 * done "in place", i.e \c src and \c dest can be the same matrix.
 */
static inline void invertBoxMatrix(const matrix src, matrix dest)
{
    assertMatrixIsBoxMatrix(src);

    double tmp = src[XX][XX] * src[YY][YY] * src[ZZ][ZZ];
    if (std::fabs(tmp) <= 100 * GMX_REAL_MIN)
    {
        GMX_THROW(RangeError("Cannot invert matrix, determinant is too close to zero"));
    }

    dest[XX][XX] = 1 / src[XX][XX];
    dest[YY][YY] = 1 / src[YY][YY];
    dest[ZZ][ZZ] = 1 / src[ZZ][ZZ];
    dest[ZZ][XX] = (src[YY][XX] * src[ZZ][YY] * dest[YY][YY] - src[ZZ][XX]) * dest[XX][XX] * dest[ZZ][ZZ];
    dest[YY][XX] = -src[YY][XX] * dest[XX][XX] * dest[YY][YY];
    dest[ZZ][YY] = -src[ZZ][YY] * dest[YY][YY] * dest[ZZ][ZZ];
    dest[XX][YY] = 0.0;
    dest[XX][ZZ] = 0.0;
    dest[YY][ZZ] = 0.0;
}

} // namespace gmx

#endif
