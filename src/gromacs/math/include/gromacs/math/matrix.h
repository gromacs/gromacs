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
/*! \libinternal
 * \file
 * \brief Declares special case of 3x3 matrix frequently used, and associated functions.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_MATRIX_H_
#define GMX_MATH_MATRIX_H_

#include <array>

#include "gromacs/math/multidimarray.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief Three-by-three matrix of ElementType.
 * \tparam ElementType type of element to be stored in matrix
 */
template<class ElementType>
using BasicMatrix3x3 = MultiDimArray<std::array<ElementType, 3 * 3>, extents<3, 3>>;

/*! \brief Three-by-three real number matrix.
 * \note will replace the C-style real[3][3] "matrix"
 */
using Matrix3x3 = BasicMatrix3x3<real>;
//! Convenience alias for a matrix view
using Matrix3x3Span = Matrix3x3::view_type;
//! Convenience alias for a const matrix view
using Matrix3x3ConstSpan = Matrix3x3::const_view_type;

//! Determinant of a 3x3 matrix
constexpr real determinant(Matrix3x3ConstSpan matrix)
{
    return (matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(2, 1) * matrix(1, 2))
            - matrix(1, 0) * (matrix(0, 1) * matrix(2, 2) - matrix(2, 1) * matrix(0, 2))
            + matrix(2, 0) * (matrix(0, 1) * matrix(1, 2) - matrix(1, 1) * matrix(0, 2)));
}

//! Calculates the trace of a 3x3 matrix view
constexpr real trace(Matrix3x3ConstSpan matrixView)
{
    return matrixView(0, 0) + matrixView(1, 1) + matrixView(2, 2);
}

/*! \brief Create a diagonal matrix of ElementType with N * M elements.
 *
 * \tparam ElementType type of matrix elements
 * \tparam N           number of rows
 * \tparam M           number of columns, defaults to number of rows if not set
 * \param  value       The value that fills the leading diagonal
 *
 * \returns a matrix with values \c value where row equals column index and null
 *          where row does not equal column index
 */
template<typename ElementType, int N, int M = N>
MultiDimArray<std::array<ElementType, N * M>, extents<N, M>> diagonalMatrix(const ElementType value)
{
    std::array<ElementType, N * M>                               matrixEntries{};
    MultiDimArray<std::array<ElementType, N * M>, extents<N, M>> matrix(matrixEntries);
    for (int i = 0; i < std::min(N, M); i++)
    {
        matrix(i, i) = value;
    }
    return matrix;
}

/*! \brief Create an identity matrix of ElementType with N * M elements.
 *
 * \tparam ElementType type of matrix elements
 * \tparam N number of rows
 * \tparam M number of columns, defaults to number of rows if not set
 *
 * \returns a matrix with values one where row equals column index and null
 *          where row does not equal column index
 */
template<typename ElementType, int N, int M = N>
MultiDimArray<std::array<ElementType, N * M>, extents<N, M>> identityMatrix()
{
    return diagonalMatrix<ElementType, N, M>(1);
}

//! Calculate the transpose of a 3x3 matrix, from its view
Matrix3x3 transpose(Matrix3x3ConstSpan matrixView);

//! Multiply matrix with vector.
void matrixVectorMultiply(Matrix3x3ConstSpan matrix, RVec* v);

//! Create new matrix type from legacy type.
static inline Matrix3x3 createMatrix3x3FromLegacyMatrix(const matrix legacyMatrix)
{
    GMX_RELEASE_ASSERT(legacyMatrix, "Need valid legacy matrix");
    Matrix3x3 newMatrix;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            newMatrix(i, j) = legacyMatrix[i][j];
        }
    }
    return newMatrix;
}

//! Fill legacy matrix from new matrix type.
static inline void fillLegacyMatrix(Matrix3x3ConstSpan newMatrix, matrix legacyMatrix)
{
    GMX_RELEASE_ASSERT(legacyMatrix, "Need valid legacy matrix");
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            legacyMatrix[i][j] = newMatrix(i, j);
        }
    }
}

//! Return the product of multiplying the vector \c v by the 3x3 matrix \c m
template<typename ElementType>
BasicVector<ElementType> multiplyVectorByMatrix(const BasicMatrix3x3<ElementType>& m, const rvec v)
{
    BasicVector<ElementType> result;
    for (int d = 0; d < DIM; ++d)
    {
        result[d] = m(d, 0) * v[0] + m(d, 1) * v[1] + m(d, 2) * v[2];
    }
    return result;
}

//! Return the sum of two 3x3 matrices \c a and \c b
template<typename ElementType>
BasicMatrix3x3<ElementType> operator+(const BasicMatrix3x3<ElementType>& a,
                                      const BasicMatrix3x3<ElementType>& b)
{
    BasicMatrix3x3<ElementType> result;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            result(i, j) = a(i, j) + b(i, j);
        }
    }
    return result;
}

//! Return the inner product of multiplication of two 3x3 matrices \c a and \c b
template<typename ElementType>
BasicMatrix3x3<ElementType> inner(const BasicMatrix3x3<ElementType>& a,
                                  const BasicMatrix3x3<ElementType>& b)
{
    BasicMatrix3x3<ElementType> result({ a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0) + a(0, 2) * b(2, 0),
                                         a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1) + a(0, 2) * b(2, 1),
                                         a(0, 0) * b(0, 2) + a(0, 1) * b(1, 2) + a(0, 2) * b(2, 2),
                                         a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0) + a(1, 2) * b(2, 0),
                                         a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1),
                                         a(1, 0) * b(0, 2) + a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2),
                                         a(2, 0) * b(0, 0) + a(2, 1) * b(1, 0) + a(2, 2) * b(2, 0),
                                         a(2, 0) * b(0, 1) + a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1),
                                         a(2, 0) * b(0, 2) + a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2) });
    return result;
}

//! Return the inner product of multiplication of two 3x3 matrices \c a and \c b
template<typename ElementType>
BasicMatrix3x3<ElementType> operator*(const BasicMatrix3x3<ElementType>& a,
                                      const BasicMatrix3x3<ElementType>& b)
{
    return inner(a, b);
}

//! Return the product of multiplying the 3x3 matrix \c m by the scalar \c s
template<typename ElementType>
BasicMatrix3x3<ElementType> operator*(const BasicMatrix3x3<ElementType>& m, const ElementType s)
{
    BasicMatrix3x3<ElementType> result;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            result(i, j) = m(i, j) * s;
        }
    }
    return result;
}

//! Return the product of multiplying the 3x3 matrix \c m by the scalar \c s
template<typename ElementType>
BasicMatrix3x3<ElementType> operator*(const ElementType s, const BasicMatrix3x3<ElementType>& m)
{
    return m * s;
}

//! Return a vector that is the diagonal of the 3x3 matrix \c m
template<typename ElementType>
BasicVector<ElementType> diagonal(const BasicMatrix3x3<ElementType>& m)
{
    return { m(XX, XX), m(YY, YY), m(ZZ, ZZ) };
}

} // namespace gmx

#endif
