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
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_math
 */

#ifndef GMX_MATH_MATRIX_H_
#define GMX_MATH_MATRIX_H_

#include <array>

#include "gromacs/math/multidimarray.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

/*! \brief Three-by-three matrix of ElementType.
 * \tparam ElementType type of element to be stored in matrix
 */
template<class ElementType>
class BasicMatrix3x3 : public MultiDimArray<std::array<ElementType, DIM * DIM>, extents<DIM, DIM>>
{
public:
    //! Default constructor
    BasicMatrix3x3() : MultiDimArray<std::array<ElementType, DIM * DIM>, extents<DIM, DIM>>({}) {}

    //! Propagate constructor from Base
    template<typename... Args>
    BasicMatrix3x3(Args&&... args) :
        MultiDimArray<std::array<ElementType, DIM * DIM>, extents<DIM, DIM>>(std::forward<Args>(args)...)
    {
    }

    //! Proxy class for row-wise access via a[i][j] and a[i] syntax
    class RowProxy
    {
        BasicMatrix3x3& mat_;
        size_t          row_;

    public:
        //! Construct proxy for given row of matrix
        RowProxy(BasicMatrix3x3& mat, size_t row) : mat_(mat), row_(row) {}

        //! Non-const element access: returns reference to element at (row, col)
        ElementType& operator[](size_t col) { return mat_(row_, col); }

        //! Const element access: returns reference to element at (row, col)
        const ElementType& operator[](size_t col) const { return mat_(row_, col); }

        //! Implicit copy as const BasicVector for convenient row-wise read access
        operator const BasicVector<ElementType>() const
        {
            return BasicVector<ElementType>{ (*this)[0], (*this)[1], (*this)[2] };
        }
    };

    //! Proxy class for row-wise const access via a[i][j] and a[i] syntax
    class ConstRowProxy
    {
        const BasicMatrix3x3& mat_;
        size_t                row_;

    public:
        //! Construct proxy for given row of matrix
        ConstRowProxy(const BasicMatrix3x3& mat, const size_t row) : mat_(mat), row_(row) {}

        //! Const element access: returns reference to element at (row, col)
        const ElementType& operator[](size_t col) const { return mat_(row_, col); }

        //! Implicit copy as const BasicVector for convenient row-wise read access
        operator const BasicVector<ElementType>() const
        {
            return BasicVector<ElementType>{ (*this)[0], (*this)[1], (*this)[2] };
        }
    };

    //! Row access: returns proxy for the given row, enabling a[i][j] syntax
    RowProxy operator[](size_t row) { return RowProxy(*this, row); }

    //! Const row access: returns const proxy for the given row, enabling a[i][j] syntax (const matrix)
    ConstRowProxy operator[](size_t row) const { return ConstRowProxy(*this, row); }

    //! Constructor for list initializer (either 1 or 9 elements in row major ordering)
    BasicMatrix3x3(std::initializer_list<ElementType> initList)
    {
        if ((initList.size() != DIM * DIM) && (initList.size() != 1))
        {
            throw std::invalid_argument(
                    "Initializer list must contain exactly either 1 or 9 elements");
        }
        if (initList.size() != 1)
        {
            auto it = initList.begin();
            for (size_t i = 0; i < DIM; ++i)
            {
                for (size_t j = 0; j < DIM; ++j)
                {
                    (*this)(i, j) = *it++;
                }
            }
        }
        else
        {
            auto initValue = *initList.begin();
            for (size_t i = 0; i < DIM; ++i)
            {
                for (size_t j = 0; j < DIM; ++j)
                {
                    (*this)(i, j) = initValue;
                }
            }
        }
    }


    //! Return result of adding \c other matrix to this one
    BasicMatrix3x3 operator+(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    //! Return result of subtracting \c other matrix from this one
    BasicMatrix3x3 operator-(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    //! Return negation of all values of this matrix
    BasicMatrix3x3 operator-() const
    {
        BasicMatrix3x3<ElementType> result;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                result(i, j) = -(*this)(i, j);
            }
        }
        return result;
    }

    //! Return result of multiplication (inner product) of this matrix by \c other matrix
    BasicMatrix3x3 operator*(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result(
                { (*this)(0, 0) * other(0, 0) + (*this)(0, 1) * other(1, 0) + (*this)(0, 2) * other(2, 0),
                  (*this)(0, 0) * other(0, 1) + (*this)(0, 1) * other(1, 1) + (*this)(0, 2) * other(2, 1),
                  (*this)(0, 0) * other(0, 2) + (*this)(0, 1) * other(1, 2) + (*this)(0, 2) * other(2, 2),
                  (*this)(1, 0) * other(0, 0) + (*this)(1, 1) * other(1, 0) + (*this)(1, 2) * other(2, 0),
                  (*this)(1, 0) * other(0, 1) + (*this)(1, 1) * other(1, 1) + (*this)(1, 2) * other(2, 1),
                  (*this)(1, 0) * other(0, 2) + (*this)(1, 1) * other(1, 2) + (*this)(1, 2) * other(2, 2),
                  (*this)(2, 0) * other(0, 0) + (*this)(2, 1) * other(1, 0) + (*this)(2, 2) * other(2, 0),
                  (*this)(2, 0) * other(0, 1) + (*this)(2, 1) * other(1, 1) + (*this)(2, 2) * other(2, 1),
                  (*this)(2, 0) * other(0, 2) + (*this)(2, 1) * other(1, 2) + (*this)(2, 2) * other(2, 2) });
        return result;
    }

    //! Return result of multiplication this matrix by constant \c scalar
    BasicMatrix3x3 operator*(const ElementType& scalar) const
    {
        BasicMatrix3x3<ElementType> result;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
        return result;
    }

    //! Return result of division this matrix by constant scalar
    BasicMatrix3x3 operator/(const ElementType& scalar) const
    {
        GMX_RELEASE_ASSERT(scalar != 0, "Division by zero.");
        BasicMatrix3x3<ElementType> result;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                result(i, j) = (*this)(i, j) / scalar;
            }
        }
        return result;
    }

    //! Return result of addition with assignment of \c other matrix to this matrix
    BasicMatrix3x3& operator+=(const BasicMatrix3x3& other)
    {
        *this = *this + other;
        return *this;
    }

    //! Return result of subtraction with assignment of \c other matrix from this matrix
    BasicMatrix3x3& operator-=(const BasicMatrix3x3& other)
    {
        *this = *this - other;
        return *this;
    }

    //! Return result of multiplication with assignment of this matrix by \c scalar
    BasicMatrix3x3& operator*=(const ElementType& scalar)
    {
        *this = *this * scalar;
        return *this;
    }

    //! Return result of division with assignment of this matrix by \c scalar
    BasicMatrix3x3& operator/=(const ElementType& scalar)
    {
        GMX_RELEASE_ASSERT(scalar != 0, "Division by zero.");
        *this = *this / scalar;
        return *this;
    }

    //! Returns result of multiplication of this matrix by a \c vector
    BasicVector<ElementType> operator*(const BasicVector<ElementType>& vector) const
    {
        BasicVector<ElementType> result = { 0, 0, 0 };
        for (auto i = 0; i < DIM; ++i)
        {
            for (auto j = 0; j < DIM; ++j)
            {
                result[i] += (*this)(i, j) * vector[j];
            }
        }
        return result;
    }
};

//! Return the product of multiplying the 3x3 matrix \c other by the \c scalar
template<typename ElementType>
BasicMatrix3x3<ElementType> operator*(const ElementType scalar, const BasicMatrix3x3<ElementType>& other)
{
    return other * scalar;
}

//! Return a vector that is the diagonal of the 3x3 \c matrix
template<typename ElementType>
BasicVector<ElementType> diagonal(const BasicMatrix3x3<ElementType>& matrix)
{
    return { matrix(XX, XX), matrix(YY, YY), matrix(ZZ, ZZ) };
}

//! Returns the transposition of \c matrix
template<typename ElementType>
BasicMatrix3x3<ElementType> transpose(BasicMatrix3x3<ElementType> matrix)
{

    return { matrix(0, 0), matrix(1, 0), matrix(2, 0), matrix(0, 1), matrix(1, 1),
             matrix(2, 1), matrix(0, 2), matrix(1, 2), matrix(2, 2) };
}

//! Returns the determinant of \c matrix
template<typename ElementType>
constexpr ElementType determinant(BasicMatrix3x3<ElementType> matrix)
{
    return { matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(2, 1) * matrix(1, 2))
             - matrix(1, 0) * (matrix(0, 1) * matrix(2, 2) - matrix(2, 1) * matrix(0, 2))
             + matrix(2, 0) * (matrix(0, 1) * matrix(1, 2) - matrix(1, 1) * matrix(0, 2)) };
}

//! Returns the trace of \c matrix
template<typename ElementType>
constexpr ElementType trace(BasicMatrix3x3<ElementType> matrix)
{
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}

//! Return the inner product of multiplication of two 3x3 matrices \c a and \c b
template<typename ElementType>
BasicMatrix3x3<ElementType> inner(const BasicMatrix3x3<ElementType>& a,
                                  const BasicMatrix3x3<ElementType>& b)
{
    return a * b;
}

/*! \brief Three-by-three real number matrix.
 * \note will replace the C-style real[3][3] "matrix"
 */
using Matrix3x3 = BasicMatrix3x3<real>;

/*! \brief Create a diagonal matrix of ElementType.
 *
 * \tparam ElementType type of matrix elements
 * \param  value       The value that fills the leading diagonal
 *
 * \returns a matrix with values \c value where row equals column index and null
 *          where row does not equal column index
 */
template<typename ElementType>
BasicMatrix3x3<ElementType> diagonalMatrix(const ElementType value)
{
    BasicMatrix3x3<ElementType> matrix = { 0 };
    for (auto i = 0; i < DIM; i++)
    {
        matrix(i, i) = value;
    }
    return matrix;
}

/*! \brief Create an identity matrix of ElementType.
 * \tparam ElementType type of matrix elements
 * \returns a matrix with values one where row equals column index and null
 *          where row does not equal column index
 */
template<typename ElementType>
BasicMatrix3x3<ElementType> identityMatrix()
{
    return diagonalMatrix(static_cast<ElementType>(1));
}

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

//! Fill existing legacy matrix from new matrix type.
static inline void fillLegacyMatrix(Matrix3x3 newMatrix, matrix legacyMatrix)
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

} // namespace gmx

#endif
