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
 * \author Anatolii Titov <Wapuk-cobaka@yandex.ru>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_math
 */

#ifndef GMX_MATH_MATRIX_H_
#define GMX_MATH_MATRIX_H_

#include <array>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

/*! \brief Three-by-three matrix of ElementType.
 * \tparam ElementType type of element to be stored in matrix
 */
template<class ElementType>
class BasicMatrix3x3
{
private:
    //! array size
    static const size_t matrixSize_ = DIM * DIM;

public:
    //! Default constructor
    constexpr BasicMatrix3x3() noexcept : storage_{ 0, 0, 0, 0, 0, 0, 0, 0, 0 } {}

    //! Constructor for a matrix with all entries equal
    constexpr BasicMatrix3x3(const ElementType& initValue) noexcept :
        storage_{ initValue, initValue, initValue, initValue, initValue,
                  initValue, initValue, initValue, initValue }
    {
    }

    //! Constructor from a 9 element std::array
    constexpr BasicMatrix3x3(const std::array<ElementType, matrixSize_>& array) noexcept :
        storage_{ array }
    {
    }

    //! Old style access operator()
    //[[deprecated("Use [][] instead, operator () will be deleted at some point")]]
    ElementType& operator()(const size_t row, const size_t col)
    {
        return storage_[(row * DIM) + col];
    }

    //! Old style access operator()
    //[[deprecated("Use [][] instead, operator () will be deleted at some point")]]
    const ElementType& operator()(const size_t row, const size_t col) const
    {
        return storage_[(row * DIM) + col];
    }

    //! Proxy class for row-wise access via a[i][j] syntax
    //! * \tparam MatrixType type of matrix (const or non-const)
    template<class MatrixType, class LocalElementType>
    class RowProxy
    {
    public:
        //! Construct proxy for given row of matrix
        RowProxy(MatrixType& mat, size_t row) : mat_(mat), row_(row) {}

        //! Non-const element access: returns reference to element at (row, col)
        LocalElementType& operator[](const size_t col) { return mat_(row_, col); }

        //! Const element access: returns reference to element at (row, col)
        const LocalElementType& operator[](const size_t col) const { return mat_(row_, col); }

        //! Implicit conversion to BasicVector for convenient row access
        operator BasicVector<ElementType>()
        {
            return BasicVector<ElementType>{ mat_(row_, XX), mat_(row_, YY), mat_(row_, ZZ) };
        }

        //! Implicit conversion to BasicVector for convenient row access
        operator BasicVector<ElementType>() const
        {
            return BasicVector<ElementType>{ mat_(row_, XX), mat_(row_, YY), mat_(row_, ZZ) };
        }

    private:
        //! stored reference to matrix
        MatrixType& mat_;
        //! store row number
        size_t row_;
    };

    //! Row access: returns proxy for the given row, enabling a[i][j] syntax
    RowProxy<BasicMatrix3x3, ElementType> operator[](const size_t row)
    {
        return RowProxy<BasicMatrix3x3, ElementType>(*this, row);
    }

    //! const Row access: returns a const proxy for the given row, enabling a[i][j] syntax
    RowProxy<const BasicMatrix3x3, const ElementType> operator[](const size_t row) const
    {
        return RowProxy<const BasicMatrix3x3, const ElementType>(*this, row);
    }

    //! Return result of adding \c other matrix to this one
    BasicMatrix3x3 operator+(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result;
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            result.storage_[i] = storage_[i] + other.storage_[i];
        }
        return result;
    }

    //! Return result of subtracting \c other matrix from this one
    BasicMatrix3x3 operator-(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result;
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            result.storage_[i] = storage_[i] - other.storage_[i];
        }
        return result;
    }

    //! Return negation of all values of this matrix
    BasicMatrix3x3 operator-() const
    {
        BasicMatrix3x3<ElementType> result;
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            result.storage_[i] = -storage_[i];
        }
        return result;
    }

    //! Return == result, true if all element pairs are equal
    bool operator==(const BasicMatrix3x3& other) const { return storage_ == other.storage_; }

    //! Return result of multiplication (inner product) of this matrix by \c other matrix
    BasicMatrix3x3 operator*(const BasicMatrix3x3& other) const
    {
        BasicMatrix3x3<ElementType> result(
                { (*this)(XX, XX) * other(XX, XX) + (*this)(XX, YY) * other(YY, XX)
                          + (*this)(XX, ZZ) * other(ZZ, XX),
                  (*this)(XX, XX) * other(XX, YY) + (*this)(XX, YY) * other(YY, YY)
                          + (*this)(XX, ZZ) * other(ZZ, YY),
                  (*this)(XX, XX) * other(XX, ZZ) + (*this)(XX, YY) * other(YY, ZZ)
                          + (*this)(XX, ZZ) * other(ZZ, ZZ),
                  (*this)(YY, XX) * other(XX, XX) + (*this)(YY, YY) * other(YY, XX)
                          + (*this)(YY, ZZ) * other(ZZ, XX),
                  (*this)(YY, XX) * other(XX, YY) + (*this)(YY, YY) * other(YY, YY)
                          + (*this)(YY, ZZ) * other(ZZ, YY),
                  (*this)(YY, XX) * other(XX, ZZ) + (*this)(YY, YY) * other(YY, ZZ)
                          + (*this)(YY, ZZ) * other(ZZ, ZZ),
                  (*this)(ZZ, XX) * other(XX, XX) + (*this)(ZZ, YY) * other(YY, XX)
                          + (*this)(ZZ, ZZ) * other(ZZ, XX),
                  (*this)(ZZ, XX) * other(XX, YY) + (*this)(ZZ, YY) * other(YY, YY)
                          + (*this)(ZZ, ZZ) * other(ZZ, YY),
                  (*this)(ZZ, XX) * other(XX, ZZ) + (*this)(ZZ, YY) * other(YY, ZZ)
                          + (*this)(ZZ, ZZ) * other(ZZ, ZZ) });
        return result;
    }

    //! Return result of multiplication this matrix by constant \c scalar
    BasicMatrix3x3 operator*(const ElementType& scalar) const
    {
        BasicMatrix3x3<ElementType> result;
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            result.storage_[i] = storage_[i] * scalar;
        }
        return result;
    }

    //! Return result of division this matrix by constant scalar
    BasicMatrix3x3 operator/(const ElementType& scalar) const
    {
        GMX_RELEASE_ASSERT(scalar != 0, "Division by zero.");
        BasicMatrix3x3<ElementType> result;
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            result.storage_[i] = storage_[i] / scalar;
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
        *this = *this / scalar;
        return *this;
    }

    //! Returns result of multiplication of this matrix by a \c vector
    BasicVector<ElementType> operator*(const BasicVector<ElementType>& vector) const
    {
        BasicVector<ElementType> result = { 0, 0, 0 };
        for (size_t i = XX; i < DIM; ++i)
        {
            for (size_t j = XX; j < DIM; ++j)
            {
                result[i] += storage_[(i * DIM) + j] * vector[j];
            }
        }
        return result;
    }

    //! Conversion to gmx::ArrayRef
    ArrayRef<ElementType> toArrayRef()
    {
        return { storage_.data(), storage_.data() + storage_.size() };
    }

    //! Conversion to const ArrayRef
    ArrayRef<const ElementType> toConstArrayRef() const
    {
        return { storage_.data(), storage_.data() + storage_.size() };
    }

    //! Returns a pointer to the start of the storage
    ElementType* data() { return storage_; }

    //! Returns a const pointer to start of the storage
    const ElementType* data() const { return storage_; }

    //! Swaps two matrixies
    void swap(BasicMatrix3x3& other) noexcept { std::swap(storage_, other.storage_); }

    //! Sets all elements to 0
    void clear() { std::fill(storage_.begin(), storage_.end(), 0); }

private:
    //! local "multi dim array" 3x3
    std::array<ElementType, matrixSize_> storage_;
};

//! Stand alone swap function based on https://en.cppreference.com/w/cpp/named_req/Swappable
template<class ElementType>
void swap(BasicMatrix3x3<ElementType>& mat1, BasicMatrix3x3<ElementType>& mat2) noexcept
{
    mat1.swap(mat2);
}

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
    return { { matrix(XX, XX),
               matrix(YY, XX),
               matrix(ZZ, XX),
               matrix(XX, YY),
               matrix(YY, YY),
               matrix(ZZ, YY),
               matrix(XX, ZZ),
               matrix(YY, ZZ),
               matrix(ZZ, ZZ) } };
}

//! Returns the determinant of \c matrix
template<typename ElementType>
constexpr ElementType determinant(BasicMatrix3x3<ElementType> matrix)
{
    return { matrix(XX, XX) * (matrix(YY, YY) * matrix(ZZ, ZZ) - matrix(ZZ, YY) * matrix(YY, ZZ))
             - matrix(YY, XX) * (matrix(XX, YY) * matrix(ZZ, ZZ) - matrix(ZZ, YY) * matrix(XX, ZZ))
             + matrix(ZZ, XX) * (matrix(XX, YY) * matrix(YY, ZZ) - matrix(YY, YY) * matrix(XX, ZZ)) };
}

//! Returns the trace of \c matrix
template<typename ElementType>
constexpr ElementType trace(BasicMatrix3x3<ElementType> matrix)
{
    return matrix(XX, XX) + matrix(YY, YY) + matrix(ZZ, ZZ);
}

//! Return the inner product of multiplication of two 3x3 matrices \c a and \c b
template<typename ElementType>
BasicMatrix3x3<ElementType> inner(const BasicMatrix3x3<ElementType>& matrixA,
                                  const BasicMatrix3x3<ElementType>& matrixB)
{
    return matrixA * matrixB;
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
    BasicMatrix3x3<ElementType> matrix;
    for (size_t i = XX; i < DIM; ++i)
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
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
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
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            legacyMatrix[i][j] = newMatrix(i, j);
        }
    }
}

} // namespace gmx

#endif
