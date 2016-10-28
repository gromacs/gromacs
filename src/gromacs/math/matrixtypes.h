/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#ifndef GMX_MATH_MATRIXTYPES_H
#define GMX_MATH_MATRIXTYPES_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus

namespace gmx
{

//! Types of matrix constraints
enum class MatrixConstraint
{
    UpperTriangle,
    LowerTriangle
};

/*! \brief
 * C++ class for 3x3 matrices.
 *
 * \tparam ValueType       Value type
 *
 * This class provides a C++ version of matrix,
 * that can be put into STL containers etc. It is more or less
 * a drop-in replacement for `matrix` and friends:
 * it can be used in most contexts that accept the equivalent C type.
 *
 * Modifying individual values directly is not allowed; please use setValue() instead.
 * The reason for that is that the descendant class ConstrainedMatrix3x3
 * would be trickier to implement if the direct modification
 * through the row array pointer would be allowed.
 *
 * \inpublicapi
 */
template <typename ValueType>
class BasicMatrix3x3
{
    public:
        //! 1D raw array as a return type for operator[]
        typedef ValueType RawArray[DIM];
        //! Underlying raw C 2D array (same as matrix in vectypes.h).
        typedef RawArray RawMatrix[DIM];

        //! Constructs default (uninitialized) matrix.
        BasicMatrix3x3() = default;
        //! Constructs a matrix from given values.
        BasicMatrix3x3(ValueType xx, ValueType xy, ValueType xz,
                       ValueType yx, ValueType yy, ValueType yz,
                       ValueType zx, ValueType zy, ValueType zz)
        {
            matrix_[XX][XX] = xx;
            matrix_[XX][YY] = xy;
            matrix_[XX][ZZ] = xz;
            matrix_[YY][XX] = yx;
            matrix_[YY][YY] = yy;
            matrix_[YY][ZZ] = yz;
            matrix_[ZZ][XX] = zx;
            matrix_[ZZ][YY] = zy;
            matrix_[ZZ][ZZ] = zz;
        }
        //! Constructs a matrix from C-style 2D array.
        BasicMatrix3x3(const RawMatrix matrix)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    matrix_[i][j] = matrix[i][j];
                }
            }
        }

        /*! \brief Indexing operator to make the class work as the raw 2D array.
         * This is only good for reading the const values.
         * For altering individual values there is an explicit setValue() function.
         */
        const ValueType *operator[](int i) const { return matrix_[i]; }
        /*! \brief
         * Changes the individual value.
         * Explicit replacement of the indexing operator for writing.
         */
        void setValue(size_t row, size_t col, const ValueType &value)
        {
            GMX_RELEASE_ASSERT(row < DIM, "Wrong row index");
            GMX_RELEASE_ASSERT(col < DIM, "Wrong column index");
            matrix_[row][col] = value;
        }

        //! Makes BasicMatrix3x3 usable in contexts where a raw C 2D array is expected.
        operator RawArray *() { return matrix_; }
        //! Makes BasicMatrix3x3 usable in contexts where a raw C 2D array is expected.
        operator const RawArray *() const { return matrix_; }

    protected:
        //! Storage
        RawMatrix matrix_;
};

/*! \brief
 * C++ class for constrained 3x3 matrices - e.g. triangular.
 *
 * \tparam ValueType       Value type
 * \tparam ConstraintType  Matrix constraint type
 *
 * This class inherits from BasicMatrix3x3.
 * Same idea applies: it can be used in most contexts that accept the matrix C type.
 * Unlike BasicMatrix3x3, it also upholds the provided constraint type,
 * calling the internal Validate() function on any assignment.
 * ConstrainedMatrix3x3 should also cast transparently into BasicMatrix3x3 of the same ValueType.
 *
 * Modifying individual values directly with an indexing operator
 * is not allowed; please use setValue() instead.
 *
 * \inpublicapi
 */
template <typename ValueType, MatrixConstraint ConstraintType>
class ConstrainedMatrix3x3 : public BasicMatrix3x3<ValueType>
{
    public:
        //! 1D raw array as a return type for operator[]
        typedef ValueType RawArray[DIM];
        //! Underlying raw C 2D array (same as matrix in vectypes.h).
        typedef RawArray RawMatrix[DIM];

        //! Constructs default (uninitialized) matrix.
        ConstrainedMatrix3x3() = default;

        /*! \brief Constructs a matrix from given values.
         * Asserts if the constraint is violated.
         */
        ConstrainedMatrix3x3(ValueType xx, ValueType xy, ValueType xz,
                             ValueType yx, ValueType yy, ValueType yz,
                             ValueType zx, ValueType zy, ValueType zz)
            : BasicMatrix3x3<ValueType>(xx, xy, xz, yx, yy, yz, zx, zy, zz)
        {
            Validate();
        }
        /*! \brief Constructs a matrix from C-style 2D array.
         * Asserts if the constraint is violated.
         */
        ConstrainedMatrix3x3(const RawMatrix matrix) : BasicMatrix3x3<ValueType>(matrix)
        {
            Validate();
        }

        /*! \brief
         * Changes the individual value.
         * Asserts if the constraint is violated.
         * Explicit replacement of the indexing operator for writing.
         */
        void setValue(size_t row, size_t col, const ValueType &value)
        {
            BasicMatrix3x3<ValueType>::setValue(row, col, value);
            Validate();
        }

        //! Makes constrained matrix able to implicitly convert to unconstrained.
        operator BasicMatrix3x3<ValueType> &() const
        {
            return (BasicMatrix3x3<ValueType> &) this->matrix_;
        }

    private:
        /*! \brief Upholds the contraint. Should be called from any assignment.
         * Asserts if the constraint is violated.
         */
        void Validate()
        {
            bool correct = true;
#ifdef __INTEL_COMPILER
#pragma warning( push )
#pragma warning( disable : 280 ) // "selector expression is constant"
#endif
            switch (ConstraintType)
            {
                case MatrixConstraint::UpperTriangle:
                    correct = (this->matrix_[YY][XX] == 0) && (this->matrix_[ZZ][XX] == 0) && (this->matrix_[ZZ][YY] == 0);
                    break;

                case MatrixConstraint::LowerTriangle:
                    correct = (this->matrix_[XX][YY] == 0) && (this->matrix_[XX][ZZ] == 0) && (this->matrix_[YY][ZZ] == 0);
                    break;

                default:
                    GMX_RELEASE_ASSERT(false, "This matrix constraint type is not implemented");
            }
            GMX_RELEASE_ASSERT(correct, "The matrix constraint was violated");
#ifdef __INTEL_COMPILER
#pragma warning( pop )
#endif
        }
};

//! Shorthand for C++ `matrix`-equivalent type.
typedef BasicMatrix3x3<real> Matrix3x3;
//! Shorthand for lower triangular 3x3 matrix type.
typedef ConstrainedMatrix3x3<real, MatrixConstraint::LowerTriangle> Matrix3x3Lower;
//! Shorthand for upper triangular 3x3 matrix type.
typedef ConstrainedMatrix3x3<real, MatrixConstraint::UpperTriangle> Matrix3x3Upper;

} // namespace gmx

#endif

#endif
