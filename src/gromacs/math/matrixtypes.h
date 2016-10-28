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
    None,
    UpperTriangle,
    LowerTriangle
};

/*! \brief
 * C++ class for 3x3 matrices.
 *
 * \tparam ValueType       Value type
 * \tparam ConstraintType  Matrix constraint type
 *
 * This class provides a C++ version of matrix,
 * that can be put into STL containers etc. It is more or less
 * a drop-in replacement for `matrix` and friends:
 * it can be used in most contexts that accept the equivalent C type.
 * Provided constraint is enforced internally.
 * If it is violated, gmx::InternalError is thrown.
 *
 * \inpublicapi
 */
template <typename ValueType, MatrixConstraint ConstraintType>
class BasicMatrix3x3
{
    public:
        //! 1D raw array as a return type for operator[]
        typedef ValueType RawArray[DIM];
        //! Underlying raw C 2D array (same as matrix in vectypes.h).
        typedef RawArray RawMatrix[DIM];

        //! Constructs default (uninitialized) matrix.
        BasicMatrix3x3() {}
        /*! \brief Constructs a matrix from given values.
         * \throws InternalError if the constraint is violated.
         */
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
            Validate();
        }
        /*! \brief Constructs a matrix from C-style 2D array.
         * \throws InternalError if the constraint is violated.
         */
        BasicMatrix3x3(const RawMatrix matrix)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    matrix_[i][j] = matrix[i][j];
                }
            }
            Validate();
        }

        //! Indexing operator to make the class work as the raw 2D array.
        ValueType *operator[](int i) { return matrix_[i]; }
        //! Indexing operator to make the class work as the raw 2D array.
        const ValueType *operator[](int i) const { return matrix_[i]; }

        //! Makes BasicMatrix3x3 usable in contexts where a raw C 2D array is expected.
        operator RawArray *() { return matrix_; }
        //! Makes BasicMatrix3x3 usable in contexts where a raw C 2D array is expected.
        operator const RawArray *() const { return matrix_; }

    private:
        //! Storage
        RawMatrix matrix_;
        /*! \brief Upholds the contraint. Should be called from any assignment.
         * \throws InternalError if the constraint is violated.
         */
        void Validate()
        {
            bool correct;
            switch (ConstraintType)
            {
                case MatrixConstraint::None:
                    correct = true;
                    break;

                case MatrixConstraint::UpperTriangle:
                    correct = (matrix_[YY][XX] == 0) && (matrix_[ZZ][XX] == 0) && (matrix_[ZZ][YY] == 0);
                    break;

                case MatrixConstraint::LowerTriangle:
                    correct = (matrix_[XX][YY] == 0) && (matrix_[XX][ZZ] == 0) && (matrix_[YY][ZZ] == 0);
                    break;

                default:
                    GMX_THROW(InternalError("This matrix constraint type is not implemented"));
            }
            if (!correct)
            {
                GMX_THROW(InternalError("The matrix constraint was violated"));
            }
        }
};

//! Shorthand for C++ `matrix`-equivalent type.
typedef BasicMatrix3x3<real, MatrixConstraint::None> Matrix3x3;
//! Shorthand for lower triangular 3x3 matrix type.
typedef BasicMatrix3x3<real, MatrixConstraint::LowerTriangle> Matrix3x3Lower;
//! Shorthand for upper triangular 3x3 matrix type.
typedef BasicMatrix3x3<real, MatrixConstraint::UpperTriangle> Matrix3x3Upper;

} // namespace gmx

#endif

#endif
