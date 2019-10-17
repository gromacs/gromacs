/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declares functions for comparing views of vector-like data.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 * \inlibraryapi
 */
#ifndef GMX_MATH_TESTARRAYREFS_H
#define GMX_MATH_TESTARRAYREFS_H

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"

namespace gmx
{
namespace test
{

//! Initialization overload for non-BasicVector
template<typename T>
void fillInputContents(ArrayRef<T> inputRef, int scaleFactor)
{
    inputRef[0] = 1;
    inputRef[1] = 2;
    inputRef[2] = 3;
    for (auto& element : inputRef)
    {
        element *= scaleFactor;
    }
}

//! Initialization overload for BasicVector
template<typename T>
void fillInputContents(ArrayRef<BasicVector<T>> inputRef, int scaleFactor)
{
    inputRef[0] = { 1, 2, 3 };
    inputRef[1] = { 4, 5, 6 };
    inputRef[2] = { 7, 8, 9 };
    for (auto& element : inputRef)
    {
        element *= scaleFactor;
    }
}

//! Dispatcher function for filling.
template<typename PaddedVectorOfT>
void fillInput(PaddedVectorOfT* input, int scaleFactor)
{
    // Use a size for the vector in tests that is prime enough to
    // expose problems where they exist.
    int sizeOfVector = 3;
    input->resizeWithPadding(sizeOfVector);
    fillInputContents(makeArrayRef(*input), scaleFactor);
    EXPECT_LE(sizeOfVector, input->paddedSize());
}

//! Comparison overload for non-BasicVector
template<typename T>
void compareViews(ArrayRef<T> input, ArrayRef<T> output)
{
    ASSERT_EQ(input.size(), output.size());
    for (index i = 0; i != input.ssize(); ++i)
    {
        EXPECT_EQ(input[i], output[i]) << "for index " << i;
    }
}

//! Comparison overload for BasicVector<T>
template<typename T>
void compareViews(ArrayRef<BasicVector<T>> input, ArrayRef<BasicVector<T>> output)
{
    ASSERT_EQ(input.size(), output.size());
    for (index i = 0; i != input.ssize(); ++i)
    {
        EXPECT_EQ(input[i][XX], output[i][XX]) << "for index " << i;
        EXPECT_EQ(input[i][YY], output[i][YY]) << "for index " << i;
        EXPECT_EQ(input[i][ZZ], output[i][ZZ]) << "for index " << i;
    }
}

} // namespace test
} // namespace gmx

#endif
