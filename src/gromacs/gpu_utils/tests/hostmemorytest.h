/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Tests for GPU host vector.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#ifndef GMX_GPU_UTILS_TESTS_HOSTMEMORYTEST_H
#define GMX_GPU_UTILS_TESTS_HOSTMEMORYTEST_H

#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "devicetransfers.h"
#include "gputest.h"

namespace gmx
{

namespace
{

/*! \internal \brief Typed test fixture for infrastructure for
 * host-side memory used for GPU transfers. */
template <typename T>
class HostMemoryTest : public test::GpuTest
{
    public:
        //! Convenience type
        using ValueType = T;
        //! Convenience type
        using ViewType = ArrayRef<ValueType>;
        //! Convenience type
        using ConstViewType = ArrayRef<const ValueType>;
        //! Prepare contents of a VectorType.
        template <typename VectorType>
        void fillInput(VectorType *input) const;
        //! Compares input and output vectors.
        void compareVectors(ConstViewType input,
                            ConstViewType output) const;
        //! Do some transfers and test the results.
        void runTest(ConstViewType input, ViewType output) const;
};

// Already documented
template <typename T> template <typename VectorType>
void HostMemoryTest<T>::fillInput(VectorType *input) const
{
    input->resize(3);
    (*input)[0] = 1;
    (*input)[1] = 2;
    (*input)[2] = 3;
}

//! Initialization specialization for RVec
template <> template <typename VectorType>
void HostMemoryTest<RVec>::fillInput(VectorType *input) const
{
    input->resize(1);
    (*input)[0] = {1, 2, 3};
}

// Already documented
template <typename T>
void HostMemoryTest<T>::compareVectors(ConstViewType input,
                                       ConstViewType output) const
{
    for (size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_EQ(input[i], output[i]) << "for index " << i;
    }
}

//! Comparison specialization for RVec
template <>
void HostMemoryTest<RVec>::compareVectors(ConstViewType input,
                                          ConstViewType output) const
{
    for (size_t i = 0; i != input.size(); ++i)
    {
        EXPECT_EQ(input[i][XX], output[i][XX]) << "for index " << i;
        EXPECT_EQ(input[i][YY], output[i][YY]) << "for index " << i;
        EXPECT_EQ(input[i][ZZ], output[i][ZZ]) << "for index " << i;
    }
}

/*! \brief Convenience function to transform a view into one with base
 * type of (non-const) char.
 *
 * This transformation is useful for using containers with C APIs
 * where the function signature is not declared const even where the
 * semantics of the usage actually are const.
 *
 * \param[in]    data   The data pointer.
 * \param[in]    size   The size of the data pointer (in T).
 * \tparam       T      The base type of the container
 * */
template <typename T>
ArrayRef<char> charArrayRefFromArray(T *data, size_t size)
{
    // Make a type like T, but without its possible const qualifier.
    using NonConstT = typename std::remove_const<T>::type;
    return arrayRefFromArray<char>(reinterpret_cast<char *>(const_cast<NonConstT *>(data)), size * sizeof(T));
}

template <typename T>
void HostMemoryTest<T>::runTest(ConstViewType input, ViewType output) const
{
    // We can't do a test that does a transfer unless we have a
    // compatible device.
    if (!this->haveValidGpus())
    {
        return;
    }

    // Convert the views of input and output to flat non-const chars,
    // so that there's no templating when we call doDeviceTransfers.
    auto inputRef  = charArrayRefFromArray(input.data(), input.size());
    auto outputRef = charArrayRefFromArray(output.data(), output.size());

    doDeviceTransfers(*this->gpuInfo_, inputRef, outputRef);
    this->compareVectors(input, output);
}

} // namespace
} // namespace

#endif
