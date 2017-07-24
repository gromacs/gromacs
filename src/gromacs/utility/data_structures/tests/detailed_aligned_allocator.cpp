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
   \brief
   Tests for the the irreg_array classes

   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/data_structures/detailed_aligned_allocator.h"

#include <vector>

#include "gromacs/utility/data_structures/flat_irreg_array_4d.h"
#include "gromacs/utility/data_structures/irreg_array_4d.h"

#include "data_struct_test_commons.h"

namespace
{

using std::endl;


class DetailedAlignedAllocatorTest : public gmx::data_struct_test::DataStructTest
{
};

/*********************************************************************/

// the following two are just adopted from the AlignedAllocatorTest
// DefaultAlignedAllocator should do exactly the same as AlignedAllocator

TEST_F(DetailedAlignedAllocatorTest, AllocatorAlign)
{
    gmx::DefaultAlignedAllocator<real>   a;
    real *                               p    = a.allocate(1000);

    // Mask for 128-byte alignment is 128-1 - these bits should be zero in p
    std::size_t                          mask = static_cast<std::size_t>(128-1);

    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a.deallocate(p, 1000);
}


TEST_F(DetailedAlignedAllocatorTest, Vector)
{
    // Mask for 128-byte alignment is 128-1 - these bits should be zero in pointers
    std::size_t mask = static_cast<std::size_t>(128-1);

    std::vector<double, gmx::DefaultAlignedAllocator<double> > v(10);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(v.data()) & mask);

    for (std::size_t i = 10000; i <= 100000; i += 10000)
    {
        v.resize(i);
        EXPECT_EQ(0, reinterpret_cast<std::size_t>(v.data()) & mask);
    }
}

TEST_F(DetailedAlignedAllocatorTest, ConstructDestructIrregArrays)
{

    typedef gmx::DefaultAlignedAllocator<double> allocator_type1;
    EXPECT_NO_THROW_GMX
    (
            {
                allocator_type1 test_alloc1;
                double *tmpPtr = test_alloc1.allocate(1);
                test_alloc1.deallocate(tmpPtr, 1);
            }
    )
    << "Construction destruction of \"DefaultAlignedAllocator\" failed with exception";

    typedef gmx::DetailedAlignedAllocator<double, 64u, 128u> allocator_type2;
    EXPECT_NO_THROW_GMX
    (
            {
                allocator_type2 test_alloc2;
                double *tmpPtr = test_alloc2.allocate(1);
                test_alloc2.deallocate(tmpPtr, 1);
            }
    )
    << "Construction destruction of \"DetailedAlignedAllocator\" failed with exception";

}

TEST_F(DetailedAlignedAllocatorTest, DataStructuresIrregArray1DAlloc)
{

    typedef gmx::DefaultAlignedAllocator<unsigned int> default_allocator_type;

    gmx::IrregArray1D<unsigned int, default_allocator_type> test_arr1(1, 6, 0u);

    EXPECT_NO_THROW_GMX(
            {
                test_arr1[1] = 1;
                test_arr1[2] = 2;
                test_arr1[3] = 3;
                test_arr1[4] = 4;
                test_arr1[5] = 5;
                test_arr1[6] = 6;
            }
            ) << "Caught exception while assigning via operator[]." << endl;

    EXPECT_TRUE(test_arr1[1] == 1 &&
                test_arr1[2] == 2 &&
                test_arr1[3] == 3 &&
                test_arr1[4] == 4 &&
                test_arr1[5] == 5 &&
                test_arr1[6] == 6)
    << "Test of assignment and retrieval via IrregArray1D::operator[] failed.";

    debug << "test_arr1: " << endl << test_arr1 << endl;
}


TEST_F(DetailedAlignedAllocatorTest, DataStructuresIrregArray4DUsage)
{

    typedef gmx::DefaultAlignedAllocator<unsigned int> default_allocator_type;

    gmx::IrregArray4D<unsigned int, default_allocator_type> test_arr1(1, 2, 1, 3, 1, 2, 1, 7, 0);

    EXPECT_NO_THROW_GMX(
            {
                test_arr1[1][1][1][1] = 1;
                test_arr1[1][2][1][1] = 2;
                test_arr1[1][3][1][1] = 3;
                test_arr1[2][1][1][1] = 4;
                test_arr1[2][2][1][1] = 5;
                test_arr1[2][3][1][1] = 6;
                test_arr1[1][1][2][1] = 7;
                test_arr1[1][2][2][1] = 8;
                test_arr1[1][3][2][1] = 9;
                test_arr1[2][1][2][1] = 10;
                test_arr1[2][2][2][1] = 11;
                test_arr1[2][3][2][1] = 12;
            }
            ) << "Caught exception while constructing an IrregArray1D" << endl;

    debug << "test_arr1: " << endl << test_arr1 << endl;


    gmx::IrregArray1D<size_t> sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::IrregArray2D<size_t>       sizes2(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 3u);
    gmx::IrregArray3D<size_t>       sizes3(sizes1.getBegin1(), sizes1.getLast1(), sizes1, sizes2, 1u);

    gmx::IrregArray4D<unsigned int, default_allocator_type> test_arr2((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);

    EXPECT_NO_THROW_GMX(
            {
                unsigned int m = 1;
                for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
                {
                    for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                    {
                        for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                            {
                                test_arr2[i][j][k][l] = m;
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

}


} // namespace
