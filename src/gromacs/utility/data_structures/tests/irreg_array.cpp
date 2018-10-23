/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "data_struct_test_commons.h"

namespace gmx
{

namespace data_struct_test
{

namespace
{

using std::endl;

/*********************************************************************/

TEST_F(IrregArrayTest, ConstructDestructIrregArrays)
{

    EXPECT_NO_THROW_GMX({ gmx::IrregArray1D<double> test_1d_irreg_array; })
    << "Construction destruction of \"IrregArray1D\" failed with exception";

    EXPECT_NO_THROW_GMX({ gmx::IrregArray2D<double> test_2d_irreg_array; })
    << "Construction destruction of \"IrregArray2D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::IrregArray3D<double> test_3d_irreg_array; })
    << "Construction destruction of \"IrregArray3D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::IrregArray4D<double> test_4d_irreg_array; })
    << "Construction destruction of \"IrregArray4D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::FlatIrregArray2D<double> test_flat_2d_irreg_array; })
    << "Construction destruction of \"IrregArray2D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::FlatIrregArray3D<double> test_flat_3d_irreg_array; })
    << "Construction destruction of \"IrregArray3D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::FlatIrregArray4D<double> test_flat_4d_irreg_array; })
    << "Construction destruction of \"IrregArray4D\" failed";

}

TEST_F(IrregArrayTest, ConstructIrregArray1DFromInitializerList)
{
    typedef double                          value_type;
    typedef gmx::IrregArray1D<value_type>   array_type;
    typedef array_type::index_type          index_type;

    bool caughtException = false;
    try
    {
        const std::initializer_list<value_type> iList = {0, 1, 2, 3, 4, 5};
        array_type testArr = iList;

        debug << "Initializer list " << iList   << std::endl;
        debug << "testArr          " << testArr << std::endl;

        value_type const * iListPtr = iList.begin();
        for (index_type i = testArr.getFirst1(); i <= testArr.getLast1(); ++i)
        {
            EXPECT_TRUE(iListPtr != iList.end()) << "The test array contains more elements than the initializer list.";
            EXPECT_EQ(testArr[i], *iListPtr) << "Value stored in the test array doesn't match input value from the initializer list.";
            iListPtr++;
        }
        EXPECT_TRUE(iListPtr == iList.end()) << "The test array contains less elements than the initializer list.";

    }
    catch (const std::exception &e)
    {
        caughtException = true;
        debug << "Exception message: " << e.what() << std::endl;
    }
    EXPECT_FALSE(caughtException) << "Construction destruction of \"IrregArray1D\" failed with exception";

}

TYPED_TEST(ParamIrregArray2DTest, Construct2DArrayFromInitializerList)
{
    typedef TypeParam                                     array_type;
    typedef typename array_type::value_type               value_type;
    typedef typename array_type::size_type                index_type;
    typedef std::initializer_list<value_type>          l1_ilist_type;
    typedef std::initializer_list<l1_ilist_type>       l2_ilist_type;
    typedef const value_type*                      l1_ilist_ptr_type;
    typedef const l1_ilist_type*                   l2_ilist_ptr_type;

    bool caughtException = false;
    try
    {
        const l2_ilist_type iList = {
            {0, 1, 2, 3, 4, 5},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            {0, 1, 2},
            {0, 1, 2, 3, 4},
            {0},
            {0, 1, 2, 3},
            {0, 1, 2, 3, 4, 5, 6, 7, 8},
            {0, 1, 2, 3, 4, 5, 6},
            {0, 1, 2, 3, 4, 5, 6, 7},
            {0, 1}
        };
        array_type          testArr = iList;

        this->debug << "Initializer list:\n" << iList   << std::endl;
        this->debug << "testArr         :\n" << testArr << std::endl;

        l2_ilist_ptr_type iListPtrL2 = iList.begin();
        for (index_type j = 0; j < testArr.getLength1(); ++j)
        {
            EXPECT_TRUE(iListPtrL2 != iList.end()) << "The test array contains more elements than the initializer list (L2).";
            l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
            for (index_type i = 0; i < testArr.getLength2(j); ++i)
            {
                EXPECT_TRUE(iListPtrL1 != iListPtrL2->end()) << "The test array contains more elements than the initializer list (L1).";
                EXPECT_EQ(testArr[j][i], *iListPtrL1) << "Value stored in the test array doesn't match input value from the initializer list.";
                iListPtrL1++;
            }
            EXPECT_TRUE(iListPtrL1 == iListPtrL2->end()) << "The test array contains less elements than the initializer list (L1).";
            iListPtrL2++;
        }
        EXPECT_TRUE(iListPtrL2 == iList.end()) << "The test array contains less elements than the initializer list (L2).";

    }
    catch (const std::exception &e)
    {
        caughtException = true;
        this->debug << "Exception message: " << e.what() << std::endl;
    }
    EXPECT_FALSE(caughtException) << "Construction destruction of \"IrregArray2D\" failed with exception";
}

TYPED_TEST(ParamIrregArray3DTest, Construct3DArrayFromInitializerList)
{
    typedef TypeParam                                     array_type;
    typedef typename array_type::value_type               value_type;
    typedef typename array_type::size_type                index_type;
    typedef std::initializer_list<value_type>          l1_ilist_type;
    typedef std::initializer_list<l1_ilist_type>       l2_ilist_type;
    typedef std::initializer_list<l2_ilist_type>       l3_ilist_type;
    typedef const value_type*                      l1_ilist_ptr_type;
    typedef const l1_ilist_type*                   l2_ilist_ptr_type;
    typedef const l2_ilist_type*                   l3_ilist_ptr_type;

    bool caughtException = false;
    try
    {
        const l3_ilist_type iList = {
            {
                {0, 1, 2, 3, 4, 5},
                {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                {0, 1, 2},
                {0, 1, 2, 3, 4},
                {0},
                {0, 1, 2, 3},
                {0, 1, 2, 3, 4, 5, 6, 7, 8},
                {0, 1, 2, 3, 4, 5, 6},
                {0, 1, 2, 3, 4, 5, 6, 7},
                {0, 1}
            },
            {
                {0},
                {0, 1, 2, 3},
                {0, 1, 2, 3, 4, 5, 6, 7, 8},
                {0, 1, 2, 3, 4},
                {0, 1}
            },
            {
                {0},
                {0, 1, 2, 3},
                {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                {0, 1, 2},
            }
        };
        array_type          testArr = iList;

        this->debug << "Initializer list:\n" << iList   << std::endl;
        this->debug << "testArr         :\n" << testArr << std::endl;

        l3_ilist_ptr_type iListPtrL3 = iList.begin();
        for (index_type k = 0; k < testArr.getLength1(); ++k)
        {
            EXPECT_TRUE(iListPtrL3 != iList.end()) << "The test array contains more elements than the initializer list (L3).";
            l2_ilist_ptr_type iListPtrL2 = iListPtrL3->begin();
            for (index_type j = 0; j < testArr.getLength2(k); ++j)
            {
                EXPECT_TRUE(iListPtrL2 != iListPtrL3->end()) << "The test array contains more elements than the initializer list (L2).";
                l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
                for (index_type i = 0; i < testArr.getLength3(k, j); ++i)
                {
                    EXPECT_TRUE(iListPtrL1 != iListPtrL2->end()) << "The test array contains more elements than the initializer list (L1).";
                    EXPECT_EQ(testArr(k, j, i), *iListPtrL1) << "Value stored in the test array doesn't match input value from the initializer list.";
                    iListPtrL1++;
                }
                EXPECT_TRUE(iListPtrL1 == iListPtrL2->end()) << "The test array contains less elements than the initializer list (L1).";
                iListPtrL2++;
            }
            EXPECT_TRUE(iListPtrL2 == iListPtrL3->end()) << "The test array contains less elements than the initializer list (L2).";
            iListPtrL3++;
        }
        EXPECT_TRUE(iListPtrL3 == iList.end()) << "The test array contains less elements than the initializer list (L3).";

    }
    catch (const std::exception &e)
    {
        caughtException = true;
        this->debug << "Exception message: " << e.what() << std::endl;
    }
    EXPECT_FALSE(caughtException) << "Construction destruction of \"IrregArray3D\" failed with exception";

}


TYPED_TEST(ParamIrregArray4DTest, Construct4DArrayFromInitializerList)
{
    typedef TypeParam                                     array_type;
    typedef typename array_type::value_type               value_type;
    typedef typename array_type::size_type                index_type;
    typedef std::initializer_list<value_type>          l1_ilist_type;
    typedef std::initializer_list<l1_ilist_type>       l2_ilist_type;
    typedef std::initializer_list<l2_ilist_type>       l3_ilist_type;
    typedef std::initializer_list<l3_ilist_type>       l4_ilist_type;
    typedef const value_type*                      l1_ilist_ptr_type;
    typedef const l1_ilist_type*                   l2_ilist_ptr_type;
    typedef const l2_ilist_type*                   l3_ilist_ptr_type;
    typedef const l3_ilist_type*                   l4_ilist_ptr_type;

    bool caughtException = false;
    try
    {
        const l4_ilist_type iList = {
            {
                {
                    {0},
                    {0, 1}
                },
                {
                    {0},
                    {0, 1, 2, 3},
                    {0, 1, 2, 3, 4, 5, 6, 7, 8},
                    {0, 1, 2, 3, 4},
                    {0, 1}
                },
                {
                    {0},
                    {0, 1, 2, 3},
                    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                    {0, 1, 2}
                }
            },
            {
                {
                    {0},
                    {0, 1, 2, 3},
                    {0, 1, 2, 3, 4, 5, 6, 7, 8},
                    {0, 1, 2, 3, 4},
                    {0, 1}
                },
                {
                    {0, 1}
                },
                {
                    {0},
                    {0, 1, 2, 3},
                    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                    {0, 1, 2}
                }
            }
        };
        array_type          testArr = iList;

        this->debug << "Initializer list:\n" << iList   << std::endl;
        this->debug << "testArr         :\n" << testArr << std::endl;

        l4_ilist_ptr_type iListPtrL4 = iList.begin();
        for (index_type l = 0; l < testArr.getLength1(); ++l)
        {
            EXPECT_TRUE(iListPtrL4 != iList.end()) << "The test array contains more elements than the initializer list (L4).";
            l3_ilist_ptr_type iListPtrL3 = iListPtrL4->begin();
            for (index_type k = 0; k < testArr.getLength2(l); ++k)
            {
                EXPECT_TRUE(iListPtrL3 != iListPtrL4->end()) << "The test array contains more elements than the initializer list (L3).";
                l2_ilist_ptr_type iListPtrL2 = iListPtrL3->begin();
                for (index_type j = 0; j < testArr.getLength3(l, k); ++j)
                {
                    EXPECT_TRUE(iListPtrL2 != iListPtrL3->end()) << "The test array contains more elements than the initializer list (L2).";
                    l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
                    for (index_type i = 0; i < testArr.getLength4(l, k, j); ++i)
                    {
                        EXPECT_TRUE(iListPtrL1 != iListPtrL2->end()) << "The test array contains more elements than the initializer list (L1).";
                        EXPECT_EQ(testArr(l, k, j, i), *iListPtrL1) << "Value stored in the test array doesn't match input value from the initializer list.";
                        iListPtrL1++;
                    }
                    EXPECT_TRUE(iListPtrL1 == iListPtrL2->end()) << "The test array contains less elements than the initializer list (L1).";
                    iListPtrL2++;
                }
                EXPECT_TRUE(iListPtrL2 == iListPtrL3->end()) << "The test array contains less elements than the initializer list (L2).";
                iListPtrL3++;
            }
            EXPECT_TRUE(iListPtrL3 == iListPtrL4->end()) << "The test array contains less elements than the initializer list (L3).";
            iListPtrL4++;
        }
        EXPECT_TRUE(iListPtrL4 == iList.end()) << "The test array contains less elements than the initializer list (L4).";
    }
    catch (const std::exception &e)
    {
        caughtException = true;
        this->debug << "Exception message: " << e.what() << std::endl;
    }
    EXPECT_FALSE(caughtException) << "Construction destruction of \"IrregArray4D\" failed with exception";

}


TEST_F(IrregArrayTest, DataStructuresIrregArray1DUsage)
{
    typedef  typename gmx::IrregArray1D<unsigned int>::index_type  index_type;
    constexpr index_type            first1 = 1;
    constexpr index_type            last1  = 6;
    constexpr unsigned int          iniVal = 0;

    gmx::IrregArray1D<unsigned int> test_arr1(first1, last1, iniVal);
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

    gmx::IrregArray1D<unsigned int> test_arr2(first1, last1, iniVal);
    EXPECT_NO_THROW_GMX(
            {
                test_arr2(1) = 1;
                test_arr2(2) = 2;
                test_arr2(3) = 3;
                test_arr2(4) = 4;
                test_arr2(5) = 5;
                test_arr2(6) = 6;
            }
            ) << "Caught exception while assigning via operator()." << endl;

    EXPECT_TRUE(test_arr2(1) == 1 &&
                test_arr2(2) == 2 &&
                test_arr2(3) == 3 &&
                test_arr2(4) == 4 &&
                test_arr2(5) == 5 &&
                test_arr2(6) == 6)
    << "Test of assignment and retrieval via IrregArray1D::operator() failed.";
}

TEST_F(IrregArrayTest, DataStructuresIrregArray1DResize)
{
    typedef  typename gmx::IrregArray1D<unsigned int>::index_type  index_type;
    constexpr index_type            first1 = 1;
    constexpr index_type            last1  = 6;
    constexpr unsigned int          iniVal = 0;
    constexpr index_type            first2 = 0;
    constexpr index_type            last2  = 6;
    constexpr unsigned int          newVal = 1;
    gmx::IrregArray1D<unsigned int> test_arr1(first1, last1, iniVal);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(first2, last2, newVal);
            }
            ) << "Caught exception while reinitializing array with new bounds.";

    EXPECT_TRUE(test_arr1[0] == newVal &&
                test_arr1[1] == newVal &&
                test_arr1[2] == newVal &&
                test_arr1[3] == newVal &&
                test_arr1[4] == newVal &&
                test_arr1[5] == newVal &&
                test_arr1[6] == newVal)
    << "Test of resize and assignment failed.";
}

TEST_F(IrregArrayTest, IrregArray1DCopy)
{
    const gmx::IrregArray1D<unsigned int> test_arr1 = {1, 2, 3, 4, 5, 6};

    gmx::IrregArray1D<unsigned int>       test_arr2(test_arr1);

    EXPECT_TRUE(test_arr2(0) == 1 &&
                test_arr2(1) == 2 &&
                test_arr2(2) == 3 &&
                test_arr2(3) == 4 &&
                test_arr2(4) == 5 &&
                test_arr2(5) == 6)
    << "Test of copy constructor failed.";

    debug << "test_arr1: " << endl << test_arr1 << endl;
    debug << "test_arr2: " << endl << test_arr2 << endl;

// clang-tidy complains about an unnecessary copy
#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX(
            {
                const gmx::IrregArray1D<unsigned int>  test_arr3;
                const gmx::IrregArray1D<unsigned int>  test_arr4(test_arr3);
                debug << "test_arr4 (should be empty) = " << test_arr4 << std::endl;
            }
            ) << "Caught exception while copy constructing an empty array.";
#endif
}

TEST_F(IrregArrayTest, IrregArray1DCast)
{
    const gmx::IrregArray1D<unsigned int>       test_arr1 = {1, 2, 3, 4, 5, 6};

    const gmx::IrregArray1D<unsigned short int> test_arr2 =
        static_cast<gmx::IrregArray1D<unsigned short int> >(test_arr1);

    EXPECT_TRUE(test_arr2(0) == 1 &&
                test_arr2(1) == 2 &&
                test_arr2(2) == 3 &&
                test_arr2(3) == 4 &&
                test_arr2(4) == 5 &&
                test_arr2(5) == 6)
    << "Test of explicit conversion operator failed.";

    debug << "test_arr1: " << endl << test_arr1 << endl;
    debug << "test_arr2: " << endl << test_arr2 << endl;
}

TEST_F(IrregArrayTest, IrregArray1DMove)
{
// clang-tidy complains about test_arrN being used after it was moved
// and about the values in the initializer list being magic constants
// since test_arr1 is not const. Both issues are unavoidable here.
#ifndef __clang_analyzer__
    gmx::IrregArray1D<unsigned int> test_arr1 = { 1, 2, 4, 6, 7, 8 };

    debug << "test_arr1 before moving: " << test_arr1 << endl;

    gmx::IrregArray1D<unsigned int> test_arr2(std::move(test_arr1));

    debug << "test_arr2 after move construction: " << test_arr2 << endl;

    EXPECT_TRUE(test_arr1.getSize() == sizeof(decltype(test_arr1))) << "Array was not actually moved";

    EXPECT_TRUE(test_arr2(0) == 1 &&
                test_arr2(1) == 2 &&
                test_arr2(2) == 4 &&
                test_arr2(3) == 6 &&
                test_arr2(4) == 7 &&
                test_arr2(5) == 8)
    << "Test of move constructor failed.";

    test_arr1 = std::move(test_arr2);

    EXPECT_TRUE(test_arr2.getSize() == sizeof(decltype(test_arr1))) << "Array was not actually moved";

    EXPECT_TRUE(test_arr1(0) == 1 &&
                test_arr1(1) == 2 &&
                test_arr1(2) == 4 &&
                test_arr1(3) == 6 &&
                test_arr1(4) == 7 &&
                test_arr1(5) == 8)
    << "Test of move constructor failed.";

    debug << "test_arr1 after move assignment: " << test_arr1 << endl;
#endif
}

TYPED_TEST(ParamIrregArray2DTest, 2DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    array_type test_arr(1, 2, 1, 3, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr[1][1] = 1;
                test_arr[1][2] = 2;
                test_arr[1][3] = 3;
                test_arr[2][1] = 4;
                test_arr[2][2] = 5;
                test_arr[2][3] = 6;
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr: " << endl << test_arr << endl;

    array_type test_arr2(1, 2, 1, 3, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr2(1, 1) = 1;
                test_arr2(1, 2) = 2;
                test_arr2(1, 3) = 3;
                test_arr2(2, 1) = 4;
                test_arr2(2, 2) = 5;
                test_arr2(2, 3) = 6;
            }
            )<< "Test of assignment to a regular array via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1, 1) == 1 &&
                test_arr2(1, 2) == 2 &&
                test_arr2(1, 3) == 3 &&
                test_arr2(2, 1) == 4 &&
                test_arr2(2, 2) == 5 &&
                test_arr2(2, 3) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1, 1) == test_arr[1][1] &&
                test_arr(1, 2) == test_arr[1][2] &&
                test_arr(1, 3) == test_arr[1][3] &&
                test_arr(2, 1) == test_arr[2][1] &&
                test_arr(2, 2) == test_arr[2][2] &&
                test_arr(2, 3) == test_arr[2][3])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    // let the array start with index 0 at variable stripe length and initialize it with value 99
    const size_1d_array_type sizes {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr value_type iniVal = 99;
    array_type           test_arr3(sizes, iniVal);

    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(j); ++j)
                    {
                        test_arr3[i][j] = static_cast<value_type>(j);
                    }
                }
            }
            ) << "Test of assignment an irregular array via operator[] failed with an exception.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    constexpr index_type begin1     = 1;
    constexpr index_type end1       = 7;
    constexpr index_type iniLength2 = 0;
    size_1d_array_type   sizes2(begin1, end1, iniLength2);
// Yes, the lengths of the array stripes are arbitrary numbers,
// just for testing but naming each of them doesn't help.
#ifndef __clang_analyzer__
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
#endif
    constexpr index_type length3 = 2;
    array_type           test_arr4(sizes2, length3);
    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(j); ++j)
                    {
                        test_arr4[i][j] = static_cast<value_type>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
                    {
                        if (test_arr4(i, j) != test_arr4[i][j])
                        {
                            matchedAll = false;
                            break;
                        }
                    }
                }
            }
            ) << "Consistency check for an irregular array between data accession via operator[] and operator() failed with an exception.";

    EXPECT_TRUE(matchedAll)
    << "Consistency check for an irregular array between accession via operator[] and operator() failed.";
}


TYPED_TEST(ParamIrregArray2DTest, 2DArrayResize)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    constexpr index_type first1A = 1;
    constexpr index_type last1A  = 2;
    constexpr index_type first2A = 1;
    constexpr index_type last2A  = 3;
    constexpr value_type iniValA = 0;

    constexpr index_type first1B = 1;
    constexpr index_type last1B  = 6;
    constexpr index_type first2B = 1;
    constexpr index_type last2B  = 6;
    constexpr value_type iniValB = 2;
    array_type           test_arr1(first1A, last1A, first2A, last2A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(first1B, last1B, first2B, last2B, iniValB);
            }
            ) << "Test of resize from regular to regular array dimensions failed with an exception.";

    EXPECT_TRUE(test_arr1.getLength1() == static_cast<size_type>(last1B - first1B + 1) &&
                test_arr1.getLength2() == static_cast<size_type>(last2B - first2B + 1))
    << "Test of array resize failed because actual array dimensions != requested array dimensions";

    // let the array start with index 1
    constexpr index_type first1C     = 1;
    constexpr index_type last1C      = 7;
    constexpr size_type  iniLength2C = 0;
    constexpr value_type iniValC     = 2;
    size_1d_array_type   sizesC(first1C, last1C, iniLength2C);
#ifndef __clang_analyzer__
    sizesC[1] = 2;
    sizesC[2] = 1;
    sizesC[3] = 4;
    sizesC[4] = 8;
    sizesC[5] = 5;
    sizesC[6] = 3;
    sizesC[7] = 9;
#endif
    // let the array start with index 0 and vary the lengths in dimension 2
    constexpr index_type     first1D = 0;
    constexpr index_type     last1D  = 7;
    constexpr value_type     iniValD = 3;
    const size_1d_array_type sizesD {
        5, 9, 1, 4, 8, 5, 3, 2
    };

    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(first1C, last1C, sizesC, iniValC);
            }
            ) << "Test of resize from regular to irregular array dimensions failed with an exception.";

    for (index_type i = test_arr1.getFirst1(); i <= test_arr1.getLast1(); ++i)
    {
        EXPECT_EQ(sizesC(i), test_arr1.getLength2(i))
        << "Test of array resize failed because actual array dimensions != requested array dimensions";
    }

    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(first1D, last1D, sizesD, iniValD);
            }
            ) << "Test of resize from irregular to irregular array dimensions failed with an exception.";

    for (index_type i = test_arr1.getFirst1(); i <= test_arr1.getLast1(); ++i)
    {
        EXPECT_EQ(sizesD(i), test_arr1.getLength2(i))
        << "Test of array resize failed because actual array dimensions != requested array dimensions";
    }

    constexpr index_type first1E = 1;
    constexpr index_type last1E  = 2;
    constexpr index_type first2E = 1;
    constexpr index_type last2E  = 3;
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(first1E, last1E, first2E, last2E);
            }
            ) << "Test of resize from irregular to regular array dimensions failed with an exception.";

    EXPECT_TRUE(test_arr1.getLength1() == static_cast<size_type>(last1E - first1E + 1) &&
                test_arr1.getLength2() == static_cast<size_type>(last2E - first2E + 1))
    << "Test of array resize failed because actual array dimensions != requested array dimensions";
}


TYPED_TEST(ParamIrregArray2DTest, 2DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    const array_type test_arr1 {
        {
            1, 2, 3
        },
        {
            4, 5, 6
        },
    };

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    array_type test_arr2(test_arr1);

    EXPECT_TRUE(test_arr1.getLength1() == test_arr2.getLength1() &&
                test_arr1.getLength2() == test_arr2.getLength2())
    << "Test of copy constructor for a regular array via operator() failed because of non-matching array dimensions";

    EXPECT_TRUE(test_arr2(0, 0) == test_arr1(0, 0) &&
                test_arr2(0, 1) == test_arr1(0, 1) &&
                test_arr2(0, 2) == test_arr1(0, 2) &&
                test_arr2(1, 0) == test_arr1(1, 0) &&
                test_arr2(1, 1) == test_arr1(1, 1) &&
                test_arr2(1, 2) == test_arr1(1, 2))
    << "Test of copy constructor for a regular array via operator() failed because of non-matching array content";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    array_type test_arr3 = test_arr1;

    EXPECT_TRUE(test_arr1.getLength1() == test_arr3.getLength1() &&
                test_arr1.getLength2() == test_arr3.getLength2())
    << "Test of copy constructor for a regular array via operator() failed because of non-matching array dimensions";

    EXPECT_TRUE(test_arr3(0, 0) == test_arr1(0, 0) &&
                test_arr3(0, 1) == test_arr1(0, 1) &&
                test_arr3(0, 2) == test_arr1(0, 2) &&
                test_arr3(1, 0) == test_arr1(1, 0) &&
                test_arr3(1, 1) == test_arr1(1, 1) &&
                test_arr3(1, 2) == test_arr1(1, 2))
    << "Test of assignment operator= for a regular array failed.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;


    // let the array start with index 1, and initialize it with value 2
    constexpr index_type begin1  = 1;
    constexpr index_type end1    = 7;
    constexpr size_type  iniLen1 = 0;
    constexpr value_type iniVal1 = 2;
    size_1d_array_type   sizes(begin1, end1, iniLen1);
#ifndef __clang_analyzer__
    sizes[1] = 2;
    sizes[2] = 1;
    sizes[3] = 4;
    sizes[4] = 8;
    sizes[5] = 5;
    sizes[6] = 3;
    sizes[7] = 9;
#endif
    array_type test_arr4(begin1, end1, sizes, iniVal1);
    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(j); ++j)
                    {
                        test_arr4[i][j] = static_cast<value_type>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    array_type test_arr5(test_arr4);

    bool       matchedAll = true;
    for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
    {
        for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
        {
            if (test_arr5(i, j) != test_arr4(i, j))
            {
                matchedAll = false;
                break;
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Test of copy constructor for an irregular array failed.";

    array_type test_arr6 = test_arr4;

    matchedAll = true;
    for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
    {
        for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
        {
            if (test_arr6(i, j) != test_arr4(i, j))
            {
                matchedAll = false;
                break;
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Test of assignment operator= for an irregular array failed.";

#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX(
            {
                gmx::IrregArray2D<unsigned int>  test_arr5;
                gmx::IrregArray2D<unsigned int>  test_arr6(test_arr5);
                this->debug << "test_arr6 (should be empty) = " << test_arr6 << std::endl;
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
#endif
}


TYPED_TEST(ParamIrregArray2DTest, 2DArrayMove)
{
// clang-tidy complains about test_arrN being used after it was moved
// and about the values in the initializer list being magic constants
// since test_arr1 is not const. Both issues are unavoidable here.
#ifndef __clang_analyzer__
    typedef TypeParam array_type;

    array_type test_arr1 =
    {
        {1},
        {2, 3},
        {4, 5, 6},
        {7, 8},
    };
    this->debug << "test_arr1 before moving:\n" << test_arr1 << endl;

    array_type test_arr2(std::move(test_arr1));

    EXPECT_TRUE(test_arr1.getSize() == sizeof(decltype(test_arr1))) << "Array was not actually moved";

    this->debug << "test_arr2 after move construction:\n" << test_arr2 << endl;

    EXPECT_TRUE(test_arr2(0, 0) == 1 &&
                test_arr2(1, 0) == 2 &&
                test_arr2(1, 1) == 3 &&
                test_arr2(2, 0) == 4 &&
                test_arr2(2, 1) == 5 &&
                test_arr2(2, 2) == 6 &&
                test_arr2(3, 0) == 7 &&
                test_arr2(3, 1) == 8
                )
    << "Test of move constructor failed.";

    test_arr1 = std::move(test_arr2);

    EXPECT_TRUE(test_arr2.getSize() == sizeof(decltype(test_arr2))) << "Array was not actually moved";

    this->debug << "test_arr1 after move assignment:\n" << test_arr1 << endl;

    EXPECT_TRUE(test_arr1(0, 0) == 1 &&
                test_arr1(1, 0) == 2 &&
                test_arr1(1, 1) == 3 &&
                test_arr1(2, 0) == 4 &&
                test_arr1(2, 1) == 5 &&
                test_arr1(2, 2) == 6 &&
                test_arr1(3, 0) == 7 &&
                test_arr1(3, 1) == 8
                )
    << "Test of move assignment operator failed.";
#endif
}


TYPED_TEST(ParamIrregArray3DTest, 3DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    constexpr index_type first1A = 1;
    constexpr index_type last1A  = 2;
    constexpr index_type first2A = 1;
    constexpr index_type last2A  = 3;
    constexpr index_type first3A = 1;
    constexpr index_type last3A  = 2;
    constexpr value_type iniValA = 0;
    array_type           test_arr1(first1A, last1A, first2A, last2A, first3A, last3A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                size_type v = 3;
                for (index_type i = first1A; i <= last1A; ++i)
                {
                    for (index_type j = first2A; j <= last2A; ++j)
                    {
                        for (index_type k = first3A; k <= last3A; ++k)
                        {
                            v++;
                            test_arr1[i][j][k] = static_cast<value_type>(v);
                        }
                    }
                }
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    {
        size_type v = 3;
        for (index_type i = first1A; i <= last1A; ++i)
        {
            for (index_type j = first2A; j <= last2A; ++j)
            {
                for (index_type k = first3A; k <= last3A; ++k)
                {
                    v++;
                    EXPECT_EQ(test_arr1(i, j, k),  static_cast<value_type>(v))
                    << "Test of assignment to a regular array and retrieval via operator() failed.";
                }
            }
        }
    }

    constexpr index_type first1B = 1;
    constexpr index_type last1B  = 2;
    constexpr index_type first2B = 1;
    constexpr index_type last2B  = 3;
    constexpr index_type first3B = 1;
    constexpr index_type last3B  = 2;
    constexpr value_type iniValB = 0;
    array_type           test_arr2(first1B, last1B, first2B, last2B, first3B, last3B, iniValB);
    EXPECT_NO_THROW_GMX(
            {
                size_type v = 0;
                for (index_type i = first1A; i <= last1A; ++i)
                {
                    for (index_type j = first2A; j <= last2A; ++j)
                    {
                        for (index_type k = first3A; k <= last3A; ++k)
                        {
                            v++;
                            test_arr2(i, j, k) = static_cast<value_type>(v);
                        }
                    }
                }
            }
            ) << "Test of assignment to a regular array via operator() failed with an exception.";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    {
        size_type v = 0;
        for (index_type i = first1A; i <= last1A; ++i)
        {
            for (index_type j = first2A; j <= last2A; ++j)
            {
                for (index_type k = first3A; k <= last3A; ++k)
                {
                    v++;
                    EXPECT_EQ(test_arr2(i, j, k),  static_cast<value_type>(v))
                    << "Test of assignment to a regular array and retrieval via operator() failed.";

                    EXPECT_EQ(test_arr2(i, j, k), test_arr2[i][j][k])
                    << "Consistency check for a regular array between accession via operator[] and operator() failed.";
                }
            }
        }
    }

    const size_1d_array_type sizes1C {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr size_type  length2C = 7;
    constexpr value_type iniValC  = 99;
    // let the array start with index 0, and initialize it with value 7
    size_2d_array_type   sizes2C(sizes1C, length2C);

    array_type           test_arr3(sizes2C, iniValC);
    EXPECT_NO_THROW_GMX(
            {
                size_type l = 1;
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
                        {
                            test_arr3(i, j, k) = static_cast<value_type>(l);
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    constexpr index_type first1D  = 1;
    constexpr index_type last1D   = 7;
    constexpr size_type  iniLen2D = 0;
    size_1d_array_type   sizes1D(first1D, last1D, iniLen2D);
#ifndef __clang_analyzer__
    sizes1D[1] = 2;
    sizes1D[2] = 1;
    sizes1D[3] = 4;
    sizes1D[4] = 8;
    sizes1D[5] = 5;
    sizes1D[6] = 3;
    sizes1D[7] = 9;
#endif
    // let the array start with index 0, and initialize it with index 1
    constexpr size_type      len2D = 2;
    const size_2d_array_type sizes2D(sizes1D, len2D);
    constexpr value_type     iniValD = 99;

    array_type               test_arr4(sizes2D, iniValD);
    EXPECT_NO_THROW_GMX(
            {
                size_type l = 1;
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
                        {
                            test_arr4(i, j, k) = static_cast<value_type>(l);
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
                        {
                            if (test_arr4(i, j, k) != test_arr4[i][j][k])
                            {
                                matchedAll = false;
                                break;
                            }
                        }
                    }
                }
            }
            ) << "Consistency check for an irregular array between data accession via operator[] and operator() failed with an exception.";

    EXPECT_TRUE(matchedAll)
    << "Consistency check for an irregular array between accession via operator[] and operator() failed.";

#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX(
            {
                gmx::IrregArray3D<unsigned int>  test_arr5;
                gmx::IrregArray3D<unsigned int>  test_arr6(test_arr5);
                this->debug << "test_arr6 (should be empty) = " << test_arr6 << std::endl;
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
#endif
}

TYPED_TEST(ParamIrregArray3DTest, 3DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    constexpr index_type first1A = 1;
    constexpr index_type last1A  = 2;
    constexpr index_type first2A = 1;
    constexpr index_type last2A  = 3;
    constexpr index_type first3A = 1;
    constexpr index_type last3A  = 2;
    constexpr value_type iniValA = 0;
    array_type           test_arr1(first1A, last1A, first2A, last2A, first3A, last3A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                size_type v = 0;
                for (index_type i = first1A; i <= last1A; ++i)
                {
                    for (index_type j = first2A; j <= last2A; ++j)
                    {
                        for (index_type k = first3A; k <= last3A; ++k)
                        {
                            v++;
                            test_arr1[i][j][k] = static_cast<value_type>(v);
                        }
                    }
                }
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    array_type test_arr2(test_arr1);

    this->debug << "test_arr1: " << endl << test_arr1 << endl;
    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    bool matchedAll = true;
    for (index_type i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
    {
        for (index_type j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
        {
            for (index_type k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
            {
                if (test_arr1(i, j, k) != test_arr2(i, j, k))
                {
                    matchedAll = false;
                    break;
                }
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Copy construction of regular array failed";


    // create an irregular array with variable length in dimension 2
    const size_1d_array_type sizes1B {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr size_type  length3B = 7;
    constexpr value_type iniValB  = 99;
    size_2d_array_type   sizes2B(sizes1B.getFirst1(), sizes1B.getLast1(), sizes1B, length3B);

    array_type           test_arr3(sizes2B, iniValB);
    EXPECT_NO_THROW_GMX(
            {
                size_type l = 1;
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
                        {
                            test_arr3[i][j][k] = static_cast<value_type>(l);
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    array_type test_arr4(test_arr3);

    matchedAll = true;
    for (ssize_t i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
    {
        for (ssize_t j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
        {
            for (ssize_t k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
            {
                if (test_arr3(i, j, k) != test_arr4(i, j, k))
                {
                    matchedAll = false;
                    break;
                }
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Copy construction of irregular array failed";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;
    this->debug << "test_arr4: " << endl << test_arr3 << endl;

// quite complaints about an unnecessary copy
#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX(
            {
                array_type  test_arr5;
                array_type  test_arr6(test_arr5);
                this->debug << "test_arr6 (should be empty) = " << test_arr6 << std::endl;
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
#endif
}


TYPED_TEST(ParamIrregArray3DTest, 3DArrayMove)
{
// quite complaints about magic values in the initializer_lists and
// variables used after moving them, both are necessary here
#ifndef __clang_analyzer__
    typedef TypeParam array_type;

    array_type test_arr1 =
    {
        {
            {1},
            {2, 3},
            {4, 5, 6},
            {7, 8},
        },
        {
            {4, 5, 6},
            {1},
            {2, 3},
            {7, 8},
        },
    };
    this->debug << "test_arr1 before moving:\n" << test_arr1 << endl;

    array_type test_arr2(std::move(test_arr1));

    EXPECT_TRUE(test_arr1.getSize() == sizeof(decltype(test_arr1))) << "Array was not actually moved";

    this->debug << "test_arr2 after move construction:\n" << test_arr2 << endl;

    EXPECT_TRUE(test_arr2(0, 0, 0) == 1 &&
                test_arr2(0, 1, 0) == 2 &&
                test_arr2(0, 1, 1) == 3 &&
                test_arr2(0, 2, 0) == 4 &&
                test_arr2(0, 2, 1) == 5 &&
                test_arr2(0, 2, 2) == 6 &&
                test_arr2(0, 3, 0) == 7 &&
                test_arr2(0, 3, 1) == 8 &&
                test_arr2(1, 0, 0) == 4 &&
                test_arr2(1, 0, 1) == 5 &&
                test_arr2(1, 0, 2) == 6 &&
                test_arr2(1, 1, 0) == 1 &&
                test_arr2(1, 2, 0) == 2 &&
                test_arr2(1, 2, 1) == 3 &&
                test_arr2(1, 3, 0) == 7 &&
                test_arr2(1, 3, 1) == 8
                )
    << "Test of move constructor failed.";

    test_arr1 = std::move(test_arr2);

    EXPECT_TRUE(test_arr2.getSize() == sizeof(decltype(test_arr2))) << "Array was not actually moved";

    this->debug << "test_arr1 after move assignment:\n" << test_arr1 << endl;

    EXPECT_TRUE(test_arr1(0, 0, 0) == 1 &&
                test_arr1(0, 1, 0) == 2 &&
                test_arr1(0, 1, 1) == 3 &&
                test_arr1(0, 2, 0) == 4 &&
                test_arr1(0, 2, 1) == 5 &&
                test_arr1(0, 2, 2) == 6 &&
                test_arr1(0, 3, 0) == 7 &&
                test_arr1(0, 3, 1) == 8 &&
                test_arr1(1, 0, 0) == 4 &&
                test_arr1(1, 0, 1) == 5 &&
                test_arr1(1, 0, 2) == 6 &&
                test_arr1(1, 1, 0) == 1 &&
                test_arr1(1, 2, 0) == 2 &&
                test_arr1(1, 2, 1) == 3 &&
                test_arr1(1, 3, 0) == 7 &&
                test_arr1(1, 3, 1) == 8
                )
    << "Test of move assignment operator failed.";
#endif          // __clang_analyzer__
}


TYPED_TEST(ParamIrregArray4DTest, 4DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;

    constexpr index_type first1A = 1;
    constexpr index_type last1A  = 2;
    constexpr index_type first2A = 1;
    constexpr index_type last2A  = 3;
    constexpr index_type first3A = 1;
    constexpr index_type last3A  = 2;
    constexpr index_type first4A = 1;
    constexpr index_type last4A  = 7;
    constexpr value_type iniValA = 0;
    array_type           test_arr1(first1A, last1A, first2A, last2A, first3A, last3A, first4A, last4A, iniValA);
    EXPECT_NO_THROW_GMX (
            {
                size_type v = 0;
                for (index_type i = test_arr1.getFirst1(); i <= test_arr1.getLast1(); ++i)
                {
                    for (index_type j = test_arr1.getFirst2(); j <= test_arr1.getLast2(); ++j)
                    {
                        for (index_type k = test_arr1.getFirst3(); k <= test_arr1.getLast3(); ++k)
                        {
                            for (index_type l = test_arr1.getFirst4(); l <= test_arr1.getLast4(); ++l)
                            {
                                v++;
                                test_arr1[i][j][k][l] = static_cast<value_type>(v);
                            }
                        }
                    }
                }
            }
            )
    << "Test of assignment to a regular array via operator[] failed with an exception.";

    {
        size_type v = 0;
        for (index_type i = test_arr1.getFirst1(); i <= test_arr1.getLast1(); ++i)
        {
            for (index_type j = test_arr1.getFirst2(); j <= test_arr1.getLast2(); ++j)
            {
                for (index_type k = test_arr1.getFirst3(); k <= test_arr1.getLast3(); ++k)
                {
                    for (index_type l = test_arr1.getFirst4(); l <= test_arr1.getLast4(); ++l)
                    {
                        v++;
                        EXPECT_EQ(test_arr1[i][j][k][l], static_cast<value_type>(v))
                        << "retrieval of assigned values from a regular 4D array failed";

                        EXPECT_EQ(test_arr1[i][j][k][l], test_arr1(i, j, k, l))
                        << "Consistency check for a regular array between accession via operator[] and operator() failed.";
                    }
                }
            }
        }
    }

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    constexpr index_type first1B = 1;
    constexpr index_type last1B  = 2;
    constexpr index_type first2B = 1;
    constexpr index_type last2B  = 3;
    constexpr index_type first3B = 1;
    constexpr index_type last3B  = 1;
    constexpr index_type first4B = 1;
    constexpr index_type last4B  = 1;
    constexpr value_type iniValB = 0;
    array_type           test_arr2(first1B, last1B, first2B, last2B, first3B, last3B, first4B, last4B, iniValB);
    EXPECT_NO_THROW_GMX (
            {
                size_type v = 3;
                for (index_type i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
                {
                    for (index_type j = test_arr2.getFirst2(); j <= test_arr2.getLast2(); ++j)
                    {
                        for (index_type k = test_arr2.getFirst3(); k <= test_arr2.getLast3(); ++k)
                        {
                            for (index_type l = test_arr2.getFirst4(); l <= test_arr2.getLast4(); ++l)
                            {
                                v++;
                                test_arr2(i, j, k, l) = static_cast<value_type>(v);
                            }
                        }
                    }
                }
            }
            )
    << "Test of assignment to a regular array via operator() failed with an exception.";

    {
        size_type v = 3;
        for (index_type i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
        {
            for (index_type j = test_arr2.getFirst2(); j <= test_arr2.getLast2(); ++j)
            {
                for (index_type k = test_arr2.getFirst3(); k <= test_arr2.getLast3(); ++k)
                {
                    for (index_type l = test_arr2.getFirst4(); l <= test_arr2.getLast4(); ++l)
                    {
                        v++;
                        EXPECT_EQ(test_arr2[i][j][k][l], static_cast<value_type>(v))
                        << "retrieval of assigned values from a regular 4D array failed";

                        EXPECT_EQ(test_arr2[i][j][k][l], test_arr2(i, j, k, l))
                        << "Consistency check for a regular array between accession via operator[] and operator() failed.";
                    }
                }
            }
        }
    }

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    constexpr size_type      length3C = 3;
    constexpr size_type      length4C = 1;
    constexpr value_type     iniValC  = 99;
    const size_1d_array_type sizes1C {
        2, 1, 4, 8, 5, 3, 9
    };
    const size_2d_array_type sizes2C(sizes1C.getFirst1(), sizes1C.getLast1(), sizes1C, length3C);
    const size_3d_array_type sizes3C(sizes1C.getFirst1(), sizes1C.getLast1(), sizes1C, sizes2C, length4C);

    array_type               test_arr3(sizes1C.getFirst1(), sizes1C.getLast1(), sizes1C, sizes2C, sizes3C, iniValC);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
                        {
                            for (index_type l = test_arr3.getFirst4(); l <= test_arr3.getLast4(i, j, k); ++l)
                            {
                                test_arr3[i][j][k][l] = static_cast<value_type>(m);
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";
    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    array_type         test_arr4(sizes3C, iniValC);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
                        {
                            for (index_type l = test_arr3.getFirst4(); l <= test_arr3.getLast4(i, j, k); ++l)
                            {
                                test_arr4[i][j][k][l] = static_cast<value_type>(m);
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    constexpr index_type first13D     = 1;
    constexpr size_type  last13D      = 7;
    constexpr size_type  iniLength13D = 1;
    constexpr value_type iniValD      = 5;
    size_1d_array_type   sizes13D(first13D, last13D, iniLength13D);
#ifndef __clang_analyzer__
    sizes13D[1] = 2;
    sizes13D[2] = 1;
    sizes13D[3] = 4;
    sizes13D[4] = 8;
    sizes13D[5] = 5;
    sizes13D[6] = 3;
    sizes13D[7] = 9;
#endif
    array_type test_arr5(first13D, last13D, sizes13D, first13D, last13D, sizes13D, iniValD);
    this->debug << "test_arr5: " << endl << test_arr5 << endl;

    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (index_type i = test_arr5.getFirst1(); i <= test_arr5.getLast1(); ++i)
                {
                    for (index_type j = test_arr5.getFirst2(); j <= test_arr5.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr5.getFirst3(); k <= test_arr5.getLast3(i, j); ++k)
                        {
                            for (index_type l = test_arr5.getFirst4(); l <= test_arr5.getLast4(i, j, k); ++l)
                            {
                                test_arr5(i, j, k, l) = static_cast<value_type>(m);
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    bool matchedAll = true;

    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (ssize_t j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr4.getFirst4(); l <= test_arr4.getLast4(i, j, k); ++l)
                            {
                                if (test_arr4(i, j, k, l) != test_arr4[i][j][k][l])
                                {
                                    matchedAll = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            ) << "Consistency check for an irregular array between data accession via operator[] and operator() failed with an exception.";

    EXPECT_TRUE(matchedAll)
    << "Consistency check for an irregular array between accession via operator[] and operator() failed.";
}


TYPED_TEST(ParamIrregArray4DTest, 4DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;


    array_type test_arr1(1, 2, 1, 3, 1, 1, 1, 1, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1(1, 1, 1, 1) = 1;
                test_arr1(1, 2, 1, 1) = 2;
                test_arr1(1, 3, 1, 1) = 3;
                test_arr1(2, 1, 1, 1) = 4;
                test_arr1(2, 2, 1, 1) = 5;
                test_arr1(2, 3, 1, 1) = 6;
            }
            ) << "Test of assignment to a regular array via operator() failed with an exception.";

    array_type test_arr2(test_arr1);

    bool       matchedAll = true;
    for (index_type i = test_arr1.getFirst1(); i <= test_arr1.getLast1(); ++i)
    {
        for (index_type j = test_arr1.getFirst2(); j <= test_arr1.getLast2(i); ++j)
        {
            for (index_type k = test_arr1.getFirst3(); k <= test_arr1.getLast3(i, j); ++k)
            {
                for (index_type l = test_arr1.getFirst4(); l <= test_arr1.getLast4(i, j, k); ++l)
                {
                    if (test_arr1(i, j, k, l) != test_arr2(i, j, k, l))
                    {
                        this->debug << "test_arr1(" << i << ", " << j << ", " << k << ", " << l << ") != "
                        << "test_arr2(" << i << ", " << j << ", " << k << ", " << l << ")"
                        << endl;
                        this->debug << "test_arr1(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr1(i, j, k, l) << endl;
                        this->debug << "test_arr2(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr2(i, j, k, l) << endl;
                        matchedAll = false;
                    }
                }
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Copy construction of regular array failed";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;
    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    // let the array start with index 0, and initialize it with value 99
    const size_1d_array_type sizes1 {
        2, 1, 4, 8, 5, 3, 9
    };
    const size_t       length3 = 3;
    size_2d_array_type sizes2(sizes1, length3);
    const size_t       length4 = 1;
    size_3d_array_type sizes3(sizes2, length4);
    const value_type   iniVal = 99;
    array_type         test_arr3(sizes1.getFirst1(), sizes1.getLast1(), sizes1, sizes2, sizes3, iniVal);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
                        {
                            for (index_type l = test_arr3.getFirst4(); l <= test_arr3.getLast4(i, j, k); ++l)
                            {
                                test_arr3[i][j][k][l] = static_cast<value_type>(m);
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    array_type test_arr4(test_arr3);

    matchedAll = true;
    for (index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
    {
        for (index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(i); ++j)
        {
            for (index_type k = test_arr3.getFirst3(); k <= test_arr3.getLast3(i, j); ++k)
            {
                for (index_type l = test_arr3.getFirst4(); l <= test_arr3.getLast4(i, j, k); ++l)
                {
                    if (test_arr3(i, j, k, l) != test_arr4(i, j, k, l))
                    {
                        this->debug << "test_arr3(" << i << ", " << j << ", " << k << ", " << l << ") != "
                        << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ")"
                        << endl;
                        this->debug << "test_arr3(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr3(i, j, k, l) << endl;
                        this->debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr4(i, j, k, l) << endl;
                        matchedAll = false;
                    }
                }
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Copy construction of regular array failed";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    matchedAll = true;
    array_type test_arr5 = test_arr4;
    for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
    {
        for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
        {
            for (index_type k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
            {
                for (index_type l = test_arr4.getFirst4(); l <= test_arr4.getLast4(i, j, k); ++l)
                {
                    if (test_arr4(i, j, k, l) != test_arr5(i, j, k, l))
                    {
                        this->debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") != "
                        << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ")"
                        << endl;
                        this->debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr4(i, j, k, l) << endl;
                        this->debug << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr5(i, j, k, l) << endl;
                        matchedAll = false;
                    }
                }
            }
        }
    }
    this->debug << "test_arr4: " << endl << test_arr4 << endl;
    this->debug << "test_arr5: " << endl << test_arr5 << endl;
    EXPECT_TRUE(matchedAll) << "Test of assignment operator failed.";

#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX(
            {
                array_type test_arr6;
                array_type test_arr7(test_arr6);
                this->debug << "test_arr7 (should be empty) = " << test_arr7 << std::endl;
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
#endif          // __clang_analyzer__
}

TYPED_TEST(ParamIrregArray4DTest, 4DArrayMove)
{
// quite complaints about magic values in the initializer_lists and
// variables used after moving them, both are necessary here
#ifndef __clang_analyzer__
    typedef TypeParam array_type;

    array_type test_arr1 =
    {
        {
            {
                {1},
                {2, 3},
                {4, 5, 6},
                {7, 8},
            },
            {
                {4, 5, 6},
                {1},
            },
        },
        {
            {
                {4, 5, 6},
                {1},
                {7, 8},
            },
            {
                {1},
                {2, 3},
            },
        },
    };
    this->debug << "test_arr1 before moving:\n" << test_arr1 << endl;

    array_type test_arr2(std::move(test_arr1));

    EXPECT_TRUE(test_arr1.getSize() == sizeof(decltype(test_arr1))) << "Array was not actually moved";

    this->debug << "test_arr2 after move construction:\n" << test_arr2 << endl;

    EXPECT_TRUE(test_arr2(0, 0, 0, 0) == 1 &&
                test_arr2(0, 0, 1, 0) == 2 &&
                test_arr2(0, 0, 1, 1) == 3 &&
                test_arr2(0, 0, 2, 0) == 4 &&
                test_arr2(0, 0, 2, 1) == 5 &&
                test_arr2(0, 0, 2, 2) == 6 &&
                test_arr2(0, 0, 3, 0) == 7 &&
                test_arr2(0, 0, 3, 1) == 8 &&
                test_arr2(0, 1, 0, 0) == 4 &&
                test_arr2(0, 1, 0, 1) == 5 &&
                test_arr2(0, 1, 0, 2) == 6 &&
                test_arr2(0, 1, 1, 0) == 1 &&
                test_arr2(1, 0, 0, 0) == 4 &&
                test_arr2(1, 0, 0, 1) == 5 &&
                test_arr2(1, 0, 0, 2) == 6 &&
                test_arr2(1, 0, 1, 0) == 1 &&
                test_arr2(1, 0, 2, 0) == 7 &&
                test_arr2(1, 0, 2, 1) == 8 &&
                test_arr2(1, 1, 0, 0) == 1 &&
                test_arr2(1, 1, 1, 0) == 2 &&
                test_arr2(1, 1, 1, 1) == 3
                )
    << "Test of move constructor failed.";

    test_arr1 = std::move(test_arr2);

    EXPECT_TRUE(test_arr2.getSize() == sizeof(decltype(test_arr2))) << "Array was not actually moved";

    this->debug << "test_arr1 after move assignment:\n" << test_arr1 << endl;

    EXPECT_TRUE(test_arr1(0, 0, 0, 0) == 1 &&
                test_arr1(0, 0, 1, 0) == 2 &&
                test_arr1(0, 0, 1, 1) == 3 &&
                test_arr1(0, 0, 2, 0) == 4 &&
                test_arr1(0, 0, 2, 1) == 5 &&
                test_arr1(0, 0, 2, 2) == 6 &&
                test_arr1(0, 0, 3, 0) == 7 &&
                test_arr1(0, 0, 3, 1) == 8 &&
                test_arr1(0, 1, 0, 0) == 4 &&
                test_arr1(0, 1, 0, 1) == 5 &&
                test_arr1(0, 1, 0, 2) == 6 &&
                test_arr1(0, 1, 1, 0) == 1 &&
                test_arr1(1, 0, 0, 0) == 4 &&
                test_arr1(1, 0, 0, 1) == 5 &&
                test_arr1(1, 0, 0, 2) == 6 &&
                test_arr1(1, 0, 1, 0) == 1 &&
                test_arr1(1, 0, 2, 0) == 7 &&
                test_arr1(1, 0, 2, 1) == 8 &&
                test_arr1(1, 1, 0, 0) == 1 &&
                test_arr1(1, 1, 1, 0) == 2 &&
                test_arr1(1, 1, 1, 1) == 3
                )
    << "Test of move assignment operator failed.";
#endif          // __clang_analyzer__
}


}         // namespace

}         // namespace data_struct_test

}         // namespace gmx
