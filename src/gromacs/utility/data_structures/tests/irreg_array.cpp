/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
    EXPECT_NO_THROW_GMX({ gmx::IrregArray2D<double> test_2d_irreg_array; })
    << "Construction destruction of \"IrregArray2D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::IrregArray3D<double> test_3d_irreg_array; })
    << "Construction destruction of \"IrregArray3D\" failed";

    EXPECT_NO_THROW_GMX({ gmx::IrregArray4D<double> test_4d_irreg_array; })
    << "Construction destruction of \"IrregArray4D\" failed";
}

TYPED_TEST(ParamIrregArray2DTest, ConstructZeroLengthArray)
{
    typedef  TypeParam                                        array_type;
    typedef  typename array_type::size_1d_array_type  size_1d_array_type;

    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0); })
    << "Construction of regular array with l1 = l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 0); })
    << "Construction of regular array with l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 1); })
    << "Construction of regular array with l1 = 0 failed";

    EXPECT_NO_THROW_GMX({ const size_1d_array_type l2(0, 0); const array_type testArr(l2); })
    << "Construction of irregular array with l1 = l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const size_1d_array_type l2(1, 0); const array_type testArr(l2); })
    << "Construction of irregular array with l2 = 0 failed";
}

TYPED_TEST(ParamIrregArray3DTest, ConstructZeroLengthArray)
{
    typedef  TypeParam                                        array_type;
    typedef  typename array_type::size_2d_array_type  size_2d_array_type;

    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0, 0); })
    << "Construction of regular array with l1 = l2 = l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 0, 0); })
    << "Construction of regular array with l2 = l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 1, 0); })
    << "Construction of regular array with l1 = l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0, 1); })
    << "Construction of regular array with l1 = l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 1, 1); })
    << "Construction of regular array with l1 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 0, 1); })
    << "Construction of regular array with l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 1, 0); })
    << "Construction of regular array with l3 = 0 failed";

    EXPECT_NO_THROW_GMX({ const size_2d_array_type l3(0, 0, 0); const array_type testArr(l3); })
    << "Construction of irregular array with l1 = l2 = l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const size_2d_array_type l3(1, 1, 0); const array_type testArr(l3); })
    << "Construction of irregular array with l3 = 0 failed";
}

TYPED_TEST(ParamIrregArray4DTest, ConstructZeroLengthArray)
{
    typedef  TypeParam                                        array_type;
    typedef  typename array_type::size_3d_array_type  size_3d_array_type;

    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0, 0, 0); })
    << "Construction of regular array with l1 = l2 = l3 = l4 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 0, 0, 0); })
    << "Construction of regular array with l2 = l3 = l4 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 1, 0, 0); })
    << "Construction of regular array with l1 = l3 = l4 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0, 1, 0); })
    << "Construction of regular array with l1 = l2 = l4 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 0, 0, 1); })
    << "Construction of regular array with l1 = l2 = l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(0, 1, 1, 1); })
    << "Construction of regular array with l1 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 0, 1, 1); })
    << "Construction of regular array with l2 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 1, 0, 1); })
    << "Construction of regular array with l3 = 0 failed";
    EXPECT_NO_THROW_GMX({ const array_type testArr(1, 1, 1, 0); })
    << "Construction of regular array with l4 = 0 failed";

    EXPECT_NO_THROW_GMX({ const size_3d_array_type l4(0, 0, 0, 0); const array_type testArr(l4); })
    << "Construction of irregular array with l1 = l2 = l3 = 0 failed";
// clang static analyzer gives a false positive nullptr dereference for the flat variant
#ifndef __clang_analyzer__
    EXPECT_NO_THROW_GMX({ const size_3d_array_type l4(1, 1, 1, 0); const array_type testArr(l4); })
    << "Construction of irregular array with l3 = 0 failed";
#endif
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
        for (index_type j = 0; j < testArr.length1(); ++j)
        {
            EXPECT_TRUE(iListPtrL2 != iList.end()) << "The test array contains more elements than the initializer list (L2).";
            l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
            for (index_type i = 0; i < testArr.length2(j); ++i)
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
        for (index_type k = 0; k < testArr.length1(); ++k)
        {
            EXPECT_TRUE(iListPtrL3 != iList.end()) << "The test array contains more elements than the initializer list (L3).";
            l2_ilist_ptr_type iListPtrL2 = iListPtrL3->begin();
            for (index_type j = 0; j < testArr.length2(k); ++j)
            {
                EXPECT_TRUE(iListPtrL2 != iListPtrL3->end()) << "The test array contains more elements than the initializer list (L2).";
                l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
                for (index_type i = 0; i < testArr.length3(k, j); ++i)
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
        for (index_type l = 0; l < testArr.length1(); ++l)
        {
            EXPECT_TRUE(iListPtrL4 != iList.end()) << "The test array contains more elements than the initializer list (L4).";
            l3_ilist_ptr_type iListPtrL3 = iListPtrL4->begin();
            for (index_type k = 0; k < testArr.length2(l); ++k)
            {
                EXPECT_TRUE(iListPtrL3 != iListPtrL4->end()) << "The test array contains more elements than the initializer list (L3).";
                l2_ilist_ptr_type iListPtrL2 = iListPtrL3->begin();
                for (index_type j = 0; j < testArr.length3(l, k); ++j)
                {
                    EXPECT_TRUE(iListPtrL2 != iListPtrL3->end()) << "The test array contains more elements than the initializer list (L2).";
                    l1_ilist_ptr_type iListPtrL1 = iListPtrL2->begin();
                    for (index_type i = 0; i < testArr.length4(l, k, j); ++i)
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


TYPED_TEST(ParamIrregArray2DTest, 2DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    array_type test_arr(2, 3, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr[0][0] = 1;
                test_arr[0][1] = 2;
                test_arr[0][2] = 3;
                test_arr[1][0] = 4;
                test_arr[1][1] = 5;
                test_arr[1][2] = 6;
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr: " << endl << test_arr << endl;

    array_type test_arr2(2, 3, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr2(0, 0) = 1;
                test_arr2(0, 1) = 2;
                test_arr2(0, 2) = 3;
                test_arr2(1, 0) = 4;
                test_arr2(1, 1) = 5;
                test_arr2(1, 2) = 6;
            }
            )<< "Test of assignment to a regular array via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(0, 0) == 1 &&
                test_arr2(0, 1) == 2 &&
                test_arr2(0, 2) == 3 &&
                test_arr2(1, 0) == 4 &&
                test_arr2(1, 1) == 5 &&
                test_arr2(1, 2) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(0, 0) == test_arr[0][0] &&
                test_arr(0, 1) == test_arr[0][1] &&
                test_arr(0, 2) == test_arr[0][2] &&
                test_arr(1, 0) == test_arr[1][0] &&
                test_arr(1, 1) == test_arr[1][1] &&
                test_arr(1, 2) == test_arr[1][2])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    // let the array start with index 0 at variable stripe length and initialize it with value 99
    const size_1d_array_type sizes {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr value_type iniVal = 99;
    array_type           test_arr3(sizes, iniVal);

    EXPECT_NO_THROW_GMX(
            {
                for (size_type i = 0; i < test_arr3.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr3.length2(j); ++j)
                    {
                        test_arr3[i][j] = static_cast<value_type>(j);
                    }
                }
            }
            ) << "Test of assignment an irregular array via operator[] failed with an exception.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    constexpr index_type length1    = 8;
    constexpr index_type iniLength2 = 0;
    size_1d_array_type   sizes2(length1, iniLength2);
// Yes, the lengths of the array stripes are arbitrary numbers,
// just for testing but naming each of them doesn't help.
#ifndef __clang_analyzer__
    sizes2[0] = 2;
    sizes2[1] = 1;
    sizes2[2] = 4;
    sizes2[3] = 8;
    sizes2[4] = 5;
    sizes2[5] = 3;
    sizes2[6] = 9;
#endif
    constexpr value_type iniVal4 = 2;
    array_type           test_arr4(sizes2, iniVal4);
    EXPECT_NO_THROW_GMX(
            {
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
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
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
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
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    constexpr size_type  length1A = 2;
    constexpr size_type  length2A = 3;
    constexpr value_type iniValA  = 0;

    constexpr size_type  length1B = 6;
    constexpr size_type  length2B = 5;
    constexpr value_type iniValB  = 2;
    array_type           test_arr1(length1A, length2A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(length1B, length2B, iniValB);
            }
            ) << "Test of resize from regular to regular array dimensions failed with an exception.";

    EXPECT_TRUE(test_arr1.length1() == static_cast<size_type>(length1B) &&
                test_arr1.length2() == static_cast<size_type>(length2B))
    << "Test of array resize failed because actual array dimensions != requested array dimensions";

    constexpr size_type  length1C    = 7;
    constexpr size_type  iniLength2C = 0;
    constexpr value_type iniValC     = 2;
    size_1d_array_type   length2C(length1C, iniLength2C);
#ifndef __clang_analyzer__
    length2C[0] = 2;
    length2C[1] = 1;
    length2C[2] = 4;
    length2C[3] = 8;
    length2C[4] = 5;
    length2C[5] = 3;
    length2C[6] = 9;
#endif
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(length2C, iniValC);
            }
            ) << "Test of resize from regular to irregular array dimensions failed with an exception.";
    EXPECT_EQ(length1C, test_arr1.length1()) << "array length in dimension 1 doesn't match request.";
    for (size_type i = 0; i < test_arr1.length1(); ++i)
    {
        EXPECT_EQ(length2C[i], test_arr1.length2(i))
        << "Test of array resize failed because actual array dimensions != requested array dimensions";
    }

    // vary the lengths in dimension 2
    constexpr value_type     iniValD  = 3;
    const size_1d_array_type length2D {
        5, 9, 1, 4, 8, 5, 3, 2
    };
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(length2D, iniValD);
            }
            ) << "Test of resize from irregular to irregular array dimensions failed with an exception.";
    EXPECT_EQ(length2D.size(), test_arr1.length1()) << "array length in dimension 1 doesn't match request.";
    for (size_type i = 0; i < test_arr1.length1(); ++i)
    {
        EXPECT_EQ(length2D[i], test_arr1.length2(i))
        << "Test of array resize failed because actual array dimensions != requested array dimensions";
    }

    constexpr size_type length1E  = 2;
    constexpr size_type length2E  = 3;
    EXPECT_NO_THROW_GMX(
            {
                test_arr1.initArray(length1E, length2E);
            }
            ) << "Test of resize from irregular to regular array dimensions failed with an exception.";

    EXPECT_TRUE(test_arr1.length1() == length1E &&
                test_arr1.length2() == length2E)
    << "Test of array resize failed because actual array dimensions != requested array dimensions";
}


TYPED_TEST(ParamIrregArray2DTest, 2DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
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

    EXPECT_TRUE(test_arr1.length1() == test_arr2.length1() &&
                test_arr1.length2() == test_arr2.length2())
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

    EXPECT_TRUE(test_arr1.length1() == test_arr3.length1() &&
                test_arr1.length2() == test_arr3.length2())
    << "Test of copy assignment for a regular array via operator() failed because of non-matching array dimensions";

    EXPECT_TRUE(test_arr3(0, 0) == test_arr1(0, 0) &&
                test_arr3(0, 1) == test_arr1(0, 1) &&
                test_arr3(0, 2) == test_arr1(0, 2) &&
                test_arr3(1, 0) == test_arr1(1, 0) &&
                test_arr3(1, 1) == test_arr1(1, 1) &&
                test_arr3(1, 2) == test_arr1(1, 2))
    << "Test of copy assignment for a regular array failed.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 2
    const value_type           iniVal1 = 4;
    const size_1d_array_type   length2 = {2, 1, 4, 8, 5, 3, 9};
    array_type                 test_arr4(length2, iniVal1);
    EXPECT_NO_THROW_GMX(
            {
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(j); ++j)
                    {
                        test_arr4(i, j) = static_cast<value_type>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    array_type test_arr5(test_arr4);
    this->debug << "test_arr5: " << endl << test_arr5 << endl;

    bool       matchedAll = true;
    for (size_type i = 0; i < test_arr4.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr4.length2(i); ++j)
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
    for (size_type i = 0; i < test_arr4.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr4.length2(i); ++j)
        {
            if (test_arr6(i, j) != test_arr4(i, j))
            {
                matchedAll = false;
                break;
            }
        }
    }

    EXPECT_TRUE(matchedAll) << "Test of copy assignement for an irregular array failed.";

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

    const void * const saved_address1 = test_arr1.data();

    array_type         test_arr2(std::move(test_arr1));

    EXPECT_TRUE(saved_address1 == test_arr2.data() && test_arr1.data() == nullptr)
    << "Array was not actually moved";

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

    EXPECT_TRUE(saved_address1 == test_arr1.data() && test_arr2.data() == nullptr)
    << "Array was not actually moved";

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
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    constexpr size_type  length1A = 2;
    constexpr size_type  length2A = 3;
    constexpr size_type  length3A = 2;
    constexpr value_type iniValA  = 0;
    array_type           test_arr1(length1A, length2A, length3A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                size_type v = 3;
                for (size_type i = 0; i < length1A; ++i)
                {
                    for (size_type j = 0; j < length2A; ++j)
                    {
                        for (size_type k = 0; k < length3A; ++k)
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
        for (size_type i = 0; i < length1A; ++i)
        {
            for (size_type j = 0; j < length2A; ++j)
            {
                for (size_type k = 0; k < length3A; ++k)
                {
                    v++;
                    EXPECT_EQ(test_arr1(i, j, k),  static_cast<value_type>(v))
                    << "Test of assignment to a regular array and retrieval via operator() failed.";
                }
            }
        }
    }

    EXPECT_NO_THROW_GMX(
            {
                size_type v = 0;
                for (size_type i = 0; i < length1A; ++i)
                {
                    for (size_type j = 0; j < length2A; ++j)
                    {
                        for (size_type k = 0; k < length3A; ++k)
                        {
                            v++;
                            test_arr1(i, j, k) = static_cast<value_type>(v);
                        }
                    }
                }
            }
            ) << "Test of assignment to a regular array via operator() failed with an exception.";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    {
        size_type v = 0;
        for (size_type i = 0; i < length1A; ++i)
        {
            for (size_type j = 0; j < length2A; ++j)
            {
                for (size_type k = 0; k < length3A; ++k)
                {
                    v++;
                    EXPECT_EQ(test_arr1(i, j, k),  static_cast<value_type>(v))
                    << "Test of assignment to a regular array and retrieval via operator() failed.";

                    EXPECT_EQ(test_arr1(i, j, k), test_arr1[i][j][k])
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
                for (size_type i = 0; i < test_arr3.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr3.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr3.length3(i, j); ++k)
                        {
                            test_arr3(i, j, k) = static_cast<value_type>(l);
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    size_1d_array_type   length2D {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr size_type      len3D = 2;
    const size_2d_array_type length3D(length2D, len3D);
    constexpr value_type     iniValD = 99;

    array_type               test_arr4(length3D, iniValD);
    EXPECT_NO_THROW_GMX(
            {
                size_type l = 1;
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
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
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
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
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    constexpr size_type  length1A = 2;
    constexpr size_type  length2A = 4;
    constexpr size_type  length3A = 7;
    constexpr value_type iniValA  = 0;
    array_type           test_arr1(length1A, length2A, length3A, iniValA);
    EXPECT_NO_THROW_GMX(
            {
                size_type v = 0;
                for (size_type i = 0; i < length1A; ++i)
                {
                    for (size_type j = 0; j < length2A; ++j)
                    {
                        for (size_type k = 0; k < length3A; ++k)
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
    for (size_type i = 0; i < test_arr2.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr2.length2(i); ++j)
        {
            for (size_type k = 0; k < test_arr2.length3(i, j); ++k)
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
    const size_1d_array_type length2B {
        2, 1, 4, 8, 5, 3, 9
    };
    constexpr size_type  length3B = 7;
    constexpr value_type iniValB  = 99;
    size_2d_array_type   sizes2B(length2B, length3B);

    array_type           test_arr3(sizes2B, iniValB);
    EXPECT_NO_THROW_GMX(
            {
                size_type l = 1;
                for (size_type i = 0; i < test_arr3.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr3.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr3.length3(i, j); ++k)
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
    for (size_type i = 0; i < test_arr4.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr4.length2(i); ++j)
        {
            for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
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

    const void * const saved_address1 = test_arr1.data();

    array_type         test_arr2(std::move(test_arr1));

    EXPECT_TRUE(saved_address1 == test_arr2.data() && test_arr1.data() == nullptr)
    << "Array was not actually moved";

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

    EXPECT_TRUE(saved_address1 == test_arr1.data() && test_arr2.data() == nullptr)
    << "Array was not actually moved";

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
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;

    constexpr size_type  length1A = 2;
    constexpr size_type  length2A = 3;
    constexpr size_type  length3A = 2;
    constexpr size_type  length4A = 7;
    constexpr value_type iniValA  = 0;
    array_type           test_arr1(length1A, length2A, length3A, length4A, iniValA);
    EXPECT_NO_THROW_GMX (
            {
                size_type v = 0;
                for (size_type i = 0; i < test_arr1.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr1.length2(); ++j)
                    {
                        for (size_type k = 0; k < test_arr1.length3(); ++k)
                        {
                            for (size_type l = 0; l < test_arr1.length4(); ++l)
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
        for (size_type i = 0; i < test_arr1.length1(); ++i)
        {
            for (size_type j = 0; j < test_arr1.length2(); ++j)
            {
                for (size_type k = 0; k < test_arr1.length3(); ++k)
                {
                    for (size_type l = 0; l < test_arr1.length4(); ++l)
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

    constexpr size_type  length1B = 2;
    constexpr size_type  length2B = 3;
    constexpr size_type  length3B = 1;
    constexpr size_type  length4B = 1;
    constexpr value_type iniValB  = 0;
    array_type           test_arr2(length1B, length2B, length3B, length4B, iniValB);
    EXPECT_NO_THROW_GMX (
            {
                size_type v = 3;
                for (size_type i = 0; i < test_arr2.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr2.length2(); ++j)
                    {
                        for (size_type k = 0; k < test_arr2.length3(); ++k)
                        {
                            for (size_type l = 0; l < test_arr2.length4(); ++l)
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
        for (size_type i = 0; i < test_arr2.length1(); ++i)
        {
            for (size_type j = 0; j < test_arr2.length2(); ++j)
            {
                for (size_type k = 0; k < test_arr2.length3(); ++k)
                {
                    for (size_type l = 0; l < test_arr2.length4(); ++l)
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

    constexpr size_type      l3C      = 3;
    constexpr size_type      l4C      = 1;
    constexpr value_type     iniValC  = 99;
    const size_1d_array_type length2C {
        2, 1, 4, 8, 5, 3, 9
    };
    const size_2d_array_type length3C(length2C, l3C);
    const size_3d_array_type length4C(length3C, l4C);

    array_type               test_arr3(length4C, iniValC);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (size_type i = 0; i < test_arr3.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr3.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr3.length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < test_arr3.length4(i, j, k); ++l)
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

    array_type test_arr4(length4C, iniValC);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < test_arr4.length4(i, j, k); ++l)
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

    bool matchedAll = true;

    EXPECT_NO_THROW_GMX(
            {
                for (size_type i = 0; i < test_arr4.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr4.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < test_arr4.length4(i, j, k); ++l)
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
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;

    const array_type test_arr1
    {
        {
            {
                {
                    1
                }
            },
            {
                {
                    2
                }
            },
            {
                {
                    3
                }
            },
        },
        {
            {
                {
                    4
                }
            },
            {
                {
                    5
                }
            },
            {
                {
                    6
                }
            },
        }
    };

    array_type     test_arr2(test_arr1);

    bool           matchedAll = true;
    for (size_type i = 0; i < test_arr1.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr1.length2(i); ++j)
        {
            for (size_type k = 0; k < test_arr1.length3(i, j); ++k)
            {
                for (size_type l = 0; l < test_arr1.length4(i, j, k); ++l)
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
    const size_1d_array_type length2 {
        2, 1, 4, 8, 5, 3, 9
    };
    const size_t       l3 = 3;
    size_2d_array_type length3(length2, l3);
    const size_t       l4 = 1;
    size_3d_array_type length4(length3, l4);
    const value_type   iniVal = 99;
    array_type         test_arr3(length4, iniVal);
    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (size_type i = 0; i < test_arr3.length1(); ++i)
                {
                    for (size_type j = 0; j < test_arr3.length2(i); ++j)
                    {
                        for (size_type k = 0; k < test_arr3.length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < test_arr3.length4(i, j, k); ++l)
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
    for (size_type i = 0; i < test_arr3.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr3.length2(i); ++j)
        {
            for (size_type k = 0; k < test_arr3.length3(i, j); ++k)
            {
                for (size_type l = 0; l < test_arr3.length4(i, j, k); ++l)
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
    for (size_type i = 0; i < test_arr4.length1(); ++i)
    {
        for (size_type j = 0; j < test_arr4.length2(i); ++j)
        {
            for (size_type k = 0; k < test_arr4.length3(i, j); ++k)
            {
                for (size_type l = 0; l < test_arr4.length4(i, j, k); ++l)
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

    const void * const saved_address1 = test_arr1.data();

    array_type         test_arr2(std::move(test_arr1));

    EXPECT_TRUE(saved_address1 == test_arr2.data() && test_arr1.data() == nullptr)
    << "Array was not actually moved";

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

    EXPECT_TRUE(saved_address1 == test_arr1.data() && test_arr2.data() == nullptr)
    << "Array was not actually moved";

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
