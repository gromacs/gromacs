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

#define GMX_LAMBDA_SITE_IS_DEVEL

namespace gmx
{

namespace data_struct_test
{

namespace
{

using std::endl;

/*********************************************************************/

/*! \brief  output stream operator to insert a string representation of a
            single-layer std::initializer_list object into an output stream, e.g., cout

    \tparam      T        data type stored in the initializer list

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   list     the object whose contents are to be inserted   */
template<typename T>
std::ostream &operator<<(std::ostream &output, const std::initializer_list<T> &list)
{
    // safe the current ostream format for restoring it after output insertion
    std::ios  state(NULL);
    state.copyfmt(output);

    for (size_t i = 0; i < list.size(); ++i)
    {
        output << std::setw(10) << std::fixed << std::setprecision(2) << *(list.begin() + static_cast<std::ptrdiff_t>(i)) << " ";
    }

    // restore the original ostream format
    output.copyfmt(state);

    return output;
}

/*! \brief  output stream operator to insert a string representation of a
            two-layer std::initializer_list object into an output stream, e.g., cout

    \tparam      T        data type stored in the initializer list

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   list     the object whose contents are to be inserted   */
template<typename T>
std::ostream &operator<<(std::ostream &output, const std::initializer_list<std::initializer_list<T> > &list)
{
    // safe the current ostream format for restoring it after output insertion
    std::ios  state(NULL);
    state.copyfmt(output);

    for (size_t j = 0; j < list.size(); ++j)
    {
        const std::initializer_list<T> * const jPtr = list.begin() + static_cast<std::ptrdiff_t>(j);
        for (size_t i = 0; i < jPtr->size(); ++i)
        {
            const T * const iPtr = jPtr->begin() + static_cast<std::ptrdiff_t>(i);
            output << std::setw(10) << std::fixed << std::setprecision(2) << (*iPtr) << " ";
        }
        output << std::endl;
    }

    // restore the original ostream format
    output.copyfmt(state);

    return output;
}

/*! \brief  output stream operator to insert a string representation of a
            three-layer std::initializer_list object into an output stream, e.g., cout

    \tparam      T        data type stored in the initializer list

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   list     the object whose contents are to be inserted   */
template<typename T>
std::ostream &operator<<(std::ostream &output, const std::initializer_list<std::initializer_list<std::initializer_list<T> > > &list)
{
    // safe the current ostream format for restoring it after output insertion
    std::ios  state(NULL);
    state.copyfmt(output);

    for (size_t k = 0; k < list.size(); ++k)
    {
        output << "------------------------------" << std::endl << k << std::endl;
        const std::initializer_list<std::initializer_list<T> > * const kPtr = list.begin() + static_cast<std::ptrdiff_t>(k);
        for (size_t j = 0; j < kPtr->size(); ++j)
        {
            const std::initializer_list<T> * const jPtr = kPtr->begin() + static_cast<std::ptrdiff_t>(j);
            for (size_t i = 0; i < jPtr->size(); ++i)
            {
                const T * const iPtr = jPtr->begin() + static_cast<std::ptrdiff_t>(i);
                output << std::setw(10) << std::fixed << std::setprecision(2) << (*iPtr) << " ";
            }
            output << std::endl;
        }
    }

    // restore the original ostream format
    output.copyfmt(state);

    return output;
}

/*! \brief  output stream operator to insert a string representation of a
            four-layer std::initializer_list object into an output stream, e.g., cout

    \tparam      T        data type stored in the initializer list

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   list     the object whose contents are to be inserted   */
template<typename T>
std::ostream &operator<<(std::ostream &output, const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > &list)
{
    // safe the current ostream format for restoring it after output insertion
    std::ios  state(NULL);
    state.copyfmt(output);

    for (size_t l = 0; l < list.size(); ++l)
    {
        const std::initializer_list<std::initializer_list<std::initializer_list<T> > >* const lPtr = list.begin() + static_cast<std::ptrdiff_t>(l);
        for (size_t k = 0; k < lPtr->size(); ++k)
        {
            output << "------------------------------" << std::endl << l << ", " << k << std::endl;
            const std::initializer_list<std::initializer_list<T> > * const kPtr = lPtr->begin() + static_cast<std::ptrdiff_t>(k);
            for (size_t j = 0; j < kPtr->size(); ++j)
            {
                const std::initializer_list<T> * const jPtr = kPtr->begin() + static_cast<std::ptrdiff_t>(j);
                for (size_t i = 0; i < jPtr->size(); ++i)
                {
                    const T * const iPtr = jPtr->begin() + static_cast<std::ptrdiff_t>(i);
                    output << std::setw(10) << std::fixed << std::setprecision(2) << (*iPtr) << " ";
                }
                output << std::endl;
            }
        }
    }

    // restore the original ostream format
    output.copyfmt(state);

    return output;
}

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
        std::initializer_list<value_type> iList = {0, 1, 2, 3, 4, 5};
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
    EXPECT_FALSE(caughtException) << "Construction destruction of \"IrregArray1D\" failed with exception";

}


TEST_F(IrregArrayTest, DataStructuresIrregArray1DUsage)
{

    gmx::IrregArray1D<unsigned int> test_arr1(1, 6, 0u);
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

    gmx::IrregArray1D<unsigned int> test_arr2(1, 6, 0u);
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

TEST_F(IrregArrayTest, CopyIrregArray1D)
{
    gmx::IrregArray1D<unsigned int> test_arr1(1, 6, 0u);
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

    gmx::IrregArray1D<unsigned int> test_arr2(test_arr1);

    EXPECT_TRUE(test_arr2(1) == 1 &&
                test_arr2(2) == 2 &&
                test_arr2(3) == 3 &&
                test_arr2(4) == 4 &&
                test_arr2(5) == 5 &&
                test_arr2(6) == 6)
    << "Test of copy constructor failed.";

    debug << "test_arr1: " << endl << test_arr1 << endl;
    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_NO_THROW_GMX(
            {
                gmx::IrregArray1D<unsigned int>  test_arr3;
                gmx::IrregArray1D<unsigned int>  test_arr4(test_arr3);
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
}

TYPED_TEST(ParamIrregArray2DTest, 2DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
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

    size_1d_array_type sizes((index_type)0, (index_type)6, (size_type)0);
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 4;
    sizes[3] = 8;
    sizes[4] = 5;
    sizes[5] = 3;
    sizes[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    array_type test_arr3((index_type)0, (index_type)6, sizes, (value_type)99);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::IrregArray2D<unsigned int>::index_type i = test_arr3.getFirst1(); i <= test_arr3.getLast1(); ++i)
                {
                    for (gmx::IrregArray2D<unsigned int>::index_type j = test_arr3.getFirst2(); j <= test_arr3.getLast2(j); ++j)
                    {
                        test_arr3[i][j] = static_cast<gmx::IrregArray2D<unsigned int>::value_type>(j);
                    }
                }
            }
            ) << "Test of assignment an irregular array via operator[] failed with an exception.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    size_1d_array_type sizes2((index_type)1, (index_type)7, (size_type)0);
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
    array_type test_arr4(1, 7, sizes2, (value_type)2);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::IrregArray2D<float>::index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (gmx::IrregArray2D<float>::index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(j); ++j)
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
                for (ssize_t i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (ssize_t j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
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


TYPED_TEST(ParamIrregArray2DTest, 2DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    array_type test_arr1(1, 2, 1, 3, (value_type)0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1[1][1] = 1;
                test_arr1[1][2] = 2;
                test_arr1[1][3] = 3;
                test_arr1[2][1] = 4;
                test_arr1[2][2] = 5;
                test_arr1[2][3] = 6;
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    array_type test_arr2(test_arr1);

    EXPECT_TRUE(test_arr1.getLength1() == test_arr2.getLength1() &&
                test_arr1.getLength2() == test_arr2.getLength2())
    << "Test of copy constructor for a regular array via operator() failed because of non-matching array dimensions";

    EXPECT_TRUE(test_arr2(1, 1) == 1 &&
                test_arr2(1, 2) == 2 &&
                test_arr2(1, 3) == 3 &&
                test_arr2(2, 1) == 4 &&
                test_arr2(2, 2) == 5 &&
                test_arr2(2, 3) == 6)
    << "Test of copy constructor for a regular array via operator() failed because of non-matching array content";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    array_type test_arr3 = test_arr1;

    EXPECT_TRUE(test_arr3(1, 1) == 1 &&
                test_arr3(1, 2) == 2 &&
                test_arr3(1, 3) == 3 &&
                test_arr3(2, 1) == 4 &&
                test_arr3(2, 2) == 5 &&
                test_arr3(2, 3) == 6)
    << "Test of assignment operator= for a regular array failed.";

    this->debug << "test_arr3: " << endl << test_arr3 << endl;


    // let the array start with index 1, and initialize it with value 1
    size_1d_array_type sizes((index_type)1, (index_type)7, (size_type)0);
    sizes[1] = 2;
    sizes[2] = 1;
    sizes[3] = 4;
    sizes[4] = 8;
    sizes[5] = 5;
    sizes[6] = 3;
    sizes[7] = 9;
    array_type test_arr4(1, 7, sizes, 2u);
    EXPECT_NO_THROW_GMX(
            {
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(j); ++j)
                    {
                        test_arr4[i][j] = static_cast<float>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    array_type               test_arr5(test_arr4);

    bool                     matchedAll = true;
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

    EXPECT_NO_THROW_GMX(
            {
                gmx::IrregArray2D<unsigned int>  test_arr5;
                gmx::IrregArray2D<unsigned int>  test_arr6(test_arr5);
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
}


TYPED_TEST(ParamIrregArray3DTest, 3DArrayUsage)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    array_type test_arr1(1, 2, 1, 3, 1, 2, (value_type)0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1[1][1][1] = 1;
                test_arr1[1][2][1] = 2;
                test_arr1[1][3][1] = 3;
                test_arr1[2][1][1] = 4;
                test_arr1[2][2][1] = 5;
                test_arr1[2][3][1] = 6;
                test_arr1[1][1][2] = 7;
                test_arr1[1][2][2] = 8;
                test_arr1[1][3][2] = 9;
                test_arr1[2][1][2] = 10;
                test_arr1[2][2][2] = 11;
                test_arr1[2][3][2] = 12;
            }
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    array_type test_arr2(1, 2, 1, 3, 1, 1, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr2(1, 1, 1) = 1;
                test_arr2(1, 2, 1) = 2;
                test_arr2(1, 3, 1) = 3;
                test_arr2(2, 1, 1) = 4;
                test_arr2(2, 2, 1) = 5;
                test_arr2(2, 3, 1) = 6;
            }
            ) << "Test of assignment to a regular array via operator() failed with an exception.";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr2(1, 1, 1) == 1 &&
                test_arr2(1, 2, 1) == 2 &&
                test_arr2(1, 3, 1) == 3 &&
                test_arr2(2, 1, 1) == 4 &&
                test_arr2(2, 2, 1) == 5 &&
                test_arr2(2, 3, 1) == 6)
    << "Test of assignment to a regular array and retrieval via operator() failed.";

    EXPECT_TRUE(test_arr1(1, 1, 1) == test_arr1[1][1][1] &&
                test_arr1(1, 2, 1) == test_arr1[1][2][1] &&
                test_arr1(1, 3, 1) == test_arr1[1][3][1] &&
                test_arr1(2, 1, 1) == test_arr1[2][1][1] &&
                test_arr1(2, 2, 1) == test_arr1[2][2][1] &&
                test_arr1(2, 3, 1) == test_arr1[2][3][1])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    size_1d_array_type sizes1((index_type)0, (index_type)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    size_2d_array_type sizes2(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 7u);

    array_type         test_arr3((index_type)0, (index_type)6, sizes1, sizes2, 99u);
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

    this->debug << "test_arr3: " << endl << test_arr3 << endl;

    sizes1.initArray((index_type)1, (index_type)7, (size_type)0);
    sizes1[1] = 2;
    sizes1[2] = 1;
    sizes1[3] = 4;
    sizes1[4] = 8;
    sizes1[5] = 5;
    sizes1[6] = 3;
    sizes1[7] = 9;
    // let the array start with index 0, and initialize it with index 1
    sizes2.initArray(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 2u);

    array_type test_arr4((ssize_t)1, (ssize_t)7, sizes1, sizes2, 99u);
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

    EXPECT_NO_THROW_GMX(
            {
                gmx::IrregArray3D<unsigned int>  test_arr5;
                gmx::IrregArray3D<unsigned int>  test_arr6(test_arr5);
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
}

TYPED_TEST(ParamIrregArray3DTest, 3DArrayCopy)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::value_type                  value_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;

    array_type test_arr1(1, 2, 1, 3, 1, 2, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr1[1][1][1] = 1;
                test_arr1[1][2][1] = 2;
                test_arr1[1][3][1] = 3;
                test_arr1[2][1][1] = 4;
                test_arr1[2][2][1] = 5;
                test_arr1[2][3][1] = 6;
                test_arr1[1][1][2] = 7;
                test_arr1[1][2][2] = 8;
                test_arr1[1][3][2] = 9;
                test_arr1[2][1][2] = 10;
                test_arr1[2][2][2] = 11;
                test_arr1[2][3][2] = 12;
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


    size_1d_array_type sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    size_2d_array_type sizes2(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 7u);

    array_type         test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, 99u);
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

    EXPECT_NO_THROW_GMX(
            {
                array_type  test_arr5;
                array_type  test_arr6(test_arr5);
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
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

    bool       caught_exception = false;

    array_type test_arr1(1, 2, 1, 3, 1, 2, 1, 7, 0);
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
            ) << "Test of assignment to a regular array via operator[] failed with an exception.";

    this->debug << "test_arr1: " << endl << test_arr1 << endl;

    array_type test_arr2(1, 2, 1, 3, 1, 1, 1, 1, 0);
    EXPECT_NO_THROW_GMX(
            {
                test_arr2(1, 1, 1, 1) = 1;
                test_arr2(1, 2, 1, 1) = 2;
                test_arr2(1, 3, 1, 1) = 3;
                test_arr2(2, 1, 1, 1) = 4;
                test_arr2(2, 2, 1, 1) = 5;
                test_arr2(2, 3, 1, 1) = 6;
            }
            ) << "Test of assignment to a regular array via operator() failed with an exception.";

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1, 1, 1, 1) == 1 &&
                test_arr2(1, 2, 1, 1) == 2 &&
                test_arr2(1, 3, 1, 1) == 3 &&
                test_arr2(2, 1, 1, 1) == 4 &&
                test_arr2(2, 2, 1, 1) == 5 &&
                test_arr2(2, 3, 1, 1) == 6)
    << "Test of assignment and retrieval for a regular array via operator() failed.";

    this->debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr1(1, 1, 1, 1) == test_arr1[1][1][1][1] &&
                test_arr1(1, 2, 1, 1) == test_arr1[1][2][1][1] &&
                test_arr1(1, 3, 1, 1) == test_arr1[1][3][1][1] &&
                test_arr1(2, 1, 1, 1) == test_arr1[2][1][1][1] &&
                test_arr1(2, 2, 1, 1) == test_arr1[2][2][1][1] &&
                test_arr1(2, 3, 1, 1) == test_arr1[2][3][1][1])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    size_1d_array_type sizes1((index_type)0, (index_type)6, (size_type)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    size_2d_array_type sizes2(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 3u);
    size_3d_array_type sizes3(sizes1.getBegin1(), sizes1.getLast1(), sizes1, sizes2, 1u);

    array_type         test_arr3((index_type)0, (index_type)6, sizes1, sizes2, sizes3, 99u);
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

    const index_type s1 = 1;
    sizes1.initArray(s1, s1 + (index_type)6, (index_type)0);
    sizes1[s1 + 0] = 2;
    sizes1[s1 + 1] = 1;
    sizes1[s1 + 2] = 4;
    sizes1[s1 + 3] = 8;
    sizes1[s1 + 4] = 5;
    sizes1[s1 + 5] = 3;
    sizes1[s1 + 6] = 9;
    array_type test_arr4(s1, s1 + (index_type)6, sizes1, s1, s1 + (index_type)6, sizes1, 5);
    this->debug << "test_arr4: " << endl << test_arr4 << endl;

    EXPECT_NO_THROW_GMX(
            {
                size_type m = 1;
                for (index_type i = test_arr4.getFirst1(); i <= test_arr4.getLast1(); ++i)
                {
                    for (index_type j = test_arr4.getFirst2(); j <= test_arr4.getLast2(i); ++j)
                    {
                        for (index_type k = test_arr4.getFirst3(); k <= test_arr4.getLast3(i, j); ++k)
                        {
                            for (index_type l = test_arr4.getFirst4(); l <= test_arr4.getLast4(i, j, k); ++l)
                            {
                                test_arr4(i, j, k, l) = static_cast<value_type>(m);
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


    size_1d_array_type sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    size_2d_array_type sizes2(sizes1.getBegin1(), sizes1.getLast1(), sizes1, 3u);
    size_3d_array_type sizes3(sizes1.getBegin1(), sizes1.getLast1(), sizes1, sizes2, 1u);

    array_type         test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);
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

    EXPECT_NO_THROW_GMX(
            {
                array_type test_arr6;
                array_type test_arr7(test_arr6);
            }
            ) << "Caught exception while copy constructing an empty array." << endl;
}


} // namespace

} // namespace data_struct_test

} // namespace gmx
