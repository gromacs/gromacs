/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#include "gromacs/utility/data_structures/flat_irreg_array_4d.h"
#include "gromacs/utility/data_structures/irreg_array_4d.h"

#include "data_struct_test_commons.h"

namespace
{

using std::endl;


class IrregArrayTest : public ::testing::Test
{
    public:
        IrregArrayTest()
            : debug(nullptr), error(nullptr)
        {
        }
        void SetUp()
        {
            // trigger debug output
            if (gmx::data_struct_test::g_debugBool)
            {
                activateDebugOut();
                activateErrorOut();
            }
        }
        void TearDown()
        {
            deactivateDebugOut();
            deactivateErrorOut();
        }

        void activateDebugOut()
        {
            debug.rdbuf(std::cout.rdbuf());
        }
        void activateErrorOut()
        {
            error.rdbuf(std::cerr.rdbuf());
        }
        void deactivateDebugOut()
        {
            debug.rdbuf(nullptr);
        }
        void deactivateErrorOut()
        {
            error.rdbuf(nullptr);
        }

        gmx::test::TestFileManager   fileManager_;
    public:
        // debugging output stream, silent by default
        std::ostream                 debug;
        // error output stream, silent by default
        std::ostream                 error;
};

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

    debug << "test_arr2: " << endl << test_arr2 << endl;

}

TEST_F(IrregArrayTest, DataStructuresIrregArray2DUsage)
{

    gmx::IrregArray2D<unsigned int> test_arr(1, 2, 1, 3, 0);
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

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::IrregArray2D<unsigned int> test_arr2(1, 2, 1, 3, 0);
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

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1, 1) == test_arr[1][1] &&
                test_arr(1, 2) == test_arr[1][2] &&
                test_arr(1, 3) == test_arr[1][3] &&
                test_arr(2, 1) == test_arr[2][1] &&
                test_arr(2, 2) == test_arr[2][2] &&
                test_arr(2, 3) == test_arr[2][3])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    gmx::IrregArray1D<size_t> sizes((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 4;
    sizes[3] = 8;
    sizes[4] = 5;
    sizes[5] = 3;
    sizes[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::IrregArray2D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes, 99u);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::IrregArray2D<unsigned int>::index_type i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (gmx::IrregArray2D<unsigned int>::index_type j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(j); ++j)
                    {
                        test_arr3[i][j] = static_cast<gmx::IrregArray2D<unsigned int>::value_type>(j);
                    }
                }
            }
            ) << "Test of assignment an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    gmx::IrregArray1D<size_t> sizes2((ssize_t)1, (ssize_t)7, (size_t)0);
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
    gmx::IrregArray2D<float> test_arr4(1, 7, sizes2, 2u);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::IrregArray2D<float>::index_type i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (gmx::IrregArray2D<float>::index_type j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(j); ++j)
                    {
                        test_arr4[i][j] = static_cast<float>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
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

TEST_F(IrregArrayTest, DataStructuresIrregArray3DUsage)
{

    gmx::IrregArray3D<unsigned int> test_arr1(1, 2, 1, 3, 1, 2, 0);
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

    debug << "test_arr1: " << endl << test_arr1 << endl;

    gmx::IrregArray3D<unsigned int> test_arr2(1, 2, 1, 3, 1, 1, 0);
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

    debug << "test_arr2: " << endl << test_arr2 << endl;

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

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    gmx::IrregArray2D<size_t>       sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 7u);

    gmx::IrregArray3D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int l = 1;
                for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i, j); ++k)
                        {
                            test_arr3[i][j][k] = l;
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    sizes1.initArray((ssize_t)1, (ssize_t)7, (size_t)0);
    sizes1[1] = 2;
    sizes1[2] = 1;
    sizes1[3] = 4;
    sizes1[4] = 8;
    sizes1[5] = 5;
    sizes1[6] = 3;
    sizes1[7] = 9;
    // let the array start with index 0, and initialize it with index 1
    sizes2.initArray(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 2u);

    gmx::IrregArray3D<unsigned int> test_arr4((ssize_t)1, (ssize_t)7, sizes1, sizes2, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int l = 1;
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            test_arr4(i, j, k) = l;
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
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

}


TEST_F(IrregArrayTest, DataStructuresIrregArray4DUsage)
{

    bool caught_exception = false;

    gmx::IrregArray4D<unsigned int> test_arr1(1, 2, 1, 3, 1, 2, 1, 7, 0);
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

    debug << "test_arr1: " << endl << test_arr1 << endl;

    gmx::IrregArray4D<unsigned int> test_arr2(1, 2, 1, 3, 1, 1, 1, 1, 0);
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

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr1(1, 1, 1, 1) == test_arr1[1][1][1][1] &&
                test_arr1(1, 2, 1, 1) == test_arr1[1][2][1][1] &&
                test_arr1(1, 3, 1, 1) == test_arr1[1][3][1][1] &&
                test_arr1(2, 1, 1, 1) == test_arr1[2][1][1][1] &&
                test_arr1(2, 2, 1, 1) == test_arr1[2][2][1][1] &&
                test_arr1(2, 3, 1, 1) == test_arr1[2][3][1][1])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::IrregArray2D<size_t>       sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 3u);
    gmx::IrregArray3D<size_t>       sizes3(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, sizes2, 1u);

    gmx::IrregArray4D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int m = 1;
                for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr3.getBegin4(); l <= test_arr3.getEnd4(i, j, k); ++l)
                            {
                                test_arr3[i][j][k][l] = m;
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    const ssize_t s1 = 1;
    sizes1.initArray(s1, s1 + (ssize_t)6, (ssize_t)0);
    sizes1[s1 + 0] = 2;
    sizes1[s1 + 1] = 1;
    sizes1[s1 + 2] = 4;
    sizes1[s1 + 3] = 8;
    sizes1[s1 + 4] = 5;
    sizes1[s1 + 5] = 3;
    sizes1[s1 + 6] = 9;
    gmx::IrregArray4D<double> test_arr4(s1, s1 + (ssize_t)6, sizes1, s1, s1 + (ssize_t)6, sizes1, 2.5);
    debug << "test_arr4: " << endl << test_arr4 << endl;

    EXPECT_NO_THROW_GMX(
            {
                unsigned int m = 1;
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
                            {
                                test_arr4(i, j, k, l) = m;
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    bool matchedAll                     = true;
    gmx::IrregArray4D<size_t> test_arr5 = static_cast<gmx::IrregArray4D<size_t> >(test_arr4);
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
                {
                    if (static_cast<size_t>(test_arr4(i, j, k, l)) != test_arr5(i, j, k, l))
                    {
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") != "
                        << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ")"
                        << endl;
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr4(i, j, k, l) << endl;
                        debug << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr5(i, j, k, l) << endl;
                        matchedAll = false;
                    }
                }
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;
    debug << "test_arr5: " << endl << test_arr5 << endl;
    EXPECT_TRUE(matchedAll)
    << "Test of explicit conversion operator for changing the content data type failed.";

    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
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

TEST_F(IrregArrayTest, DataStructuresFlatIrregArray2DUsage)
{

    gmx::FlatIrregArray2D<unsigned int> test_arr(1, 2, 1, 3, 0);
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

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::FlatIrregArray2D<unsigned int> test_arr2(1, 2, 1, 3, 0);
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
    << "Test of assignment and retrieval via FlatIrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1, 1) == test_arr[1][1] &&
                test_arr(1, 2) == test_arr[1][2] &&
                test_arr(1, 3) == test_arr[1][3] &&
                test_arr(2, 1) == test_arr[2][1] &&
                test_arr(2, 2) == test_arr[2][2] &&
                test_arr(2, 3) == test_arr[2][3])
    << "Consistency check for FlatIrregArray2D::operator[] and FlatIrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 4;
    sizes[3] = 8;
    sizes[4] = 5;
    sizes[5] = 3;
    sizes[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::FlatIrregArray2D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes, 99u);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::FlatIrregArray2D<unsigned int>::index_type i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (gmx::FlatIrregArray2D<unsigned int>::index_type j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(j); ++j)
                    {
                        test_arr3[i][j] = static_cast<gmx::FlatIrregArray2D<unsigned int>::value_type>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    gmx::IrregArray1D<size_t> sizes2((ssize_t)1, (ssize_t)7, (size_t)0);
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
    gmx::FlatIrregArray2D<float> test_arr4(1, 7, sizes2, 2u);
    EXPECT_NO_THROW_GMX(
            {
                for (gmx::FlatIrregArray2D<float>::index_type i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (gmx::FlatIrregArray2D<float>::index_type j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(j); ++j)
                    {
                        test_arr4[i][j] = static_cast<float>(j);
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";
    debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
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
    << "Consistency check for an irregular array for operator[] and operator() failed.";

}

TEST_F(IrregArrayTest, DataStructuresFlatIrregArray3DUsage)
{

    gmx::FlatIrregArray3D<unsigned int> test_arr1(1, 2, 1, 3, 1, 2, 0);
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

    debug << "test_arr1: " << endl << test_arr1 << endl;

    gmx::FlatIrregArray3D<unsigned int> test_arr2(1, 2, 1, 3, 1, 1, 0);
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

    debug << "test_arr2: " << endl << test_arr2 << endl;

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
    << "Consistency check for a regular array for operator[] and operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    gmx::FlatIrregArray2D<size_t>       sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 7u);

    gmx::FlatIrregArray3D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int l = 1;
                for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i, j); ++k)
                        {
                            test_arr3[i][j][k] = l;
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    sizes1.initArray((ssize_t)1, (ssize_t)7, (size_t)0);
    sizes1[1] = 2;
    sizes1[2] = 1;
    sizes1[3] = 4;
    sizes1[4] = 8;
    sizes1[5] = 5;
    sizes1[6] = 3;
    sizes1[7] = 9;
    // let the array start with index 0, and initialize it with index 1
    sizes2.initArray(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 2u);

    gmx::FlatIrregArray3D<unsigned int> test_arr4((ssize_t)1, (ssize_t)7, sizes1, sizes2, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int l = 1;
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            test_arr4(i, j, k) = l;
                            l++;
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    debug << "test_arr4: " << endl << test_arr4 << endl;

    bool matchedAll = true;
    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
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
    << "Consistency check for an irregular regular array for operator[] and operator() failed.";

}


TEST_F(IrregArrayTest, DataStructuresFlatIrregArray4DUsage)
{

    bool caught_exception = false;

    gmx::FlatIrregArray4D<unsigned int> test_arr1(1, 2, 1, 3, 1, 2, 1, 7, 0);
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

    debug << "test_arr1: " << endl << test_arr1 << endl;

    gmx::FlatIrregArray4D<unsigned int> test_arr2(1, 2, 1, 3, 1, 1, 1, 1, 0);
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

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr1(1, 1, 1, 1) == test_arr1[1][1][1][1] &&
                test_arr1(1, 2, 1, 1) == test_arr1[1][2][1][1] &&
                test_arr1(1, 3, 1, 1) == test_arr1[1][3][1][1] &&
                test_arr1(2, 1, 1, 1) == test_arr1[2][1][1][1] &&
                test_arr1(2, 2, 1, 1) == test_arr1[2][2][1][1] &&
                test_arr1(2, 3, 1, 1) == test_arr1[2][3][1][1])
    << "Consistency check for a regular array between accession via operator[] and operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0, (ssize_t)6, (size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::FlatIrregArray2D<size_t>       sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 3u);
    gmx::FlatIrregArray3D<size_t>       sizes3(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, sizes2, 1u);

    gmx::FlatIrregArray4D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);
    EXPECT_NO_THROW_GMX(
            {
                unsigned int m = 1;
                for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr3.getBegin4(); l <= test_arr3.getEnd4(i, j, k); ++l)
                            {
                                test_arr3[i][j][k][l] = m;
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator[] failed with an exception.";

    debug << "test_arr3: " << endl << test_arr3 << endl;

    const ssize_t s1 = 1;
    sizes1.initArray(s1, s1 + (ssize_t)6, (ssize_t)0);
    sizes1[s1 + 0] = 2;
    sizes1[s1 + 1] = 1;
    sizes1[s1 + 2] = 4;
    sizes1[s1 + 3] = 8;
    sizes1[s1 + 4] = 5;
    sizes1[s1 + 5] = 3;
    sizes1[s1 + 6] = 9;
    gmx::FlatIrregArray4D<double> test_arr4(s1, s1 + (ssize_t)6, sizes1, s1, s1 + (ssize_t)6, sizes1, 2.5);
    debug << "test_arr4: " << endl << test_arr4 << endl;

    EXPECT_NO_THROW_GMX(
            {
                unsigned int m = 1;
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
                            {
                                test_arr4(i, j, k, l) = m;
                                m++;
                            }
                        }
                    }
                }
            }
            ) << "Test of assignment to an irregular array via operator() failed with an exception.";

    bool matchedAll = true;
    gmx::FlatIrregArray4D<size_t> test_arr5 = static_cast<gmx::FlatIrregArray4D<size_t> >(test_arr4);
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
                {
                    if (static_cast<size_t>(test_arr4(i, j, k, l)) != test_arr5(i, j, k, l))
                    {
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") != "
                        << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ")"
                        << endl;
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr4(i, j, k, l) << endl;
                        debug << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ") = "
                        << test_arr5(i, j, k, l) << endl;
                        matchedAll = false;
                    }
                }
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;
    debug << "test_arr5: " << endl << test_arr5 << endl;
    EXPECT_TRUE(matchedAll)
    << "Test of explicit conversion operator for changing the content data type failed.";

    EXPECT_NO_THROW_GMX(
            {
                for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
                {
                    for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
                    {
                        for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i, j); ++k)
                        {
                            for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i, j, k); ++l)
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


} // namespace
