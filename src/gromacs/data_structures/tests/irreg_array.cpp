/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015 by the GROMACS development team, led by
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

   \ingroup module_data_structures
 */

#if __cplusplus >= 201103L
 #define DO_TIMING
 #ifdef DO_TIMING
  // for timing
  #include <chrono>
 #endif
#endif
#include <stdexcept>
#include <iostream> // ::std::cout
#include <ios>      // ::std::left, ::std::rightusing std::cout;
#include <ostream>  // ::std::ostream (<<)
#include <iomanip>  // ::std::setw, ::std::setfill
#include <fstream>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

// define if boost::serialization is also to be compiled in and tested, add corresponding compiler flags -DHAVE_BOOST_SERIALIZATION and for linking
//#define HAVE_BOOST_SERIALIZATION
#include "gromacs/data_structures/IrregArray4D.h"
#include "gromacs/data_structures/FlatIrregArray4D.h"

#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

//#define DEBUG

namespace
{

using namespace std;

class IrregArrayTest : public ::testing::Test
{
    public:
        IrregArrayTest()
        {
        }
        gmx::test::TestFileManager      fileManager_;
};

//! dummy stream as sink for unwanted debugging output, this is guaranteed to work/give a no-effect iostream by the C++ standard according to a post on stackoverflow.com
#if __cplusplus >= 201103L
std::ostream  cnull(nullptr);
#else
std::ostream  cnull(0);
#endif
//std::wostream wcnull(nullptr);
#ifndef DEBUG
//! By default, the debugpt points to the dummy stream cnull silencing debugging information. If macro DEBUG is defined, the debugpt points to cout (stdout) and debugging information will be printed.
std::ostream* debugpt = &cnull;
#else
//! By default, the debugpt points to the dummy stream cnull silencing debugging information. If macro DEBUG is defined, the debugpt points to cout (stdout) and debugging information will be printed.
std::ostream* debugpt = &cout;
#endif
//! used for stearing debug output to cout (stdout) or the dummy stream cnull
#define debug (*debugpt)

/*********************************************************************/

TEST_F(IrregArrayTest, ConstructDestruct)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test 1: creation and destruction FMM interface-related objects"          << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test creation and destruction of data structures:"       << endl;

    bool caught_exception = false;
    try
    {
        gmx::IrregArray1D<double> test_d1_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"IrregArray1D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"IrregArray1D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray1D\" failed";

    caught_exception = false;
    try
    {
        gmx::IrregArray2D<double> test_d2_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"IrregArray2D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"IrregArray2D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray2D\" failed";

    caught_exception = false;
    try
    {
        gmx::IrregArray3D<double> test_d3_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"IrregArray3D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"IrregArray3D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray3D\" failed";

    caught_exception = false;
    try
    {
        gmx::IrregArray4D<double> test_d4_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"IrregArray4D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"IrregArray4D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray4D\" failed";

    caught_exception = false;
    try
    {
        gmx::FlatIrregArray2D<double> test_d2_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"FlatIrregArray2D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"FlatIrregArray2D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"FlatIrregArray2D\" failed";

    caught_exception = false;
    try
    {
        gmx::FlatIrregArray3D<double> test_d3_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"FlatIrregArray3D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"FlatIrregArray3D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"FlatIrregArray3D\" failed";

    caught_exception = false;
    try
    {
        gmx::FlatIrregArray4D<double> test_d4_irreg_array;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception in constructing/destructing a \"FlatIrregArray4D\" object." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception in constructing/destructing a \"FlatIrregArray4D\" object." << endl;
    }
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"FlatIrregArray4D\" failed";

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresIrregArray1DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure IrregArray1D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::IrregArray1D<unsigned int> test_arr(1,6,0u);
    try
    {
        test_arr[1] = 1;
        test_arr[2] = 2;
        test_arr[3] = 3;
        test_arr[4] = 4;
        test_arr[5] = 5;
        test_arr[6] = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::IrregArray1D<unsigned int> test_arr2(1,6,0u);
    try
    {
        test_arr2(1) = 1;
        test_arr2(2) = 2;
        test_arr2(3) = 3;
        test_arr2(4) = 4;
        test_arr2(5) = 5;
        test_arr2(6) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1) == 1 &&
                test_arr2(2) == 2 &&
                test_arr2(3) == 3 &&
                test_arr2(4) == 4 &&
                test_arr2(5) == 5 &&
                test_arr2(6) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1) == test_arr[1] &&
                test_arr(2) == test_arr[2] &&
                test_arr(3) == test_arr[3] &&
                test_arr(4) == test_arr[4] &&
                test_arr(5) == test_arr[5] &&
                test_arr(6) == test_arr[6])
    << "Consistency check for IrregArray2D::operator[] and IrregArray2D::operator() failed.";

    // let the array start with index 1, and initialize it with index 2
    gmx::IrregArray1D<size_t> sizes2((ssize_t)2,(ssize_t)7,(size_t)0);

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresIrregArray2DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure IrregArray2D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::IrregArray2D<unsigned int> test_arr(1,2,1,3,0);
    try
    {
        test_arr[1][1] = 1;
        test_arr[1][2] = 2;
        test_arr[1][3] = 3;
        test_arr[2][1] = 4;
        test_arr[2][2] = 5;
        test_arr[2][3] = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::IrregArray2D<unsigned int> test_arr2(1,2,1,3,0);
    try
    {
        test_arr2(1,1) = 1;
        test_arr2(1,2) = 2;
        test_arr2(1,3) = 3;
        test_arr2(2,1) = 4;
        test_arr2(2,2) = 5;
        test_arr2(2,3) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1) == 1 &&
                test_arr2(1,2) == 2 &&
                test_arr2(1,3) == 3 &&
                test_arr2(2,1) == 4 &&
                test_arr2(2,2) == 5 &&
                test_arr2(2,3) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1) == test_arr[1][1] &&
                test_arr(1,2) == test_arr[1][2] &&
                test_arr(1,3) == test_arr[1][3] &&
                test_arr(2,1) == test_arr[2][1] &&
                test_arr(2,2) == test_arr[2][2] &&
                test_arr(2,3) == test_arr[2][3])
    << "Consistency check for IrregArray2D::operator[] and IrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 4;
    sizes[3] = 8;
    sizes[4] = 5;
    sizes[5] = 3;
    sizes[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::IrregArray2D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes,99u);
    debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    gmx::IrregArray1D<size_t> sizes2((ssize_t)1,(ssize_t)7,(size_t)0);
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
    gmx::IrregArray2D<unsigned int> test_arr4(1, 7, sizes2, 2u);
    debug << "test_arr4: " << endl << test_arr4 << endl;

    gmx::IrregArray2D<double> test_arr5(1, 7, sizes2, 2.5f);
    debug << "test_arr5: " << endl << test_arr5 << endl;

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresIrregArray3DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure IrregArray3D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::IrregArray3D<unsigned int> test_arr(1,2,1,3,1,2,0);
    try
    {
        test_arr[1][1][1] = 1;
        test_arr[1][2][1] = 2;
        test_arr[1][3][1] = 3;
        test_arr[2][1][1] = 4;
        test_arr[2][2][1] = 5;
        test_arr[2][3][1] = 6;
        test_arr[1][1][2] = 7;
        test_arr[1][2][2] = 8;
        test_arr[1][3][2] = 9;
        test_arr[2][1][2] = 10;
        test_arr[2][2][2] = 11;
        test_arr[2][3][2] = 12;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::IrregArray3D<unsigned int> test_arr2(1,2,1,3,1,1,0);
    try
    {
        test_arr2(1,1,1) = 1;
        test_arr2(1,2,1) = 2;
        test_arr2(1,3,1) = 3;
        test_arr2(2,1,1) = 4;
        test_arr2(2,2,1) = 5;
        test_arr2(2,3,1) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1,1) == 1 &&
                test_arr2(1,2,1) == 2 &&
                test_arr2(1,3,1) == 3 &&
                test_arr2(2,1,1) == 4 &&
                test_arr2(2,2,1) == 5 &&
                test_arr2(2,3,1) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1,1) == test_arr[1][1][1] &&
                test_arr(1,2,1) == test_arr[1][2][1] &&
                test_arr(1,3,1) == test_arr[1][3][1] &&
                test_arr(2,1,1) == test_arr[2][1][1] &&
                test_arr(2,2,1) == test_arr[2][2][1] &&
                test_arr(2,3,1) == test_arr[2][3][1])
    << "Consistency check for IrregArray2D::operator[] and IrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    gmx::IrregArray2D<size_t> sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 7u);

    gmx::IrregArray3D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, 99u);
    unsigned int l = 1;
    for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i,j); ++k)
            {
                test_arr3(i, j, k) = l;
                l++;
            }
        }
    }
    debug << "test_arr3: " << endl << test_arr3 << endl;

    sizes1.initArray((ssize_t)1,(ssize_t)7,(size_t)0);
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
    l = 1;
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                test_arr4(i, j, k) = l;
                l++;
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;

    gmx::IrregArray3D<double> test_arr5((ssize_t)1, (ssize_t)7, sizes1, sizes2, 2.5f);
    debug << "test_arr5: " << endl << test_arr5 << endl;

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresIrregArray4DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure IrregArray4D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::IrregArray4D<unsigned int> test_arr(1,2,1,3,1,2,1,7,0);
    try
    {
        test_arr[1][1][1][1] = 1;
        test_arr[1][2][1][1] = 2;
        test_arr[1][3][1][1] = 3;
        test_arr[2][1][1][1] = 4;
        test_arr[2][2][1][1] = 5;
        test_arr[2][3][1][1] = 6;
        test_arr[1][1][2][1] = 7;
        test_arr[1][2][2][1] = 8;
        test_arr[1][3][2][1] = 9;
        test_arr[2][1][2][1] = 10;
        test_arr[2][2][2][1] = 11;
        test_arr[2][3][2][1] = 12;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::IrregArray4D<unsigned int> test_arr2(1,2,1,3,1,1,1,1,0);
    try
    {
        test_arr2(1,1,1,1) = 1;
        test_arr2(1,2,1,1) = 2;
        test_arr2(1,3,1,1) = 3;
        test_arr2(2,1,1,1) = 4;
        test_arr2(2,2,1,1) = 5;
        test_arr2(2,3,1,1) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1,1,1) == 1 &&
                test_arr2(1,2,1,1) == 2 &&
                test_arr2(1,3,1,1) == 3 &&
                test_arr2(2,1,1,1) == 4 &&
                test_arr2(2,2,1,1) == 5 &&
                test_arr2(2,3,1,1) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1,1,1) == test_arr[1][1][1][1] &&
                test_arr(1,2,1,1) == test_arr[1][2][1][1] &&
                test_arr(1,3,1,1) == test_arr[1][3][1][1] &&
                test_arr(2,1,1,1) == test_arr[2][1][1][1] &&
                test_arr(2,2,1,1) == test_arr[2][2][1][1] &&
                test_arr(2,3,1,1) == test_arr[2][3][1][1])
    << "Consistency check for IrregArray2D::operator[] and IrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::IrregArray2D<size_t> sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 3u);
    gmx::IrregArray3D<size_t> sizes3(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, sizes2, 1u);

    gmx::IrregArray4D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);
    unsigned int m = 1;
    for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr3.getBegin4(); l <= test_arr3.getEnd4(i,j,k); ++l)
                {
                    test_arr3(i, j, k, l) = m;
                    m++;
                }
            }
        }
    }
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

    bool all_equal = true;
    m = 1;
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i,j,k); ++l)
                {
                    test_arr4(i, j, k, l) = m;
                    m++;
                }
            }
        }
    }
    gmx::IrregArray4D<size_t> test_arr5 = static_cast<gmx::IrregArray4D<size_t> >(test_arr4);
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i,j,k); ++l)
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
                        all_equal = false;
                    }
                }
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;
    debug << "test_arr5: " << endl << test_arr5 << endl;
    EXPECT_TRUE(all_equal)
    << "Test of explicit conversion operator for FlatIrregArray4D for changing the content data type failed.";

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresFlatIrregArray2DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure FlatIrregArray2D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::FlatIrregArray2D<unsigned int> test_arr(1,2,1,3,0);
    try
    {
        test_arr[1][1] = 1;
        test_arr[1][2] = 2;
        test_arr[1][3] = 3;
        test_arr[2][1] = 4;
        test_arr[2][2] = 5;
        test_arr[2][3] = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::FlatIrregArray2D<unsigned int> test_arr2(1,2,1,3,0);
    try
    {
        test_arr2(1,1) = 1;
        test_arr2(1,2) = 2;
        test_arr2(1,3) = 3;
        test_arr2(2,1) = 4;
        test_arr2(2,2) = 5;
        test_arr2(2,3) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1) == 1 &&
                test_arr2(1,2) == 2 &&
                test_arr2(1,3) == 3 &&
                test_arr2(2,1) == 4 &&
                test_arr2(2,2) == 5 &&
                test_arr2(2,3) == 6)
    << "Test of assignment and retrieval via FlatIrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1) == test_arr[1][1] &&
                test_arr(1,2) == test_arr[1][2] &&
                test_arr(1,3) == test_arr[1][3] &&
                test_arr(2,1) == test_arr[2][1] &&
                test_arr(2,2) == test_arr[2][2] &&
                test_arr(2,3) == test_arr[2][3])
    << "Consistency check for FlatIrregArray2D::operator[] and FlatIrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes[0] = 2;
    sizes[1] = 1;
    sizes[2] = 4;
    sizes[3] = 8;
    sizes[4] = 5;
    sizes[5] = 3;
    sizes[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::FlatIrregArray2D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes, 99u);
    debug << "test_arr3: " << endl << test_arr3 << endl;

    // let the array start with index 1, and initialize it with value 1
    gmx::IrregArray1D<size_t> sizes2((ssize_t)1,(ssize_t)7,(size_t)0);
    sizes2[1] = 2;
    sizes2[2] = 1;
    sizes2[3] = 4;
    sizes2[4] = 8;
    sizes2[5] = 5;
    sizes2[6] = 3;
    sizes2[7] = 9;
    gmx::FlatIrregArray2D<unsigned int> test_arr4(1, 7, sizes2, 2u);
    debug << "test_arr4: " << endl << test_arr4 << endl;

    gmx::FlatIrregArray2D<double> test_arr5(1, 7, sizes2, 2.5f);
    debug << "test_arr5: " << endl << test_arr5 << endl;

    debug << "-----------------------------------------------------<<" << endl;
}

TEST_F(IrregArrayTest, DataStructuresFlatIrregArray3DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure IrregArray3D"        << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::FlatIrregArray3D<unsigned int> test_arr(1,2,1,3,1,2,0);
    try
    {
        test_arr[1][1][1] = 1;
        test_arr[1][2][1] = 2;
        test_arr[1][3][1] = 3;
        test_arr[2][1][1] = 4;
        test_arr[2][2][1] = 5;
        test_arr[2][3][1] = 6;
        test_arr[1][1][2] = 7;
        test_arr[1][2][2] = 8;
        test_arr[1][3][2] = 9;
        test_arr[2][1][2] = 10;
        test_arr[2][2][2] = 11;
        test_arr[2][3][2] = 12;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::FlatIrregArray3D<unsigned int> test_arr2(1,2,1,3,1,1,0);
    try
    {
        test_arr2(1,1,1) = 1;
        test_arr2(1,2,1) = 2;
        test_arr2(1,3,1) = 3;
        test_arr2(2,1,1) = 4;
        test_arr2(2,2,1) = 5;
        test_arr2(2,3,1) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1,1) == 1 &&
                test_arr2(1,2,1) == 2 &&
                test_arr2(1,3,1) == 3 &&
                test_arr2(2,1,1) == 4 &&
                test_arr2(2,2,1) == 5 &&
                test_arr2(2,3,1) == 6)
    << "Test of assignment and retrieval via IrregArray2D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1,1) == test_arr[1][1][1] &&
                test_arr(1,2,1) == test_arr[1][2][1] &&
                test_arr(1,3,1) == test_arr[1][3][1] &&
                test_arr(2,1,1) == test_arr[2][1][1] &&
                test_arr(2,2,1) == test_arr[2][2][1] &&
                test_arr(2,3,1) == test_arr[2][3][1])
    << "Consistency check for IrregArray2D::operator[] and IrregArray2D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 7
    gmx::FlatIrregArray2D<size_t> sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 7u);

    gmx::FlatIrregArray3D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, 99u);
    unsigned int l = 1;
    for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i,j); ++k)
            {
                test_arr3(i, j, k) = l;
                l++;
            }
        }
    }
    debug << "test_arr3: " << endl << test_arr3 << endl;

    sizes1.initArray((ssize_t)1,(ssize_t)7,(size_t)0);
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
    l = 1;
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                test_arr4(i, j, k) = l;
                l++;
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;

    gmx::FlatIrregArray3D<double> test_arr5((ssize_t)1, (ssize_t)7, sizes1, sizes2, 2.5f);
    debug << "test_arr5: " << endl << test_arr5 << endl;

    debug << "-----------------------------------------------------<<" << endl;
}


TEST_F(IrregArrayTest, DataStructuresFlatIrregArray4DUsage)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: usage of FMM interface-related data structure FlatIrregArray4D"    << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    bool caught_exception = false;

    gmx::FlatIrregArray4D<unsigned int> test_arr(1,2,1,3,1,2,1,7,0);
    try
    {
        test_arr[1][1][1][1] = 1;
        test_arr[1][2][1][1] = 2;
        test_arr[1][3][1][1] = 3;
        test_arr[2][1][1][1] = 4;
        test_arr[2][2][1][1] = 5;
        test_arr[2][3][1][1] = 6;
        test_arr[1][1][2][1] = 7;
        test_arr[1][2][2][1] = 8;
        test_arr[1][3][2][1] = 9;
        test_arr[2][1][2][1] = 10;
        test_arr[2][2][2][1] = 11;
        test_arr[2][3][2][1] = 12;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator[] failed with an exception.";

    debug << "test_arr: " << endl << test_arr << endl;

    gmx::FlatIrregArray4D<unsigned int> test_arr2(1,2,1,3,1,1,1,1,0);
    try
    {
        test_arr2(1,1,1,1) = 1;
        test_arr2(1,2,1,1) = 2;
        test_arr2(1,3,1,1) = 3;
        test_arr2(2,1,1,1) = 4;
        test_arr2(2,2,1,1) = 5;
        test_arr2(2,3,1,1) = 6;
    }
    catch (std::exception const &e)
    {
        caught_exception = true;
        cout << "Caught std::exception while assigning values via operator[]." << endl;
        cout << "Exception message: " << e.what() << endl;
    }
    catch (...)
    {
        caught_exception = true;
        cout << "Caught non-standard exception while assigning via operator[]." << endl;
    };

    EXPECT_FALSE(caught_exception)
    << "Test of assignment via operator() failed with an exception.";

    EXPECT_TRUE(test_arr2(1,1,1,1) == 1 &&
                test_arr2(1,2,1,1) == 2 &&
                test_arr2(1,3,1,1) == 3 &&
                test_arr2(2,1,1,1) == 4 &&
                test_arr2(2,2,1,1) == 5 &&
                test_arr2(2,3,1,1) == 6)
    << "Test of assignment and retrieval via FlatIrregArray4D::operator() failed.";

    debug << "test_arr2: " << endl << test_arr2 << endl;

    EXPECT_TRUE(test_arr(1,1,1,1) == test_arr[1][1][1][1] &&
                test_arr(1,2,1,1) == test_arr[1][2][1][1] &&
                test_arr(1,3,1,1) == test_arr[1][3][1][1] &&
                test_arr(2,1,1,1) == test_arr[2][1][1][1] &&
                test_arr(2,2,1,1) == test_arr[2][2][1][1] &&
                test_arr(2,3,1,1) == test_arr[2][3][1][1])
    << "Consistency check for FlatIrregArray4D::operator[] and FlatIrregArray4D::operator() failed.";

    gmx::IrregArray1D<size_t> sizes1((ssize_t)0,(ssize_t)6,(size_t)0);
    sizes1[0] = 2;
    sizes1[1] = 1;
    sizes1[2] = 4;
    sizes1[3] = 8;
    sizes1[4] = 5;
    sizes1[5] = 3;
    sizes1[6] = 9;
    // let the array start with index 0, and initialize it with value 99
    gmx::FlatIrregArray2D<size_t> sizes2(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, 3u);
    gmx::FlatIrregArray3D<size_t> sizes3(sizes1.getBegin1(), sizes1.getEnd1(), sizes1, sizes2, 1u);

    gmx::FlatIrregArray4D<unsigned int> test_arr3((ssize_t)0, (ssize_t)6, sizes1, sizes2, sizes3, 99u);
    size_t m = 1;
    for (ssize_t i = test_arr3.getBegin1(); i <= test_arr3.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr3.getBegin2(); j <= test_arr3.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr3.getBegin3(); k <= test_arr3.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr3.getBegin4(); l <= test_arr3.getEnd4(i,j,k); ++l)
                {
                    test_arr3(i, j, k, l) = m;
                    m++;
                }
            }
        }
    }
    debug << "test_arr3: " << endl << test_arr3 << endl;

    // copy constructor and implicit conversion operator
    // implicitly also tests the lower dimensional arrays
    const ssize_t s1 = 1;
    sizes1.initArray(s1, s1 + (ssize_t)6, (ssize_t)0);
    sizes1[s1 + 0] = 2;
    sizes1[s1 + 1] = 1;
    sizes1[s1 + 2] = 4;
    sizes1[s1 + 3] = 8;
    sizes1[s1 + 4] = 5;
    sizes1[s1 + 5] = 3;
    sizes1[s1 + 6] = 9;
    gmx::FlatIrregArray4D<size_t> test_arr4(s1, s1 + (ssize_t)6, sizes1, s1, s1 + (ssize_t)6, sizes1, 0);

    bool all_equal = true;
    m = 1;
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i,j,k); ++l)
                {
                    test_arr4(i, j, k, l) = m;
                    m++;
                }
            }
        }
    }
    gmx::FlatIrregArray4D<double> test_arr5 = static_cast<gmx::FlatIrregArray4D<double> >(test_arr4);
    for (ssize_t i = test_arr4.getBegin1(); i <= test_arr4.getEnd1(); ++i)
    {
        for (ssize_t j = test_arr4.getBegin2(); j <= test_arr4.getEnd2(i); ++j)
        {
            for (ssize_t k = test_arr4.getBegin3(); k <= test_arr4.getEnd3(i,j); ++k)
            {
                for (ssize_t l = test_arr4.getBegin4(); l <= test_arr4.getEnd4(i,j,k); ++l)
                {
                    if (static_cast<double>(test_arr4(i, j, k, l)) != test_arr5(i, j, k, l))
                    {
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") != "
                              << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ")"
                              << endl;
                        debug << "test_arr4(" << i << ", " << j << ", " << k << ", " << l << ") = "
                              << test_arr4(i, j, k, l) << endl;
                        debug << "test_arr5(" << i << ", " << j << ", " << k << ", " << l << ") = "
                              << test_arr5(i, j, k, l) << endl;
                        all_equal = false;
                    }
                }
            }
        }
    }
    debug << "test_arr4: " << endl << test_arr4 << endl;
    debug << "test_arr5: " << endl << test_arr5 << endl;
    EXPECT_TRUE(all_equal)
    << "Test of explicit conversion operator for FlatIrregArray4D for changing the content data type failed.";

    debug << "-----------------------------------------------------<<" << endl;
}


} // namespace

#ifdef DEBUG
#undef DEBUG
#endif
