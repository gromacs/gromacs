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
   Tests of different access variants of the irregular array data structures

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
#include <vector>

// define if boost::serialization is also to be compiled in and tested, add corresponding compiler flags -DHAVE_BOOST_SERIALIZATION and for linking
//#define HAVE_BOOST_SERIALIZATION
#include "gromacs/data_structures/AlignedAllocator.h"
#include "gromacs/data_structures/FlatIrregArray4D.h"
#include "gromacs/data_structures/IrregArray4D.h"
#include "gromacs/data_structures/IrregArray2D.h"

#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

//#define DEBUG

namespace
{

using namespace std;

class AlignedAllocatorTest : public ::testing::Test
{
    public:
        AlignedAllocatorTest()
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

TEST_F(AlignedAllocatorTest, VectorConstructDestruct)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: AlignedAllocator std::vector, IrregArray test"                     << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    const size_t                    dim1    = 100;
    const size_t                    dim2    = 10;
    const gmx::IrregArray1D<size_t> sizes(1, dim1, dim2);

    bool caught_exception = false;
    try
    {
        gmx::IrregArray1D<double, AlignedAllocator<double, 32u> > test_d1_irreg_array;
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
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray1D\" with alignment 32 failed";

    caught_exception = false;
    try
    {
        gmx::IrregArray2D<double, AlignedAllocator<double, 32u> > test_d2_irreg_array(1, dim1, sizes, 1.0f);
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
    EXPECT_FALSE(caught_exception) << "Construction destruction of \"IrregArray1D\" with alignment 32 failed";

    debug << "-----------------------------------------------------<<" << endl;
}


TEST_F(AlignedAllocatorTest, IrregArray2DAccessOperatorPerformance)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running test: access operators IrregArray2D"                                     << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    const size_t                    ncycles = 1;
    const size_t                    dim1    = 100;
    const size_t                    dim2    = 10;
    const gmx::IrregArray1D<size_t> sizes(1, dim1, dim2);
    //! test with 32 byte alignment as suitable for AVX
    gmx::IrregArray2D<double, AlignedAllocator<double, 32> > test_arr(1, dim1, 1, dim2, 1.0f);
    gmx::IrregArray2D<double, AlignedAllocator<double, 32> > test_arr2(1, dim1, sizes, 1.0f);

    // for size comparison
    debug << "Size of the regular special  case of IrregArray2D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of IrregArray2D " << test_arr2.getSize()
    << " for " << test_arr2.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    << static_cast<double>(test_arr2.getSize() - test_arr2.getNelements() * sizeof(double))
    / static_cast<double>(test_arr2.getNelements()) << endl;

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr.getBegin1(); i <= test_arr.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr.getBegin2(); j <= test_arr.getEnd2(); ++j)
                {
                    sum += test_arr.getArray()[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr.getBegin1(); i <= test_arr.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr.getBegin2(); j <= test_arr.getEnd2(); ++j)
                {
                    sum += test_arr(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr.getBegin1(); i <= test_arr.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr.getBegin2(); j <= test_arr.getEnd2(); ++j)
                {
                    sum += test_arr[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr2.getBegin1(); i <= test_arr2.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr2.getBegin2(); j <= test_arr2.getEnd2(i); ++j)
                {
                    sum += test_arr2.getArray()[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr2.getBegin1(); i <= test_arr2.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr2.getBegin2(); j <= test_arr2.getEnd2(i); ++j)
                {
                    sum += test_arr2(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (ssize_t i = test_arr2.getBegin1(); i <= test_arr2.getEnd1(); ++i)
            {
                for (ssize_t j = test_arr2.getBegin2(); j <= test_arr2.getEnd2(i); ++j)
                {
                    sum += test_arr2[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
//              << chrono::duration_cast<chrono::microseconds>(stop - start).count()
        << chrono::duration<double, std::milli>(stop - start).count() << " ms"
//              << chrono::duration<double,std::micro>(stop - start).count() << " microseconds"
//              << chrono::duration<double,std::nano>(stop - start).count() << " nanoseconds"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

    debug << "-----------------------------------------------------<<" << endl;
}


} // namespace

#ifdef DEBUG
#undef DEBUG
#endif
