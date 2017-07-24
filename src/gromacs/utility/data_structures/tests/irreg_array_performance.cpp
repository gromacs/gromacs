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
   Tests of different access variants of the irregular array data structures

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

#define DO_TIMING
#ifdef DO_TIMING
#include <chrono>
#endif

namespace
{

using std::endl;

class IrregArrayPerformanceTest : public gmx::data_struct_test::DataStructTest
{
};

/*********************************************************************/

// NOTE:
// precomputing the start and end indices of the array stripe in the last dimension is much faster
// but here I want to compare the relative cost of the different access variants
// results: operator[] is nearly twofold slower than operator() -- is there a better implementation?
//          the ratio between book-keeping and storage data is heavily deplastent on the length of the
//           array stripes in the highest dimension,
//          For realistic examples, 3.2 bytes book-keeping data per stored array element was
//           the worst case tested (for a 1000 x 5 x 1000 x 5 IrregArray4D).

TEST_F(IrregArrayPerformanceTest, IrregArray2DAccessOperatorPerformance)
{

    const size_t                    ncycles = 1;
    const size_t                    dim1    = 100;
    const size_t                    dim2    = 10;
    const gmx::IrregArray1D<size_t> sizes(1, dim1, dim2);
    gmx::IrregArray2D<double>       test_arr(1, dim1, 1, dim2, 1.0);
    gmx::IrregArray2D<double>       test_arr2(1, dim1, sizes, 1.0);

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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    sum += test_arr.getArray()[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    sum += test_arr(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    sum += test_arr[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    sum += test_arr2.getArray()[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    sum += test_arr2(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    sum += test_arr2[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}

TEST_F(IrregArrayPerformanceTest, FlatIrregArray2DAccessOperatorPerformance)
{

    const size_t                    ncycles = 1;
    const size_t                    dim1    = 100;
    const size_t                    dim2    = 10;
    const gmx::IrregArray1D<size_t> sizes(1, dim1, dim2);
    gmx::FlatIrregArray2D<double>   test_arr(1, dim1, 1, dim2, 1.0);
    gmx::FlatIrregArray2D<double>   test_arr2(1, dim1, sizes, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of FlatIrregArray2D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of FlatIrregArray2D " << test_arr2.getSize()
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
// with optimized index precomputation twice as fast, but not a fair comparison, because the trick
// could also be used for optimized determination of bounds for the other access variants
//                const size_t jstart = test_arr.getArrayIndex(i, test_arr.getBegin2());
//                const size_t jlast  = jstart + test_arr.getLength2() - 1;
//                for (size_t j = jstart; j <= jlast; ++j)
                for (size_t j = test_arr.getArrayIndex(i, test_arr.getBegin2()); j <= test_arr.getArrayIndex(i, test_arr.getBegin2()) + test_arr.getLength2() - 1; ++j)
                {
                    sum += test_arr.getArray()[j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    sum += test_arr(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    sum += test_arr[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
//                const size_t jstart = test_arr2.getArrayIndex(i, test_arr2.getBegin2());
//                const size_t jlast  = jstart + test_arr2.getLength2(i) - 1;
//                for (size_t j = jstart; j <= jlast; ++j)
                for (size_t j = test_arr.getArrayIndex(i, test_arr.getBegin2()); j <= test_arr.getArrayIndex(i, test_arr.getBegin2()) + test_arr.getLength2(i) - 1; ++j)
                {
                    sum += test_arr2.getArray()[j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    sum += test_arr2(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    sum += test_arr2[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}


TEST_F(IrregArrayPerformanceTest, IrregArray3DAccessOperatorPerformance)
{
    const size_t                    ncycles = 1;
    const size_t                    dim1    = 30;
    const size_t                    dim2    = 5;
    const size_t                    dim3    = 30;
    const gmx::IrregArray1D<size_t> sizes1(1, dim1, dim2);
    const gmx::IrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    gmx::IrregArray3D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1.0);
    gmx::IrregArray3D<double>       test_arr2(1, dim1, sizes1, sizes2, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of IrregArray3D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of IrregArray3D " << test_arr2.getSize()
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        sum += test_arr.getArray()[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug  << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        sum += test_arr(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        sum += test_arr[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        sum += test_arr2.getArray()[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        sum += test_arr2(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        sum += test_arr2[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}

TEST_F(IrregArrayPerformanceTest, FlatIrregArray3DAccessOperatorPerformance)
{

    const size_t                        ncycles = 1;
    const size_t                        dim1    = 30;
    const size_t                        dim2    = 5;
    const size_t                        dim3    = 30;
    const gmx::IrregArray1D<size_t>     sizes1(1, dim1, dim2);
    const gmx::FlatIrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    gmx::FlatIrregArray3D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1.0);
    gmx::FlatIrregArray3D<double>       test_arr2(1, dim1, sizes1, sizes2, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of FlatIrregArray3D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of FlatIrregArray3D " << test_arr2.getSize()
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (size_t k = test_arr.getArrayIndex(i, j, test_arr.getBegin3()); k <= test_arr.getArrayIndex(i, j, test_arr.getBegin3()) + test_arr.getLength3() - 1; ++k)
                    {
                        sum += test_arr.getArray()[k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        sum += test_arr(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        sum += test_arr[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (size_t k = test_arr2.getArrayIndex(i, j, test_arr2.getBegin3()); k <= test_arr2.getArrayIndex(i, j, test_arr2.getBegin3()) + test_arr2.getLength3(i, j) - 1; ++k)
                    {
                        sum += test_arr2.getArray()[k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        sum += test_arr2(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        sum += test_arr2[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}

TEST_F(IrregArrayPerformanceTest, IrregArray4DAccessOperatorPerformance)
{
    const size_t                    ncycles = 1;
    const size_t                    dim1    = 30;
    const size_t                    dim2    = 5;
    const size_t                    dim3    = 30;
    const size_t                    dim4    = 5;
    const gmx::IrregArray1D<size_t> sizes1(1, dim1, dim2);
    const gmx::IrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    const gmx::IrregArray3D<size_t> sizes3(1, dim1, 1, dim2, 1, dim3, dim4);
    gmx::IrregArray4D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1, dim4, 1.0);
    gmx::IrregArray4D<double>       test_arr2(1, dim1, sizes1, sizes2, sizes3, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of IrregArray4D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of IrregArray4D " << test_arr2.getSize()
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (ssize_t l = test_arr.getFirst4(); l <= test_arr.getLast4(); ++l)
                        {
                            sum += test_arr.getArray()[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (ssize_t l = test_arr.getFirst4(); l <= test_arr.getLast4(); ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (ssize_t l = test_arr.getFirst4(); l <= test_arr.getLast4(); ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                        {
                            sum += test_arr2.getArray()[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }
}

TEST_F(IrregArrayPerformanceTest, FlatIrregArray4DAccessOperatorPerformance)
{
    const size_t                        ncycles = 1;
    const size_t                        dim1    = 30;
    const size_t                        dim2    = 5;
    const size_t                        dim3    = 30;
    const size_t                        dim4    = 5;
    const gmx::IrregArray1D<size_t>     sizes1(1, dim1, dim2);
    const gmx::FlatIrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    const gmx::FlatIrregArray3D<size_t> sizes3(1, dim1, 1, dim2, 1, dim3, dim4);
    gmx::FlatIrregArray4D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1, dim4, 1.0);
    gmx::FlatIrregArray4D<double>       test_arr2(1, dim1, sizes1, sizes2, sizes3, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of FlatIrregArray4D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of FlatIrregArray4D " << test_arr2.getSize()
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (size_t l = test_arr.getArrayIndex(i, j, k, test_arr.getBegin4()); l <= test_arr.getArrayIndex(i, j, k, test_arr.getBegin4()) + test_arr.getLength4() - 1; ++l)
                        {
                            sum += test_arr.getArray()[l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (ssize_t l = test_arr.getFirst4(); l <= test_arr.getLast4(); ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug  << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr.getFirst1(); i <= test_arr.getLast1(); ++i)
            {
                for (ssize_t j = test_arr.getFirst2(); j <= test_arr.getLast2(); ++j)
                {
                    for (ssize_t k = test_arr.getFirst3(); k <= test_arr.getLast3(); ++k)
                    {
                        for (ssize_t l = test_arr.getFirst4(); l <= test_arr.getLast4(); ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug  << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (size_t l = test_arr.getArrayIndex(i, j, k, test_arr.getBegin4()); l <= test_arr.getArrayIndex(i, j, k, test_arr.getBegin4()) + test_arr.getLength4(i, j, k) - 1; ++l)
                        {
                            sum += test_arr.getArray()[l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            for (ssize_t i = test_arr2.getFirst1(); i <= test_arr2.getLast1(); ++i)
            {
                for (ssize_t j = test_arr2.getFirst2(); j <= test_arr2.getLast2(i); ++j)
                {
                    for (ssize_t k = test_arr2.getFirst3(); k <= test_arr2.getLast3(i, j); ++k)
                    {
                        for (ssize_t l = test_arr2.getFirst4(); l <= test_arr2.getLast4(i, j, k); ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}

TEST_F(IrregArrayPerformanceTest, IrregArray4DAccessOperatorPerformanceCachedIndices)
{
    const size_t                    ncycles = 1;
    const size_t                    dim1    = 30;
    const size_t                    dim2    = 5;
    const size_t                    dim3    = 30;
    const size_t                    dim4    = 5;
    const gmx::IrregArray1D<size_t> sizes1(1, dim1, dim2);
    const gmx::IrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    const gmx::IrregArray3D<size_t> sizes3(1, dim1, 1, dim2, 1, dim3, dim4);
    gmx::IrregArray4D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1, dim4, 1.0);
    gmx::IrregArray4D<double>       test_arr2(1, dim1, sizes1, sizes2, sizes3, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of IrregArray4D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of IrregArray4D " << test_arr2.getSize()
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
            const ssize_t first1 = test_arr.getFirst1();
            const ssize_t last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr.getFirst2();
                const ssize_t last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr.getFirst3();
                    const ssize_t last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr.getFirst4();
                        const ssize_t last4  = test_arr.getLast4();
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr.getArray()[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t first1 = test_arr.getFirst1();
            const ssize_t last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr.getFirst2();
                const ssize_t last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr.getFirst3();
                    const ssize_t last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr.getFirst4();
                        const ssize_t last4  = test_arr.getLast4();
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t first1 = test_arr.getFirst1();
            const ssize_t last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr.getFirst2();
                const ssize_t last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr.getFirst3();
                    const ssize_t last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr.getFirst4();
                        const ssize_t last4  = test_arr.getLast4();
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t first1 = test_arr2.getFirst1();
            const ssize_t last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr2.getFirst2();
                const ssize_t last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr2.getFirst3();
                    const ssize_t last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr2.getFirst4();
                        const ssize_t last4  = test_arr2.getLast4(i, j, k);
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2.getArray()[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t first1 = test_arr2.getFirst1();
            const ssize_t last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr2.getFirst2();
                const ssize_t last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr2.getFirst3();
                    const ssize_t last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr2.getFirst4();
                        const ssize_t last4  = test_arr2.getLast4(i, j, k);
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug  << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t first1 = test_arr2.getFirst1();
            const ssize_t last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t first2 = test_arr2.getFirst2();
                const ssize_t last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t first3 = test_arr2.getFirst3();
                    const ssize_t last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t first4 = test_arr2.getFirst4();
                        const ssize_t last4  = test_arr2.getLast4(i, j, k);
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}

TEST_F(IrregArrayPerformanceTest, FlatIrregArray4DAccessOperatorPerformanceCachedIndices)
{
    const size_t                        ncycles = 1;
    const size_t                        dim1    = 30;
    const size_t                        dim2    = 5;
    const size_t                        dim3    = 30;
    const size_t                        dim4    = 5;
    const gmx::IrregArray1D<size_t>     sizes1(1, dim1, dim2);
    const gmx::FlatIrregArray2D<size_t> sizes2(1, dim1, 1, dim2, dim3);
    const gmx::FlatIrregArray3D<size_t> sizes3(1, dim1, 1, dim2, 1, dim3, dim4);
    gmx::FlatIrregArray4D<double>       test_arr(1, dim1, 1, dim2, 1, dim3, 1, dim4, 1.0);
    gmx::FlatIrregArray4D<double>       test_arr2(1, dim1, sizes1, sizes2, sizes3, 1.0);

    // for size comparison
    debug << "Size of the regular special  case of FlatIrregArray4D " << test_arr.getSize()
    << " for " << test_arr.getNelements() << " array elements, "
    << "storage overhead per stored element relative to pure storage data: "
    <<  static_cast<double>(test_arr.getSize() - test_arr.getNelements() * sizeof(double))
    / static_cast<double>(test_arr.getNelements()) << endl;
    debug << "Size of the irregular normal case of FlatIrregArray4D " << test_arr2.getSize()
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
            const ssize_t   first1 = test_arr.getFirst1();
            const ssize_t   last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr.getFirst2();
                const ssize_t   last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr.getFirst3();
                    const ssize_t   last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const size_t   first4 = test_arr.getArrayIndex(i, j, k, test_arr.getBegin4());
                        const size_t   last4  = first4+ test_arr.getLength4() - 1;
                        for (size_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr.getArray()[l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t   first1 = test_arr.getFirst1();
            const ssize_t   last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr.getFirst2();
                const ssize_t   last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr.getFirst3();
                    const ssize_t   last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t   first4 = test_arr.getFirst4();
                        const ssize_t   last4  = test_arr.getLast4();
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t   first1 = test_arr.getFirst1();
            const ssize_t   last1  = test_arr.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr.getFirst2();
                const ssize_t   last2  = test_arr.getLast2();
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr.getFirst3();
                    const ssize_t   last3  = test_arr.getLast3();
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t   first4 = test_arr.getFirst4();
                        const ssize_t   last4  = test_arr.getLast4();
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t   first1 = test_arr2.getFirst1();
            const ssize_t   last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr2.getFirst2();
                const ssize_t   last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr2.getFirst3();
                    const ssize_t   last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const size_t   first4 = test_arr2.getArrayIndex(i, j, k, test_arr2.getBegin4());
                        const size_t   last4  = first4+ test_arr2.getLength4(i, j, k) - 1;
                        for (size_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2.getArray()[l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t   first1 = test_arr2.getFirst1();
            const ssize_t   last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr2.getFirst2();
                const ssize_t   last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr2.getFirst3();
                    const ssize_t   last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t   first4 = test_arr2.getFirst4();
                        const ssize_t   last4  = test_arr2.getLast4(i, j, k);
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
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
            const ssize_t   first1 = test_arr2.getFirst1();
            const ssize_t   last1  = test_arr2.getLast1();
            for (ssize_t i = first1; i <= last1; ++i)
            {
                const ssize_t   first2 = test_arr2.getFirst2();
                const ssize_t   last2  = test_arr2.getLast2(i);
                for (ssize_t j = first2; j <= last2; ++j)
                {
                    const ssize_t   first3 = test_arr2.getFirst3();
                    const ssize_t   last3  = test_arr2.getLast3(i, j);
                    for (ssize_t k = first3; k <= last3; ++k)
                    {
                        const ssize_t   first4 = test_arr2.getFirst4();
                        const ssize_t   last4  = test_arr2.getLast4(i, j, k);
                        for (ssize_t l = first4; l <= last4; ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        debug << "sum2 " << sum << endl;
    }

}


} // namespace

#ifdef DO_TIMING
#undef DO_TIMING
#endif
