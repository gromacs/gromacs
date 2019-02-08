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
   Definitions commonly used by the data structure unit tests.

   \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \ingroup module_math
 */

#ifndef GMX_MATH_IRREG_ARRAY_MATH_TESTS_COMMONS_H
#define GMX_MATH_IRREG_ARRAY_MATH_TESTS_COMMONS_H

#include <fstream>
#include <iomanip>  // ::std::setw, ::std::setfill
#include <ios>      // ::std::left, ::std::right
#include <iostream> // ::debug
#include <ostream>  // ::std::ostream (<<)
#include <stdexcept>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/irreg_array_math.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{

namespace irreg_array_math_test
{

//! floating point data types to be considered in the typed tests
typedef ::testing::Types<float, double, long double> TestFloatTypes;
//typedef ::testing::Types<float> TestFloatTypes;

//! flag triggering debug output on user-request via commandline switch [no]debug
extern bool g_debugBool;
//! flag triggering error output on user-request via commandline switch [no]error
extern bool g_errorBool;


class IrregArrayMathTest : public ::testing::Test
{
    public:
        IrregArrayMathTest()
            : debug(nullptr), error(nullptr)
        {
        }
        void SetUp() override
        {
            // trigger debug output
            if (g_debugBool)
            {
                activateDebugOut();
            }
            if (g_errorBool)
            {
                activateErrorOut();
            }

        }

        void TearDown() override
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

    protected:
        gmx::test::TestFileManager   fileManager_;
        // debugging output stream, silent by default
        std::ostream                 debug;
        // error output stream, silent by default
        std::ostream                 error;
};

/*! \brief  output stream operator to insert a string representation of
            a vector into an output stream, e.g., cout

    \tparam   TReal        floating point data type

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   vec      the object whose contents are to be inserted
 */
template<typename TReal = double, class Alloc = std::allocator<TReal> >
std::ostream &operator<<(std::ostream &output, const std::vector<TReal, Alloc> &vec)
{
    if (vec.empty())
    {
        return output;
    }
    constexpr size_t extraFloatWidth = 4;
    for (typename std::vector<TReal, Alloc>::size_type i = 0; i < vec.size(); ++i)
    {
        output << std::setw(extraFloatWidth + std::numeric_limits<TReal>::digits) << std::fixed << std::setprecision(std::numeric_limits<TReal>::digits) << vec[i];
        if (i < (vec.size() - 1))
        {
            output << "  ";
        }
    }
    return output;
}


}      // namespace irreg_array_math_test

}      // namespace gmx

#endif // end header
