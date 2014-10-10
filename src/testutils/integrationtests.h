/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares test fixture for integration tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_INTEGRATION_TESTS_MODULETEST_H
#define GMX_INTEGRATION_TESTS_MODULETEST_H

#include <gtest/gtest.h>

#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

/*! \libinternal \brief
 * Test fixture for integration tests.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class IntegrationTestFixture : public ::testing::Test
{
    protected:
        IntegrationTestFixture();
        virtual ~IntegrationTestFixture();

        /*! \brief Accepts a string as input, writes it to a temporary
         * file and then reopens stdin to read the contents of that
         * string.
         *
         * \throws FileIOError  when the freopen() fails
         */
        void redirectStringToStdin(const char* theString);

        /*! \brief Discards stdout while running a test
         *
         * \todo Implement this when the output routines are
         * sufficiently modular to permit it to work. */
        void redirectStdoutToDevNull();
        /*! \brief Discards stderr while running a test
         *
         * \todo Implement this when the output routines are
         * sufficiently modular to permit it to work. */
        void redirectStderrToDevNull();

        /* TEST_F() constructs derived classes, and those classes and
         * their member objects might need to access implementation
         * details, so we cannot use the private access specifer
         * here. */
    public:

        /*! \brief Object that manages finding input files, writing
         * temporary output files and cleaning up files.
         */
        ::gmx::test::TestFileManager fileManager_;
};

} // namespace test
} // namespace gmx

#endif
