/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \brief
 * Declares test fixture for integration tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_INTEGRATION_TESTS_MODULETEST_H
#define GMX_INTEGRATION_TESTS_MODULETEST_H

#include <gtest/gtest.h>

#include "gromacs/utility/common.h"
#include "testutils/testfilemanager.h"
#include "testutils/programcaller.h"

namespace gmx
{

namespace test
{

/*! \internal \brief
 * Test fixture for integration tests.
 *
 * // Any method in this class may throw std::bad_alloc if out of
 * memory.
 *
 * \ingroup module_integration_tests
 */
class IntegrationTestFixture : public ::testing::Test
{
    protected:
        IntegrationTestFixture();
        virtual ~IntegrationTestFixture();

        /*! \brief
         * Sets up the test
         */
        //virtual void SetUp();
        /*! \brief Copies the value of a string option to be used when running the program.
         */
        void setProgramOption(std::string const& programName, const char *optionName, std::string const& value);
        /*! \brief Copies the value of an integer option to be used when running the program.
         */
        void setProgramOption(std::string const& programName, const char *optionName, int const& value);
        /*! \brief Copies the value of a double option to be used when running the program.
         */
        void setProgramOption(std::string const& programName, const char *optionName, double const& value);
        // TODO helper functions for other option types (enum, bool)
        /*! \brief Defines the named option for use with the program.
         *
         * Also copies the value to be used when running the program.
         *
         * The OptionType must implement the AbstractOption interface.
         */
        template <class OptionType> void
            setProgramOption(std::string const& programName, const char *optionName, std::string const& value);
        /*! \brief Sets the value of the integer option to the program.
         *
         * TODO cater for boolean options, since the output has to be "-(no)?optionname" when the visitor gets there later.
         */
        //        void setBoolOption(std::string const& programName, const char *optionName, gmx_bool value);
        /*! \brief Calls a GROMACS program
         *
         * \param[in] programName  Name of the GROMACS program
         */
        int runProgram(std::string const& programName);
        /*! \brief Accepts a string as input, writes it to a temporary
         * file and then reopens stdin to read the contents of that
         * string.
         *
         * \throws FileIOError  when the freopen() fails
         */
        void redirectStringToStdin(const char* theString);
        /*! \brief Factory method to build an object that can manage
         * the options for and call to a GROMACS program via its
         * cmain function.
         *
         * No options are defined for use with this program.
         *
         * \param[in] programName  Name of the GROMACS program
         * \param[in] cmain        The cmain function of that GROMACS program
         */
        void createProgramCaller(std::string const& programName, ProgramCaller::cmain_func_ptr cmain);

        /* TEST_F() constructs derived classes, and those classes
         * might need to customize the argument lists and/or files
         * used by the programs called, so these instance variables
         * cannot be private. */

        /*! /brief Object that manages finding input files, writing
         * temporary output files and cleaning up files.
         */
        ::gmx::test::TestFileManager fileManager;

    protected:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};
 
} // namespace test
} // namespace gmx

#endif
