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
 * Declares test fixture for mdrun tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_TESTS_MODULETEST_H
#define GMX_MDRUN_TESTS_MODULETEST_H

#include <gtest/gtest.h>

#include "gromacs/utility/common.h"
#include "testutils/integrationtests.h"

namespace gmx
{

namespace test
{

/*! \internal \brief
 * Test fixture for mdrun.
 *
 * Specifying the execution context (such as numbers of threads and
 * processors) is normally sensible to specify from the test harness
 * (i.e. when CMake/CTest/the user runs a test executable), because
 * only there is information about the hardware available. The default
 * values for such context provided in test fixtures for mdrun should
 * mirror the defaults for mdrun, but currently mdrun.c hard-codes
 * those in a gmx_hw_opt_t.
 *
 * Derived fixture classes (or individual test cases) that might have
 * specific requirements should assert that behaviour, rather than
 * hard-code the requirements. A test that (for example) can't run
 * with more than one thread should report that as a diagnostic, so the
 * person running the test (or designing the test harness) can get
 * feedback on what tests need what conditions without having to read
 * the code of lots of tests.
 *
 * The setup phase creates various temporary files for input and
 * output that are common for mdrun tests. Individual tests should
 * create any extra filenames similarly, so that the test users's CWD
 * does not get littered with files left over from all manner of
 * tests.
 *
 * // Any method in this class may throw std::bad_alloc if out of
 * memory.
 *
 * \ingroup module_mdrun
 */
class MdrunTestFixture : public IntegrationTestFixture
{
    protected:
        MdrunTestFixture();
        virtual ~MdrunTestFixture();

        /*! \brief
         * Sets up the test
         */
        //virtual void SetUp();
        //! \brief Use an empty .mdp file as input to grompp
        void useEmptyMdpFile();
        //! \brief Use an empty .mdp file as input to grompp
        void useStringAsMdpFile(const char *mdpString);
        //! \brief Use a standard .top and .gro file as input to grompp
        void useTopAndGroFromDatabase(const char *name);
        /*! \brief Configures GROMACS programs with default options
         * generally suitable for testing of mdrun.
         *
         * Other options can be configured later with \p setProgramOption.
         */
        void configureProgramWithDefaultOptions(std::string const& programName);
        /*! \brief Calls grompp to prepare for the mdrun test
         */
        int callGrompp();
        /*! \brief Calls mdrun with default arguments for testing
         */
        int callMdrun();
        /*! /brief Number of tMPI threads for child mdrun call */
        int numThreads;
        /*! /brief Number of OpenMP threads for child mdrun call */
        int numOpenMPThreads;
        /*! /brief Input filenames for grompp or mdrun
         */
        std::string mdpFileName;
        std::string topFileName;
        std::string groFileName;
        std::string trrFileName;
        /*! /brief Filenames for output files.
         *
         * These are created in a temporary directory and deleted when
         * the test is complete.
         */
        std::string mdpOutputFileName;
        std::string tprFileName;
        std::string logFileName;
        std::string edrFileName;

    private:
};
 
} // namespace test
} // namespace gmx

#endif
