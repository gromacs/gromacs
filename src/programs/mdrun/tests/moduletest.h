/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \libinternal
 * \defgroup module_mdrun_integration_tests Integration test utilities
 * \ingroup group_mdrun
 *
 * \brief Functionality for testing mdrun as a whole
 */
#ifndef GMX_MDRUN_TESTS_MODULETEST_H
#define GMX_MDRUN_TESTS_MODULETEST_H

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/integrationtests.h"

namespace gmx
{

namespace test
{

/*! \libinternal \brief Helper object for running grompp and mdrun in
 * integration tests of mdrun functionality
 *
 * Objects of this class are intended to be owned by
 * IntegrationTestFixture objects, and an IntegrationTestFixture
 * object might own more than one SimulationRunner.
 *
 * The setup phase creates various temporary files for input and
 * output that are common for mdrun tests, using the file manager
 * object of the fixture that owns this object. Individual tests
 * should create any extra filenames similarly, so that the test
 * users's current working directory does not get littered with files
 * left over from tests.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class SimulationRunner
{
    public:
        /*! \brief Constructor, which establishes the fixture that
         * will own each object */
        explicit SimulationRunner(IntegrationTestFixture *fixture_);

        //! Use an empty .mdp file as input to grompp
        void useEmptyMdpFile();
        //! Use a given string as input to grompp
        void useStringAsMdpFile(const char *mdpString);
        //! Use a given string as input to grompp
        void useStringAsMdpFile(const std::string &mdpString);
        //! Use a string as -n input to grompp
        void useStringAsNdxFile(const char *ndxString);
        //! Use a standard .top and .gro file as input to grompp
        void useTopGroAndNdxFromDatabase(const char *name);
        //! Use a standard .gro file as input to grompp
        void useGroFromDatabase(const char *name);
        //! Calls grompp (on rank 0) to prepare for the mdrun test
        int callGrompp();
        //! Calls grompp (on this rank) to prepare for the mdrun test
        int callGromppOnThisRank();
        //! Calls mdrun for testing with a customized command line
        int callMdrun(const CommandLine &callerRef);
        /*! \brief Convenience wrapper for calling mdrun for testing
         * with default command line */
        int callMdrun();

    private:
        //! Provides access to the test fixture, e.g. for the TestFileManager
        IntegrationTestFixture *fixture_;
    public:
        //@{
        /*! \name Names for frequently used grompp and mdrun output files
         *
         * These strings can be set to point to files present in the
         * source tree, or to temporary files created for the test
         * fixture. In the latter case,
         * IntegrationTestFixture::fileManager_ should be used to fill
         * these strings with paths to files, so that they are created
         * in a temporary directory and (by default behaviour of
         * TestFileManager) deleted when the test is complete.
         */
        std::string topFileName_;
        std::string groFileName_;
        std::string fullPrecisionTrajectoryFileName_;
        std::string reducedPrecisionTrajectoryFileName_;
        std::string groOutputFileName_;
        std::string ndxFileName_;
        std::string mdpInputFileName_;
        std::string mdpOutputFileName_;
        std::string tprFileName_;
        std::string logFileName_;
        std::string edrFileName_;
        std::string cptFileName_;
        std::string swapFileName_;
        int         nsteps_;
        //@}
};

/*! \libinternal \brief Declares test fixture base class for
 * integration tests of mdrun functionality
 *
 * Derived fixture classes (or individual test cases) that might have
 * specific requirements should assert that behaviour, rather than
 * hard-code the requirements. A test that (for example) can't run
 * with more than one thread should report that as a diagnostic, so the
 * person running the test (or designing the test harness) can get
 * feedback on what tests need what conditions without having to read
 * the code of lots of tests.
 *
 * Specifying the execution context (such as numbers of threads and
 * processors) is normally sensible to specify from the test harness
 * (i.e. when CMake/CTest/the user runs a test executable), because
 * only there is information about the hardware available. The default
 * values for such context provided in test fixtures for mdrun should
 * mirror the defaults for mdrun, but currently mdrun.c hard-codes
 * those in a gmx_hw_opt_t.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunTestFixtureBase : public IntegrationTestFixture
{
    public:
        MdrunTestFixtureBase();
        virtual ~MdrunTestFixtureBase();
};

/*! \libinternal \brief Declares test fixture class for integration
 * tests of mdrun functionality that use a single call of mdrun
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunTestFixture : public IntegrationTestFixture
{
    public:
        MdrunTestFixture();
        virtual ~MdrunTestFixture();

        //! Helper object to manage the preparation for and call of mdrun
        SimulationRunner runner_;
};

/*! \libinternal \brief
 * Parameterized test fixture for mdrun integration tests
 */
class ParameterizedMdrunTestFixture : public gmx::test::MdrunTestFixture,
                                      public ::testing::WithParamInterface<const char *>
{
};

} // namespace test
} // namespace gmx

#endif
