/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \defgroup module_mdrun_integration_tests Integration test utilities
 * \ingroup group_mdrun
 *
 * \brief Functionality for testing mdrun as a whole
 */
/*! \internal \file
 * \brief
 * Declares test fixtures for general mdrun functionality.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_MDRUN_TESTS_MODULETEST_H
#define GMX_MDRUN_TESTS_MODULETEST_H

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

struct gmx_hw_info_t;

namespace gmx
{
namespace test
{

/*! \internal
 * \brief How the mdp file of the SimulationRunner is defined
 */
enum class SimulationRunnerMdpSource
{
    //! The default behavior. Will result in an empty mdp file.
    Undefined,
    //! Mdp options are set via string using SimulationRunner::useStringAsMdpFile
    String,
    //! Mdp options are read from a file set in SimulationRunner::useTopGroAndMdpFromFepTestDatabase
    File,
    //! Signals the last enum entry
    Count
};

/*! \internal
 * \brief Helper object for running grompp and mdrun in
 * integration tests of mdrun functionality
 *
 * Objects of this class must be owned by objects descended from
 * MdrunTestFixtureBase, which sets up necessary infrastructure for
 * it. Such an object may own more than one SimulationRunner.
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
 * By default, the convenience methods callGrompp() and callMdrun()
 * just prepare and run a default call to mdrun. If there is a need to
 * customize the command-line for grompp or mdrun (e.g. to invoke
 * -maxwarn n, or -reprod), then make a CommandLine object with the
 * appropriate flags and pass that into the routines that accept such.
 *
 * \ingroup module_mdrun_integration_tests
 */
class SimulationRunner
{
public:
    //! Initializes a runner with given manager for temporary files.
    explicit SimulationRunner(TestFileManager* fileManager);

    //! Use an empty .mdp file as input to grompp
    void useEmptyMdpFile();
    //! Use a given string as input to grompp
    void useStringAsMdpFile(const char* mdpString);
    //! Use a given string as input to grompp
    void useStringAsMdpFile(const std::string& mdpString);
    //! Use a string as -n input to grompp
    void useStringAsNdxFile(const char* ndxString) const;
    //! Use a standard .top and .g96 file as input to grompp
    void useTopG96AndNdxFromDatabase(const std::string& name);
    //! Use a standard .top and .gro file as input to grompp
    void useTopGroAndNdxFromDatabase(const std::string& name);
    //! Use a standard .gro file as input to grompp
    void useGroFromDatabase(const char* name);
    //! Use a standard .ndx as input to grompp
    void useNdxFromDatabase(const std::string& name);
    //! Use .top, .gro, and .mdp from FEP test database
    void useTopGroAndMdpFromFepTestDatabase(const std::string& name);
    //! Set a maxmum number of acceptable warnings.
    void setMaxWarn(int maxwarn);
    //! Calls grompp (on rank 0, with a customized command line) to prepare for the mdrun test
    int callGrompp(const CommandLine& callerRef);
    //! Convenience wrapper for a default call to \c callGrompp
    int callGrompp();
    //! Calls grompp (on this rank, with a customized command line) to prepare for the mdrun test
    int callGromppOnThisRank(const CommandLine& callerRef);
    //! Convenience wrapper for a default call to \c callGromppOnThisRank
    int callGromppOnThisRank();
    //! Calls nmeig for testing
    int callNmeig() const;
    //! Calls mdrun for testing with a customized command line
    int callMdrun(const CommandLine& callerRef);
    /*! \brief Convenience wrapper for calling mdrun for testing
     * with default command line */
    int callMdrun();
    //! Calls convert-tpr on this rank to set a new number of steps in the tpr.
    int changeTprNsteps(int nsteps) const;

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
    std::string mdpFileName_;
    std::string fullPrecisionTrajectoryFileName_;
    std::string reducedPrecisionTrajectoryFileName_;
    std::string groOutputFileName_;
    std::string cptOutputFileName_;
    std::string ndxFileName_;
    std::string mdpOutputFileName_;
    std::string tprFileName_;
    std::string logFileName_;
    std::string edrFileName_;
    std::string mtxFileName_;
    std::string swapFileName_;
    std::string dhdlFileName_;
    int         nsteps_;
    int         maxwarn_;
    //@}
    //! How the mdp options are defined
    SimulationRunnerMdpSource mdpSource_;
    //! What will be written into a temporary mdp file before the grompp call
    std::string mdpInputContents_;

private:
    //! The file manager used to manage I/O
    TestFileManager& fileManager_;

    GMX_DISALLOW_COPY_AND_ASSIGN(SimulationRunner);
};

/*! \internal
 * \brief Declares test fixture base class for
 * integration tests of mdrun functionality
 *
 * Heavyweight resources are set up here and shared
 * across all tests in the test case fixture, e.g.
 * the MPI communicator for the tests and the hardware
 * detected that is available to it.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunTestFixtureBase : public ::testing::Test
{
public:
    //! Per-test-case setup for lengthy processes that need run only once.
    static void SetUpTestSuite();
    //! Per-test-case tear down
    static void TearDownTestSuite();

    MdrunTestFixtureBase();
    ~MdrunTestFixtureBase() override;

    //! Communicator over which the test fixture works
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static MPI_Comm s_communicator;
    /*! \brief Hardware information object
     *
     * Detected within \c communicator_ and available to re-use
     * over all tests in the test case of this text fixture. */
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static std::unique_ptr<gmx_hw_info_t> s_hwinfo;
};

/*! \internal
 * \brief Declares test fixture class for integration
 * tests of mdrun functionality that use a single call of mdrun
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunTestFixture : public MdrunTestFixtureBase
{
public:
    MdrunTestFixture();
    ~MdrunTestFixture() override;

    //! Manages temporary files during the test.
    TestFileManager fileManager_;
    //! Helper object to manage the preparation for and call of mdrun
    SimulationRunner runner_;
};

/*! \internal
 * \brief
 * Parameterized test fixture for mdrun integration tests
 */
class ParameterizedMdrunTestFixture :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<const char*>
{
};

/*! \internal
 * \brief
 * Returns the number of OpenMP threads to use.
 *
 * \returns the number specified using \c -ntomp option, or the default.
 */
int getNumberOfTestOpenMPThreads();

} // namespace test
} // namespace gmx

#endif
