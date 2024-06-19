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
#ifndef GMX_MDRUN_TESTS_MULTISIMTEST_H
#define GMX_MDRUN_TESTS_MULTISIMTEST_H

/*! \internal \file
 * \brief
 * Declares test fixture for the mdrun multi-simulation functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */

#include <memory>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"

enum class IntegrationAlgorithm : int;
enum class PressureCoupling : int;
enum class TemperatureCoupling : int;

namespace gmx
{
namespace test
{

//! Convenience typedef
typedef std::unique_ptr<CommandLine> CommandLinePointer;

/*! \internal
 * \brief Test flag setting the number of ranks per simulation
 */
class NumRanksPerSimulation final
{
public:
    //! Explicit constructor
    explicit NumRanksPerSimulation(int value) : value_(value) {}
    //! Implicit conversion to int
    operator int() const { return value_; }

private:
    //! Internal state
    int value_;
};

//! The parameters of the MultiSimTest class
typedef std::tuple<NumRanksPerSimulation, IntegrationAlgorithm, TemperatureCoupling, PressureCoupling> MultiSimTestParams;

/*! \internal
 * \brief Test fixture for multi-sim functionality.
 *
 * This is intended to be re-used also for tests of functionality that
 * are derived from multi-sim, e.g. REMD.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MultiSimTest : public MdrunTestFixtureBase, public ::testing::WithParamInterface<MultiSimTestParams>
{
public:
    MultiSimTest();

    /*! \brief Check whether the MPI setup is valid
     *
     * Excludes MPI setups which are not supported by multi-sim
     */
    bool mpiSetupValid() const;

    /*! \brief Organize the .mdp file for this rank
     *
     * For testing multi-simulation, this .mdp file is more
     * complicated than it needs to be, but it does little harm,
     * and doing it this way allows this function to be re-used
     * for testing replica-exchange.
     *
     * The mdp options, specifically the temperature and pressure
     * coupling, allow parameterization to work with
     * T, P or (later) lambda as the control variable, by passing a
     * string with a value for the respective mdp option such that
     * different paths in init_replica_exchange() are followed.
     *
     * \param runner      The simulation runner that uses the
     *                    mdp file that is organized.
     * \param integrator  Value for the mdp option `integrator`
     * \param tcoupl      Value for the mdp option `tcoupl`
     * \param pcoupl      Value for the mdp option `pcoupl`
     * \param numSteps    Number of MD steps to perform.
     * \param doRegression  Whether the mdp file will be used for
     *                      regression tests, request use of
     *                      reproducible parameters
     */
    void organizeMdpFile(SimulationRunner*    runner,
                         IntegrationAlgorithm integrator,
                         TemperatureCoupling  tcoupl,
                         PressureCoupling     pcoupl,
                         int                  numSteps,
                         bool                 doRegression) const;

    /*! \brief Run grompp on the ranks that need to run it
     *
     * Adds an MPI barrier to ensure that all ranks have written before returning
     *
     * \param runner        The simulation runner to run grompp on
     * \param numSteps      Number of MD steps to perform
     * \param doRegression  Whether to write trajectories during the simulation
     * \param maxWarnings   Number of grompp warning tolerated
     */
    void runGrompp(SimulationRunner* runner, int numSteps = 2, bool doRegression = false, int maxWarnings = 0) const;
    //! Test that a basic simulation works
    void runExitsNormallyTest();
    //! Test that mdrun -maxh and restart works
    void runMaxhTest();
    //! Number of MPI ranks
    int size_;
    //! MPI rank of this process
    int rank_;
    //! Number of ranks per simulation
    int numRanksPerSimulation_;
    //! The simulation this rank belongs to (equal to `int( rank_ / numRanksPerSimulation_)`)
    int simulationNumber_;
    //! Object for building the mdrun command line
    CommandLinePointer mdrunCaller_;
    //! Manages temporary files during the test.
    TestFileManager fileManager_;
};

} // namespace test
} // namespace gmx

#endif
