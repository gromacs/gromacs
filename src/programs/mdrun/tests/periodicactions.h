/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

/*! \internal \file
 * \brief Interfaces of related classes for tests to verify that a simulator that only does some
 * actions periodically produces the expected results.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "programs/mdrun/tests/comparison_helpers.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{

/*! \brief Mdp parameters that determine the manner of simulation
 * propagation. */
using PropagationParameters = MdpFieldValues;

/*! \brief Mdp parameters that should only affect the observations,
 *  not the simulation propagation. */
using PeriodicOutputParameters = MdpFieldValues;

//! \internal \brief Helper type of output file names for the reference mdrun call
struct ReferenceFileNames
{
    //! Name of energy file
    std::string edrFileName_;
};

//! Function type to produce sets of .mdp parameters for testing periodic output
using OutputParameterGeneratorFunction = std::vector<PeriodicOutputParameters> (*)();

/*! \internal
 *  \brief Test fixture base for comparing a simulator with one that does everything every step
 *
 * This test ensures that two simulator code paths called via
 * different mdp options yield identical energy trajectories,
 * up to some (arbitrary) precision.
 *
 * These tests are useful to check that periodic actions implemented
 * in simulators are correct, and that different code paths expected
 * to yield identical results are equivalent.
 */
class PeriodicActionsTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<PropagationParameters, OutputParameterGeneratorFunction>>
{
public:
    // PeriodicActionsTest();
    //! Run mdrun with given output parameters
    void doMdrun(const PeriodicOutputParameters& output);
    //! Generate reference data from mdrun writing everything every step.
    void prepareReferenceData();
    //! Names for the output files from the reference mdrun call
    ReferenceFileNames referenceFileNames_ = { fileManager_.getTemporaryFilePath("reference.edr").string() };
    //! Functor for energy comparison
    EnergyComparison energyComparison_{ EnergyComparison::defaultEnergyTermsToCompare(),
                                        MaxNumFrames::compareAllFrames() };
    //! Names of energies compared by energyComparison_
    std::vector<std::string> namesOfEnergiesToMatch_ = energyComparison_.getEnergyNames();
};

/*! \brief Return vector of mdp parameter sets to test
 *
 * These are constructed to observe the mdp parameter choices that
 * only affect when output is written agree with those that were
 * written from a reference run where output was done every step. The
 * numbers are chosen in the context of the defaults in
 * prepareDefaultMdpFieldValues().
 *
 * \todo Test nstlog, nstdhdl, nstxout-compressed */
std::vector<PeriodicOutputParameters> outputParameters();

//! Returns sets of simple simulation propagators
std::vector<PropagationParameters> simplePropagationParameters();

/*! \brief Returns sets of simulation propagators including coupling
 *
 * These are chosen to cover the commonly used space of propagation
 * algorithms togther with the perdiods between their actions. The
 * periods tested are chosen to be mutually co-prime and distinct from
 * the pair search and user output period (4), so that over the
 * duration of a short simulation many kinds of simulation step
 * behavior are tested. */
std::vector<PropagationParameters> propagationParametersWithCoupling();

/*! \brief Returns sets of simulation propagators including coupling
 *
 * These are chosen to cover the commonly used space of propagation
 * algorithms on systems with constraints. */
std::vector<PropagationParameters> propagationParametersWithConstraints();

} // namespace test
} // namespace gmx
