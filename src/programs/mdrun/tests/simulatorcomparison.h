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
 * \brief
 * Helper classes for tests that compare the results of equivalent
 * simulation runs. Currently used for the rerun and the simulator
 * tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */

#ifndef GMX_PROGRAMS_MDRUN_TESTS_SIMULATORCOMPARISON_H
#define GMX_PROGRAMS_MDRUN_TESTS_SIMULATORCOMPARISON_H

#include <string>

namespace gmx
{
namespace test
{
/*!
 * \brief Run and compare a simulator run with and without an environment variable
 *
 * Run grompp, and repeat mdrun with and without the environment variable set.
 * Compare energies (via EnergyComparator) and trajectories.
 */
template<typename... Args>
void executeSimulatorComparisonTest(const std::string& environmentVariable, Args&&... args);

/*!
 * \brief Run and compare a simulator run to its rerun
 *
 * Run grompp, run mdrun and rerun the resulting trajectory.
 * Compare energies (via EnergyComparator) and trajectories.
 */
template<typename... Args>
void executeRerunTest(Args&&... args);

} // namespace test
} // namespace gmx

// Including this here avoid having to put everything in the header file,
// or to explicitly declare the templates (which would render the parameter
// pack useless)
#include "simulatorcomparison.cpp"

#endif // GMX_PROGRAMS_MDRUN_TESTS_SIMULATORCOMPARISON_H
