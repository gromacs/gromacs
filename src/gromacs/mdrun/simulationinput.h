/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \file
 * \brief Public interface for SimulationInput facilities.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_mdrun
 */

#ifndef GMX_MDRUN_SIMULATIONINPUT_H
#define GMX_MDRUN_SIMULATIONINPUT_H

#include <memory>

namespace gmx
{

// Forward declarations for types from other modules that are opaque to the public API.
class LegacyMdrunOptions;

/*!
 * \brief Prescription for molecular simulation.
 *
 * Represents the complete and unique information needed to generate a simulation
 * trajectory segment. SimulationInput objects are opaque to the public API.
 * Clients can acquire SimulationInput objects through makeSimulationInput().
 *
 * See also https://gitlab.com/gromacs/gromacs/-/issues/3379 for design and development road map.
 */
class SimulationInput;

/*!
 * \brief Owning handle to a SimulationInput object.
 *
 * SimulationInput has no public API. Acquire a SimulationInputHolder with makeSimulationInput().
 * See https://gitlab.com/gromacs/gromacs/-/issues/3374
 */
class SimulationInputHolder
{
public:
    SimulationInputHolder();
    ~SimulationInputHolder();

    SimulationInputHolder(SimulationInputHolder&&) noexcept = default;
    SimulationInputHolder& operator=(SimulationInputHolder&&) noexcept = default;

    std::unique_ptr<SimulationInput> object_;
};

namespace detail
{

/*! \brief Get a SimulationInput.
 *
 * \internal
 * Design notes: SimulationInput creation will warrant a builder protocol, and
 * this helper can evolve into a director to apply the contents of LegacyMdrunOptions,
 * while such an operation is still relevant.
 *
 * Example:
 *     // After preparing a LegacyMdrunOptions and calling handleRestart()...
 *     SimulationInputBuilder builder;
 *     auto simulationInputHandle = makeSimulationInput(options, &builder);
 *
 *     // In addition to MdrunnerBuilder::addFiles(),
 *     mdrunnerBuilder.addInput(simulationInputHandle.get());
 *
 */
SimulationInputHolder makeSimulationInput(const char* tprFilename, const char* cpiFilename);

} // end namespace detail

} // end namespace gmx

#endif // GMX_MDRUN_SIMULATIONINPUT_H
