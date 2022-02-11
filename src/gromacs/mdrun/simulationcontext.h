/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#ifndef GROMACS_SIMULATIONCONTEXT_H
#define GROMACS_SIMULATIONCONTEXT_H

/*! \libinternal
 * \file
 * \brief Provide ways for client code to own simulation resources.
 *
 * For `gmx mdrun` to be implemented as a client program, public API needs to
 * provide a way to create and manipulate the SimulationContext.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_mdrun
 * \inlibraryapi
 */

#include <memory>
#include <string>

#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_multisim_t;

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \cond libapi
 * \libinternal
 * \brief Simulation environment and configuration.
 *
 * SimulationContext allows a simulation module (\relates gmx::mdrun) to retrieve
 * runtime parameters and resources from client code. The client retains ownership
 * of the context and its resources, with exceptions as noted.
 *
 * A client must share ownership of some resources and transfer ownership of
 * other resources to create or configure the context. The simulation may in
 * turn consume or borrow some resources from the context. This is a new
 * framework that will evolve in the contexts of
 * https://gitlab.com/gromacs/gromacs/-/issues/2375
 * https://gitlab.com/gromacs/gromacs/-/issues/2587
 *
 * The public interface of SimulationContext is not yet well-specified.
 * Client code can create an instance with gmx::createSimulationContext()
 *
 * \warning The SimulationContext does not yet encapsulate the resources allocated
 * to the simulator. See https://gitlab.com/gromacs/gromacs/-/issues/3650
 *
 * \warning On thread-MPI code paths, the SimulationContext *does not* represent
 * constraints on the number of simulation ranks, nor does it represent initialized
 * communication resources.
 *
 * \todo Clarify and strengthen the invariant represented by SimulationContext.
 *
 * \todo This class should also handle aspects of simulation
 * environment such as working directory and environment variables.
 *
 * \ingroup module_mdrun
 * \inlibraryapi
 *
 * \internal
 * This is a minimal placeholder for a more complete implementation.
 * Interfaces for different API levels are not yet final.
 * \todo Impose sensible access restrictions.
 * Either the SimulationContext should be passed to the Mdrunner as logically constant or
 * a separate handle class can provide access to resources that have been
 * allocated by (negotiated with) the client for the current simulation
 * (or simulation segment).
 * Non-const accessors to shared resources need to be associated with update
 * signals that the simulation components (modules and runner) can subscribe to.
 *
 * Also reference https://gitlab.com/gromacs/gromacs/-/issues/2587
 */
class SimulationContext final
{
public:
    /*!
     * \brief Object must be initialized with non-default constructor.
     */
    SimulationContext() = delete;
    /*!
     * \brief Initialize from borrowed communicator.
     *
     * \param communicator Communicator for all collaborating simulation processes.
     * \param multiSimDirectoryNames  Names of any directories used with -multidir
     *
     * Caller is responsible for keeping *communicator* valid for the life of SimulationContext,
     * and then freeing the communicator, if appropriate.
     *
     * With an external MPI library (non-thread-MPI chosen when configuring with CMake),
     * the client promises that MPI has been initialized (such as by calling gmx::init()).
     * This communicator is "borrowed" (not duplicated) from the caller.
     * Additional communicators may be split from the provided communicator
     * during the life of the SimulationContext or its consumers.
     *
     * With thread-MPI, *communicator* must be MPI_COMM_NULL.
     * The communicator is set up later
     * during the process of spawning the threads that will be the MPI
     * ranks. (Multi-simulation is not supported with thread-MPI.)
     *
     * \todo Refine distribution of environment management.
     *       There should be a context level at which only the current simulator directory matters,
     *       and a level above which encapsulates multisim details in a specialized type.
     */
    explicit SimulationContext(MPI_Comm communicator, ArrayRef<const std::string> multiSimDirectoryNames);

    /*!
     * \brief MPI resources for the entire simulation context.
     *
     * With an external MPI library (non-thread-MPI chosen when configuring with CMake),
     * gmx::init() has called MPI_Init and the provided communicator is valid to use.
     * The communicator is "borrowed" (not duplicated) from the caller.
     *
     * With thread-MPI, the communicator is set up later
     * during the process of spawning the threads that will be the MPI
     * ranks.
     */
    MPI_Comm libraryWorldCommunicator_ = MPI_COMM_NULL;

    /*!
     * \brief MPI communicator object for this simulation.
     *
     * With an external MPI library (non-thread-MPI chosen when configuring with CMake),
     * gmx::init() has called MPI_Init and the provided communicator is valid to use.
     * The communicator is "borrowed" (not duplicated) from the caller.
     *
     * With thread-MPI, the communicator is set up later
     * during the process of spawning the threads that will be the MPI
     * ranks.
     */
    MPI_Comm simulationCommunicator_ = MPI_COMM_NULL;

    /*!
     * \brief Multi-sim handler (if required by e.g. gmx mdrun
     * -multidir; only supported with real MPI)
     *
     * If a multi-simulation is in use, then its
     * communicator(s) are found in multiSimulation_. This
     * communicator is that of all ranks from all simulations, and
     * will later be split into one for each simulation.
     * TODO: Perhaps (for simplicity) communicator splitting
     *       can be undertaken during multi-sim setup (acquisition of the multisim resource).
     *
     * Multi-simulation is not supported with thread-MPI.
     */
    std::unique_ptr<gmx_multisim_t> multiSimulation_;
};
//! \endcond

} // end namespace gmx

#endif // GROMACS_SIMULATIONCONTEXT_H
