/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GROMACS_SIMULATIONCONTEXT_FORWARD_H
#define GROMACS_SIMULATIONCONTEXT_FORWARD_H

/*!
 * \file
 * \brief Provide ways for client code to own simulation resources.
 *
 * For `gmx mdrun` to be implemented as a client program, public API needs to
 * provide a way to create and manipulate the SimulationContext.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_mdrun
 * \inpublicapi
 */

#include <cstdio>

#include <memory>

struct gmx_hw_opt_t;
struct t_filenm;
struct t_commrec;
struct gmx_output_env_t;

namespace gmx
{

class MdFilenames;

/*!
 * \brief Simulation environment and configuration.
 *
 * SimulationContext allows a simulation module (\relates gmx::mdrun) to retrieve
 * runtime parameters and resources from client code. The client retains ownership
 * of the context and its resources, with exceptions as noted.
 *
 * The public interface of SimulationContext is not yet specified. Client code
 * can create an instance with gmx::createSimulationContext()
 *
 * \ingroup module_mdrun
 * \libinternal
 *
 * Placeholder for more complete implementation. Interfaces for different API
 * levels are not yet final.
 *
 * \inlibraryapi
 * \internal
 * Either the Context should be passed to the Runner as constant, or non-const
 * accessors need to be associated with update signals that the simulation
 * components (modules and runner) can subscribe to.
 */
class SimulationContext;

/*!
 * \brief Get ownership of a new SimulationContext object.
 *
 * A client must share ownership of some resources and transfer ownership of
 * other resources to create or configure the context. The simulation may in
 * turn consume or borrow some resources from the context.
 *
 * This call signature sets all parameters processed by the command-line client
 * code. Additional call signatures can allow implementation-specified default
 * values or newer initialization interfaces.
 *
 * \param simulationCommunicator Handle to communication data structure.
 * \param hardwareOptions Parallelism-related user options.
 * \param filenames Filenames and properties from command-line argument values.
 * \param outputEnvironment Output context for writing text files.
 * \param logFileHandle Handle to file used for logging.
 *
 */
std::unique_ptr<SimulationContext>
createSimulationContext(t_commrec        ** simulationCommunicator,
                        const gmx_hw_opt_t &hardwareOptions,
                        const MdFilenames  &filenames,
                        gmx_output_env_t ** outputEnvironment,
                        FILE             ** logFileHandle);

}      // end namespace gmx

#endif //GROMACS_SIMULATIONCONTEXT_FORWARD_H
