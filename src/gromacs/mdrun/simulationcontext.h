/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#include "gromacs/gmxlib/network.h"

namespace gmx
{

/*! \cond libapi
 * \libinternal
 * \brief Simulation environment and configuration.
 *
 * SimulationContext allows a simulation module (\relates gmx::mdrun) to retrieve
 * runtime parameters and resources from client code. The client retains ownership
 * of the context and its resources, with exceptions as noted.
 *
 * The public interface of SimulationContext is not yet well-specified.
 * Client code can create an instance with gmx::createSimulationContext()
 *
 * \todo This class should also handle aspects of simulation
 * environment such as working directory and environment variables.
 *
 * \ingroup module_mdrun
 * \inlibraryapi
 *
 * \internal
 * This is a minimal placeholder for a more complete implementation.
 * Interfaces for different API levels are not yet final, but also depend on
 * additional development of t_commrec and other resources.
 * \todo Impose sensible access restrictions.
 * Either the SimulationContext should be passed to the Mdrunner as logically constant or
 * a separate handle class can provide access to resources that have been
 * allocated by (negotiated with) the client for the current simulation
 * (or simulation segment).
 * Non-const accessors to shared resources need to be associated with update
 * signals that the simulation components (modules and runner) can subscribe to.
 *
 * Also reference https://redmine.gromacs.org/issues/2587
 */
class SimulationContext final
{
    public:
        /*!
         * \brief Object must be initialized with non-default constructor.
         */
        SimulationContext() = delete;

        /*!
         * \brief Initializate with borrowed values.
         *
         * \param communicationRecord non-owning communication record handle.
         *
         * Client code is responsible for cleaning up communicationRecord after
         * SimulationContext is destroyed.
         *
         * \internal
         * SimulationContext should be the owner of these objects and these implementation
         * details are subject to change as ownership semantics are clarified in future
         * development.
         */
        explicit SimulationContext(CommrecHandle communicationRecord);
        //! Move constructor
        SimulationContext(SimulationContext && source) noexcept;

        /*! \brief Communicator handle. */
        CommrecHandle communicationRecord_;
};
//! \endcond

/*! \cond libapi
 * \brief Get ownership of a new SimulationContext object.
 *
 * A client must share ownership of some resources and transfer ownership of
 * other resources to create or configure the context. The simulation may in
 * turn consume or borrow some resources from the context. This is a new
 * framework that will evolve in the contexts of
 * https://redmine.gromacs.org/issues/2375
 * https://redmine.gromacs.org/issues/2587
 * https://redmine.gromacs.org/issues/2605
 *
 * Ownership is returned as a smart pointer because ownership or reference may need to
 * be transferred through higher level code or maintained as a heap object,
 * and there is not yet a public interface, wrapper, or separate handle object.
 *
 * This call signature sets all parameters processed by the command-line client
 * code. Additional call signatures can allow implementation-specified default
 * values or newer initialization interfaces.
 *
 * \param simulationCommunicator Handle to communication data structure.
 *
 * \ingroup module_mdrun
 */
SimulationContext createSimulationContext(CommrecHandle simulationCommunicator);
//! \endcond

}      // end namespace gmx

#endif //GROMACS_SIMULATIONCONTEXT_H
