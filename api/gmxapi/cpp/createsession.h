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

#ifndef GMXAPI_LIBRARY_SESSION_H
#define GMXAPI_LIBRARY_SESSION_H
/*! \file
 * \brief Library internal details for Session API class(es).
 *
 * \ingroup gmxapi
 */

#include "gromacs/mdrunutility/logging.h"

// Above are headers for dependencies.
// Following are public headers for the current module.
#include "gmxapi/context.h"
#include "gmxapi/session.h"

namespace gmx
{
class MdrunnerBuilder;
class SimulationContext;
} // end namespace gmx

namespace gmxapi
{

/*!
 * \brief Factory free function for creating new Session objects.
 *
 * \param context Shared ownership of a Context implementation instance.
 * \param runnerBuilder MD simulation builder to take ownership of.
 * \param simulationContext Take ownership of the simulation resources.
 * \param logFilehandle Take ownership of filehandle for MD logging
 * \param multiSim Take ownership of resources for Mdrunner multi-sim.
 *
 * \todo Log file management will be updated soon.
 *
 * \return Ownership of new Session implementation instance.
 */
std::shared_ptr<Session> createSession(std::shared_ptr<ContextImpl> context,
                                       gmx::MdrunnerBuilder&&       runnerBuilder,
                                       gmx::SimulationContext&&     simulationContext,
                                       gmx::LogFilePtr              logFilehandle);


} // end namespace gmxapi

#endif // GMXAPI_LIBRARY_SESSION_H
