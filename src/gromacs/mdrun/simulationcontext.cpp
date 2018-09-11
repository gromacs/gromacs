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

#include "gmxpre.h"

#include "simulationcontext.h"

#include <cassert>

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/hardware/hw_info.h"

#include "runner.h"

namespace gmx
{
//! \cond libapi
SimulationContext::SimulationContext(t_commrec         * communicatonRecord,
                                     const gmx_hw_opt_t &hardwareOptions,
                                     const MdFilenames  &filenames,
                                     gmx_output_env_t  * outputEnvironment,
                                     FILE             ** logFile) :
    communicationRecord_(communicatonRecord),
    hardwareOptions_(hardwareOptions),
    filenames_(filenames),
    outputEnvironment_(outputEnvironment),
    logFile_(logFile)
{
    // SimulationContext receives a valid, initialized communications record.
    // \todo Initialization of the communications record is an implementation detail
    // dependent on the build configuration of the GROMACS installation and on
    // the run time computing context, not the client code that owns this Context.
    assert(communicationRecord_);
    // \todo Responsibility for owning, opening and closing the log file should be consolidated.
    assert(logFile_);
    // \todo Clarify ownership and lifetime management for gmx_output_env_t
    assert(outputEnvironment_);
    // \todo determine an invariant to check or confirm that all gmx_hw_opt_t objects are valid
    // \todo determine an invariant to checkout or establish that all MdFilenames objects are valid
}

std::unique_ptr<SimulationContext>
createSimulationContext(t_commrec         * simulationCommunicator,
                        const gmx_hw_opt_t &hardwareOptions,
                        const MdFilenames  &filenames,
                        gmx_output_env_t  * outputEnvironment,
                        FILE             ** logFileHandle)
{
    auto context = gmx::compat::make_unique<SimulationContext>(simulationCommunicator,
                                                               hardwareOptions,
                                                               filenames,
                                                               outputEnvironment,
                                                               logFileHandle);
    return context;
}
//! \endcond
} // end namespace gmx
