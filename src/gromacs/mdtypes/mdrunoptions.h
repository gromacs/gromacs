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
/*! \libinternal \file
 *
 * \brief This file declares helper functionality for legacy option handling for mdrun
 *
 * It is likely that much of this content will move closer to the
 * functionality that supports the respective features. For example,
 * modules that change behaviour according to whether it is a rerun
 * could register themselves with the rerun module and get notified at
 * setup time to set their own boolean, rather than rely on a central
 * glob of mdrun options being passed around.
 *
 * \ingroup module_mdtypes
 * \inlibraryapi
 */
#ifndef GMX_MDTYPES_MDRUNOPTIONS_H
#define GMX_MDTYPES_MDRUNOPTIONS_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

//! Enumeration for mdrun appending behavior
enum class AppendingBehavior
{
    //! Append only if user command-line and file input is correct
    Auto,
    //! Must append
    Appending,
    //! Must not append
    NoAppending
};

//! \internal \brief Options for writing checkpoint files
struct CheckpointOptions
{
    //! True means keep all checkpoint file and add the step number to the name
    gmx_bool keepAndNumberCheckpointFiles = FALSE;
    //! The period in minutes for writing checkpoint files
    real period = 15;
};

//! \internal \brief Options for timing (parts of) mdrun
struct TimingOptions
{
    //! Reset timers at the start of this MD step, -1 means do not reset
    int resetStep = -1;
    //! If true, reset timers half-way the run
    gmx_bool resetHalfway = FALSE;
};

//! \internal \brief Options for IMD
struct ImdOptions
{
    //! IMD listening port
    int port = 8888;
    //! If true, pause the simulation while no IMD client is connected
    gmx_bool wait = FALSE;
    //! If true, allow termination of the simulation from IMD client
    gmx_bool terminatable = FALSE;
    //! If true, allow COM pulling in the simulation from IMD client
    gmx_bool pull = FALSE;
};

//! \internal \brief Collection of all options of mdrun that are not processed separately
struct MdrunOptions
{
    //! Re-compute energies, and possibly forces, for frames from an input trajectory
    gmx_bool rerun = FALSE;
    //! Re-construct virtual sites durin a rerun simulation
    gmx_bool rerunConstructVsites = FALSE;
    //! Try to make the simulation binary reproducible
    gmx_bool reproducible = FALSE;
    //! Write confout.gro at the end of the run
    gmx_bool writeConfout = TRUE;
    //! User option for appending.
    AppendingBehavior appendingBehavior = AppendingBehavior::Auto;
    //! Options for checkpointing th simulation
    CheckpointOptions checkpointOptions;
    //! Number of steps to run, -2 is use inputrec, -1 is infinite
    int64_t numStepsCommandline = -2;
    //! Maximum duration of this simulation in wall-clock hours, -1 is no limit
    real maximumHoursToRun = -1;
    //! Options for timing the run
    TimingOptions timingOptions;
    //! If true and supported, will tune the PP-PME load balance
    gmx_bool tunePme = TRUE;
    //! True if the user explicitly set the -ntomp command line option
    gmx_bool ntompOptionIsSet = FALSE;
    //! Options for IMD
    ImdOptions imdOptions;
    //! Increase the verbosity level in the logging and/or stdout/stderr
    gmx_bool verbose = FALSE;
    //! If verbose=true, print remaining runtime at this step interval
    int verboseStepPrintInterval = 100;
};

} // end namespace gmx

#endif
