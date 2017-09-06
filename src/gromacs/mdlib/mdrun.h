/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief This file declares types and functions for initializing an MD run
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_MDRUN_H
#define GMX_MDLIB_MDRUN_H

#include "gromacs/timing/wallcycle.h"

struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;
class t_state;

//! \internal \brief Options and settings for continuing from checkpoint
struct ContinuationOptions
{
    //! \brief Constructor
    ContinuationOptions() :
        appendFiles(false),
        appendFilesOptionSet(false),
        startedFromCheckpoint(false),
        haveReadEkin(false)
    {
    }

    //! True if we are continuing from a checkpoint and should append output files
    bool appendFiles;
    //! True if the -append option was explicitly set by the user (either to true of false
    bool appendFilesOptionSet;
    //! True if we started from a checkpoint file
    bool startedFromCheckpoint;
    //! True if we read the kinetic energy from checkpoint file
    bool haveReadEkin;
};

//! \internal \brief Options for writing checkpoint files
struct CheckpointOptions
{
    //! \brief Constructor
    CheckpointOptions() :
        keepAndNumberCheckpointFiles(FALSE),
        period(15)
    {
    }

    //! True means keep all checkpoint file and add the step number to the name
    gmx_bool keepAndNumberCheckpointFiles;
    //! The period in minutes for writing checkpoint files
    real     period;
};

//! \internal \brief Options for timing (parts of) mdrun
struct TimingOptions
{
    //! \brief Constructor
    TimingOptions() :
        resetStep(-1),
        resetHalfway(FALSE)
    {
    }

    //! Reset timers at the start of this MD step, -1 means do not reset
    int      resetStep;
    //! If true, reset timers half-way the run
    gmx_bool resetHalfway;
};

//! \internal \brief Options for IMD
struct ImdOptions
{
    //! Constructor
    ImdOptions() :
        port(8888),
        wait(FALSE),
        terminatable(FALSE),
        pull(FALSE)
    {
    }

    //! IMD listening port
    int      port;
    //! If true, pause the simulation while no IMD client is connected
    gmx_bool wait;
    //! If true, allow termination of the simulation from IMD client
    gmx_bool terminatable;
    //! If true, allow COM pulling in the simulation from IMD client
    gmx_bool pull;
};

//! \internal \brief Collection of all options of mdrun that are not processed separately
struct MdrunOptions
{
    //! \brief Constructor
    MdrunOptions() :
        rerun(FALSE),
        rerunConstructVsites(FALSE),
        globalCommunicationInterval(-1),
        reproducible(FALSE),
        writeConfout(TRUE),
        continuationOptions(),
        checkpointOptions(),
        numStepsCommandline(-2),
        maximumHoursToRun(-1),
        timingOptions(),
        tunePme(TRUE),
        ntompOptionIsSet(FALSE),
        imdOptions(),
        verbose(FALSE),
        verboseStepPrintInterval(100)
    {
    }

    //! Re-compute energies, and possibly forces, for frames from an input tracjectory
    gmx_bool            rerun;
    //! Re-construct virual sites durin a rerun simulation
    gmx_bool            rerunConstructVsites;
    //! Request to do global communication at this interval in steps, -1 is determine from inputrec
    int                 globalCommunicationInterval;
    //! Try to make the simulation binary reproducible
    gmx_bool            reproducible;
    //! Write confout.gro at the end of the run
    gmx_bool            writeConfout;
    //! Options for continuing a simulation from a checkpoint file
    ContinuationOptions continuationOptions;
    //! Options for checkpointing th simulation
    CheckpointOptions   checkpointOptions;
    //! Number of steps to run, -2 is use inputrec, -1 is infinite
    gmx_int64_t         numStepsCommandline;
    //! Maximum duration of this simulation in wall-clock hours, -1 is no limit
    real                maximumHoursToRun;
    //! Options for timing the run
    TimingOptions       timingOptions;
    //! If true and supported, will tune the PP-PME load balance
    gmx_bool            tunePme;
    //! True if the user explicitly set the -ntomp command line option
    gmx_bool            ntompOptionIsSet;
    //! Options for IMD
    ImdOptions          imdOptions;
    //! Increase the verbosity level in the logging and/or stdout/stderr
    gmx_bool            verbose;
    //! If verbose=true, print remaining runtime at this step interval
    int                 verboseStepPrintInterval;
};

//! \brief Allocate and initialize node-local state entries
void set_state_entries(t_state *state, const t_inputrec *ir);

//! \brief Broadcast inputrec and mtop and allocate node-specific settings
void init_parallel(t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop);

//! \brief Broadcasts the, non-dynamic, state from the master to all ranks in cr->mpi_comm_mygroup
//
// This is intended to be used with MPI parallelization without
// domain decompostion (currently with NM and TPI).
void broadcastStateWithoutDynamics(const t_commrec *cr, t_state *state);

#endif
