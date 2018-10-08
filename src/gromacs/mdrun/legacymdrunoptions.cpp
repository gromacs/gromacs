/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 *
 * \brief This file declares helper functionality for legacy option handling for mdrun
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "legacymdrunoptions.h"

#include "config.h"

#include <cstring>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdrun/multisim.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

/*! \brief Return whether the command-line parameter that
 *  will trigger a multi-simulation is set */
static bool is_multisim_option_set(int argc, const char *const argv[])
{
    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "-multidir") == 0)
        {
            return true;
        }
    }
    return false;
}

int LegacyMdrunOptions::updateFromCommandLine(int argc, char **argv, ArrayRef<const char *> desc)
{
    unsigned long     PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multidir, the working directory still needs to be
    // changed, so we can't check for the existence of files during
    // parsing.  It isn't useful to do any completion based on file
    // system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    if (!parse_common_args(&argc, argv, PCA_Flags,
                           static_cast<int>(filenames.size()), filenames.data(), asize(pa), pa,
                           static_cast<int>(desc.size()), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }

    // Handle the options that permits the user to either declare
    // which compatible GPUs are availble for use, or to select a GPU
    // task assignment. Either could be in an environment variable (so
    // that there is a way to customize it, when using MPI in
    // heterogeneous contexts).
    {
        // TODO Argument parsing can't handle std::string. We should
        // fix that by changing the parsing, once more of the roles of
        // handling, validating and implementing defaults for user
        // command-line options have been seperated.
        hw_opt.gpuIdsAvailable       = gpuIdsAvailable;
        hw_opt.userGpuTaskAssignment = userGpuTaskAssignment;

        const char *env = getenv("GMX_GPU_ID");
        if (env != nullptr)
        {
            if (!hw_opt.gpuIdsAvailable.empty())
            {
                gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
            }
            hw_opt.gpuIdsAvailable = env;
        }

        env = getenv("GMX_GPUTASKS");
        if (env != nullptr)
        {
            if (!hw_opt.userGpuTaskAssignment.empty())
            {
                gmx_fatal(FARGS, "GMX_GPUTASKS and -gputasks can not be used at the same time");
            }
            hw_opt.userGpuTaskAssignment = env;
        }

        if (!hw_opt.gpuIdsAvailable.empty() && !hw_opt.userGpuTaskAssignment.empty())
        {
            gmx_fatal(FARGS, "-gpu_id and -gputasks cannot be used at the same time");
        }
    }

    hw_opt.thread_affinity = nenum(thread_aff_opt_choices);

    // now check for a multi-simulation
    ArrayRef<const std::string> multidir = opt2fnsIfOptionSet("-multidir",
                                                              static_cast<int>(filenames.size()),
                                                              filenames.data());

    if (replExParams.exchangeInterval != 0 && multidir.size() < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (use option -multidir)");
    }

    if (replExParams.numExchanges < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    ms = init_multisystem(MPI_COMM_WORLD, multidir);

    /* Prepare the intra-simulation communication */
    // TODO consolidate this with init_commrec, after changing the
    // relative ordering of init_commrec and init_multisystem
#if GMX_MPI
    if (ms != nullptr)
    {
        cr->nnodes = cr->nnodes / ms->nsim;
        MPI_Comm_split(MPI_COMM_WORLD, ms->sim, cr->sim_nodeid, &cr->mpi_comm_mysim);
        cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
        MPI_Comm_rank(cr->mpi_comm_mysim, &cr->sim_nodeid);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
    }
#endif

    if (!opt2bSet("-cpi",
                  static_cast<int>(filenames.size()), filenames.data()))
    {
        // If we are not starting from a checkpoint we never allow files to be appended
        // to, since that has caused a ton of strange behaviour and bugs in the past.
        if (opt2parg_bSet("-append", asize(pa), pa))
        {
            // If the user explicitly used the -append option, explain that it is not possible.
            gmx_fatal(FARGS, "GROMACS can only append to files when restarting from a checkpoint.");
        }
        else
        {
            // If the user did not say anything explicit, just disable appending.
            bTryToAppendFiles = FALSE;
        }
    }

    ContinuationOptions &continuationOptions = mdrunOptions.continuationOptions;

    continuationOptions.appendFilesOptionSet = opt2parg_bSet("-append", asize(pa), pa);

    handleRestart(cr, ms, bTryToAppendFiles,
                  static_cast<int>(filenames.size()),
                  filenames.data(),
                  &continuationOptions.appendFiles,
                  &continuationOptions.startedFromCheckpoint);

    mdrunOptions.rerun            = opt2bSet("-rerun",
                                             static_cast<int>(filenames.size()),
                                             filenames.data());
    mdrunOptions.ntompOptionIsSet = opt2parg_bSet("-ntomp", asize(pa), pa);

    domdecOptions.rankOrder    = static_cast<DdRankOrder>(nenum(ddrank_opt_choices));
    domdecOptions.dlbOption    = static_cast<DlbOption>(nenum(dddlb_opt_choices));
    domdecOptions.numCells[XX] = roundToInt(realddxyz[XX]);
    domdecOptions.numCells[YY] = roundToInt(realddxyz[YY]);
    domdecOptions.numCells[ZZ] = roundToInt(realddxyz[ZZ]);

    return 1;
}

LegacyMdrunOptions::~LegacyMdrunOptions()
{
    if (GMX_LIB_MPI)
    {
        done_commrec(cr);
    }
    done_multisim(ms);
}

} // namespace gmx
