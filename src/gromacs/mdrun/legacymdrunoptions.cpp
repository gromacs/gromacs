/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include <cstdlib>
#include <cstring>

#include <filesystem>

#include "gromacs/fileio/oenv.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

/*! \brief Return whether the command-line parameter that
 *  will trigger a multi-simulation is set */
static bool is_multisim_option_set(int argc, const char* const argv[])
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

int LegacyMdrunOptions::updateFromCommandLine(int argc, char** argv, ArrayRef<const char*> desc)
{
    unsigned long PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multidir, the working directory still needs to be
    // changed, so we can't check for the existence of files during
    // parsing.  It isn't useful to do any completion based on file
    // system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    if (!parse_common_args(&argc,
                           argv,
                           PCA_Flags,
                           gmx::ssize(filenames),
                           filenames.data(),
                           asize(pa),
                           pa,
                           ssize(desc),
                           desc.data(),
                           0,
                           nullptr,
                           &oenv))
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
        // command-line options have been separated.
        hw_opt.devicesSelectedByUser = devicesSelectedByUser;
        hw_opt.userGpuTaskAssignment = userGpuTaskAssignment;

        const char* env = getenv("GMX_GPU_ID");
        if (env != nullptr)
        {
            if (!hw_opt.devicesSelectedByUser.empty())
            {
                gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
            }
            hw_opt.devicesSelectedByUser = env;
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

        if (!hw_opt.devicesSelectedByUser.empty() && !hw_opt.userGpuTaskAssignment.empty())
        {
            gmx_fatal(FARGS, "-gpu_id and -gputasks cannot be used at the same time");
        }
    }

    hw_opt.threadAffinity = static_cast<ThreadAffinity>(nenum(thread_aff_opt_choices));

    if (!opt2parg_bSet("-append", asize(pa), pa))
    {
        mdrunOptions.appendingBehavior = AppendingBehavior::Auto;
    }
    else
    {
        if (opt2parg_bool("-append", asize(pa), pa))
        {
            mdrunOptions.appendingBehavior = AppendingBehavior::Appending;
        }
        else
        {
            mdrunOptions.appendingBehavior = AppendingBehavior::NoAppending;
        }
    }

    mdrunOptions.rerun            = opt2bSet("-rerun", gmx::ssize(filenames), filenames.data());
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
    output_env_done(oenv);
}

} // namespace gmx
