/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief This file defines high-level functions for managing doing
 * mdrun restarts.
 *
 * This includes reading checkpoints and deciding whether and how to
 * append output files.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */

#include "gmxpre.h"

#include "handlerestart.h"

#include <string>

#include "gromacs/fileio/fileiohandlerinterface.h"
#include "gromacs/utility/mpihandlerinterface.h"
#include "gromacs/utility/stringutil.h"

#include "handlerestart-impl.h"

namespace gmx
{

/* This routine cannot print tons of data, since it is called before
   the log file is opened. */
RestartInformation
handleRestart(FileIOHandlerInterface *inputHandler,
              MpiHandlerInterface    *mpiHandler,
              CommandLineFilenames    mdrunFilenames,
              bool                    bTryToAppendFiles)
{
    // Construct and use object TODO
    HandleRestart::Impl impl(inputHandler, mpiHandler, constArrayRefFromArrayRef(mdrunFilenames));
    //    RestartInformation info = impl.getRestartInformation(bTryToAppendFiles);

    HandleRestart::DataFromCheckpoint dataFromCheckpoint =
        impl.getDataFromCheckpoint(mpiHandler->isSimMaster(), mpiHandler->isMultiMaster());
    RestartInformation                info;
    info.bAppendFiles_ = impl.canAppend(dataFromCheckpoint,
                                        bTryToAppendFiles);

    /* Communicate results of checks */
    if (mpiHandler->hasMultipleRanks())
    {
        mpiHandler->broadcast(sizeof(dataFromCheckpoint.partNumber_),
                              &dataFromCheckpoint.partNumber_);

        if (dataFromCheckpoint.partNumber_ > 0 && bTryToAppendFiles)
        {
            /* In this case, we need to coordinate whether we'll do
               appending. Otherwise, the default for info.bAppendFiles_
               will do. */
            mpiHandler->broadcast(sizeof(info.bAppendFiles_),
                                  &info.bAppendFiles_);
        }
    }
    // Finish filling info
    info.simulationPartNumber_ = dataFromCheckpoint.partNumber_;
    info.bStartFromCpt_        = info.simulationPartNumber_ > 0;


    /* Log file is not yet available, so if there's a problem we can
     * only write to stderr. */
    mpiHandler->checkAcrossMultiSim(stderr, dataFromCheckpoint.partNumber_,
                                    "simulation part", TRUE);

    /* Add suffixes to the output file names that index the simulation
       part, if appropriate. */
    if (info.bStartFromCpt_ && info.bAppendFiles_)
    {
        /* Rename all output files (except checkpoint files) using the
           suffix containing the new part number. The new part number
           is the one from the checkpoint file, plus one. */
        std::string suffix = formatString(".part%04d", info.simulationPartNumber_ + 1);
        add_suffix_to_output_names(mdrunFilenames.data(), mdrunFilenames.size(), suffix.c_str());
        if (mpiHandler->isMultiMaster())
        {
            fprintf(stdout, "Checkpoint file is from part %d, new output files will be suffixed '%s'.\n",
                    info.simulationPartNumber_, suffix.c_str());
        }
    }

    return info;
}

} // namespace
