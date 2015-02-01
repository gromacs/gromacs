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
/*! \defgroup module_mdrunutility Implementation of mdrun utility functionality
 * \ingroup group_mdrun
 *
 * \brief This module contains code that implements general
 * infrastructure for mdrun that does not suit any other module.
 */
/*! \libinternal \file
 * \brief This file declares a function and class used by mdrun to
 * manage the details of doing restarts (ie. reading checkpoints,
 * whether it can or should append output files).
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 * \inlibraryapi
 */

#ifndef GMX_MDRUNUTILITY_HANDLERESTART_H
#define GMX_MDRUNUTILITY_HANDLERESTART_H

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/types/commrec.h"

namespace gmx
{

/*! \internal
 * \brief POD class that reports how mdrun will manage restarting
 * and appending. */
class RestartInformation
{
    public:
        //! Whether mdrun will start from the -cpi file
        bool bWillStartFromCpt_;
        //! Whether mdrun will append to files
        bool bWillAppendFiles_;
};

/*! \brief Handle mdrun restart from checkpoint
 *
 * Use an input checkpoint file only if it is valid (and consistent
 * among all simulations) and the names of output files used in the
 * run that generated it are consistent with those found on disk and
 * present on the mdrun command line. If a checkpoint file is found,
 * and all other conditions are consistent, use it for the restart,
 * otherwise give a fatal error. If a checkpoint file is not found,
 * proceed normally.
 *
 * Append to output files only if an input checkpoint file is used,
 * the previous job part didn't number the output files with part
 * numbers, and the user requested appending.
 *
 * Give a fatal error if members of a multi-simulation are not
 * starting from the same part.
 *
 * Does communication to coordinate behaviour between all ranks of a
 * simulation, and/or simulations.
 *
 * \param[in]    nfile              Size of fnm struct
 * \param[inout] fnm                Filename parameters to mdrun
 * \param[in]    cr                 Communication structure
 * \param[in]    bTryToAppendFiles  Whether mdrun -append was used (including by default)
 *
 * \throws std::bad_alloc          if out of memory
 * \throws APIError                if constructor was passed empty \p fnm
 * \throws InternalError           if the contents of the list of output filenames in the checkpoint file do not conform to expectations
 * \throws InconsistentInputError  if the output files on disk are an empty or only partial match for the names in the checkpoint, or if checkpoint file is missing and a log file exists
 * \throws FileIOError             if checkpoint file exists but cannot be read
 *
 * \returns A RestartInfo object whose members describe how mdrun will
 * operate. If mdrun will be doing a restart from checkpoint without
 * appending, \p fnm is updated with suffix strings for part numbers.
 */
RestartInformation
handleRestart(const int  nfile,
              t_filenm   fnm[],
              t_commrec *cr,
              bool       bTryToAppendFiles);

} // namespace

#endif
