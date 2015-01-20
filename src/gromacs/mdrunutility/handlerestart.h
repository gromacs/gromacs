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
 *
 * \brief This file declares functions for mdrun to call to manage the
 * details of doing a restart (ie. reading checkpoints, appending
 * output files).
 *
 * \todo There may be other code in runner.cpp etc. that can usefully
 * live here
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */

#ifndef GMX_MDRUNUTILITY_HANDLERESTART_H
#define GMX_MDRUNUTILITY_HANDLERESTART_H

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/basedefinitions.h"

/*! \brief Handle startup of mdrun, particularly regarding -cpi and -append
 *
 * If there is a checkpoint file, then prepare to start from that
 * state. If there is also appending, then organize the file naming
 * accordingly, if possible. Issues various fatal errors if the input
 * conditions are inconsistent or broken. \p fnm is updated with
 * suffix strings for part numbers if we are doing a restart from
 * checkpoint and are not appending.
 *
 * Does communication to coordinate behaviour between all ranks of a
 * simulation, and/or simulations.
 *
 * \param[in]    cr                 Communication structure
 * \param[in]    bTryToAppendFiles  Whether mdrun -append was used
 * \param[in]    NFILE              Size of fnm struct
 * \param[inout] fnm                Filename parameters to mdrun
 * \param[out]   bDoAppendFiles     Whether mdrun will append to files
 * \param[out]   bStartFromCpt      Whether mdrun will start from the -cpi file
 */
void handleRestart(t_commrec *cr,
                   gmx_bool   bTryToAppendFiles,
                   const int  NFILE,
                   t_filenm   fnm[],
                   gmx_bool  *bDoAppendFiles,
                   gmx_bool  *bStartFromCpt);

#endif
