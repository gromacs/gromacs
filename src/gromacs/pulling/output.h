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

/*! \libinternal \file
 *
 * \brief
 * This file declares functions for pull output writing.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_OUTPUT_H
#define GMX_PULLING_OUTPUT_H

#include <cstdint>

#include "gromacs/utility/basedefinitions.h"

struct pull_t;
struct gmx_output_env_t;
struct ObservablesHistory;
struct t_filenm;

namespace gmx
{
enum class StartingBehavior;
}

/*! \brief Set up and open the pull output files, when requested.
 *
 * NOTE: This should only be called on the main rank and only when
 *       doing dynamics (e.g. not with energy minimization).
 *
 * \param pull        The pull work data struct
 * \param nfile       Number of files.
 * \param fnm         Standard filename struct.
 * \param oenv        Output options.
 * \param startingBehavior  Describes whether this is a restart appending to output files
 */
void init_pull_output_files(pull_t*                 pull,
                            int                     nfile,
                            const t_filenm          fnm[],
                            const gmx_output_env_t* oenv,
                            gmx::StartingBehavior   startingBehavior);

/*! \brief Print the pull output (x and/or f)
 *
 * \param pull     The pull data structure.
 * \param step     Time step number.
 * \param time     Time.
 */
void pull_print_output(pull_t* pull, int64_t step, double time);

/*! \brief Allocate and initialize pull work history (for average pull output) and set it in a pull work struct
 *
 * \param pull                The pull work struct
 * \param observablesHistory  Container of history data, e.g., pull history.
 */
void initPullHistory(pull_t* pull, ObservablesHistory* observablesHistory);

#endif
