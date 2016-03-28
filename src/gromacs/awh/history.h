/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * \brief
 * Contains datatypes and function declarations needed by AWH to
 * have its data checkpointed.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 */

#ifndef GMX_AWH_HISTORY_H
#define GMX_AWH_HISTORY_H

#include "gromacs/utility/basedefinitions.h"

struct t_awhbias;
struct awhbiashistory_t;
struct t_commrec;


/*! \brief
 * Initialize an AWH history struct trivially.
 *
 * This would be called if there is not an AWH working struct
 * with values to initialize with.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 */
void init_awhbiashistory(awhbiashistory_t *awhbiashist);

/*! \brief
 * Initialize an AWH history with the given AWH state.
 *
 * This function would be called at the start of a new simulation.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 * \param[in] awhbias          AWH state to initialize with.
 */
void init_awhbiashistory_from_state(awhbiashistory_t *awhbiashist, const t_awhbias *awhbias);

/*! \brief
 * Initialize an AWH history when restarting a from checkpoint.
 *
 * This function assumes that the master rank has read a checkpoint
 * and initialized its AWH history.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 * \param[in]     cr           Communicator needed for broadcasting.
 */
void init_awhbiashistory_from_checkpoint(awhbiashistory_t *awhbiashist, const t_commrec *cr);

/*! \brief
 * Restores the AWH bias state from the AWH history.
 *
 * \param[in] awhbiashist      AWH bias history to read.
 * \param[in,out] awhbias      AWH state to set.
 */
void restore_awhbias_state_from_history(const awhbiashistory_t       *awhbiashist,
                                        t_awhbias                    *awhbias);

#endif
