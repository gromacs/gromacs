/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * \brief
 * Contains datatypes and function declarations needed by AWH to
 * have its data checkpointed.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_HISTORY_H
#define GMX_AWH_HISTORY_H

#include "gromacs/utility/basedefinitions.h"

struct AwhBiasCollection;
struct awh_history_t;
struct t_commrec;


/*! \brief
 * Allocate and initialize an AWH history with the given AWH state.
 *
 * This function would be called at the start of a new simulation.
 * Note that no data only constant data will be initialized here.
 * History data is set by update_awh_history.
 *
 * \param[in,out] awh_history  AWH bias history to initialize.
 * \param[in] awh          AWH state to initialize with.
 */
void init_awh_history_from_state(awh_history_t *awh_history, const AwhBiasCollection *awh);

/*! \brief
 * Initialize an AWH history when restarting a from checkpoint.
 *
 * This function assumes that the master rank has read a checkpoint
 * and initialized its AWH history.
 *
 * \param[in,out] awh_history  AWH bias history to initialize.
 * \param[in]     cr           Communicator needed for broadcasting.
 */
void init_awh_history_from_checkpoint(awh_history_t    *awh_history,
                                      const t_commrec  *cr);

/*! \brief
 * Restores the AWH  state from history.
 *
 * \param[in] awh_history      AWH bias history to read.
 * \param[in,out] awh          AWH state to set.
 */
void restore_awh_state_from_history(const awh_history_t  *awh_history,
                                    AwhBiasCollection    *awh);

#endif
