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

/*! \internal \file
 *
 * \brief
 * Contains datatypes and function declarations needed by AWH to
 * have its force correlation data checkpointed.
 *
 * \author Viveca Lindahl
 * \ingroup module_awh
 */

#ifndef GMX_AWH_CORRELATION_HISTORY_H
#define GMX_AWH_CORRELATION_HISTORY_H

struct correlation_grid_t;
struct correlation_grid_history_t;
struct t_commrec;

/*! \brief
 * Allocate a correlation grid history with the same structure as the given correlation grid.
 *
 * This function would be called at the start of a new simulation.
 * Note that no data only constant data will be initialized here.
 * History data is set by update_correlation_grid_history.
 *
 * \param[in,out] corrgrid      Correlation grid state to initialize with.
 * \returns the correlation grid history struct.
 */
correlation_grid_history_t *init_correlation_grid_history_from_state(const correlation_grid_t *corrgrid);

/*! \brief
 * Allocate a correlation grid history when restarting from a checkpoint.
 *
 * This function assumes that the master rank has read a checkpoint
 * and allocated and initialized its correlation grid history.
 *
 * \param[in,out] corrgrid_hist_in  Correlation grid history for master rank.
 * \param[in]     cr                Communicator needed for broadcasting.
 * \returns the correlation grid history struct for non-master and NULL for master rank.
 */
correlation_grid_history_t *init_correlation_grid_history_from_checkpoint(correlation_grid_history_t *corrgrid_hist_in, const t_commrec *cr);

/*! \brief
 * Restores the correlation grid state from the correlation grid history.
 *
 * \param[in] corrgrid_hist    Correlation grid history to read.
 * \param[in,out] corrgrid     Correlation grid state to set.
 */
void restore_correlation_grid_state_from_history(const correlation_grid_history_t *corrgrid_hist, correlation_grid_t *corrgrid);

/*! \brief
 * Update the correlation grid history for checkpointing.
 *
 * \param[in,out] corrgrid_hist   Correlation grid history to set.
 * \param[in]     corrgrid        Correlation grid state to read.
 */
void update_correlation_grid_history(correlation_grid_history_t *corrgrid_hist, const correlation_grid_t *corrgrid);

#endif
