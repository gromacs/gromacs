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
 * have its force correlation data checkpointed.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_CORRELATIONHISTORY_H
#define GMX_AWH_CORRELATIONHISTORY_H

struct t_commrec;

namespace gmx
{
class CorrelationGrid;
struct CorrelationGridHistory;

/*! \brief
 * Allocate a correlation grid history with the same structure as the given correlation grid.
 *
 * This function would be called at the start of a new simulation.
 * Note that only sizes and memory are initialized here.
 * History data is set by \ref updateCorrelationGridHistory.
 *
 * \param[in,out] corrGrid      Correlation grid state to initialize with.
 * \returns the correlation grid history struct.
 */
CorrelationGridHistory initCorrelationGridHistoryFromState(const CorrelationGrid &corrGrid);

/*! \brief
 * Restores the correlation grid state from the correlation grid history.
 *
 * \param[in] corrGridHist  Correlation grid history to read.
 * \param[in,out] corrGrid  Correlation grid state to set.
 */
void restoreCorrelationGridStateFromHistory(const CorrelationGridHistory &corrGridHist, CorrelationGrid *corrGrid);

/*! \brief
 * Update the correlation grid history for checkpointing.
 *
 * \param[in,out] corrGridHist  Correlation grid history to set.
 * \param[in]     corrGrid      Correlation grid state to read.
 */
void updateCorrelationGridHistory(CorrelationGridHistory *corrGridHist, const CorrelationGrid &corrGrid);

}      // namespace gmx

#endif /* GMX_AWH_CORRELATIONHISTORY_H */
