/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * This file contains the definition of a container for history data
 * for simulation observables.
 *
 * The container is used for storing the simulation state data that needs
 * to be written to / read from checkpoint file. This struct should only
 * contain pure observable data. Microstate data should be in t_state.
 * The state of the mdrun machinery is also stored elsewhere.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDLIB_OBSERVABLESHISTORY_H
#define GMX_MDLIB_OBSERVABLESHISTORY_H

#include <memory>

class energyhistory_t;
class PullHistory;
struct edsamhistory_t;
struct swaphistory_t;

/*! \libinternal \brief Observables history, for writing/reading to/from checkpoint file
 */
struct ObservablesHistory
{
    //! History for energy observables, used for output only
    std::unique_ptr<energyhistory_t> energyHistory;

    //! History for pulling observables, used for output only
    std::unique_ptr<PullHistory> pullHistory;

    //! Essential dynamics and flooding history
    std::unique_ptr<edsamhistory_t> edsamHistory;

    //! Ion/water position swapping history
    std::unique_ptr<swaphistory_t> swapHistory;

    ObservablesHistory();

    ~ObservablesHistory();
};

#endif
