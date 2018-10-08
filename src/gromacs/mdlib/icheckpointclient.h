/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Declares and ICheckpointClient interface.
 *
 * Any module that contributes to checkpointing needs to implement this
 * interface.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */

#ifndef GROMACS_ICHECKPOINTCLIENT_H
#define GROMACS_ICHECKPOINTCLIENT_H

#include <map>

#include "gromacs/utility/real.h"

namespace gmx
{
template <typename T> class ArrayRef;

/*! \brief Checkpointing keyword
 *
 * This enum allows modules to register for reading or writing a specific set of data. The
 * enum is saved along the checkpointed data. Only one client can read or write for each
 * enum entry. If a checkpoint containing a specific enum is read, it is passed to the
 * module which registered to that enum (if no module registered, an error is thrown).
 * This also uniquely defines the order in which checkpoints are written or read, so changing
 * the order requires a bump in checkpointing version.
 */
enum class CheckpointKeyword
{
    state, ekinstate, swaphist,
    enerhist, dfhist, edsamhist,
    correlationGridHistory, awhBiasHistory, awhHistory,
    count
};

class ICheckpointClient
{
    protected:
        virtual ~ICheckpointClient() = default;

    public:
        virtual CheckpointKeyword getKeyword() = 0;
        virtual int getVersion()            = 0;

        virtual ArrayRef<int> getIntView() = 0;
        virtual ArrayRef<int64_t> getInt64View() = 0;
        virtual ArrayRef<real> getRealView() = 0;
        virtual ArrayRef<double> getDoubleView() = 0;

        virtual void notifyRead()  = 0;
        virtual void notifyWrite() = 0;
};

}      // namespace gmx

#endif //GROMACS_ICHECKPOINTCLIENT_H
