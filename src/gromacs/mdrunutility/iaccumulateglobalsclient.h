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
 * \brief Declares an interface intended to be implemented by
 * components that wish to accumulate simulation variables across PP
 * ranks during an MD step.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_IACCUMULATEGLOBALSCLIENT_H
#define GMX_MDRUNUTILITY_IACCUMULATEGLOBALSCLIENT_H

#include "gromacs/utility/arrayref.h"

namespace gmx
{

class AccumulateGlobals;

/*! \libinternal
 * \brief Interface for classes of objects that need to register with
 * the AccumulateGlobalsBuilder that they have interest in
 * accumulating values across PP ranks.
 *
 * All the registered client objects of this type will be queried by
 * \c getNumGlobalsRequired() for the number of values they wish to
 * accumulate. That number could vary when different modules or their
 * options are active. The AccumulateGlobalsBuilder will allocate
 * enough contiguous memory for all clients, which is thus efficient
 * for in-place MPI reduction.
 *
 * In principle, one component could have responsibility for computing
 * the inputs for globals and another component for using them after
 * reduction, but so far the design does not cater for that
 * possibility.
 *
 * In debug mode, the clients are notified after communication
 * completes, so that correctness can be verified.
 *
 * \todo Once all/enough modules have been converted to use this
 * framework, then AccumulateGlobals should take the responsibility
 * for doing the reduction.
 */
class IAccumulateGlobalsClient
{
    protected:
        virtual ~IAccumulateGlobalsClient()                   = 0;
    public:
        //! Return the number of values to reduce required by this module.
        virtual int getNumGlobalsRequired() const             = 0;
        //! Called to notify this module where to write and later read the values for reduction.
        virtual void setViewForGlobals(AccumulateGlobals *accumulateGlobals,
                                       ArrayRef<double>   view) = 0;
        //! Called (in debug mode) after MPI reduction is complete.
        virtual void notifyAfterCommunication()               = 0;
};

} // namespace gmx

#endif
