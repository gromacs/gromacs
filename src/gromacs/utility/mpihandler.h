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
/*! \libinternal \file
 * \brief This file declares a class used for functionality that
 * requires MPI communication, so that such dependency may be
 * injected into modules.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */

#ifndef GMX_UTILITY_MPIHANDLER_H
#define GMX_UTILITY_MPIHANDLER_H

#include <stdio.h>

#include "gromacs/utility/mpihandlerinterface.h"
#include "gromacs/utility/uniqueptr.h"

//! Forward declaration
struct t_commrec;

namespace gmx
{

/*! \internal
 * \brief Wrapper used for dependency injection of methods used by
 * Impl for MPI-related behaviour of t_commrec.
 *
 * \todo Move other implementations of such functionality into this
 * class. */
class MpiHandler : public MpiHandlerInterface
{
    public:
        //! Constructor
        MpiHandler(const t_commrec *cr);
        //! Virtual destructor required
        virtual ~MpiHandler();
        //! Getter
        virtual bool hasMultipleRanks() const;
        //! Getter
        virtual bool isMaster() const;
        //! Getter
        virtual bool isSimMaster() const;
        //! Getter
        virtual bool isMultiSim() const;
        //! Getter
        virtual bool isMultiMaster() const;
        //! \copydoc MpiHandlerInterface::broadcast
        // TODO name this better
        virtual void broadcast(int byteSize, void *data) const;
        //! \copydoc MpiHandlerInterface::checkAcrossMultiSim
        virtual void checkAcrossMultiSim(FILE *fp, int theInteger, const char *description, bool bQuiet) const;

    private:
        const t_commrec *cr_;
};

//! Convenience typedef
typedef gmx_unique_ptr<MpiHandler>::type MpiHandlerPointer;

} // namespace

#endif
