/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Defines the interface for external pull potentials and registers them.
 *
 */

#ifndef GMX_PULLING_PULL_POTENTIAL_H
#define GMX_PULLING_PULL_POTENTIAL_H

#include "gmxpre.h"

#include <map>
#include <memory>
#include <string>

namespace gmx
{

/*! \internal \brief
 * Abstract base class for pull potential function types.
 *
 * \internal Derived classes are registered in gmx::pull_potentials.
 */
class PullPotentialFunction
{

    public:

        virtual ~PullPotentialFunction();

        /*! \internal \brief
         * Function type for an external pull potential function
         *
         * \param[in]  coord_ind The pull coordinate index
         * \param[in]  x         The pull coordinate value
         * \param[out] V         The pull potential
         * \param[out] f         The pull force
         */
        virtual void do_potential(int coord_ind, double x, double *V, double *f) = 0;

    protected:
        /*! \internal \brief
         * Private constructor ensures porper creation via the registered create routine
         */
        PullPotentialFunction();
};

/*! \internal \brief
 * Container for the registered pull potentials.
 */
struct registry{

    /*! \internal \brief
     * Here, the external pull potentials will be registered, e.g. {{awh::name, awh::create},{another::name, another::create}}
     */
    static const std::map<std::string, gmx::PullPotentialFunction *(*)() > pull_potentials;
};

} // namespace gmx

#endif
