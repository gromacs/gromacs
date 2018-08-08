/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

/*! \libinternal
 * \ingroup group_mdrun
 * \brief
 * Implements a part of a hybrid Monte Carlo integrator for md-vv.
 *
 * \author Sebastian Wingbermuehle
 */

/*! \libinternal \file
 *
 * \brief
 * Declares the interfaces through which the Metropolis step class is seen by
 * the classes AcceptOrRewind and HybridMCMDVelocities.
 *
 * \author Sebastian Wingbermuehle
 * \inlibraryapi
 * \ingroup module_hybridMCMD
 */

#ifndef GMX_METROPOLIS_INTERFACES_H
#define GMX_METROPOLIS_INTERFACES_H

#include "gromacs/utility/futil.h"

struct gmx_enerdata_t;

namespace gmx
{

class IMetropolisStep
{
    public:
        /*! \brief Decide to accept or reject the new configuration */
        virtual bool accept(const int64_t step, const gmx_enerdata_t *enerd) = 0;

        /*! \brief Destructor */
        virtual ~IMetropolisStep() {}
};

class IMetropolisStepVelocities
{
    public:
        /*! \brief Set the initial kinetic energy to be used in the Metropolis criterion */
        virtual void setInitialKineticEnergy(const double initialKineticEnergy) = 0;

        /*! \brief Destructor */
        virtual ~IMetropolisStepVelocities() {}
};

} // namespace

#endif
