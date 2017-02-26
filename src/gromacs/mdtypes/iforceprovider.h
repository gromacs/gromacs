/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * Declares gmx::IForceProvider.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IFORCEPROVIDER_H
#define GMX_MDTYPES_IFORCEPROVIDER_H

#include "gromacs/math/paddedvector.h"

struct t_commrec;
struct t_forcerec;
struct t_mdatoms;

/*! \libinternal \brief
 * Interface for a component that provides forces during MD.
 *
 * This is typically part of a larger structure/class managing its own
 * data, such that it has the information on what to do stored locally.
 *
 * The interface is not very generic, as it has been written purely based on
 * extraction of existing functions related to electric field handling.
 * This needs to be generalized when more modules are moved to use the
 * interface.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
struct IForceProvider
{
    public:
        /*! \brief
         * Sets relevant options in the forcerec structure.
         *
         * \param[inout] fr The forcerec structure
         *
         * \todo
         * This should be replaced by a method that returns a set of
         * flags/other options (either here, or where the IForceProvider
         * instance is returned), and forcerec should be initialized based on
         * that.
         */
        virtual void initForcerec(t_forcerec *fr) = 0;

        /*! \brief
         * Computes forces.
         *
         * \param[in]    cr      Communication record for parallel operations
         * \param[in]    mdatoms Atom information
         * \param[inout] force   The forces
         * \param[in]    t       The actual time in the simulation (ps)
         */
        virtual void calculateForces(const t_commrec  *cr,
                                     const t_mdatoms  *mdatoms,
                                     PaddedRVecVector *force,
                                     double            t) = 0;

    protected:
        ~IForceProvider() {}
};

#endif
