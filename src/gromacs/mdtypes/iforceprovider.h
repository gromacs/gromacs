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
 * Declares gmx::IForceProvider and ForceProviders.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IFORCEPROVIDER_H
#define GMX_MDTYPES_IFORCEPROVIDER_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"

struct t_commrec;
struct t_forcerec;
struct t_mdatoms;

namespace gmx
{

template <typename T>
class ArrayRef;
class ForceWithVirial;

/*! \libinternal \brief
 * Interface for a component that provides forces during MD.
 *
 * Modules implementing IMDModule generally implement this internally, and use
 * IMDModule::initForceProviders() to register their implementation in
 * ForceProviders.
 *
 * The interface most likely requires additional generalization for use in
 * other modules than the current electric field implementation.
 *
 * The forces that are produced by force providers are not taken into account
 * in the calculation of the virial. When applicable, the provider should
 * compute its own virial contribution.
 * \todo Extend this interface with a virial container and flag if the virial is needed here
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IForceProvider
{
    public:
        /*! \brief
         * Computes forces.
         *
         * \param[in]    cr               Communication record for parallel operations
         * \param[in]    mdatoms          Atom information
         * \param[in]    box              The box
         * \param[in]    t                The actual time in the simulation (ps)
         * \param[in]    x                The coordinates
         * \param[inout] forceWithVirial  The forces and virial
         */
        virtual void calculateForces(const t_commrec          *cr,
                                     const t_mdatoms          *mdatoms,
                                     const matrix              box,
                                     double                    t,
                                     const rvec               *x,
                                     gmx::ForceWithVirial     *forceWithVirial) = 0;

    protected:
        ~IForceProvider() {}
};

} // namespace gmx

/*! \libinternal \brief
 * Evaluates forces from a collection of gmx::IForceProvider.
 *
 * This class is a `struct` outside the `gmx` namespace to make it possible to
 * forward-declare it in forcerec.h, which still needs to compile when included
 * from the C group kernels.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
struct ForceProviders
{
    public:
        ForceProviders();
        ~ForceProviders();

        /*! \brief
         * Adds a provider.
         */
        void addForceProvider(gmx::IForceProvider *provider);

        //! Whether there are modules added.
        bool hasForceProvider() const;

        //! Computes forces.
        void calculateForces(const t_commrec          *cr,
                             const t_mdatoms          *mdatoms,
                             const matrix              box,
                             double                    t,
                             const rvec               *x,
                             gmx::ForceWithVirial     *forceWithVirial) const;

    private:
        class Impl;

        gmx::PrivateImplPointer<Impl> impl_;
};

#endif
