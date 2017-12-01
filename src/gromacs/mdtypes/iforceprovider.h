/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"

struct gmx_enerdata_t;
struct t_commrec;
struct t_forcerec;
struct t_mdatoms;

namespace gmx
{

template <typename T>
class ArrayRef;
class ForceWithVirial;


/*! \libinternal \brief
 * Helper struct for force provider input data
 *
 * This class assembles all input data that any of the force providers need.
 * For a particular force provider, one just needs to call the set...() methods
 * for the data this force provider expects. If in the future a new force
 * provider needs another input variable, one just adds it to the
 * public variables list and implements the corresponding set...() method.
 */
class ForceProviderInput
{
    public:
        ForceProviderInput()
            : cr(nullptr), mdatoms(nullptr), t(0.0), x(nullptr)
        {  };

        //! \brief Stores a pointer to the commrec struct in ForceProviderInput
        ForceProviderInput &setCommrec  (const t_commrec *crIn)      { cr      = crIn;      return *this; }

        //! \brief Stores a pointer to the mdatoms struct in ForceProviderInput
        ForceProviderInput &setMdatoms  (const t_mdatoms *mdatomsIn) { mdatoms = mdatomsIn; return *this; }

        //! \brief Stores a pointer to the box in ForceProviderInput
        ForceProviderInput &setBox      (const matrix     boxIn)     { copy_mat(boxIn, box); return *this; }

        //! \brief Copies the current time of the simulation into ForceProviderInput
        ForceProviderInput &setTime     (double           tIn)       { t       = tIn;       return *this; }

        //! \brief Stores a pointer to the atomic positions in ForceProviderInput
        ForceProviderInput &setPositions(const real      *xIn)       { x       = xIn;       return *this; }

        const t_commrec      *cr;         //!< communication record structure
        const t_mdatoms      *mdatoms;    //!< structure containing atomic data
        matrix                box;        //!< the simulation box
        double                t;          //!< the current time (ps) in the simulation
        const real           *x;          //!< the atomic positions
};

/*! \libinternal \brief
 * Helper struct assembling the output data of a force provider
 *
 * Same as for the ForceProviderInput class, but these variables can be written as well.
 */
class ForceProviderOutput
{
    public:
        //! \brief As the reason for a force provider is to calculate forces, a pointer to these is requested in the constructor
        ForceProviderOutput(ForceWithVirial *forceWithVirialIn)
            : forceWithVirial(forceWithVirialIn), enerd(nullptr)
        { };

        //! \brief Stores a pointer to the energy data struct in ForceProviderInput
        ForceProviderOutput &setEnerdata(gmx_enerdata_t *enerdIn) { enerd = enerdIn; return *this; }

        ForceWithVirial *forceWithVirial; //!< Container for force and virial
        gmx_enerdata_t  *enerd;           //!< structure containing energy data
};


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
         * \param[in]    forceProviderInput    struct that collects input data for the force providers
         * \param[inout] forceProviderOutput   struct that collects output data of the force providers
         */
        virtual void calculateForces(const ForceProviderInput &forceProviderInput,
                                     ForceProviderOutput      &forceProviderOutput) = 0;

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
        void calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                             gmx::ForceProviderOutput      &forceProviderOutput) const;

    private:
        class Impl;

        gmx::PrivateImplPointer<Impl> impl_;
};

#endif
