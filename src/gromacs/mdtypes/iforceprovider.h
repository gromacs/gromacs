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
#include "gromacs/utility/basedefinitions.h"
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
 * Helper struct that bundles force provider input data
 *
 * This is a short-lived container that bundles up all necessary input data for the
 * force providers. Its only purpose is to allow calling forceProviders->calculateForces()
 * with just two arguments, one being the container for the input data,
 * the other the container for the output data.
 * For a particular force provider, one just needs to call the bundle...() methods
 * for the data this force provider expects. If in the future a new force
 * provider needs another input variable, one just adds it to the
 * public variables list and implements the corresponding bundle...() method.
 *
 * Both ForceProviderInput as well as ForceProviderOutput only package existing
 * data structs together for handing it over to calcluateForces(). Apart from the
 * POD entries they own nothing.
 */
class ForceProviderInput
{
    public:
        //! \brief Constructor takes no arguments, as each provider requires different input
        ForceProviderInput() {  };

        //! \brief Stores a pointer to the commrec struct in ForceProviderInput
        ForceProviderInput &bundleCommrecPtr  (const t_commrec *crIn)      { cr      = crIn;      return *this; }

        //! \brief Stores a pointer to the mdatoms struct in ForceProviderInput
        ForceProviderInput &bundleMdatomsPtr  (const t_mdatoms *mdatomsIn) { mdatoms = mdatomsIn; return *this; }

        //! \brief Stores a pointer to the box in ForceProviderInput
        ForceProviderInput &bundleBox         (const matrix     boxIn)     { copy_mat(boxIn, box); return *this; }

        //! \brief Copies the current time of the simulation into ForceProviderInput
        ForceProviderInput &bundleTime        (double           tIn)       { t       = tIn;       return *this; }

        //! \brief Stores a pointer to the atomic positions in ForceProviderInput
        ForceProviderInput &bundlePositionsPtr(const rvec      *xIn)       { x       = xIn;       return *this; }

        const t_commrec      *cr      = nullptr;                           //!< communication record structure
        const t_mdatoms      *mdatoms = nullptr;                           //!< structure containing atomic data
        const rvec           *x       = nullptr;                           //!< the atomic positions
        gmx_int64_t           step    = 0;                                 //!< the simulation time step
        double                t       = 0.0;                               //!< the current time (ps) in the simulation
        matrix                box     = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; //!< the simulation box
};

/*! \libinternal \brief
 * Helper struct bundling the output data of a force provider
 *
 * Same as for the ForceProviderInput class, but these variables can be written as well.
 * Since the whole reason of a force provider is to calculate forces, a pointer to the
 * forces struct is explicitly requested in the constructor. Other data that is only
 * written by specific force providers is optional and therefore provided as needed by
 * the bundle...() methods.
 */
class ForceProviderOutput
{
    public:
        //! \brief Explicitly request a storage location for the computed forces in the constructor!
        explicit ForceProviderOutput(ForceWithVirial *forceWithVirialIn)
            : forceWithVirial(forceWithVirialIn)
        { };

        //! \brief Stores a pointer to the energy data struct in ForceProviderInput
        ForceProviderOutput &bundleEnerdataPtr(gmx_enerdata_t *enerdIn) { enerd = enerdIn; return *this; }

        ForceWithVirial *forceWithVirial;                                  //!< Container for force and virial
        gmx_enerdata_t  *enerd        = nullptr;                           //!< structure containing energy data
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
                                     ForceProviderOutput      *forceProviderOutput) = 0;

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
                             gmx::ForceProviderOutput      *forceProviderOutput) const;

    private:
        class Impl;

        gmx::PrivateImplPointer<Impl> impl_;
};

#endif
