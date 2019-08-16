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
/*! \internal \file
 * \brief Declares GPU implementation class for CUDA bonded
 * interactions.
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_GPUBONDED_IMPL_H
#define GMX_LISTED_FORCES_GPUBONDED_IMPL_H

#include "gromacs/gpu_utils/gpu_vec.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/topology/idef.h"

struct gmx_ffparams_t;
struct t_forcerec;

namespace gmx
{

/*! \internal \brief Version of InteractionList that supports pinning */
struct HostInteractionList
{
    /*! \brief Returns the total number of elements in iatoms */
    int size() const
    {
        return iatoms.size();
    }

    //! List of interactions, see \c HostInteractionLists
    HostVector<int> iatoms = {{}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported)};
};

/*! \internal \brief Implements GPU bondeds */
class GpuBonded::Impl
{
    public:
        //! Constructor
        Impl(const gmx_ffparams_t &ffparams,
             void                 *streamPtr);
        /*! \brief Destructor, non-default needed for freeing
         * device-side buffers */
        ~Impl();
        /*! \brief Update lists of interactions from idef suitable for the GPU,
         * using the data structures prepared for PP work.
         *
         * Intended to be called after each neighbour search
         * stage. Copies the bonded interactions assigned to the GPU
         * to device data structures, and updates device buffers that
         * may have been updated after search. */
        void updateInteractionListsAndDeviceBuffers(ArrayRef<const int>  nbnxnAtomOrder,
                                                    const t_idef        &idef,
                                                    void                *xqDevice,
                                                    void                *forceDevice,
                                                    void                *fshiftDevice);

        /*! \brief Launches bonded kernels on a GPU */
        template <bool calcVir, bool calcEner>
        void
        launchKernels(const t_forcerec *fr,
                      const matrix      box);
        /*! \brief Returns whether there are bonded interactions
         * assigned to the GPU */
        bool haveInteractions() const;
        /*! \brief Launches the transfer of computed bonded energies. */
        void launchEnergyTransfer();
        /*! \brief Waits on the energy transfer, and accumulates bonded energies to \c enerd. */
        void accumulateEnergyTerms(gmx_enerdata_t *enerd);
        /*! \brief Clears the device side energy buffer */
        void clearEnergies();
    private:
        /*! \brief The interaction lists
         *
         * \todo This is potentially several pinned allocations, which
         * could contribute to exhausting such pages. */
        std::array<HostInteractionList, F_NRE> iLists;
        //! Tells whether there are any interaction in iLists.
        bool                                   haveInteractions_;
        //! Interaction lists on the device.
        t_ilist                                iListsDevice[F_NRE];

        //! Bonded parameters for device-side use.
        t_iparams            *forceparamsDevice = nullptr;
        //! Position-charge vector on the device.
        const float4         *xqDevice = nullptr;
        //! Force vector on the device.
        fvec                 *forceDevice = nullptr;
        //! Shift force vector on the device.
        fvec                 *fshiftDevice = nullptr;
        //! \brief Host-side virial buffer
        HostVector <float>    vtot = {{}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported)};
        //! \brief Device-side total virial
        float                *vtotDevice   = nullptr;

        //! \brief Bonded GPU stream, not owned by this module
        CommandStream         stream;
};

}   // namespace gmx

#endif
