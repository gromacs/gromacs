/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * GMX compatible interface to NBLIB listed forces
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTED_FORCES_GPU_INTERFACE_H
#define NBLIB_LISTED_FORCES_GPU_INTERFACE_H

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/timing/wallcycle.h"

#include "nblib/listed_forces/conversionscommon.h"
#include "nblib/listed_forces/gpu_param.cuh"
#include "nblib/pbc.hpp"

namespace gmx
{

/*! \internal \brief Implements GPU bondeds */
class ListedForcesNblibGpuImpl
{
public:
    //! Constructor
    ListedForcesNblibGpuImpl(const gmx_ffparams_t& ffparams,
                             float                 electrostaticsScaleFactor,
                             const DeviceContext&  deviceContext,
                             const DeviceStream&   deviceStream,
                             gmx_wallcycle*        wcycle);
    /*! \brief Destructor, non-default needed for freeing
     * device-side buffers */
    ~ListedForcesNblibGpuImpl();
    /*! \brief Update lists of interactions from idef suitable for the GPU,
     * using the data structures prepared for PP work.
     *
     * Intended to be called after each neighbour search
     * stage. Copies the bonded interactions assigned to the GPU
     * to device data structures, and updates device buffers that
     * may have been updated after search. */
    void updateInteractionListsAndDeviceBuffers(ArrayRef<const int>           nbnxnAtomOrder,
                                                const InteractionDefinitions& idef,
                                                void*                         xqDevice,
                                                DeviceBuffer<RVec>            forceDevice,
                                                DeviceBuffer<RVec>            fshiftDevice);
    /*! \brief
     * Update PBC data.
     *
     * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \param[in] pbcType The type of the periodic boundary.
     * \param[in] box     The periodic boundary box matrix.
     * \param[in] canMoleculeSpanPbc  Whether one molecule can have atoms in different PBC cells.
     */
    void setPbc(PbcType pbcType, const matrix box, bool canMoleculeSpanPbc);

    /*! \brief Launches bonded kernel on a GPU */
    template<bool calcVir, bool calcEner>
    void launchKernel();
    /*! \brief Returns whether there are bonded interactions
     * assigned to the GPU */
    bool haveInteractions() const;
    /*! \brief Launches the transfer of computed bonded energies. */
    void launchEnergyTransfer();
    /*! \brief Waits on the energy transfer, and accumulates bonded energies to \c enerd. */
    void waitAccumulateEnergyTerms(gmx_enerdata_t* enerd);
    /*! \brief Clears the device side energy buffer */
    void clearEnergies();

private:
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! \brief Bonded GPU stream, not owned by this module
    const DeviceStream& deviceStream_;

    nblib::PbcHolderAiuc            pbc_;
    nblib::ListedInteractionDataGpu interactions_;

    //! \brief user supplied buffers, only pointers stored here
    const util::array<real, 4>* d_xyzq_;
    util::array<real, 3>*       d_forces_;
    util::array<real, 3>*       d_shiftForces_;

    thrust::device_vector<real> d_potentials_;
    thrust::host_vector<real>   h_potentials_;

    //! GPU kernel launch configuration
    KernelLaunchConfig kernelLaunchConfig_;

    //! \brief Pointer to wallcycle structure.
    gmx_wallcycle* wcycle_;
};


} // namespace gmx

#endif // NBLIB_LISTED_FORCES_GPU_INTERFACE_H
