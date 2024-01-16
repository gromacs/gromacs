/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares GPU implementation class for GPU bonded
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
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_GPU_IMPL_H
#define GMX_LISTED_FORCES_LISTED_FORCES_GPU_IMPL_H

#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/pbcutil/pbc_aiuc.h"

#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/syclutils.h"
#endif

struct gmx_ffparams_t;
class DeviceContext;

namespace gmx
{

/*! \internal \brief Version of InteractionList that supports pinning */
struct HostInteractionList
{
    /*! \brief Returns the total number of elements in iatoms */
    int size() const { return iatoms.size(); }

    //! List of interactions, see \c HostInteractionLists
    HostVector<int> iatoms = { {}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported) };
};

/* \brief Bonded parameters and GPU pointers
 *
 * This is used to accumulate all the parameters and pointers so they can be passed
 * to the GPU as a single structure.
 *
 */
struct BondedGpuKernelParameters
{
    //! Periodic boundary data
    PbcAiuc pbcAiuc;
    //! Scale factor
    float electrostaticsScaleFactor;
    //! The bonded types on GPU
    int fTypesOnGpu[numFTypesOnGpu];
    //! The number of bonds for every function type
    int numFTypeBonds[numFTypesOnGpu];
    //! The start index in the range of each interaction type
    int fTypeRangeStart[numFTypesOnGpu];
    //! The end index in the range of each interaction type
    int fTypeRangeEnd[numFTypesOnGpu];
    BondedGpuKernelParameters()
    {
        matrix boxDummy = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        setPbcAiuc(0, boxDummy, &pbcAiuc);
        electrostaticsScaleFactor = 1.0F;
    }
};
struct BondedGpuKernelBuffers
{
    //! Force parameters (on GPU)
    DeviceBuffer<t_iparams> d_forceParams = nullptr;
    //! Total Energy (on GPU)
    DeviceBuffer<float> d_vTot = nullptr;
    //! Interaction list atoms (on GPU)
    DeviceBuffer<t_iatom> d_iatoms[numFTypesOnGpu];
};

/*! \internal \brief Implements GPU bondeds */
class ListedForcesGpu::Impl
{
public:
    //! Constructor
    Impl(const gmx_ffparams_t&    ffparams,
         float                    electrostaticsScaleFactor,
         int                      numEnergyGroupsForListedForces,
         const DeviceInformation& deviceInfo,
         const DeviceContext&     deviceContext,
         const DeviceStream&      deviceStream,
         gmx_wallcycle*           wcycle);
    //! \brief Destructor, non-default needed for freeing device-side buffers
    ~Impl();

    /*! \brief Update flag that tells whether there are bonded interactions suitable for the GPU.
     *
     * Intended to be called early during search steps so domainWork flags can be populated.
     */
    void updateHaveInteractions(const InteractionDefinitions& idef);

    /*! \brief Update lists of interactions from idef suitable for the GPU,
     * using the data structures prepared for PP work.
     *
     * Intended to be called after each neighbour search
     * stage. Copies the bonded interactions assigned to the GPU
     * to device data structures, and updates device buffers that
     * may have been updated after search. */
    void updateInteractionListsAndDeviceBuffers(ArrayRef<const int>           nbnxnAtomOrder,
                                                const InteractionDefinitions& idef,
                                                DeviceBuffer<Float4>          d_xqPtr,
                                                DeviceBuffer<RVec>            d_fPtr,
                                                DeviceBuffer<RVec>            d_fShiftPtr);
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
    /*! \brief The interaction lists
     *
     * \todo This is potentially several pinned allocations, which
     * could contribute to exhausting such pages. */
    std::array<HostInteractionList, F_NRE> iLists_;

    //! Tells whether there are any interaction in iLists.
    bool haveInteractions_ = false;
    //! Interaction lists on the device.
    std::array<DeviceBuffer<t_iatom>, F_NRE> d_iAtoms_      = {};
    std::array<int, F_NRE>                   d_iAtomsAlloc_ = {};
    //! Bonded parameters for device-side use.
    DeviceBuffer<t_iparams> d_forceParams_ = nullptr;
    //! Position-charge vector on the device.
    DeviceBuffer<Float4> d_xq_ = nullptr;
    //! Force vector on the device.
    DeviceBuffer<Float3> d_f_ = nullptr;
    //! Shift force vector on the device.
    DeviceBuffer<Float3> d_fShift_ = nullptr;
    //! \brief Host-side virial buffer
    HostVector<float> vTot_ = { {}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported) };
    //! \brief Device-side total virial
    DeviceBuffer<float> d_vTot_ = nullptr;

    //! GPU context object
    const DeviceContext& deviceContext_;
    //! \brief Bonded GPU stream, not owned by this module
    const DeviceStream& deviceStream_;

    //! Parameters, passed to the GPU kernel
    BondedGpuKernelParameters kernelParams_;
    //! Buffers, used in the GPU kernel
    BondedGpuKernelBuffers kernelBuffers_;

    //! Device sub-group/warp size
    int deviceSubGroupSize_;

    //! GPU kernel launch configuration
    KernelLaunchConfig kernelLaunchConfig_;

    //! \brief Pointer to wallcycle structure.
    gmx_wallcycle* wcycle_;
};

} // namespace gmx

#endif // GMX_LISTED_FORCES_LISTED_FORCES_GPU_IMPL_H
