/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 *
 * \brief Declares the VsiteGpu class for GPU virtual site construction and force spreading.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_VSITE_GPU_H
#define GMX_MDLIB_VSITE_GPU_H

#include <array>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/pbcutil/pbc_aiuc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"

class InteractionDefinitions;
class DeviceContext;
class DeviceStream;
struct gmx_mtop_t;

namespace gmx
{

//! Number of currently supported virtual site types (1, 2, 3)
constexpr int c_numVSiteTypesOnGpu = 3;

/*! The supported virtual site types,
 * NOTE: these must be in the order that allows hierarchical construction
 *       (i.e. higher-order vsites can be constructed from lower-order vsites, but not vice versa)
 */
constexpr std::array<InteractionFunction, c_numVSiteTypesOnGpu> c_vsiteTypesOnGpu = {
    InteractionFunction::VirtualSite1,
    InteractionFunction::VirtualSite2,
    InteractionFunction::VirtualSite3
};

/*! \brief Returns whether the position of a vsite of \p vsiteType is a linear combination,
 * with weights summing to one, of the positions of its constructing atoms
 *
 * Only such vsites keep their correct positions when all coordinates are scaled by
 * the same matrix, as done for pressure coupling. Other types, e.g. those with a fixed
 * distance or constructed using a cross product, would need to be constructed again.
 */
constexpr bool isLinearVsiteType(InteractionFunction vsiteType)
{
    return vsiteType == InteractionFunction::VirtualSite1 || vsiteType == InteractionFunction::VirtualSite2
           || vsiteType == InteractionFunction::VirtualSite3
           || vsiteType == InteractionFunction::VirtualSiteN;
}

//! Returns whether all vsite types supported on the GPU are linear, see \c isLinearVsiteType()
constexpr bool allVsiteTypesOnGpuAreLinear()
{
    for (const InteractionFunction vsiteType : c_vsiteTypesOnGpu)
    {
        if (!isLinearVsiteType(vsiteType))
        {
            return false;
        }
    }
    return true;
}

/* With GPU update, the vsites are constructed at the end of the update, before
 * the GPU scales the coordinates for pressure coupling (see LegacySimulator::do_md()).
 * This relies on all supported vsite types being linear.
 */
static_assert(allVsiteTypesOnGpuAreLinear(),
              "GPU vsite types that are not linear need to be constructed again after "
              "coordinate scaling for pressure coupling");

/*! \brief Returns the index of \p vsiteType within \c vsiteTypesOnGpu
 *
 * This index is used to address the per-type entries of \c VsiteGpuKernelBuffers.
 */
constexpr int vsiteTypeIndex(InteractionFunction vsiteType)
{
    for (int i = 0; i < c_numVSiteTypesOnGpu; i++)
    {
        if (c_vsiteTypesOnGpu[i] == vsiteType)
        {
            return i;
        }
    }
    return -1;
}

//! \internal \brief Version of InteractionList that supports pinning
struct HostVSiteList
{
    /*! \brief Returns the total number of elements in iatoms */
    int size() const { return iatoms.size(); }

    //! List of interactions
    HostVector<int> iatoms;
};

/*! \internal
 * \brief Vsite parameters and thread ranges
 *
 * This accumulates the parameters that are not device pointers, so they can be passed
 * to the GPU as a single structure. The per-type entries are indexed by the position of
 * the vsite type in \c c_vsiteTypesOnGpu, see \c vsiteTypeIndex().
 *
 * The thread ranges are only used by the fused kernels, which handle all vsite types in
 * a single launch. They are contiguous and padded to whole sub-groups (warps), so that all
 * threads of a sub-group work on the same vsite type.
 */
struct VsiteGpuKernelParameters
{
    //! Periodic boundary data, set by VsiteGpu::construct() before every construction
    PbcAiuc pbcAiuc = {};
    //! The number of vsites of every supported type
    int numVsites[c_numVSiteTypesOnGpu] = {};
    //! The first thread index working on each vsite type
    int typeRangeStart[c_numVSiteTypesOnGpu] = {};
    //! The last thread index (inclusive) of each type, typeRangeStart - 1 when there are no vsites
    int typeRangeEnd[c_numVSiteTypesOnGpu] = {};
};

/*! \internal
 * \brief Vsite device pointers
 *
 * This accumulates the device buffers so they can be passed to the GPU as a single
 * structure. The entries are indexed as those of \c VsiteGpuKernelParameters.
 *
 * The buffer contents have strides as follows:
 * - VSite1: iatoms stride = 2: [vsite, i]
 * - VSite2: iatoms stride = 3: [vsite, i, j], params stride = 1: [a]
 * - VSite3: iatoms stride = 4: [vsite, i, j, k], params stride = 2: [a, b]
 */
struct VsiteGpuKernelBuffers
{
    //! Vsite atom indices (on GPU)
    DeviceBuffer<t_iatom> d_iAtoms[c_numVSiteTypesOnGpu] = {};
    //! Vsite construction parameters (on GPU), empty for VSite1
    DeviceBuffer<float> d_params[c_numVSiteTypesOnGpu] = {};
};

/*! \internal
 * \brief Manages GPU buffers and kernel launches for virtual site operations.
 *
 * ## Lifecycle
 *
 * 1. Construct once, passing the device context and the stream that owns the
 *    update step (DeviceStreamType::UpdateAndConstraints).
 * 2. Call set() after every neighbour-search step to upload the topology.
 * 3. Call construct() after constraints to set vsite positions. Before the first
 *    update step, the vsite positions should be initialized on the host.
 * 4. Call spreadForces() after force calculation to accumulate forces back to
 *    the parent atoms.
 *
 * ## Host/device buffer ownership
 *
 * This class owns the device copies of the per-type atom-index lists and parameter
 * arrays, held in kernelBuffers_.  They are allocated on demand in set() and freed
 * in the destructor.
 *
 * The coordinate and force buffers (d_x, d_f) are owned by UpdateConstrainGpu
 * and passed in by value to construct() and spreadForces().
 */
class VsiteGpu
{
public:
    /*! \brief Construct VsiteGpu.
     *
     * \param[in] mtop           Global topology, used to detect vsites constructed from vsites.
     * \param[in] deviceContext  GPU device context (used for buffer allocation).
     * \param[in] deviceStream   GPU stream to launch all vsite kernels on.
     *                           Should be the same stream as the update/constrain step
     *                           so that construction follows constraints without extra
     *                           synchronization.
     */
    VsiteGpu(const gmx_mtop_t& mtop, const DeviceContext& deviceContext, const DeviceStream& deviceStream);

    ~VsiteGpu();

    /*! \brief Upload vsite topology data to the GPU.
     *
     * Extracts the flat iatoms arrays and compact parameter arrays for the
     * supported vsite types from \p idef, then copies them to device buffers.
     * Must be called after every neighbour-search step (atoms may have been
     * renumbered).
     *
     * \param[in] idef  System interaction definitions.
     */
    void set(const InteractionDefinitions& idef);

    /*! \brief Run the vsite construction kernels.
     *
     * Sets the positions of all type 1, 2 and 3 virtual sites in \p d_x based on
     * the current positions of their constructing atoms. Must be called after
     * constraints have been applied.
     *
     * Uses the fused kernel, unless vsites are constructed from other vsites, in
     * which case the types are constructed one after the other, in the order of
     * \c vsiteTypesOnGpu. The latter can be forced by setting the environment
     * variable GMX_DISABLE_FUSED_VSITE_KERNELS.
     *
     * \param[in,out] d_x     Device coordinate buffer (all atoms, including vsites).
     * \param[in]     pbcAiuc PBC info needed for constructing vsites near the box boundaries.
     */
    void construct(DeviceBuffer<Float3> d_x, const PbcAiuc& pbcAiuc);

    /*! \brief Run the vsite force spreading kernels.
     *
     * Distributes the forces accumulated on virtual sites back to their
     * constructing atoms via atomicAdd.
     *
     * Uses the fused kernel, unless vsites are constructed from other vsites, in
     * which case the forces of the types are spread one after the other, in the
     * reverse order of \c vsiteTypesOnGpu. The latter can be forced by setting the
     * environment variable GMX_DISABLE_FUSED_VSITE_KERNELS.
     *
     * \param[in,out] d_f  Device force buffer (all atoms).
     */
    void spreadForces(DeviceBuffer<Float3> d_f) const;

private:
    // NOTE: The temporary stub implementation for HIP and SYCL builds (vsite_gpu_stubs.cpp)
    // leaves the private members unused. Silence the resulting clang warning for the whole class.
    // This (and the stub) can go once HIP/SYCL gain a real device path.
    CLANG_DIAGNOSTIC_IGNORE("-Wunused-private-field")

    //! GPU device context (needed for buffer allocation).
    const DeviceContext& deviceContext_;
    //! GPU stream for all kernel launches.
    const DeviceStream& deviceStream_;
    //! Device sub-group/warp size, the granularity of the fused kernel thread ranges.
    int deviceSubGroupSize_ = 0;

    /*! \brief Sets the per-type thread ranges of the fused kernels.
     *
     * Assigns a contiguous, sub-group padded range of threads to every supported vsite type,
     * based on the vsite counts in kernelParams_.
     */
    void setFusedKernelThreadRanges();
    /*! \brief Whether the fused kernels, which handle all vsite types in one launch, are used.
     *
     * These can not be used when vsites are constructed from other vsites, as those have to
     * be constructed (and their forces spread) in order of their type. Setting the environment
     * variable GMX_DISABLE_FUSED_VSITE_KERNELS forces the per-type kernels.
     */
    bool useFusedKernels_ = false;

    /*! \brief Pinned host buffers for vsite atom indices
     *
     * One buffer per supported vsite type, with the same strides and contents as the
     * corresponding device buffers in \c kernelBuffers_.
     */
    std::array<HostVSiteList, c_numVSiteTypesOnGpu> h_vsiteLists_;
    //! Pinned host buffers for vsite parameters (empty for vsite 1)
    std::array<HostVector<float>, c_numVSiteTypesOnGpu> h_params_ = {};

    //! Parameters, passed to the GPU kernels
    VsiteGpuKernelParameters kernelParams_;
    //! Device buffers, used in the GPU kernels
    VsiteGpuKernelBuffers kernelBuffers_;

    //! Allocated sizes for the device iatom buffers (elements, not vsites)
    std::array<int, c_numVSiteTypesOnGpu> numIatomsAlloc_ = {};
    //! Allocated sizes for the device parameter buffers (elements, not vsites)
    std::array<int, c_numVSiteTypesOnGpu> numParamsAlloc_ = {};
};

CLANG_DIAGNOSTIC_RESET

} // namespace gmx

#endif // GMX_MDLIB_VSITE_GPU_H
