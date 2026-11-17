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
 * \brief Implements VsiteGpu: host-side buffer management and kernel dispatch
 *        for virtual site construction and force spreading.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "vsite_gpu.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/mdlib/vsite_gpu_internal.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/functions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

//! \brief Returns whether any vsite in \p mtop is a constructing atom of another vsite.
bool haveDependentVsites(const gmx_mtop_t& mtop)
{
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        for (const InteractionFunction type : c_vsiteTypesOnGpu)
        {
            const InteractionList& il = molt.ilist[type];
            // stride: number of atoms per interaction + 1 parameter index at the start
            const int stride = interaction_function[type].nratoms + 1;
            for (int i = 0; i < il.size(); i += stride)
            {
                // i is the parameter index, i + 1 the vsite, the rest the constructing atoms
                for (int j = i + 2; j < i + stride; j++)
                {
                    if (molt.atoms.atom[il.iatoms[j]].ptype == ParticleType::VSite)
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

} // namespace

VsiteGpu::VsiteGpu(const gmx_mtop_t& mtop, const DeviceContext& deviceContext, const DeviceStream& deviceStream) :
    deviceContext_(deviceContext), deviceStream_(deviceStream)
{
    // Dependent vsites have to be constructed in order of their type, which the fused kernels
    // can not do. They can also be disabled manually, for debugging or testing purposes.
    const bool fusedKernelsDisabled = (std::getenv("GMX_DISABLE_FUSED_VSITE_KERNELS") != nullptr);
    useFusedKernels_                = !haveDependentVsites(mtop) && !fusedKernelsDisabled;
    deviceSubGroupSize_             = chooseSubGroupSizeForDevice(deviceContext_.deviceInfo(), 32);

    // Set up pinned host buffers now that a DeviceContext is available
    const HostAllocationPolicy pinnedPolicy{ deviceContext_, PinningPolicy::PinnedIfSupported };
    for (int t = 0; t < c_numVSiteTypesOnGpu; t++)
    {
        h_vsiteLists_[t].iatoms = HostVector<int>(pinnedPolicy);
        h_params_[t]            = HostVector<float>(pinnedPolicy);
    }
}

VsiteGpu::~VsiteGpu()
{
    // wait for all tasks to complete before freeing buffers
    try
    {
        deviceStream_.synchronize();

        for (int t = 0; t < c_numVSiteTypesOnGpu; t++)
        {
            if (kernelBuffers_.d_iAtoms[t])
            {
                freeDeviceBuffer(&kernelBuffers_.d_iAtoms[t]);
                kernelBuffers_.d_iAtoms[t] = nullptr;
            }
            if (kernelBuffers_.d_params[t])
            {
                freeDeviceBuffer(&kernelBuffers_.d_params[t]);
                kernelBuffers_.d_params[t] = nullptr;
            }
        }
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of VsiteGpu: %s\n", e.what());
    }
}

void VsiteGpu::set(const InteractionDefinitions& idef)
{
    // copy vsite data to host pinned buffers and then to device buffers, for all supported vsite types
    for (int t = 0; t < c_numVSiteTypesOnGpu; t++)
    {
        const InteractionFunction type     = c_vsiteTypesOnGpu[t];
        const auto&               il       = idef.il[type];
        auto&                     h_list   = h_vsiteLists_[t];
        auto&                     h_params = h_params_[t];
        // stride: number of atoms per interaction + 1 parameter index at the start
        const int stride = interaction_function[type].nratoms + 1;

        if (il.empty())
        {
            kernelParams_.numVsites[t] = 0;
            h_list.iatoms.clear();
            h_params.clear();
            continue;
        }
        // build host-side array for iatoms, stripping the leading parameter index from each entry in il
        kernelParams_.numVsites[t] = il.size() / stride;
        h_list.iatoms.resize(kernelParams_.numVsites[t] * (stride - 1));
        int k = 0;
        for (int i = 0; i < kernelParams_.numVsites[t]; i++)
        {
            for (int j = 1; j < stride; j++) // j=0 is paramIdx, skip it
            {
                h_list.iatoms[k++] = il.iatoms[stride * i + j];
            }
        }
        // reallocate device buffer
        int newSize;
        reallocateDeviceBuffer(
                &kernelBuffers_.d_iAtoms[t], h_list.size(), &newSize, &numIatomsAlloc_[t], deviceContext_);
        // copy data to device
        copyToDeviceBuffer(&kernelBuffers_.d_iAtoms[t],
                           h_list.iatoms.data(),
                           0,
                           newSize,
                           deviceStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);

        // build host-side array for parameters
        if (type == InteractionFunction::VirtualSite1)
        {
            // no parameters for VS1
            h_params.clear();
        }
        else if (type == InteractionFunction::VirtualSite2)
        {
            // 1 parameter per vsite, stride 4 in iatoms: [paramIdx, vsite, i, j]
            h_params.resize(kernelParams_.numVsites[t]);
            for (int i = 0; i < kernelParams_.numVsites[t]; i++)
            {
                int paramIdx = il.iatoms[stride * i];
                h_params[i]  = static_cast<float>(idef.iparams[paramIdx].vsite.a);
            }
        }
        else if (type == InteractionFunction::VirtualSite3)
        {
            // 2 parameters per vsite, stride 5 in iatoms: [paramIdx, vsite, i, j, k]
            h_params.resize(kernelParams_.numVsites[t] * 2);
            for (int i = 0; i < kernelParams_.numVsites[t]; i++)
            {
                int paramIdx        = il.iatoms[stride * i];
                h_params[2 * i]     = static_cast<float>(idef.iparams[paramIdx].vsite.a);
                h_params[2 * i + 1] = static_cast<float>(idef.iparams[paramIdx].vsite.b);
            }
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Unexpected vsite type in VsiteGpu::set");
        }
        // reallocate device buffer for parameters
        reallocateDeviceBuffer(
                &kernelBuffers_.d_params[t], h_params.size(), &newSize, &numParamsAlloc_[t], deviceContext_);
        // copy parameters to device
        copyToDeviceBuffer(&kernelBuffers_.d_params[t],
                           h_params.data(),
                           0,
                           newSize,
                           deviceStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
    }
    // The vsite counts changed, so the thread ranges of the fused kernels need to be updated
    setFusedKernelThreadRanges();
}

void VsiteGpu::setFusedKernelThreadRanges()
{
    int rangeStart = 0;
    for (int t = 0; t < c_numVSiteTypesOnGpu; t++)
    {
        kernelParams_.typeRangeStart[t] = rangeStart;
        // Pad each range to whole sub-groups, so no sub-group works on more than one vsite type
        kernelParams_.typeRangeEnd[t] =
                rangeStart + divideRoundUp(kernelParams_.numVsites[t], deviceSubGroupSize_) * deviceSubGroupSize_
                - 1;
        rangeStart = kernelParams_.typeRangeEnd[t] + 1;
    }
}

void VsiteGpu::construct(DeviceBuffer<Float3> d_x, const PbcAiuc& pbcAiuc)
{
    kernelParams_.pbcAiuc = pbcAiuc;
    if (useFusedKernels_)
    {
        launchConstructVSitesFusedKernel(kernelParams_, kernelBuffers_, d_x, deviceStream_);
    }
    else
    {
        // Construct in order of increasing vsite type, so that higher-order vsites can be
        // constructed from lower-order ones. The launchers are no-ops for types without vsites.
        launchConstructVSiteKernel<InteractionFunction::VirtualSite1>(
                kernelParams_, kernelBuffers_, d_x, deviceStream_);
        launchConstructVSiteKernel<InteractionFunction::VirtualSite2>(
                kernelParams_, kernelBuffers_, d_x, deviceStream_);
        launchConstructVSiteKernel<InteractionFunction::VirtualSite3>(
                kernelParams_, kernelBuffers_, d_x, deviceStream_);
    }
}

void VsiteGpu::spreadForces(DeviceBuffer<Float3> d_f) const
{
    if (useFusedKernels_)
    {
        launchSpreadVSiteForcesFusedKernel(kernelParams_, kernelBuffers_, d_f, deviceStream_);
    }
    else
    {
        // Spread forces in reverse order of construction, so that the forces on the constructing
        // atoms are complete before they are spread further. The launchers are no-ops for types
        // without vsites.
        launchSpreadVSiteForceKernel<InteractionFunction::VirtualSite3>(
                kernelParams_, kernelBuffers_, d_f, deviceStream_);
        launchSpreadVSiteForceKernel<InteractionFunction::VirtualSite2>(
                kernelParams_, kernelBuffers_, d_f, deviceStream_);
        launchSpreadVSiteForceKernel<InteractionFunction::VirtualSite1>(
                kernelParams_, kernelBuffers_, d_f, deviceStream_);
    }
}

} // namespace gmx
