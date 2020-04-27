/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *
 * \brief Implements GPU bonded lists for non-GPU builds
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include "config.h"

#include <string>

#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Returns whether there are any interactions in ilists suitable for a GPU.
static bool someInteractionsCanRunOnGpu(const InteractionLists& ilists)
{
    for (int fType : fTypesOnGpu)
    {
        if (!ilists[fType].iatoms.empty())
        {
            // Perturbation is not implemented in the GPU bonded
            // kernels. If all the interactions were actually
            // perturbed, then that will be detected later on each
            // domain, and work will never run on the GPU. This is
            // very unlikely to occur, and has little run-time cost,
            // so we don't complicate the code by catering for it
            // here.
            return true;
        }
    }
    return false;
}

//! Returns whether there are any bonded interactions in the global topology suitable for a GPU.
static bool bondedInteractionsCanRunOnGpu(const gmx_mtop_t& mtop)
{
    // Check the regular molecule types
    for (const auto& moltype : mtop.moltype)
    {
        if (someInteractionsCanRunOnGpu(moltype.ilist))
        {
            return true;
        }
    }
    // Check the inter-molecular interactions.
    if (mtop.intermolecular_ilist)
    {
        if (someInteractionsCanRunOnGpu(*mtop.intermolecular_ilist))
        {
            return true;
        }
    }
    return false;
}

/*! \brief Help build a descriptive message in \c error if there are
 * \c errorReasons why bondeds on a GPU are not supported.
 *
 * \returns Whether the lack of errorReasons indicate there is support. */
static bool addMessageIfNotSupported(ArrayRef<const std::string> errorReasons, std::string* error)
{
    bool isSupported = errorReasons.empty();
    if (!isSupported && error)
    {
        *error = "Bonded interactions cannot run on GPUs: ";
        *error += joinStrings(errorReasons, "; ") + ".";
    }
    return isSupported;
}

bool buildSupportsGpuBondeds(std::string* error)
{
    std::vector<std::string> errorReasons;

    if (GMX_DOUBLE)
    {
        errorReasons.emplace_back("not supported with double precision");
    }
    if (GMX_GPU == GMX_GPU_OPENCL)
    {
        errorReasons.emplace_back("not supported with OpenCL build of GROMACS");
    }
    else if (GMX_GPU == GMX_GPU_NONE)
    {
        errorReasons.emplace_back("not supported with CPU-only build of GROMACS");
    }
    return addMessageIfNotSupported(errorReasons, error);
}

bool inputSupportsGpuBondeds(const t_inputrec& ir, const gmx_mtop_t& mtop, std::string* error)
{
    std::vector<std::string> errorReasons;

    if (!bondedInteractionsCanRunOnGpu(mtop))
    {
        errorReasons.emplace_back("No supported bonded interactions are present");
    }
    if (!EI_DYNAMICS(ir.eI))
    {
        errorReasons.emplace_back("not a dynamical integrator");
    }
    if (EI_MIMIC(ir.eI))
    {
        errorReasons.emplace_back("MiMiC");
    }
    if (ir.opts.ngener > 1)
    {
        errorReasons.emplace_back("Cannot run with multiple energy groups");
    }
    return addMessageIfNotSupported(errorReasons, error);
}

#if GMX_GPU != GMX_GPU_CUDA

class GpuBonded::Impl
{
};

GpuBonded::GpuBonded(const gmx_ffparams_t& /* ffparams */,
                     const float /* electrostaticsScaleFactor */,
                     const DeviceContext& /* deviceContext */,
                     const DeviceStream& /* deviceStream */,
                     gmx_wallcycle* /* wcycle */) :
    impl_(nullptr)
{
}

GpuBonded::~GpuBonded() = default;

void GpuBonded::updateInteractionListsAndDeviceBuffers(ArrayRef<const int> /* nbnxnAtomOrder */,
                                                       const InteractionDefinitions& /* idef */,
                                                       void* /* xqDevice */,
                                                       DeviceBuffer<RVec> /* forceDevice */,
                                                       DeviceBuffer<RVec> /* fshiftDevice */)
{
}

void GpuBonded::setPbc(PbcType /* pbcType */, const matrix /* box */, bool /* canMoleculeSpanPbc */)
{
}

bool GpuBonded::haveInteractions() const
{
    return !impl_;
}

void GpuBonded::launchKernel(const gmx::StepWorkload& /* stepWork */) {}

void GpuBonded::setPbcAndlaunchKernel(PbcType /* pbcType */,
                                      const matrix /* box */,
                                      bool /* canMoleculeSpanPbc */,
                                      const gmx::StepWorkload& /* stepWork */)
{
}

void GpuBonded::launchEnergyTransfer() {}

void GpuBonded::waitAccumulateEnergyTerms(gmx_enerdata_t* /* enerd */) {}

void GpuBonded::clearEnergies() {}

#endif /* GMX_GPU != GMX_GPU_CUDA */

} // namespace gmx
