/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief May be used to implement force reduction interfaces for non-GPU builds.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "config.h"

#include <cstdint>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

#include "gpuforcereduction.h"

class DeviceContext;
class DeviceStream;
class GpuEventSynchronizer;
struct gmx_wallcycle;

#if !HAVE_GPU_FORCE_REDUCTION

namespace gmx
{

class GpuForceReduction::Impl
{
};

GpuForceReduction::GpuForceReduction(const DeviceContext& /* deviceContext */,
                                     const DeviceStream& /* deviceStream */,
                                     gmx_wallcycle* /*wcycle*/) :
    impl_(nullptr)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::reinit(DeviceBuffer<RVec> /*baseForcePtr*/,
                               const int /*numAtoms*/,
                               ArrayRef<const int> /*cell*/,
                               const int /*atomStart*/,
                               const bool /*accumulate*/,
                               GpuEventSynchronizer* /*completionMarker*/)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::registerNbnxmForce(DeviceBuffer<RVec> /* forcePtr */)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::registerRvecForce(DeviceBuffer<gmx::RVec> /* forcePtr */)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::registerForcesReadyNvshmemFlags(DeviceBuffer<uint64_t> /* forceSyncObjPtr */)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::addDependency(GpuEventSynchronizer* const /* dependency */)
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void GpuForceReduction::execute()
{
    GMX_RELEASE_ASSERT(false, "A CPU stub has been called instead of the correct implementation.");
}

GpuForceReduction::~GpuForceReduction() = default;

} // namespace gmx

#endif /* !HAVE_GPU_FORCE_REDUCTION */
