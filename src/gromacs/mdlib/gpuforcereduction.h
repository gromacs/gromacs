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
 * \brief Declares the GPU Force Reduction
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_GPUFORCEREDUCTION_H
#define GMX_MDLIB_GPUFORCEREDUCTION_H

#include "config.h"

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fixedcapacityvector.h"

class GpuEventSynchronizer;
class DeviceStream;
class DeviceContext;

namespace gmx
{

#define HAVE_GPU_FORCE_REDUCTION (GMX_GPU_CUDA || GMX_GPU_SYCL)

/*! \internal
 * \brief Manages the force reduction directly in GPU memory
 *
 * Manages the reduction of multiple GPU force buffers into a single
 * GPU force buffer. The reduction involves at least one (input/output)
 * Rvec-format buffer and one (input) Nbat-format buffer, where the
 * Nbat->Rvec conversion is handled internally. One additional (input)
 * Rvec-format buffer is supported as optional.
 */
class GpuForceReduction
{

public:
    /*! \brief Creates GPU force reduction object
     *
     * \param [in] deviceContext GPU device context
     * \param [in] deviceStream  Stream to use for reduction
     * \param [in] wcycle        Wall-clock cycle counter
     */
    GpuForceReduction(const DeviceContext& deviceContext,
                      const DeviceStream&  deviceStream,
                      gmx_wallcycle*       wcycle);
    ~GpuForceReduction();

    /*! \brief Register a nbnxm-format force to be reduced
     *
     * \param [in] forcePtr  Pointer to force to be reduced
     */
    void registerNbnxmForce(DeviceBuffer<RVec> forcePtr);

    /*! \brief Register a rvec-format force to be reduced
     *
     * \param [in] forcePtr  Pointer to force to be reduced
     */
    void registerRvecForce(DeviceBuffer<RVec> forcePtr);

    /*! \brief Add a dependency for this force reduction
     *
     * \param [in] dependency   Dependency for this reduction
     */
    void addDependency(GpuEventSynchronizer* dependency);

    /*! \brief Reinitialize the GPU force reduction
     *
     * \param [in] baseForcePtr     Pointer to force to be used as a base
     * \param [in] numAtoms         The number of atoms
     * \param [in] cell             Pointer to the cell array
     * \param [in] atomStart        The start atom for the reduction
     * \param [in] accumulate       Whether reduction should be accumulated
     * \param [in] completionMarker Event to be marked when launch of reduction is complete
     */
    void reinit(DeviceBuffer<RVec>    baseForcePtr,
                int                   numAtoms,
                ArrayRef<const int>   cell,
                int                   atomStart,
                bool                  accumulate,
                GpuEventSynchronizer* completionMarker = nullptr);

    /*! \brief Execute the force reduction */
    void execute();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
