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
#ifndef GMX_MDLIB_GPUFORCEREDUCTION_IMPL_H
#define GMX_MDLIB_GPUFORCEREDUCTION_IMPL_H

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/math/vectypes.h"

#include "gpuforcereduction.h"

namespace gmx
{

//! \internal
//! \brief structure to hold cell information for any nbat-format forces
struct cellInfo
{
    //! cell index mapping for any nbat-format forces
    const int* cell = nullptr;
    //! device copy of cell index mapping for any nbat-format forces
    DeviceBuffer<int> d_cell;
    //! number of atoms in cell array
    int cellSize = -1;
    //! number of atoms allocated in cell array
    int cellSizeAlloc = -1;
};

class GpuForceReduction::Impl
{

public:
    /*! \brief Creates GPU force reduction object
     *
     * \param [in] deviceStream  Stream to use for reduction
     * \param [in] deviceContext GPU device context
     * \param [in] wcycle        The wallclock counter
     */
    Impl(const DeviceContext& deviceContext, const DeviceStream& deviceStream, gmx_wallcycle* wcycle);

    ~Impl();
    /*! \brief Register a nbnxm-format force to be reduced
     *
     * \param [in] forcePtr  Pointer to force to be reduced
     */
    void registerNbnxmForce(DeviceBuffer<Float3> forcePtr);

    /*! \brief Register a rvec-format force to be reduced
     *
     * \param [in] forcePtr  Pointer to force to be reduced
     */
    void registerRvecForce(DeviceBuffer<Float3> forcePtr);

    /*! \brief Register a force synchronization NVSHMEM object
     *
     * \param [in] syncObj  Pointer to force sync object
     */
    void registerForcesReadyNvshmemFlags(DeviceBuffer<uint64_t> syncObj);


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
    void reinit(DeviceBuffer<Float3>  baseForcePtr,
                int                   numAtoms,
                ArrayRef<const int>   cell,
                int                   atomStart,
                bool                  accumulate,
                GpuEventSynchronizer* completionMarker = nullptr);

    /*! \brief Execute the force reduction */
    void execute();

private:
    //! force to be used as a base for this reduction
    DeviceBuffer<Float3> baseForce_;
    //! starting atom
    int atomStart_ = 0;
    //! number of atoms
    int numAtoms_ = 0;
    //! whether reduction is accumulated into base force buffer
    bool accumulate_ = true;
    //! cell information for any nbat-format forces
    struct cellInfo cellInfo_;
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! list of dependencies
    gmx::FixedCapacityVector<GpuEventSynchronizer*, 3> dependencyList_;
    //! stream to be used for this reduction
    const DeviceStream& deviceStream_;
    //! Nbnxm force to be added in this reduction
    DeviceBuffer<RVec> nbnxmForceToAdd_;
    //! Rvec-format force to be added in this reduction
    DeviceBuffer<RVec> rvecForceToAdd_;
    //! nvshmem sync object used in forces reduction kernel
    DeviceBuffer<uint64_t> forcesReadyNvshmemFlags;
    //! nvshmem sync object tracker used in forces reduction kernel
    uint64_t forcesReadyNvshmemFlagsCounter;

    //! event to be marked when reduction launch has been completed
    GpuEventSynchronizer* completionMarker_ = nullptr;
    //! The wallclock counter
    gmx_wallcycle* wcycle_ = nullptr;
};

} // namespace gmx

#endif
