/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
#include "gromacs/math/vectypes.h"

#include "gpuforcereduction.h"

namespace gmx
{

//! structure to hold cell information for any nbat-format forces
struct cellInfo
{
    //! cell index mapping for any nbat-format forces
    const int* cell = nullptr;
    //! device copy of cell index mapping for any nbat-format forces
    int* d_cell = nullptr;
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
    Impl(const DeviceContext& deviceContext, const DeviceStream& deviceStreami, gmx_wallcycle* wcycle);
    ~Impl();

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
    void addDependency(GpuEventSynchronizer* const dependency);

    /*! \brief Reinitialize the GPU force reduction
     *
     * \param [in] baseForcePtr     Pointer to force to be used as a base
     * \param [in] numAtoms         The number of atoms
     * \param [in] cell             Pointer to the cell array
     * \param [in] atomStart        The start atom for the reduction
     * \param [in] accumulate       Whether reduction should be accumulated
     * \param [in] completionMarker Event to be marked when launch of reduction is complete
     */
    void reinit(float3*               baseForcePtr,
                const int             numAtoms,
                ArrayRef<const int>   cell,
                const int             atomStart,
                const bool            accumulate,
                GpuEventSynchronizer* completionMarker = nullptr);

    /*! \brief Execute the force reduction */
    void execute();

private:
    //! force to be used as a base for this reduction
    float3* baseForce_ = nullptr;
    //! starting atom
    int atomStart_ = 0;
    //! number of atoms
    int numAtoms_ = 0;
    //! whether reduction is accumulated into base force buffer
    int accumulate_ = true;
    //! cell information for any nbat-format forces
    struct cellInfo cellInfo_;
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! list of dependencies
    gmx::FixedCapacityVector<GpuEventSynchronizer*, 3> dependencyList_;
    //! stream to be used for this reduction
    const DeviceStream& deviceStream_;
    //! Nbnxm force to be added in this reduction
    DeviceBuffer<RVec> nbnxmForceToAdd_ = nullptr;
    //! Rvec-format force to be added in this reduction
    DeviceBuffer<RVec> rvecForceToAdd_ = nullptr;
    //! event to be marked when redcution launch has been completed
    GpuEventSynchronizer* completionMarker_ = nullptr;
    //! The wallclock counter
    gmx_wallcycle* wcycle_ = nullptr;
};

} // namespace gmx

#endif
