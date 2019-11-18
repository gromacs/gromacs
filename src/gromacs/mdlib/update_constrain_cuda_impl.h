/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declares CUDA implementation class for update and constraints.
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_UPDATE_CONSTRAIN_CUDA_IMPL_H
#define GMX_MDLIB_UPDATE_CONSTRAIN_CUDA_IMPL_H

#include "gmxpre.h"

#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/mdlib/leapfrog_cuda.cuh"
#include "gromacs/mdlib/lincs_cuda.cuh"
#include "gromacs/mdlib/settle_cuda.cuh"
#include "gromacs/mdlib/update_constrain_cuda.h"
#include "gromacs/mdtypes/inputrec.h"

namespace gmx
{

/*! \internal \brief Class with interfaces and data for CUDA version of Update-Constraint. */
class UpdateConstrainCuda::Impl
{

public:
    /*! \brief Create Update-Constrain object.
     *
     * The constructor is given a non-nullptr \p commandStream, in which all the update and constrain
     * routines are executed. \p xUpdatedOnDevice should mark the completion of all kernels that modify
     * coordinates. The event is maintained outside this class and also passed to all (if any) consumers
     * of the updated coordinates. The \p xUpdatedOnDevice also can not be a nullptr because the
     * markEvent(...) method is called unconditionally.
     *
     * \param[in] ir                Input record data: LINCS takes number of iterations and order of
     *                              projection from it.
     * \param[in] mtop              Topology of the system: SETTLE gets the masses for O and H atoms
     *                              and target O-H and H-H distances from this object.
     * \param[in] commandStream     GPU stream to use. Can be nullptr.
     * \param[in] xUpdatedOnDevice  The event synchronizer to use to mark that update is done on the GPU.
     */
    Impl(const t_inputrec& ir, const gmx_mtop_t& mtop, const void* commandStream, GpuEventSynchronizer* xUpdatedOnDevice);

    ~Impl();

    /*! \brief Integrate
     *
     * Integrates the equation of motion using Leap-Frog algorithm and applies
     * LINCS and SETTLE constraints.
     * If computeVirial is true, constraints virial is written at the provided pointer.
     * doTempCouple should be true if:
     *   1. The temperature coupling is enabled.
     *   2. This is the temperature coupling step.
     * Parameters virial/lambdas can be nullptr if computeVirial/doTempCouple are false.
     *
     * \param[in]  fReadyOnDevice         Event synchronizer indicating that the forces are ready in
     * the device memory. \param[in]  dt                     Timestep. \param[in]  updateVelocities
     * If the velocities should be constrained. \param[in]  computeVirial          If virial should
     * be updated. \param[out] virial                 Place to save virial tensor. \param[in]
     * doTempCouple           If the temperature coupling should be performed. \param[in]  tcstat
     * Temperature coupling data. \param[in]  doPressureCouple       If the temperature coupling
     * should be applied. \param[in]  dtPressureCouple       Period between pressure coupling steps
     * \param[in]  velocityScalingMatrix  Parrinello-Rahman velocity scaling matrix
     */
    void integrate(GpuEventSynchronizer*             fReadyOnDevice,
                   real                              dt,
                   bool                              updateVelocities,
                   bool                              computeVirial,
                   tensor                            virial,
                   bool                              doTempCouple,
                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                   bool                              doPressureCouple,
                   float                             dtPressureCouple,
                   const matrix                      velocityScalingMatrix);

    /*! \brief Set the pointers and update data-structures (e.g. after NB search step).
     *
     * \param[in,out]  d_x            Device buffer with coordinates.
     * \param[in,out]  d_v            Device buffer with velocities.
     * \param[in]      d_f            Device buffer with forces.
     * \param[in] idef                System topology
     * \param[in] md                  Atoms data.
     * \param[in] numTempScaleValues  Number of temperature scaling groups. Set zero for no temperature coupling.
     */
    void set(DeviceBuffer<float>       d_x,
             DeviceBuffer<float>       d_v,
             const DeviceBuffer<float> d_f,
             const t_idef&             idef,
             const t_mdatoms&          md,
             const int                 numTempScaleValues);

    /*! \brief
     * Update PBC data.
     *
     * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \param[in] pbc The PBC data in t_pbc format.
     */
    void setPbc(const t_pbc* pbc);

    /*! \brief Return the synchronizer associated with the event indicated that the coordinates are ready on the device.
     */
    GpuEventSynchronizer* getCoordinatesReadySync();

private:
    //! CUDA stream
    CommandStream commandStream_ = nullptr;

    //! Periodic boundary data
    PbcAiuc pbcAiuc_;

    //! Number of atoms
    int numAtoms_;

    //! Local copy of the pointer to the device positions buffer
    float3* d_x_;
    //! Local copy of the pointer to the device velocities buffer
    float3* d_v_;
    //! Local copy of the pointer to the device forces buffer
    float3* d_f_;

    //! Device buffer for intermediate positions (maintained internally)
    float3* d_xp_;
    //! Number of elements in shifted coordinates buffer
    int numXp_ = -1;
    //! Allocation size for the shifted coordinates buffer
    int numXpAlloc_ = -1;


    //! 1/mass for all atoms (GPU)
    real* d_inverseMasses_;
    //! Number of elements in reciprocal masses buffer
    int numInverseMasses_ = -1;
    //! Allocation size for the reciprocal masses buffer
    int numInverseMassesAlloc_ = -1;

    //! Leap-Frog integrator
    std::unique_ptr<LeapFrogCuda> integrator_;
    //! LINCS CUDA object to use for non-water constraints
    std::unique_ptr<LincsCuda> lincsCuda_;
    //! SETTLE CUDA object for water constrains
    std::unique_ptr<SettleCuda> settleCuda_;

    //! An pointer to the event to indicate when the update of coordinates is complete
    GpuEventSynchronizer* coordinatesReady_;
};

} // namespace gmx

#endif
