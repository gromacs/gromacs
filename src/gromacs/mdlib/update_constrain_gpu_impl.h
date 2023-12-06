/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Declares GPU implementation class for update and constraints.
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_UPDATE_CONSTRAIN_GPU_IMPL_H
#define GMX_MDLIB_UPDATE_CONSTRAIN_GPU_IMPL_H

#include "gmxpre.h"

#include "config.h"

#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/math/matrix.h"
#include "gromacs/mdlib/leapfrog_gpu.h"
#include "gromacs/mdlib/lincs_gpu.h"
#include "gromacs/mdlib/settle_gpu.h"
#include "gromacs/mdlib/update_constrain_gpu.h"
#include "gromacs/mdtypes/inputrec.h"

class GpuEventSynchronizer;

namespace gmx
{

/*! \internal \brief Class with interfaces and data for GPU version of Update-Constraint. */
class UpdateConstrainGpu::Impl
{

public:
    /*! \brief Create Update-Constrain object.
     *
     * The constructor is given a non-nullptr \p deviceStream, in which all the update and constrain
     * routines are executed.
     *
     * \param[in] ir                  Input record data: LINCS takes number of iterations and order of
     *                                projection from it.
     * \param[in] mtop                Topology of the system: SETTLE gets the masses for O and H atoms
     *                                and target O-H and H-H distances from this object.
     * \param[in] numTempScaleValues  Number of temperature scaling groups. Set zero for no temperature coupling.
     * \param[in] deviceContext       GPU device context.
     * \param[in] deviceStream        GPU stream to use.
     * \param[in] wcycle              The wallclock counter
     */
    Impl(const t_inputrec&    ir,
         const gmx_mtop_t&    mtop,
         int                  numTempScaleValues,
         const DeviceContext& deviceContext,
         const DeviceStream&  deviceStream,
         gmx_wallcycle*       wcycle);

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
     * \param[in]  fReadyOnDevice           Event synchronizer indicating that the forces are ready in
     *                                      the device memory.
     * \param[in]  dt                       Timestep.
     * \param[in]  updateVelocities         If the velocities should be constrained.
     * \param[in]  computeVirial            If virial should be updated.
     * \param[out] virial                   Place to save virial tensor.
     * \param[in]  doTemperatureScaling     If velocities should be scaled for temperature coupling.
     * \param[in]  tcstat                   Temperature coupling data.
     * \param[in]  doParrinelloRahman       If current step is a Parrinello-Rahman pressure coupling step.
     * \param[in]  dtPressureCouple         Period between pressure coupling steps.
     * \param[in]  prVelocityScalingMatrix  Parrinello-Rahman velocity scaling matrix.
     */
    void integrate(GpuEventSynchronizer*             fReadyOnDevice,
                   real                              dt,
                   bool                              updateVelocities,
                   bool                              computeVirial,
                   tensor                            virial,
                   bool                              doTemperatureScaling,
                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                   bool                              doParrinelloRahman,
                   float                             dtPressureCouple,
                   const gmx::Matrix3x3&             prVelocityScalingMatrix);

    /*! \brief Scale coordinates on the GPU for the pressure coupling.
     *
     * After pressure coupling step, the box size may change. Hence, the coordinates should be
     * scaled so that all the particles fit in the new box.
     *
     * \param[in] scalingMatrix Coordinates scaling matrix.
     */
    void scaleCoordinates(const Matrix3x3& scalingMatrix);

    /*! \brief Scale velocities on the GPU for the pressure coupling.
     *
     * After pressure coupling step, the box size may change. In the C-Rescale algorithm, velocities should be scaled.
     *
     * \param[in] scalingMatrix Velocities scaling matrix.
     */
    void scaleVelocities(const Matrix3x3& scalingMatrix);

    /*! \brief Set the pointers and update data-structures (e.g. after NB search step).
     *
     * \param[in,out]  d_x            Device buffer with coordinates.
     * \param[in,out]  d_v            Device buffer with velocities.
     * \param[in]      d_f            Device buffer with forces.
     * \param[in] idef                System topology
     * \param[in] md                  Atoms data.
     */
    void set(DeviceBuffer<Float3>          d_x,
             DeviceBuffer<Float3>          d_v,
             DeviceBuffer<Float3>          d_f,
             const InteractionDefinitions& idef,
             const t_mdatoms&              md);

    /*! \brief
     * Update PBC data.
     *
     * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \param[in] pbcType The type of the periodic boundary.
     * \param[in] box     The periodic boundary box matrix.
     */
    void setPbc(PbcType pbcType, const matrix box);

    /*! \brief Return the synchronizer associated with the event indicated that the coordinates are ready on the device.
     */
    GpuEventSynchronizer* xUpdatedOnDeviceEvent();

    /*! \brief
     * Returns whether the maximum number of coupled constraints is supported
     * by the GPU LINCS code.
     *
     * \param[in] mtop The molecular topology
     */
    static bool isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop);

private:
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! GPU stream
    const DeviceStream& deviceStream_;

    //! Periodic boundary data
    PbcAiuc pbcAiuc_;

    //! Number of atoms
    int numAtoms_;

    //! Local copy of the pointer to the device positions buffer
    DeviceBuffer<Float3> d_x_;
    //! Local copy of the pointer to the device velocities buffer
    DeviceBuffer<Float3> d_v_;
    //! Local copy of the pointer to the device forces buffer
    DeviceBuffer<Float3> d_f_;

    //! Device buffer for storing a copy of previous coordinates (maintained internally)
    DeviceBuffer<Float3> d_x0_;
    //! Number of elements in shifted coordinates buffer
    int numXp_ = -1;
    //! Allocation size for the shifted coordinates buffer
    int numXpAlloc_ = -1;


    //! 1/mass for all atoms (GPU)
    DeviceBuffer<real> d_inverseMasses_;
    //! Number of elements in reciprocal masses buffer
    int numInverseMasses_ = -1;
    //! Allocation size for the reciprocal masses buffer
    int numInverseMassesAlloc_ = -1;

    //! Leap-Frog integrator
    std::unique_ptr<LeapFrogGpu> integrator_;
    //! LINCS GPU object to use for non-water constraints
    std::unique_ptr<LincsGpu> lincsGpu_;
    //! SETTLE GPU object for water constrains
    std::unique_ptr<SettleGpu> settleGpu_;

    //! The event to indicate when the update of coordinates is complete
    GpuEventSynchronizer xUpdatedOnDeviceEvent_;
    //! The wallclock counter
    gmx_wallcycle* wcycle_ = nullptr;
};

} // namespace gmx

#endif // GMX_MDLIB_UPDATE_CONSTRAIN_GPU_IMPL_H
