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
 * \brief Declaration of interfaces for GPU state data propagator object.
 *
 * This object stores and manages positions, velocities and forces for
 * all particles in the system on the GPU.
 *
 * \todo Add cycle counters.
 * \todo Add synchronization points.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_H
#define GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

#include "locality.h"

class GpuEventSynchronizer;
struct gmx_wallcycle;

namespace gmx
{

class StatePropagatorDataGpu
{
public:
    /*! \brief Constructor
     *
     * The buffers are reallocated only at the reinit call, the padding is
     * used there for the coordinates buffer. It is needed for PME and added at
     * the end of the buffer. It is assumed that if the rank has PME duties on the
     * GPU, all coordinates are copied to the GPU and hence, for this rank, the
     * coordinates buffer is not split into local and non-local ranges. For other
     * ranks, the padding size is zero. This works because only one rank ever does
     * PME work on the GPU, and if that rank also does PP work that is the only
     * rank. So all coordinates are always transferred.
     *
     * In OpenCL, only pmeStream is used since it is the only stream created in
     * PME context. The local and non-local streams are only needed when buffer
     * ops are offloaded. This feature is currently not available in OpenCL and
     * hence these streams are not set in these builds.
     *
     * \note In CUDA, the update stream is created in the constructor as a temporary
     *       solution, in place until the stream manager is introduced.
     *       Note that this makes it impossible to construct this object in CUDA
     *       builds executing on a host without any CUDA-capable device available.
     *
     * \note In CUDA, \p deviceContext is unused, hence always nullptr;
     *       all stream arguments can also be nullptr in runs where the
     *       respective streams are not required.
     *       In OpenCL, \p deviceContext needs to be a valid device context.
     *       In OpenCL runs StatePropagatorDataGpu is currently only used
     *       with PME offload, and only on ranks with PME duty. Hence, the
     *       \p pmeStream argument needs to be a valid OpenCL queue object
     *       which must have been created in \p deviceContext.
     *
     * \todo Make a \p CommandStream visible in the CPU parts of the code so we
     *       will not have to pass a void*.
     * \todo Make a \p DeviceContext object visible in CPU parts of the code so we
     *       will not have to pass a void*.
     *
     *  \param[in] pmeStream       Device PME stream, nullptr allowed.
     *  \param[in] localStream     Device NBNXM local stream, nullptr allowed.
     *  \param[in] nonLocalStream  Device NBNXM non-local stream, nullptr allowed.
     *  \param[in] deviceContext   Device context, nullptr allowed.
     *  \param[in] transferKind    H2D/D2H transfer call behavior (synchronous or not).
     *  \param[in] paddingSize     Padding size for coordinates buffer.
     *  \param[in] wcycle          Wall cycle counter data.
     */
    StatePropagatorDataGpu(const void*        pmeStream,
                           const void*        localStream,
                           const void*        nonLocalStream,
                           const void*        deviceContext,
                           GpuApiCallBehavior transferKind,
                           int                paddingSize,
                           gmx_wallcycle*     wcycle);

    /*! \brief Constructor to use in PME-only rank and in tests.
     *
     *  This constructor should be used if only a coordinate device buffer should be managed
     *  using a single stream. Any operation on force or velocity buffer as well as copy of
     *  non-local coordinates will exit with assertion failure. Note, that the pmeStream can
     *  not be a nullptr and the constructor will exit with an assertion failure.
     *
     *  \todo Currently, unsupported copy operations are blocked by assertion that the stream
     *        not nullptr. This should be improved.
     *
     *  \param[in] pmeStream       Device PME stream, nullptr is not allowed.
     *  \param[in] deviceContext   Device context, nullptr allowed for non-OpenCL builds.
     *  \param[in] transferKind    H2D/D2H transfer call behavior (synchronous or not).
     *  \param[in] paddingSize     Padding size for coordinates buffer.
     *  \param[in] wcycle          Wall cycle counter data.
     */
    StatePropagatorDataGpu(const void*        pmeStream,
                           const void*        deviceContext,
                           GpuApiCallBehavior transferKind,
                           int                paddingSize,
                           gmx_wallcycle*     wcycle);

    //! Move constructor
    StatePropagatorDataGpu(StatePropagatorDataGpu&& other) noexcept;
    //! Move assignment
    StatePropagatorDataGpu& operator=(StatePropagatorDataGpu&& other) noexcept;
    //! Destructor
    ~StatePropagatorDataGpu();

    /*! \brief Set the ranges for local and non-local atoms and reallocates buffers.
     *
     * Reallocates coordinate, velocities and force buffers on the device.
     *
     * \note
     * The coordinates buffer is (re)allocated, when required by PME, with a padding,
     * the size of which is set by the constructor. The padding region clearing kernel
     * is scheduled in the \p pmeStream_ (unlike the coordinates H2D) as only the PME
     * task uses this padding area.
     *
     * \note
     * The force buffer is cleared if its size increases, so that previously unused
     * memory is cleared before forces are accumulated.
     *
     *  \param[in] numAtomsLocal  Number of atoms in local domain.
     *  \param[in] numAtomsAll    Total number of atoms to handle.
     */
    void reinit(int numAtomsLocal, int numAtomsAll);

    /*! \brief Returns the range of atoms to be copied based on the copy type (all, local or non-local).
     *
     * \todo There are at least three versions of the function with this functionality in the code:
     *       this one and two more in NBNXM. These should be unified in a shape of a general function
     *       in DD.
     *
     * \param[in]  atomLocality    If all, local or non-local ranges are needed.
     *
     * \returns Tuple, containing the index of the first atom in the range and the total number of atoms in the range.
     */
    std::tuple<int, int> getAtomRangesFromAtomLocality(AtomLocality atomLocality);


    /*! \brief Get the positions buffer on the GPU.
     *
     *  \returns GPU positions buffer.
     */
    DeviceBuffer<float> getCoordinates();

    /*! \brief Copy positions to the GPU memory.
     *
     *  \param[in] h_x           Positions in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyCoordinatesToGpu(gmx::ArrayRef<const gmx::RVec> h_x, AtomLocality atomLocality);

    /*! \brief Get the event synchronizer of the coordinates ready for the consumption on the device.
     *
     * Returns the event synchronizer which indicates that the coordinates are ready for the
     * consumption on the device. Takes into account that the producer may be different.
     *
     * If the update is offloaded, and the current step is not a DD/search step, the returned
     * synchronizer indicates the completion of GPU update-constraint kernels. Otherwise, on search
     * steps and if update is not offloaded, the coordinates are provided by the H2D copy and the
     * returned synchronizer indicates that the copy is complete.
     *
     *  \param[in] atomLocality    Locality of the particles to wait for.
     *  \param[in] simulationWork  The simulation lifetime flags.
     *  \param[in] stepWork        The step lifetime flags.
     *
     *  \returns  The event to synchronize the stream that consumes coordinates on device.
     */
    GpuEventSynchronizer* getCoordinatesReadyOnDeviceEvent(AtomLocality              atomLocality,
                                                           const SimulationWorkload& simulationWork,
                                                           const StepWorkload&       stepWork);

    /*! \brief Blocking wait until coordinates are copied to the device.
     *
     * Synchronizes the stream in which the copy was executed.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitCoordinatesCopiedToDevice(AtomLocality atomLocality);

    /*! \brief Getter for the event synchronizer for the update is done on th GPU
     *
     *  \returns  The event to synchronize the stream coordinates wre updated on device.
     */
    GpuEventSynchronizer* xUpdatedOnDevice();

    /*! \brief Copy positions from the GPU memory.
     *
     *  \param[in] h_x           Positions buffer in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec> h_x, AtomLocality atomLocality);

    /*! \brief Wait until coordinates are available on the host.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitCoordinatesReadyOnHost(AtomLocality atomLocality);


    /*! \brief Get the velocities buffer on the GPU.
     *
     *  \returns GPU velocities buffer.
     */
    DeviceBuffer<float> getVelocities();

    /*! \brief Copy velocities to the GPU memory.
     *
     *  \param[in] h_v           Velocities in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyVelocitiesToGpu(gmx::ArrayRef<const gmx::RVec> h_v, AtomLocality atomLocality);

    /*! \brief Get the event synchronizer for the H2D velocities copy.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     *
     *  \returns  The event to synchronize the stream that consumes velocities on device.
     */
    GpuEventSynchronizer* getVelocitiesReadyOnDeviceEvent(AtomLocality atomLocality);

    /*! \brief Copy velocities from the GPU memory.
     *
     *  \param[in] h_v           Velocities buffer in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec> h_v, AtomLocality atomLocality);

    /*! \brief Wait until velocities are available on the host.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitVelocitiesReadyOnHost(AtomLocality atomLocality);


    /*! \brief Get the force buffer on the GPU.
     *
     *  \returns GPU force buffer.
     */
    DeviceBuffer<float> getForces();

    /*! \brief Copy forces to the GPU memory.
     *
     *  \param[in] h_f           Forces in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyForcesToGpu(gmx::ArrayRef<const gmx::RVec> h_f, AtomLocality atomLocality);

    /*! \brief Get the event synchronizer for the forces ready on device.
     *
     *  Returns either of the event synchronizers, depending on the offload scenario
     *  for the current simulation timestep:
     *  1. The forces are copied to the device (when GPU buffer ops are off)
     *  2. The forces are reduced on the device (GPU buffer ops are on)
     *
     *  \todo Pass step workload instead of the useGpuFBufferOps boolean.
     *
     *  \param[in] atomLocality      Locality of the particles to wait for.
     *  \param[in] useGpuFBufferOps  If the force buffer ops are offloaded to the GPU.
     *
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getForcesReadyOnDeviceEvent(AtomLocality atomLocality, bool useGpuFBufferOps);

    /*! \brief Getter for the event synchronizer for the forces are reduced on the GPU.
     *
     *  \returns  The event to mark when forces are reduced on the GPU.
     */
    GpuEventSynchronizer* fReducedOnDevice();

    /*! \brief Copy forces from the GPU memory.
     *
     *  \param[in] h_f           Forces buffer in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyForcesFromGpu(gmx::ArrayRef<gmx::RVec> h_f, AtomLocality atomLocality);

    /*! \brief Wait until forces are available on the host.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitForcesReadyOnHost(AtomLocality atomLocality);

    /*! \brief Getter for the update stream.
     *
     *  \todo This is temporary here, until the management of this stream is taken over.
     *
     *  \returns The device command stream to use in update-constraints.
     */
    void* getUpdateStream();

    /*! \brief Getter for the number of local atoms.
     *
     *  \returns The number of local atoms.
     */
    int numAtomsLocal();

    /*! \brief Getter for the total number of atoms.
     *
     *  \returns The total number of atoms.
     */
    int numAtomsAll();

private:
    class Impl;
    gmx::PrivateImplPointer<Impl> impl_;
    GMX_DISALLOW_COPY_AND_ASSIGN(StatePropagatorDataGpu);
};

} // namespace gmx

#endif // GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_H
