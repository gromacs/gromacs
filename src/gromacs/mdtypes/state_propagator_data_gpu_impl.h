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
 * \brief Declaration of low-level functions and fields of GPU state propagator object.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H
#define GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H

#include "gmxpre.h"

#include "config.h"

#include <memory>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/utility/enumerationhelpers.h"

struct gmx_wallcycle;

namespace gmx
{

class StatePropagatorDataGpu::Impl
{
public:
    Impl();


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
     *  \param[in] deviceStreamManager         Object that owns the DeviceContext and DeviceStreams.
     *  \param[in] transferKind                H2D/D2H transfer call behavior (synchronous or not).
     *  \param[in] allocationBlockSizeDivisor  Determines the padding size for coordinates buffer.
     *  \param[in] wcycle                      Wall cycle counter data.
     */
    Impl(const DeviceStreamManager& deviceStreamManager,
         GpuApiCallBehavior         transferKind,
         int                        allocationBlockSizeDivisor,
         gmx_wallcycle*             wcycle);

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
     *  \param[in] deviceContext   Device context.
     *  \param[in] transferKind    H2D/D2H transfer call behavior (synchronous or not).
     *  \param[in] allocationBlockSizeDivisor  Determines the padding size for coordinates buffer.
     *  \param[in] wcycle          Wall cycle counter data.
     */
    Impl(const DeviceStream*  pmeStream,
         const DeviceContext& deviceContext,
         GpuApiCallBehavior   transferKind,
         int                  allocationBlockSizeDivisor,
         gmx_wallcycle*       wcycle);

    ~Impl();


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
    std::tuple<int, int> getAtomRangesFromAtomLocality(AtomLocality atomLocality) const;


    /*! \brief Get the positions buffer on the GPU.
     *
     *  \returns GPU positions buffer.
     */
    DeviceBuffer<RVec> getCoordinates();

    /*! \brief Copy positions to the GPU memory.
     *
     *  \param[in] h_x           Positions in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     *  \param[in] expectedConsumptionCount How many times will the event be used later
     *                                      via \ref getCoordinatesReadyOnDeviceEvent.
     */
    void copyCoordinatesToGpu(gmx::ArrayRef<const gmx::RVec> h_x,
                              AtomLocality                   atomLocality,
                              int                            expectedConsumptionCount);

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
     *  \param[in] gpuCoordinateHaloLaunched Event recorded when GPU coordinate halo has been launched.
     *
     *  \returns  The event to synchronize the stream that consumes coordinates on device.
     */
    GpuEventSynchronizer* getCoordinatesReadyOnDeviceEvent(AtomLocality              atomLocality,
                                                           const SimulationWorkload& simulationWork,
                                                           const StepWorkload&       stepWork,
                                                           GpuEventSynchronizer* gpuCoordinateHaloLaunched = nullptr);

    /*! \brief Wait until coordinates are available on the device. */
    void waitCoordinatesUpdatedOnDevice() { xUpdatedOnDeviceEvent_->waitForEvent(); }

    /*! \brief Blocking wait until coordinates are copied to the device.
     *
     * Synchronizes the stream in which the copy was executed.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitCoordinatesCopiedToDevice(AtomLocality atomLocality);

    /*! \brief Consume the event for copying coordinates to the device.
     *
     * Used for manual event consumption. Does nothing except changing the internal event counter.
     *
     *  \param[in] atomLocality  Locality of the particles.
     */
    void consumeCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality);

    /*! \brief Reset the event for copying coordinates to the device.
     *
     * Used for manual event consumption. Does nothing except resetting the event.
     *
     *  \param[in] atomLocality  Locality of the particles.
     */
    void resetCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality);

    /*! \brief Setter for the event synchronizer for the update is done on th GPU
     *
     *  \param[in] xUpdatedOnDeviceEvent  The event to synchronize the stream coordinates wre updated on device.
     */
    void setXUpdatedOnDeviceEvent(GpuEventSynchronizer* xUpdatedOnDeviceEvent);

    /*! \brief Set the expected consumption count for the event associated with GPU update.
     *
     *  \param[in] expectedConsumptionCount  New value.
     */
    void setXUpdatedOnDeviceEventExpectedConsumptionCount(int expectedConsumptionCount);

    /*! \brief Set the expected consumption count for the event associated with GPU forces computation.
     *
     *  \param[in] atomLocality  Locality of the particles.
     *  \param[in] expectedConsumptionCount  New value.
     */
    void setFReadyOnDeviceEventExpectedConsumptionCount(AtomLocality atomLocality,
                                                        int          expectedConsumptionCount);

    /*! \brief Copy positions from the GPU memory, with an optional explicit dependency.
     *
     *  \param[in] h_x           Positions buffer in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     *  \param[in] dependency    Dependency event for this operation.
     */
    void copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec> h_x,
                                AtomLocality             atomLocality,
                                GpuEventSynchronizer*    dependency = nullptr);

    /*! \brief Wait until coordinates are available on the host.
     *
     *  \param[in] atomLocality  Locality of the particles to wait for.
     */
    void waitCoordinatesReadyOnHost(AtomLocality atomLocality);


    /*! \brief Get the velocities buffer on the GPU.
     *
     *  \returns GPU velocities buffer.
     */
    DeviceBuffer<RVec> getVelocities();

    /*! \brief Copy velocities to the GPU memory.
     *
     *  \param[in] h_v           Velocities in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyVelocitiesToGpu(gmx::ArrayRef<const gmx::RVec> h_v, AtomLocality atomLocality);

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
    DeviceBuffer<RVec> getForces();

    /*! \brief Copy forces to the GPU memory.
     *
     *  \param[in] h_f           Forces in the host memory.
     *  \param[in] atomLocality  Locality of the particles to copy.
     */
    void copyForcesToGpu(gmx::ArrayRef<const gmx::RVec> h_f, AtomLocality atomLocality);

    /*! \brief Clear forces in the GPU memory.
     *
     *  \param[in] atomLocality  Locality of the particles to clear.
     *  \param[in] dependency    Dependency event for this operation.
     */
    void clearForcesOnGpu(AtomLocality atomLocality, GpuEventSynchronizer* dependency);

    /*! \brief Get the event synchronizer for the forces ready on device.
     *
     *  Returns either of the event synchronizers, depending on the offload scenario
     *  for the current simulation timestep:
     *  1. The forces are copied to the device (when GPU buffer ops are off)
     *  2. The forces are reduced on the device (GPU buffer ops are on)
     *
     *  \param[in] stepWork        Step workload flags
     *  \param[in] simulationWork  Simulation workload flags
     *
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getLocalForcesReadyOnDeviceEvent(StepWorkload       stepWork,
                                                           SimulationWorkload simulationWork);

    /*! \brief Getter for the event synchronizer for when forces are reduced on the GPU.
     *
     *  \param[in] atomLocality      Locality of the particles to wait for.
     *  \returns                     The event to mark when forces are reduced on the GPU.
     */
    GpuEventSynchronizer* fReducedOnDevice(AtomLocality atomLocality);

    //! \brief Consume the event for when the forces are reduced on device.
    void consumeForcesReducedOnDeviceEvent(AtomLocality atomLocality);

    /*! \brief Getter for the event synchronizer for the forces are ready for GPU update.
     *
     *  \param[in] atomLocality      Locality of the particles to wait for.
     *  \returns                     The event to mark when forces are ready for GPU update.
     */
    GpuEventSynchronizer* fReadyOnDevice(AtomLocality atomLocality);

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
    const DeviceStream* getUpdateStream();

    /*! \brief Getter for the number of local atoms.
     *
     *  \returns The number of local atoms.
     */
    int numAtomsLocal() const;

    /*! \brief Getter for the total number of atoms.
     *
     *  \returns The total number of atoms.
     */
    int numAtomsAll() const;

private:
    //! GPU PME stream.
    const DeviceStream* pmeStream_;
    //! GPU NBNXM local stream.
    const DeviceStream* localStream_;
    //! GPU NBNXM non-local stream.
    const DeviceStream* nonLocalStream_;
    //! GPU Update-constraints stream.
    const DeviceStream* updateStream_;

    // Streams to use for coordinates H2D and D2H copies (one event for each atom locality)
    EnumerationArray<AtomLocality, const DeviceStream*> xCopyStreams_ = { { nullptr } };
    // Streams to use for velocities H2D and D2H copies (one event for each atom locality)
    EnumerationArray<AtomLocality, const DeviceStream*> vCopyStreams_ = { { nullptr } };
    // Streams to use for forces H2D and D2H copies (one event for each atom locality)
    EnumerationArray<AtomLocality, const DeviceStream*> fCopyStreams_ = { { nullptr } };
    // Streams internal to this module
    std::unique_ptr<DeviceStream> copyInStream_;
    std::unique_ptr<DeviceStream> memsetStream_;

    /*! \brief An array of events that indicate H2D copy is complete (one event for each atom locality)
     *
     * \todo Reconsider naming. It should be xCopiedToDevice or xH2DCopyComplete, etc.
     */
    EnumerationArray<AtomLocality, GpuEventSynchronizer> xReadyOnDevice_;
    //! A pointer to an event that the coordinates are ready after update-constraints execution
    GpuEventSynchronizer* xUpdatedOnDeviceEvent_ = nullptr;
    //! An array of events that indicate D2H copy of coordinates is complete (one event for each atom locality)
    EnumerationArray<AtomLocality, GpuEventSynchronizer> xReadyOnHost_;

    //! An array of events that indicate D2H copy of velocities is complete (one event for each atom locality)
    EnumerationArray<AtomLocality, GpuEventSynchronizer> vReadyOnHost_;

    //! An array of events that indicate H2D copy of forces is complete (one event for each atom locality)
    EnumerationArray<AtomLocality, GpuEventSynchronizer> fReadyOnDevice_;
    //! An array of events that indicate the forces were reduced on the GPU (one event for each atom locality)
    EnumerationArray<AtomLocality, GpuEventSynchronizer> fReducedOnDevice_;
    //! An array of events that indicate D2H copy of forces is complete (one event for each atom locality)
    EnumerationArray<AtomLocality, GpuEventSynchronizer> fReadyOnHost_;

    //! Device context
    const DeviceContext& deviceContext_;
    //! Default device calls behavior
    GpuApiCallBehavior transferKind_ = GpuApiCallBehavior::Async;
    //! Required minimum divisor of the allocation size of the coordinates buffer
    int allocationBlockSizeDivisor_ = 0;

    //! Number of local atoms
    int numAtomsLocal_ = -1;
    //! Total number of atoms
    int numAtomsAll_ = -1;

    //! Device positions buffer
    DeviceBuffer<RVec> d_x_;
    //! Number of particles saved in the positions buffer
    int d_xSize_ = -1;
    //! Allocation size for the positions buffer
    int d_xCapacity_ = -1;

    //! Device velocities buffer
    DeviceBuffer<RVec> d_v_;
    //! Number of particles saved in the velocities buffer
    int d_vSize_ = -1;
    //! Allocation size for the velocities buffer
    int d_vCapacity_ = -1;

    //! Device force buffer
    DeviceBuffer<RVec> d_f_;
    //! Number of particles saved in the force buffer
    int d_fSize_ = -1;
    //! Allocation size for the force buffer
    int d_fCapacity_ = -1;

    //! \brief Pointer to wallcycle structure.
    gmx_wallcycle* wcycle_;

    /*! \brief Performs the copy of data from host to device buffer.
     *
     * \todo Template on locality.
     *
     *  \param[out] d_data         Device-side buffer.
     *  \param[in]  h_data         Host-side buffer.
     *  \param[in]  dataSize       Device-side data allocation size.
     *  \param[in]  atomLocality   If all, local or non-local ranges should be copied.
     *  \param[in]  deviceStream   GPU stream to execute copy in.
     */
    void copyToDevice(DeviceBuffer<RVec>             d_data,
                      gmx::ArrayRef<const gmx::RVec> h_data,
                      int                            dataSize,
                      AtomLocality                   atomLocality,
                      const DeviceStream&            deviceStream);

    /*! \brief Performs the copy of data from device to host buffer.
     *
     *  \param[out] h_data         Host-side buffer.
     *  \param[in]  d_data         Device-side buffer.
     *  \param[in]  dataSize       Device-side data allocation size.
     *  \param[in]  atomLocality   If all, local or non-local ranges should be copied.
     *  \param[in]  deviceStream   GPU stream to execute copy in.
     */
    void copyFromDevice(gmx::ArrayRef<gmx::RVec> h_data,
                        DeviceBuffer<RVec>       d_data,
                        int                      dataSize,
                        AtomLocality             atomLocality,
                        const DeviceStream&      deviceStream);

    /*! \brief Performs the clearing of data in device buffer.
     *
     * \todo Template on locality.
     *
     *  \param[out] d_data         Device-side buffer.
     *  \param[in]  dataSize       Device-side data allocation size.
     *  \param[in]  atomLocality   If all, local or non-local ranges should be cleared.
     *  \param[in]  deviceStream   GPU stream to execute copy in.
     */
    void clearOnDevice(DeviceBuffer<RVec>  d_data,
                       int                 dataSize,
                       AtomLocality        atomLocality,
                       const DeviceStream& deviceStream) const;
};

} // namespace gmx

#endif // GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H
