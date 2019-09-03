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
 * \brief Declaration of low-level functions and fields of GPU state propagator object.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H
#define GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H

#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/utility/classhelpers.h"

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
         * \note \p commandStream and \p gpuContext are allowed to be nullptr if
         *       StatePropagatorDataGpu is not used in the OpenCL run (e.g. if PME
         *       does not run on the GPU).
         *
         * \todo Make CommandStream visible in the CPU parts of the code so we
         *       will not have to pass a void*.
         * \todo Make a Context object visible in CPU parts of the code so we
         *       will not have to pass a void*.
         *
         *  \param[in] commandStream  GPU stream, nullptr allowed.
         *  \param[in] gpuContext     GPU context, nullptr allowed.
         *  \param[in] transferKind   H2D/D2H transfer call behavior (synchronous or not).
         *  \param[in] paddingSize    Padding size for coordinates buffer.
         */
        Impl(const void        *commandStream,
             const void        *gpuContext,
             GpuApiCallBehavior transferKind,
             int                paddingSize);

        ~Impl();


        /*! \brief Set the ranges for local and non-local atoms and reallocates buffers.
         *
         * The coordinates buffer is reallocated with the padding added at the end. The
         * size of padding is set by the constructor.
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
        std::tuple<int, int> getAtomRangesFromAtomLocality(AtomLocality  atomLocality);


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
        void copyCoordinatesToGpu(gmx::ArrayRef<const gmx::RVec>  h_x,
                                  AtomLocality                    atomLocality);

        /*! \brief Copy positions from the GPU memory.
         *
         *  \param[in] h_x           Positions buffer in the host memory.
         *  \param[in] atomLocality  Locality of the particles to copy.
         */
        void copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec>  h_x,
                                    AtomLocality              atomLocality);


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
        void copyVelocitiesToGpu(gmx::ArrayRef<const gmx::RVec>  h_v,
                                 AtomLocality                    atomLocality);

        /*! \brief Copy velocities from the GPU memory.
         *
         *  \param[in] h_v           Velocities buffer in the host memory.
         *  \param[in] atomLocality  Locality of the particles to copy.
         */
        void copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec>  h_v,
                                   AtomLocality              atomLocality);


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
        void copyForcesToGpu(gmx::ArrayRef<const gmx::RVec>  h_f,
                             AtomLocality                    atomLocality);

        /*! \brief Copy forces from the GPU memory.
         *
         *  \param[in] h_f           Forces buffer in the host memory.
         *  \param[in] atomLocality  Locality of the particles to copy.
         */
        void copyForcesFromGpu(gmx::ArrayRef<gmx::RVec>  h_f,
                               AtomLocality              atomLocality);

        /*! \brief Synchronize the underlying GPU stream
         */
        void synchronizeStream();

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

        /*! \brief GPU stream.
         * \todo The stream should be set to non-nullptr once the synchronization points are restored
         */
        CommandStream        commandStream_              = nullptr;
        /*! \brief GPU context (for OpenCL builds)
         * \todo Make a Context class usable in CPU code
         */
        Context              gpuContext_                 = nullptr;
        //! Default GPU calls behavior
        GpuApiCallBehavior   transferKind_               = GpuApiCallBehavior::Async;
        //! Padding size for the coordinates buffer
        int                  paddingSize_                = 0;

        //! Number of local atoms
        int                  numAtomsLocal_              = -1;
        //! Total number of atoms
        int                  numAtomsAll_                = -1;

        //! Device positions buffer
        DeviceBuffer<float>  d_x_;
        //! Number of particles saved in the positions buffer
        int                  d_xSize_                    = -1;
        //! Allocation size for the positions buffer
        int                  d_xCapacity_                = -1;

        //! Device velocities buffer
        DeviceBuffer<float>  d_v_;
        //! Number of particles saved in the velocities buffer
        int                  d_vSize_                    = -1;
        //! Allocation size for the velocities buffer
        int                  d_vCapacity_                = -1;

        //! Device force buffer
        DeviceBuffer<float>  d_f_;
        //! Number of particles saved in the force buffer
        int                  d_fSize_                    = -1;
        //! Allocation size for the force buffer
        int                  d_fCapacity_                = -1;

        /*! \brief Performs the copy of data from host to device buffer.
         *
         * \todo Template on locality.
         *
         * \param[in,out]  d_data        Device-side buffer.
         * \param[in,out]  h_data        Host-side buffer.
         * \param[in]      dataSize      Device-side data allocation size.
         * \param[in]      atomLocality  If all, local or non-local ranges should be copied.
         */
        void copyToDevice(DeviceBuffer<float>                   d_data,
                          const gmx::ArrayRef<const gmx::RVec>  h_data,
                          int                                   dataSize,
                          AtomLocality                          atomLocality);

        /*! \brief Performs the copy of data from device to host buffer.
         *
         * \param[in,out]  h_data        Host-side buffer.
         * \param[in,out]  d_data        Device-side buffer.
         * \param[in]      dataSize      Device-side data allocation size.
         * \param[in]      atomLocality  If all, local or non-local ranges should be copied.
         */
        void copyFromDevice(gmx::ArrayRef<gmx::RVec>  h_data,
                            DeviceBuffer<float>       d_data,
                            int                       dataSize,
                            AtomLocality              atomLocality);
};

}      // namespace gmx

#endif // GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_IMPL_H
