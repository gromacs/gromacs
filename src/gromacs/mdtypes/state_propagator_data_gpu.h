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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class StatePropagatorDataGpu
{
    public:

        /*! \brief Atom locality indicator: local, non-local, all.
         *
         * \todo This should be managed by a separate object, since the localities
         *       are used here and in buffer ops.
         */
        enum class AtomLocality : int
        {
            Local    = 0, //!< Local atoms
            NonLocal = 1, //!< Non-local atoms
            All      = 2, //!< Both local and non-local atoms
            Count    = 3  //!< The number of atom locality types
        };

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
         * \todo Make \p CommandStream visible in the CPU parts of the code so we
         *       will not have to pass a void*.
         * \todo Make \p Context visible in CPU parts of the code so we will not
         *       have to pass a void*.
         *
         *  \param[in] commandStream  GPU stream, nullptr allowed.
         *  \param[in] gpuContext     GPU context, nullptr allowed.
         *  \param[in] transferKind   H2D/D2H transfer call behavior (synchronous or not).
         *  \param[in] paddingSize    Padding size for coordinates buffer.
         */
        StatePropagatorDataGpu(const void        *commandStream,
                               const void        *gpuContext,
                               GpuApiCallBehavior transferKind,
                               int                paddingSize);

        ~StatePropagatorDataGpu();

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
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;

};

}      // namespace gmx

#endif // GMX_MDTYPES_STATE_PROPAGATOR_DATA_GPU_H
