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
        /*! \brief Create Update-Constrain object
         *
         * \param[in] ir        Input record data: LINCS takes number of iterations and order of
         *                      projection from it.
         * \param[in] mtop      Topology of the system: SETTLE gets the masses for O and H atoms
         *                      and target O-H and H-H distances from this object.
         */
        Impl(const t_inputrec     &ir,
             const gmx_mtop_t     &mtop);

        ~Impl();

        /*! \brief Integrate
         *
         * Integrates the equation of motion using Leap-Frog algorithm and applies
         * LINCS and SETTLE constraints.
         * Updates d_xp_ and d_v_ fields of this object.
         * If computeVirial is true, constraints virial is written at the provided pointer.
         * doTempCouple should be true if:
         *   1. The temperature coupling is enabled.
         *   2. This is the temperature coupling step.
         * Parameters virial/lambdas can be nullptr if computeVirial/doTempCouple are false.
         *
         * \param[in]  dt                Timestep
         * \param[in]  updateVelocities  If the velocities should be constrained.
         * \param[in]  computeVirial     If virial should be updated.
         * \param[out] virial            Place to save virial tensor.
         * \param[in]  doTempCouple      If the temperature coupling should be performed.
         * \param[in]  tcstat            Temperature coupling data.
         */
        void integrate(const real                        dt,
                       const bool                        updateVelocities,
                       const bool                        computeVirial,
                       tensor                            virial,
                       const bool                        doTempCouple,
                       gmx::ArrayRef<const t_grp_tcstat> tcstat);

        /*! \brief
         * Update data-structures (e.g. after NB search step).
         *
         * \param[in] idef                 System topology
         * \param[in] md                   Atoms data.
         * \param[in] numTempScaleValues   Number of temperature scaling groups. Set zero for no temperature coupling.
         */
        void set(const t_idef    &idef,
                 const t_mdatoms &md,
                 const int        numTempScaleValues);

        /*! \brief
         * Update PBC data.
         *
         * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
         *
         * \param[in] pbc The PBC data in t_pbc format.
         */
        void setPbc(const t_pbc *pbc);

        /*! \brief
         * Copy coordinates from CPU to GPU.
         *
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_x  CPU pointer where coordinates should be copied from.
         */
        void copyCoordinatesToGpu(const rvec *h_x);

        /*! \brief
         * Copy velocities from CPU to GPU.
         *
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_v  CPU pointer where velocities should be copied from.
         */
        void copyVelocitiesToGpu(const rvec *h_v);

        /*! \brief
         * Copy forces from CPU to GPU.
         *
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_f  CPU pointer where forces should be copied from.
         */
        void copyForcesToGpu(const rvec *h_f);

        /*! \brief
         * Copy coordinates from GPU to CPU.
         *
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[out] h_xp CPU pointer where coordinates should be copied to.
         */
        void copyCoordinatesFromGpu(rvec *h_xp);

        /*! \brief
         * Copy velocities from GPU to CPU.
         *
         * The velocities are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_v  Pointer to velocities data.
         */
        void copyVelocitiesFromGpu(rvec *h_v);

        /*! \brief
         * Copy forces from GPU to CPU.
         *
         * The forces are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_f  Pointer to forces data.
         */
        void copyForcesFromGpu(rvec *h_f);

        /*! \brief
         * Set the internal GPU-memory x, xprime and v pointers.
         *
         * Data is not copied. The data are assumed to be in float3/fvec format
         * (float3 is used internally, but the data layout should be identical).
         *
         * \param[in] d_x   Pointer to the coordinates for the input (on GPU)
         * \param[in] d_xp  Pointer to the coordinates for the output (on GPU)
         * \param[in] d_v   Pointer to the velocities (on GPU)
         * \param[in] d_f   Pointer to the forces (on GPU)
         */
        void setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f);

    private:

        //! CUDA stream
        cudaStream_t        stream_;

        //! Periodic boundary data
        PbcAiuc             pbcAiuc_;

        //! Number of atoms
        int                 numAtoms_;

        //! Coordinates before the timestep (on GPU)
        float3             *d_x_;
        //! Number of elements in coordinates buffer
        int                 numX_                  = -1;
        //! Allocation size for the coordinates buffer
        int                 numXAlloc_             = -1;

        //! Coordinates after the timestep (on GPU).
        float3             *d_xp_;
        //! Number of elements in shifted coordinates buffer
        int                 numXp_                 = -1;
        //! Allocation size for the shifted coordinates buffer
        int                 numXpAlloc_            = -1;

        //! Velocities of atoms (on GPU)
        float3             *d_v_;
        //! Number of elements in velocities buffer
        int                 numV_                  = -1;
        //! Allocation size for the velocities buffer
        int                 numVAlloc_             = -1;

        //! Forces, exerted by atoms (on GPU)
        float3             *d_f_;
        //! Number of elements in forces buffer
        int                 numF_                  = -1;
        //! Allocation size for the forces buffer
        int                 numFAlloc_             = -1;

        //! 1/mass for all atoms (GPU)
        real               *d_inverseMasses_;
        //! Number of elements in reciprocal masses buffer
        int                 numInverseMasses_      = -1;
        //! Allocation size for the reciprocal masses buffer
        int                 numInverseMassesAlloc_ = -1;

        //! Leap-Frog integrator
        std::unique_ptr<LeapFrogCuda>        integrator_;
        //! LINCS CUDA object to use for non-water constraints
        std::unique_ptr<LincsCuda>           lincsCuda_;
        //! SETTLE CUDA object for water constrains
        std::unique_ptr<SettleCuda>          settleCuda_;

};

} // namespace gmx

#endif
