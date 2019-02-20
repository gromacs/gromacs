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
 * \brief Declares CUDA implementation class for LINCS
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_MDLIB_LINCS_CUDA_IMPL_H
#define GMX_MDLIB_LINCS_CUDA_IMPL_H


#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs_cuda.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/topology/idef.h"

namespace gmx
{

/*! \internal \brief Class with interfaces and data for CUDA version of LINCS. */
class LincsCuda::Impl
{

    public:
        /*! \brief Constructor.
         *
         * \param[in] nAtom    Number of atoms
         * \param[in] nIter    Number of iteration for the correction of the projection.
         * \param[in] nOrder   Order of the matrix inversion algorithm.
         */
        Impl(int nAtom,
             int nIter,
             int nOrder);
        /*! \brief Destructor.*/
        ~Impl();

        /*! \brief Apply LINCS.
         *
         * Applies LINCS to coordinates and velocities, stored on GPU.
         * Data at pointers xPrime and v (class fields) change in the GPU
         * memory. The results are not automatically copied back to the CPU
         * memory. Method uses this class data structures which should be
         * updated when needed using update method.
         *
         * \param[in] updateVelocities  If the velocities should be constrained.
         * \param[in] invdt             Inversed timestep (to scale Lagrange
         *                              multipliers when velocities are updated)
         * \param[in] computeVirial     If virial should be updated.
         * \param[in,out] virialScaled  Scaled virial tensor to be updated.
         */
        void apply(bool    updateVelocities,
                   real    invdt,
                   bool    computeVirial,
                   tensor  virialScaled);

        /*! \brief
         * Update data-structures (e.g. after NB search step).
         *
         * Updates the constraints data and copies it to the GPU. Should be
         * called if the particles were sorted, redistributed between domains, etc.
         * This version uses common data formats so it can be called from anywhere
         * in the code. Does not recycle the data preparation routines from the CPU
         * version. Works only with simple case when all the constraints in idef are
         * are handled by a single GPU. Triangles are not handled as special case.
         *
         * Information about constraints is taken from:
         *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indexes (i1, i2)
         *     idef.iparams[T].constr.dA --- target length for constraint of type T
         * From t_mdatom, the code takes:
         *     md.invmass  --- array of inverse square root of masses for each atom in the system.
         *
         * \param[in] idef  Local topology data to get information on constraints from.
         * \param[in] md    Atoms data to get atom masses from.
         */
        void set(const t_idef    &idef,
                 const t_mdatoms &md);

        /*! \brief
         * Update PBC data.
         *
         * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
         *
         * \param[in] *pbc The PBC data in t_pbc format.
         */
        void setPbc(const t_pbc *pbc);

        /*! \brief
         * Copy coordinates from provided CPU location to GPU.
         *
         * Copies the coordinates before the integration step (x) and coordinates
         * after the integration step (xp) from the provided CPU location to GPU.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *x  CPU pointer where coordinates should be copied from.
         * \param[in] *xp CPU pointer where coordinates should be copied from.
         */
        void copyCoordinatesToGpu(const rvec *x, const rvec *xp);

        /*! \brief
         * Copy velocities from provided CPU location to GPU.
         *
         * Nothing is done if the argument provided is a nullptr.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *v  CPU pointer where velocities should be copied from.
         */
        void copyVelocitiesToGpu(const rvec *v);

        /*! \brief
         * Copy coordinates from GPU to provided CPU location.
         *
         * Copies the constrained coordinates to the provided location. The coordinates
         * are assumed to be in float3/fvec format (single precision).
         *
         * \param[out] *xp CPU pointer where coordinates should be copied to.
         */
        void copyCoordinatesFromGpu(rvec *xp);

        /*! \brief
         * Copy velocities from GPU to provided CPU location.
         *
         * The velocities are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *v  Pointer to velocities data.
         */
        void copyVelocitiesFromGpu(rvec *v);

        /*! \brief
         * Set the internal GPU-memory x, xprime and v pointers.
         *
         * Data is not copied. The data are assumed to be in float3/fvec format
         * (float3 is used internally, but the data layout should be identical).
         *
         * \param[in] *xDevice  Pointer to the coordinates before integrator update (on GPU)
         * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)
         * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)
         */
        void setXVPointers(rvec *xDevice, rvec *xpDevice, rvec *vDevice);

    private:

        //! CUDA stream
        cudaStream_t        stream_;

        //! Periodic boundary data
        PbcAiuc             pbcAiuc_;

        //! Number of atoms
        int                 nAtom_;
        //! Order of expansion when inverting the matrix
        int                 nOrder_;
        //! Number of iterations used to correct the projection
        int                 nIter_;

        //! Coordinates before the timestep (in the GPU memory)
        float3             *xDevice_;
        //! Coordinates after the timestep, before constraining (on GPU)
        float3             *xpDevice_;
        //! Velocities of atoms (on GPU)
        float3             *vDevice_;
        //! 1/mass for all atoms (GPU)
        float              *inverseMassesDevice_;

        //! Scaled virial tensor (6 floats: [XX, XY, XZ, YY, YZ, ZZ], CPU)
        std::vector<float>  virialScaledHost_;
        //! Scaled virial tensor (6 floats: [XX, XY, XZ, YY, YZ, ZZ], GPU)
        float              *virialScaledDevice_;

        /*! \brief Maximum total of constraints assigned to this class so far.
         *
         * If the new assigned number is larger, the GPU data arrays are reallocated.
         */
        int                 maxConstraintsNumberSoFar_;

        /*! \brief Total number of threads.
         *
         *  This covers all constraints and gaps in the ends of the thread blocks
         *  that are necessary to avoid inter-block synchronizations.
         */
        int                 nConstraintsThreads_;

        //! List of constrained atoms (CPU memory)
        std::vector<int2>   constraintsHost_;
        //! List of constrained atoms (GPU memory)
        int2               *constraintsDevice_;

        //! Equilibrium distances for the constraints (CPU)
        std::vector<float>  constraintsTargetLengthsHost_;
        //! Equilibrium distances for the constraints (GPU)
        float              *constraintsTargetLengthsDevice_;

        //! Number of constraints, coupled with the current one (CPU)
        std::vector<int>    coupledConstraintsCountsHost_;
        //! Number of constraints, coupled with the current one (GPU)
        int                *coupledConstraintsCountsDevice_;

        //! List of coupled with the current one (CPU)
        std::vector<int>    coupledConstraintsIdxesHost_;
        //! List of coupled with the current one (GPU)
        int                *coupledConstraintsIdxesDevice_;

        //! Elements of the coupling matrix.
        float              *matrixADevice_;

        //! Mass factors (CPU)
        std::vector<float>  massFactorsHost_;
        //! Mass factors (GPU)
        float              *massFactorsDevice_;

};

} // namespace gmx

#endif
