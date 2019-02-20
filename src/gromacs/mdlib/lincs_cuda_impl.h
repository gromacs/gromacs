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
#ifndef GMX_MDLIB_LINCS_CUDA_IMPL_H
#define GMX_MDLIB_LINCS_CUDA_IMPL_H


#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs_cuda.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"

/*! \internal \brief Class with interfaces and data for CUDA version of LINCS.
 *
 * The class provides major interfaces to constrain bonds using LINCS on GPU.
 * Current implementation is developed for H_Bond constraints. Cant handle constraints triangles.
 *
 */
class LincsCuda::Impl
{

    public:
        /*! \brief Constructor.
         *
         * \param[in] nAtom    Number of atoms
         * \param[in] nIter    Number of iteration for the correction of the projection.
         * \param[in] nOrder   Order of the matrix inversion alorithm.
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
         * \param[in] invdt             Inversed timestep (to scale lagrange
         *                              multipliers when velocities are updated)
         * \param[in] bCalcVir          If virial should be updated.
         * \param[in,out] virialScaled  Scaled virial tensor to be updated.
         */
        void apply(bool       updateVelocities,
                   real       invdt,
                   gmx_bool   bCalcVir,
                   tensor     virialScaled);

        /*! \brief
         * Update data-structures (e.g. after NB search step).
         *
         * Updates the constraints data and copies it to the GPU. Should be
         * called if the particles were sorted, redistributed betwen domains, etc.
         * This version uses common data fromats so it can be called from anywhere
         * in the code. Does not recycle the data preparation routines from the CPU
         * version. Works only with simple case when all the constraints in idef are
         * are handled by a single GPU. Triangles are not handled as special case.
         *
         * Information about constraints is taken from:
         *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indeces (i1, i2)
         *     idef.iparams[T].constr.dA --- target length for constraint of type T
         * From t_mdatom, the code takes:
         *     md.invmass  --- array of inverse square root of masses for each atokm in the system.
         *
         * \param[in] idef  Local topology data to get information on constraints from.
         * \param[in] md    Atoms data to get atom masses from.
         */
        void set(const t_idef         &idef,
                 const t_mdatoms      &md);

        /*! \brief
         * Update PBC data.
         *
         * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
         *
         * \param[in] *pbc The PBC data in t_pbc format.
         */
        void setPbc(t_pbc *pbc);

        /*! \brief
         * Copy coordinates and velocities from provided CPU location to GPU.
         *
         * Copies the coordinates before the integration step (x), coordinates
         * after the integration step (xp) and velocities (v) from the provided
         * CPU location to GPU. The data are assumed to be in float3/fvec format
         * (single precision).
         *
         * \param[in] *x  CPU pointer where coordinates should be copied from.
         * \param[in] *xp CPU pointer where coordinates should be copied from.
         * \param[in] *v  CPU pointer where velocities should be copied from.
         */
        void copyCoordinatesToGpu(const rvec * x, const rvec * xp, const rvec * v);

        /*! \brief
         * Copy coordinates from GPU to provided CPU location.
         *
         * Copies the constrained coordinates to the provided location. The coordinates
         * are assumed to be in float3/fvec format (single precision).
         *
         * \param[out] *xp CPU pointer where coordinates should be copied to.
         */
        void copyCoordinatesFromGpu(rvec * xp);

        /*! \brief
         * Copy velocities from GPU to provided CPU location.
         *
         * The velocities are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *v  Pointer to velocities data.
         */
        void copyVelocitiesFromGpu(rvec * v);

        /*! \brief
         * Set the internal GPU-memory x, xprime and v pointers.
         *
         * Data is not copied. The data are assumed to be in float3/fvec format
         * (float3 is used internaly, but the data layout should be identical).
         *
         * \param[in] *xDevice  Pointer to the coordinates before integrator update (on GPU)
         * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)
         * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)
         */
        void setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice);

    private:

        cudaStream_t     stream;                           //!< CUDA stream

        PbcAiuc          pbcAiuc;                          //!< Periodic boundary data

        int              nAtom;                            //!< Number of atoms

        int              nOrder;                           //!< Order of expansion when inversing the matrix
        int              nIter;                            //!< Number of iterations used to correct the projection

        float3          *xDevice;                          //!< Coordinates before the timestep (in the GPU memory)
        float3          *xpDevice;                         //!< Coordinates after the timestep, before constraining (on GPU).
                                                           //   These will be changed upon constraining.
        float3          *vDevice;                          //!< Velocities of atoms (on GPU)

        real            *invmassDevice;                    //!< 1/mass for all atoms (GPU)
        float           *mlambdaDevice;                    //!< Scaled lagrange multipliers (GPU)

        real            *virialScaledHost;                 //!< Scaled virial tensor (9 reals, GPU)
        real            *virialScaledDevice;               //!< Scaled virial tensor (9 reals, GPU)

        int              maxConstraintsNumberSoFar;        //!< Maximum total of constraints assigned to this class so far.
                                                           //   If the new assigned number is larger, the GPU data arrays are reallocated.
        int              nConstraintsThreads;              //!< Total number of threads, which covers all constraints and gaps in the
                                                           //   ends of the thread blocks that are nessesary to avoid inter-block
                                                           //   syncronizations.

        std::vector<int2>  constraintsHost;                //!< List of constrained atoms (CPU memory)
        int2              *constraintsDevice;              //!< List of constrained atoms (GPU memory)

        std::vector<float> constraintsR0Host;              //!< Equilibrium distances for the constraints (CPU)
        float             *constraintsR0Device;            //!< Equilibrium distances for the constraints (GPU)

        std::vector<int>   coupledConstraintsCountsHost;   //!< Number of constraints, coupled with the current one (CPU)
        int               *coupledConstraintsCountsDevice; //!< Number of constraints, coupled with the current one (GPU)

        std::vector<int>   coupledConstraintsIdxesHost;    //!< List of coupled with the current one (CPU)
        int               *coupledConstraintsIdxesDevice;  //!< List of coupled with the current one (GPU)

        float             *matrixADevice;                  //!< Elements of the coupling matrix.

        std::vector<float> massFactorsHost;                //!< Mass factors (CPU)
        float            * massFactorsDevice;              //!< Mass factors (GPU)

};

#endif
