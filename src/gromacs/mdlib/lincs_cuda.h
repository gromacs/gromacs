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
/*! \libinternal \file
 *
 * \brief Declaration of high-level functions of CUDA implementation of LINCS.
 *
 * \todo This should only list interfaces needed for libgromacs clients.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LINCS_CUDA_H
#define GMX_MDLIB_LINCS_CUDA_H


#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class LincsCuda
{

    public:
        /*! \brief Constructor.
         *
         * \param[in] numAtoms         Number of atoms
         * \param[in] numIterations    Number of iteration for the correction of the projection.
         * \param[in] expansionOrder   Order of the matrix inversion algorithm.
         */
        LincsCuda(int numAtoms,
                  int numIterations,
                  int expansionOrder);
        ~LincsCuda();

        /*! \brief Apply LINCS.
         *
         * Applies LINCS to coordinates and velocities, Method uses this class data structures
         * which should be updated when needed using update method.
         *
         * \param[in] updateVelocities  If the velocities should be constrained.
         * \param[in] invdt             Reciprocal timestep (to scale Lagrange
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
         * Updates the constraints data. Should be called if the particles were sorted,
         * redistributed between domains, etc.
         *
         * Information about constraints should be taken from:
         *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indexes (i1, i2)
         *     idef.iparams[T].constr.dA --- target length for constraint of type T
         * From t_mdatom, the code should take:
         *     md.invmass  --- array of inverse square root of masses for each atom in the system.
         *
         * \param[in] idef  Local topology data to get information on constraints from.
         * \param[in] md    Atoms data to get atom masses from.
         */
        void set(const t_idef         &idef,
                 const t_mdatoms      &md);

        /*! \brief
         * Update PBC data.
         *
         * \param[in] pbc  The PBC data in t_pbc format.
         */
        void setPbc(const t_pbc *pbc);

        /*! \brief
         * Copy coordinates from provided CPU location to GPU.
         *
         * Copies the coordinates before the integration step (x) and coordinates
         * after the integration step (xp) from the provided CPU location to GPU.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_x   CPU pointer where coordinates should be copied from.
         * \param[in] h_xp  CPU pointer where coordinates should be copied from.
         */
        void copyCoordinatesToGpu(const rvec *h_x, const rvec *h_xp);

        /*! \brief
         * Copy velocities from provided CPU location to GPU.
         *
         * Nothing is done if the argument provided is a nullptr.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_v  CPU pointer where velocities should be copied from.
         */
        void copyVelocitiesToGpu(const rvec *h_v);

        /*! \brief
         * Copy coordinates from GPU to provided CPU location.
         *
         * Copies the constrained coordinates to the provided location. The coordinates
         * are assumed to be in float3/fvec format (single precision).
         *
         * \param[out] h_xp  CPU pointer where coordinates should be copied to.
         */
        void copyCoordinatesFromGpu(rvec *h_xp);

        /*! \brief
         * Copy velocities from GPU to provided CPU location.
         *
         * The velocities are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] h_v  Pointer to velocities data.
         */
        void copyVelocitiesFromGpu(rvec *h_v);

        /*! \brief
         * Set the internal GPU-memory x, xprime and v pointers.
         *
         * Data is not copied. The data are assumed to be in float3/fvec format
         * (float3 is used internally, but the data layout should be identical).
         *
         * \param[in] d_x  Pointer to the coordinates before integrator update (on GPU)
         * \param[in] d_xp Pointer to the coordinates after integrator update, before update (on GPU)
         * \param[in] d_v  Pointer to the velocities before integrator update (on GPU)
         */
        void setXVPointers(rvec *d_x, rvec *d_xp, rvec *d_v);

    private:
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;

};

} // namespace gmx

#endif
