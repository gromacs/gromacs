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
         * \param[in] numIterations    Number of iteration for the correction of the projection.
         * \param[in] expansionOrder   Order of the matrix inversion algorithm.
         */
        LincsCuda(int numIterations,
                  int expansionOrder);
        ~LincsCuda();

        /*! \brief Apply LINCS to the coordinates/velocities stored in CPU memory.
         *
         * This method should not be used in any code-path, where performance is of any value.
         * Only suitable for test and will be removed in future patch sets.
         * Allocates GPU memory, copies data from CPU, applies LINCS to coordinates and,
         * if requested, to velocities, copies the results back, frees GPU memory.
         * Method uses this class data structures which should be filled with set() and setPbc()
         * methods.
         *
         * \todo Remove this method
         *
         * \param[in]     numAtoms          Number of atoms
         * \param[in]     h_x               Coordinates before timestep (in CPU memory)
         * \param[in,out] h_xp              Coordinates after timestep (in CPU memory). The
         *                                  resulting constrained coordinates will be saved here.
         * \param[in]     updateVelocities  If the velocities should be updated.
         * \param[in,out] h_v               Velocities to update (in CPU memory, can be nullptr
         *                                  if not updated)
         * \param[in]     invdt             Reciprocal timestep (to scale Lagrange
         *                                  multipliers when velocities are updated)
         * \param[in]     computeVirial     If virial should be updated.
         * \param[in,out] virialScaled      Scaled virial tensor to be updated.
         */
        void copyApplyCopy(int         numAtoms,
                           const rvec *h_x,
                           rvec       *h_xp,
                           bool        updateVelocities,
                           rvec       *h_v,
                           real        invdt,
                           bool        computeVirial,
                           tensor      virialScaled);

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

        /*! \brief Class with hardware-specific interfaces and implementations.*/
        class Impl;

    private:
        gmx::PrivateImplPointer<Impl> impl_;

};

} // namespace gmx

#endif
