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
 * \brief Declaration of high-level functions of CUDA implementation of SETTLE.
 *
 * \todo This should only list interfaces needed for libgromacs clients.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_SETTLE_CUDA_H
#define GMX_MDLIB_SETTLE_CUDA_H

#include "gromacs/math/invertmatrix.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

// TODO: Rename to SettleGpu
class SettleCuda
{

    public:

        /*! \brief Structure containing parameters for settles.
         *
         * Contains masses of atoms, distances between them and their pre-computed
         * derivatives (to avoid recomputing them for each water molecule).
         */
        struct SettleParameters;

        /*! \brief Create SETTLE object
         *
         *  Extracts masses for oxygen and hydrogen as well as the O-H and H-H target distances
         *  from the topology data (mtop), check their values for consistency and calls the
         *  following constructor.
         *
         * \param[in] mtop      Topology of the system to get the masses for O and H atoms and
         *                      target O-H and H-H distances. These values are also checked for
         *                      consistency.
         */
        SettleCuda(const gmx_mtop_t &mtop);

        /*! \brief Create SETTLE object
         *
         * \param[in] mO        Mass of the oxygen atom.
         * \param[in] mH        Mass of the hydrogen atom.
         * \param[in] dOH       Target distance for O-H bonds.
         * \param[in] dHH       Target for the distance between two hydrogen atoms.
         */
        SettleCuda(real mO,  real mH,
                   real dOH, real dHH);

        ~SettleCuda();

        /*! \brief Apply SETTLE to the coordinates/velocities stored in CPU memory.
         *
         * This method should not be used in any code-path, where performance is of any value.
         * Only suitable for test and will be removed in future patch sets.
         * Allocates GPU memory, copies data from CPU, applies SETTLE to coordinates and,
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
         * Updates the constraints data and copies it to the GPU. Should be
         * called if the particles were sorted, redistributed between domains, etc.
         * Does not recycle the data preparation routines from the CPU version.
         * All three atoms from single water molecule should be handled by the same GPU.
         *
         * SETTLEs atom ID's are taken from idef.il[F_SETTLE].iatoms.
         *
         * \param[in] idef    System topology
         * \param[in] md      Atoms data. Can be used to update masses if needed (not used now).
         */
        void set(const t_idef     &idef,
                 const t_mdatoms  &md);

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
