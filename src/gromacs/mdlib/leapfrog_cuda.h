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
 * \brief Declaration of high-level functions of CUDA implementation of Leap-Frog.
 *
 * \todo This should only list interfaces needed for libgromacs clients (e.g.
 *       management of coordinates, velocities and forces should not be here).
 * \todo Reconsider naming towards using "gpu" suffix instead of "cuda".
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LEAPFROG_CUDA_H
#define GMX_MDLIB_LEAPFROG_CUDA_H

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class LeapFrogCuda
{

    public:

        LeapFrogCuda();
        ~LeapFrogCuda();

        /*! \brief Integrate
         *
         * Copies data from CPU to GPU, integrates the equation of motion
         * using Leap-Frog algorithm, copies the result back. Should only
         * be used for testing.
         *
         * \todo This is temporary solution and will be removed in the
         *       following revisions.
         *
         * \param[in] numAtoms  Number of atoms.
         * \param[in]     h_x   Initial coordinates.
         * \param[out]    h_xp  Place to save the resulting coordinates to.
         * \param[in,out] h_v   Velocities (will be updated).
         * \param[in]     h_f   Forces.
         * \param[in]     dt    Timestep.
         */
        void copyIntegrateCopy(int         numAtoms,
                               const rvec *h_x,
                               rvec       *h_xp,
                               rvec       *h_v,
                               const rvec *h_f,
                               real        dt);

        /*! \brief
         * Update PBC data.
         *
         * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
         *
         * \param[in] pbc The PBC data in t_pbc format.
         */
        void setPbc(const t_pbc *pbc);

        /*! \brief Set the integrator
         *
         * Copies inverse masses from CPU to GPU.
         *
         * \param[in] md    MD atoms, from which inverse masses are taken.
         */
        void set(const t_mdatoms &md);

        /*! \brief Class with hardware-specific interfaces and implementations.*/
        class Impl;

    private:
        gmx::PrivateImplPointer<Impl> impl_;

};

} //namespace gmx

#endif
