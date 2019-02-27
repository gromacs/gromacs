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
#ifndef GMX_MDLIB_SETTLE_CUDA_H
#define GMX_MDLIB_SETTLE_CUDA_H

#include "gromacs/math/invertmatrix.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class SettleCuda
{

    public:

        struct SettleParameters;

        SettleCuda(int               nAtom,
                   const gmx_mtop_t &mtop);

        SettleCuda(int  nAtom,
                   real mO,  real mH,
                   real dOH, real dHH);

        ~SettleCuda();

        void apply(bool       updateVelocities,
                   real       invdt,
                   bool       computeVirial,
                   tensor     virialScaled);

        void set(const t_idef     &idef,
                 const t_mdatoms  &md);

        void setPbc(const t_pbc *pbc);

        void copyCoordinatesToGpu(const rvec * x, const rvec * xp);

        void copyVelocitiesToGpu(const rvec * v);

        void copyCoordinatesFromGpu(rvec * xp);

        void copyVelocitiesFromGpu(rvec * v);

        void setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice);

    private:
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;



};

} // namespace gmx

#endif
