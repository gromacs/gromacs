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
#ifndef GMX_MDLIB_LINCS_CUDA_H
#define GMX_MDLIB_LINCS_CUDA_H


#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/classhelpers.h"

/*! \brief Class with interfaces and data for CUDA version of LINCS.
 *
 * The class provides major interfaces to constrain bonds using LINCS on GPU.
 * Current implementation is developed for H_Bond constraints. Cant handle constraints triangles.
 *
 *
 */
class LincsCuda
{

    public:
        /*! \brief Constructor.
         *  Initializes objects
         */
        LincsCuda(int N,
                  int nOrder,
                  int nIter);
        ~LincsCuda();

        void apply(bool       updateVelocities,
                   real       invdt,
                   gmx_bool   bCalcVir,
                   tensor     virialScaled);
        void setPbc(t_pbc *pbc);
        void set(const t_idef         &idef,
                 const t_mdatoms      &md);
        /*void set(std::vector<gmx::AtomPair>                        atoms,
                    std::vector<int>                                  blnr,
                    std::vector<int>                                  blbnb,
                    std::vector<real>                                 blmf,
                    std::vector<real, gmx::AlignedAllocator<real>>    bllen,
                    const real                                       *invmass);*/
        void copyCoordinatesToGpu(const rvec * x, const rvec * xp, const rvec * v);
        void copyCoordinatesFromGpu(rvec * xp);
        void copyVelocitiesFromGpu(rvec * v);
        void setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice);

    private:
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;

};

#endif
