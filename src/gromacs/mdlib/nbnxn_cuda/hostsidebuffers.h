/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Declare HostSideBuffers.
 *
 * \author Mark Abraham<mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_NBNXN_CUDA_HOSTSIDEBUFFERS_H
#define GMX_MDLIB_NBNXN_CUDA_HOSTSIDEBUFFERS_H

#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal
 * \brief Handles allocation of GPU host-side output buffers, and
 * provides a POD-like view to them for use by CUDA non-bonded
 * kernels.
 *
 * \todo Clarify role and scope, and thus name - should this be all
 * host-side buffers, or just output buffers, or something else?
 *
 * \todo Understand why particle forces are handled separately.
 *
 * \todo Once CUDA code can accept C++11 constructors in this
 * declaration, use a PrivateImplPointer like the rest of the C++
 * code.
 *
 * \todo Consider whether this class could also take responsibility
 * for managing the device-to-host transfer.
 */
class HostSideBuffers
{
    public:
        /*! \brief Constructor.
         *
         * \throws std::bad_alloc if any allocation fails. */
        HostSideBuffers();
        //! Destructor.
        ~HostSideBuffers();

    private:
        //! Private implementation class.
        class Impl;
        //! Pointer to implementation object.
        Impl *impl_;
    public:
        //! Host-side location for computed VDW energy.
        ArrayRef<float> vdwEnergy_;
        //! Host-side location for computed electrostatic energy.
        ArrayRef<float> electrostaticEnergy_;
        //! Host-side location for computed shift forces.
        ArrayRef<float> shiftForces_;
};

} // namespace

#endif
