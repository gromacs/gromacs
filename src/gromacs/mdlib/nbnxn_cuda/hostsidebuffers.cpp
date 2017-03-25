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
/*! \internal \file
 * \brief Define HostSideBuffers for managing allocation and
 * access to GPU host-side buffers.
 *
 * \author Mark Abraham<mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "hostsidebuffers.h"

#include "gromacs/gpu_utils/gpuhostallocator.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

//! Convenience alias.
template <typename T>
using GpuHostVector = std::vector<T, GpuHostAllocator<T> >;

/*! \libinternal
 * \brief Manage memory allocation of GPU host-side buffers
 *
 * \todo Consider consolidating these buffers into a single page for
 * more efficient usage, including fewer API calls and branches.
 *
 * \todo CUDA unit test to assert that a vector of float maps to
 * a (3x smaller) vector of float3 correctly.
 */
class HostSideBuffers::Impl
{
    public:
        //! Host-side storage for computed VDW energy.
        GpuHostVector<float> vdwEnergy_;
        //! Host-side storage for computed electrostatic energy.
        GpuHostVector<float> electrostaticEnergy_;
        /*! \brief Host-side storage for computed shift forces.
         *
         * \todo Consider using float3 for this, for simplicity and
         * maintainability, once CUDA can compile as C++11. But does
         * that make sense also for whatever OpenCL does?
         */
        GpuHostVector<float> shiftForces_;
        /*! \brief Constructor.
         *
         * \throws std::bad_alloc if any allocation fails. */
        Impl() : vdwEnergy_(1),
                 electrostaticEnergy_(1),
                 shiftForces_(SHIFTS * 3)
        {
        }
};

//! Helper function for constructing the views on the storage.
template <typename T>
ArrayRef<T> arrayRefFromGpuHostVector(GpuHostVector<T> &v)
{
    return arrayRefFromPointers(v.data(), v.data() + v.size());
}

HostSideBuffers::HostSideBuffers() : impl_(new Impl),
                                     vdwEnergy_(arrayRefFromGpuHostVector<float>(impl_->vdwEnergy_)),
                                     electrostaticEnergy_(arrayRefFromGpuHostVector<float>(impl_->electrostaticEnergy_)),
                                     shiftForces_(arrayRefFromGpuHostVector<float>(impl_->shiftForces_))
{
}

HostSideBuffers::~HostSideBuffers()
{
    delete impl_;
}

} // namespace
