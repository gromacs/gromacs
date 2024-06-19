/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements the ForceBuffers class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdtypes
 */

#include "gmxpre.h"

#include "forcebuffers.h"

#include <algorithm>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/allocator.h"

namespace gmx
{

ForceBuffers::ForceBuffers() :
    force_({}), forceMtsCombined_({}), view_({}, {}, false), useForceMtsCombined_(false)
{
}

ForceBuffers::ForceBuffers(const bool useForceMtsCombined, const PinningPolicy pinningPolicy) :
    force_({}, { pinningPolicy }),
    forceMtsCombined_({}),
    view_({}, {}, useForceMtsCombined),
    useForceMtsCombined_(useForceMtsCombined)
{
}

ForceBuffers::~ForceBuffers() = default;

ForceBuffers& ForceBuffers::operator=(ForceBuffers const& o)
{
    auto oForce = o.view().force();
    resize(oForce.size());
    std::copy(oForce.begin(), oForce.end(), view().force().begin());

    return *this;
}

PinningPolicy ForceBuffers::pinningPolicy() const
{
    return force_.get_allocator().pinningPolicy();
}

void ForceBuffers::resize(int numAtoms)
{
    force_.resizeWithPadding(numAtoms);
    if (useForceMtsCombined_)
    {
        forceMtsCombined_.resizeWithPadding(numAtoms);
    }
    view_ = ForceBuffersView(
            force_.arrayRefWithPadding(), forceMtsCombined_.arrayRefWithPadding(), useForceMtsCombined_);
}

} // namespace gmx
