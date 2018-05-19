/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Helper classes for coordinate frame modification.
 *
 * \author
 * \ingroup module_coordinatedata
 */

#include "gmxpre.h"

#include "outputselector.h"

#include <algorithm>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

void
OutputSelector::initFileOptions(IOptionsContainer * /*options*/)
{
}

void
OutputSelector::checkOptions()
{
}

void
OutputSelector::modifyFrame(const t_trxframe &input)
{
    int              natoms = sel_->atomCount();

    setFrame(input);
    setNatoms(natoms);

    rvec *xmem = nullptr;
    rvec *vmem = nullptr;
    rvec *fmem = nullptr;
    snew(xmem, natoms);
    setFrameCoordinates(xmem);
    if (getbVel())
    {
        snew(vmem, natoms);
        setFrameVelocities(vmem);
    }
    if (getbForce())
    {
        snew(fmem, natoms);
        setFrameForces(fmem);
    }

    for (int i = 0; i < natoms; i++)
    {
        int pos = sel_->position(i).refId();
        copy_rvec(input.x[pos], xmem[i]);
        if (getbVel())
        {
            copy_rvec(input.v[pos], vmem[i]);
        }
        if (getbForce())
        {
            copy_rvec(input.f[pos], fmem[i]);
        }
    }
}

} // namespace gmx
