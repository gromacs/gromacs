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
/*! \internal \file
 * \brief
 * Implements classes in topologyinformation.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "topologyinformation.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

TopologyInformation::TopologyInformation()
    : mtop_(nullptr), top_(nullptr), bTop_(false), xtop_(nullptr), boxtop_(), ePBC_(-1)
{
}


TopologyInformation::~TopologyInformation()
{
    done_top_mtop(top_, mtop_.get());
    sfree(top_);
    sfree(xtop_);
}


t_topology *TopologyInformation::topology() const
{
    if (top_ == nullptr && mtop_ != nullptr)
    {
        snew(top_, 1);
        *top_ = gmx_mtop_t_to_t_topology(mtop_.get(), false);
    }
    return top_;
}

void TopologyInformation::fillFromInputFile(const std::string &filename)
{
    mtop_ = gmx::compat::make_unique<gmx_mtop_t>();
    readConfAndTopology(filename.c_str(), &bTop_, mtop_.get(),
                        &ePBC_, &xtop_, nullptr,
                        boxtop_);
    // TODO: Only load this here if the tool actually needs it; selections
    // take care of themselves.
    for (gmx_moltype_t &moltype : mtop_->moltype)
    {
        if (!moltype.atoms.haveMass)
        {
            // Try to read masses from database, be silent about missing masses
            atomsSetMassesBasedOnNames(&moltype.atoms, FALSE);
        }
    }
}

void
TopologyInformation::getTopologyConf(rvec **x, matrix box) const
{
    if (box)
    {
        copy_mat(const_cast<rvec *>(boxtop_), box);
    }
    if (x)
    {
        if (!xtop_)
        {
            *x = nullptr;
            GMX_THROW(APIError("Topology coordinates requested without setting efUseTopX"));
        }
        *x = xtop_;
    }
}

} // namespace gmx
