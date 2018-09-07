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
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

TopologyInformation::TopologyInformation()
    : mtop_(compat::make_unique<gmx_mtop_t>()),
      hasLoadedMtop_(false),
      expandedTopology_(nullptr),
      atoms_ (nullptr),
      bTop_(false), ePBC_(-1)
{
}


TopologyInformation::~TopologyInformation()
{
}

void TopologyInformation::fillFromInputFile(const std::string &filename)
{
    mtop_ = gmx::compat::make_unique<gmx_mtop_t>();
    // TODO When filename is not a .tpr, then using readConfAndAtoms
    // would be efficient for not doing multiple conversions for
    // makeAtomsData. However we'd also need to be able to copy the
    // t_atoms that we'd keep, which we currently can't do.
    // TODO Once there are fewer callers of the file-reading
    // functionality, make them read directly into std::vector.
    rvec *x, *v;
    readConfAndTopology(filename.c_str(), &bTop_, mtop_.get(),
                        &ePBC_, &x, &v,
                        boxtop_);
    xtop_.assign(x, x + mtop_->natoms);
    vtop_.assign(v, v + mtop_->natoms);
    sfree(x);
    sfree(v);
    hasLoadedMtop_ = true;
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

const gmx_localtop_t *TopologyInformation::expandedTopology() const
{
    // Do lazy initialization
    if (expandedTopology_ == nullptr && hasTopology())
    {
        expandedTopology_.reset(gmx_mtop_generate_local_top(mtop_.get(), false));
    }

    return expandedTopology_.get();
}

namespace
{

//! Helps implement lazy initialization.
AtomsDataPtr makeAtoms(const TopologyInformation &top_)
{
    AtomsDataPtr atoms(new t_atoms);
    if (top_.hasTopology())
    {
        *atoms = gmx_mtop_global_atoms(top_.mtop());
    }
    else
    {
        init_atom(atoms.get());
    }
    return atoms;
}

}   // namespace

const t_atoms *TopologyInformation::atoms() const
{
    // Do lazy initialization
    if (atoms_ == nullptr)
    {
        atoms_ = makeAtoms(*this);
    }

    return atoms_.get();
}

AtomsDataPtr TopologyInformation::copyAtoms() const
{
    // Note that we do not return atoms_, so that regardless of
    // whether the user has already used it, or will use it in the
    // future, any transformation operations on the data structure
    // returned here cannot have unintended effects.
    return makeAtoms(*this);
}

ArrayRef<const RVec>
TopologyInformation::x() const
{
    if (xtop_.empty())
    {
        GMX_THROW(APIError("Topology coordinates requested without setting efUseTopX"));
    }
    return xtop_;
}

ArrayRef<const RVec>
TopologyInformation::v() const
{
    if (vtop_.empty())
    {
        GMX_THROW(APIError("Topology coordinates requested without setting efUseTopV"));
    }
    return vtop_;
}

void
TopologyInformation::getBox(matrix box) const
{
    GMX_RELEASE_ASSERT(box != nullptr, "Must have valid box to fill");
    copy_mat(const_cast<rvec *>(boxtop_), box);
}

const char *
TopologyInformation::name() const
{
    if (hasTopology() && mtop_->name)
    {
        return *mtop_->name;
    }
    return nullptr;
}

gmx_rmpbc_t gmx_rmpbc_init(const gmx::TopologyInformation &topInfo)
{
    GMX_RELEASE_ASSERT(topInfo.hasTopology(), "Cannot remove PBC without a topology");

    return gmx_rmpbc_init(&topInfo.expandedTopology()->idef, topInfo.ePBC(), topInfo.mtop()->natoms);
}

} // namespace gmx
