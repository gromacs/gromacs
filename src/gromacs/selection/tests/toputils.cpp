/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements test helper routines from toputils.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "toputils.h"

#include <cstring>

#include <algorithm>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

TopologyManager::TopologyManager()
    : mtop_(nullptr), frame_(nullptr)
{
}

TopologyManager::~TopologyManager()
{
    if (mtop_ != nullptr)
    {
        done_mtop(mtop_);
        sfree(mtop_);
    }

    if (frame_ != nullptr)
    {
        sfree(frame_->x);
        sfree(frame_->v);
        sfree(frame_->f);
        sfree(frame_->index);
        sfree(frame_);
    }

    for (char *atomtype : atomtypes_)
    {
        sfree(atomtype);
    }
}

void TopologyManager::requestFrame()
{
    GMX_RELEASE_ASSERT(mtop_ == nullptr,
                       "Frame must be requested before initializing topology");
    if (frame_ == nullptr)
    {
        snew(frame_, 1);
    }
}

void TopologyManager::requestVelocities()
{
    GMX_RELEASE_ASSERT(frame_ != nullptr,
                       "Velocities requested before requesting a frame");
    frame_->bV = TRUE;
    if (frame_->natoms > 0)
    {
        snew(frame_->v, frame_->natoms);
    }
}

void TopologyManager::requestForces()
{
    GMX_RELEASE_ASSERT(frame_ != nullptr,
                       "Forces requested before requesting a frame");
    frame_->bF = TRUE;
    if (frame_->natoms > 0)
    {
        snew(frame_->f, frame_->natoms);
    }
}

void TopologyManager::loadTopology(const char *filename)
{
    bool    fullTopology;
    int     ePBC;
    rvec   *xtop = nullptr;
    matrix  box;

    GMX_RELEASE_ASSERT(mtop_ == nullptr, "Topology initialized more than once");
    snew(mtop_, 1);
    readConfAndTopology(
            gmx::test::TestFileManager::getInputFilePath(filename).c_str(),
            &fullTopology, mtop_, &ePBC, frame_ != nullptr ? &xtop : nullptr,
            nullptr, box);

    if (frame_ != nullptr)
    {
        frame_->natoms = mtop_->natoms;
        frame_->bX     = TRUE;
        snew(frame_->x, frame_->natoms);
        std::memcpy(frame_->x, xtop, sizeof(*frame_->x) * frame_->natoms);
        frame_->bBox   = TRUE;
        copy_mat(box, frame_->box);
    }

    sfree(xtop);
}

void TopologyManager::initAtoms(int count)
{
    GMX_RELEASE_ASSERT(mtop_ == nullptr, "Topology initialized more than once");
    snew(mtop_, 1);
    mtop_->nmoltype = 1;
    snew(mtop_->moltype, 1);
    init_t_atoms(&mtop_->moltype[0].atoms, count, FALSE);
    mtop_->nmolblock = 1;
    snew(mtop_->molblock, 1);
    mtop_->molblock[0].type            = 0;
    mtop_->molblock[0].nmol            = 1;
    mtop_->molblock[0].natoms_mol      = count;
    mtop_->natoms                      = count;
    mtop_->maxres_renum                = 0;
    gmx_mtop_finalize(mtop_);
    GMX_RELEASE_ASSERT(mtop_->maxres_renum == 0, "maxres_renum in mtop can be modified by an env.var., that is not supported in this test");
    t_atoms &atoms = this->atoms();
    for (int i = 0; i < count; ++i)
    {
        atoms.atom[i].m = (i % 3 == 0 ? 2.0 : 1.0);
    }
    atoms.haveMass = TRUE;
    if (frame_ != nullptr)
    {
        frame_->natoms = count;
        frame_->bX     = TRUE;
        snew(frame_->x, count);
        if (frame_->bV)
        {
            snew(frame_->v, count);
        }
        if (frame_->bF)
        {
            snew(frame_->f, count);
        }
    }
}

void TopologyManager::initAtomTypes(const ArrayRef<const char *const> &types)
{
    GMX_RELEASE_ASSERT(mtop_ != nullptr, "Topology not initialized");
    atomtypes_.reserve(types.size());
    for (const char *type : types)
    {
        atomtypes_.push_back(gmx_strdup(type));
    }
    t_atoms &atoms = this->atoms();
    snew(atoms.atomtype, atoms.nr);
    size_t   j = 0;
    for (int i = 0; i < atoms.nr; ++i, ++j)
    {
        if (j == types.size())
        {
            j = 0;
        }
        atoms.atomtype[i] = &atomtypes_[j];
    }
    atoms.haveType = TRUE;
}

void TopologyManager::initUniformResidues(int residueSize)
{
    GMX_RELEASE_ASSERT(mtop_ != nullptr, "Topology not initialized");
    t_atoms &atoms        = this->atoms();
    int      residueIndex = -1;
    for (int i = 0; i < atoms.nr; ++i)
    {
        if (i % residueSize == 0)
        {
            ++residueIndex;
        }
        atoms.atom[i].resind = residueIndex;
    }
}

void TopologyManager::initUniformMolecules(int moleculeSize)
{
    GMX_RELEASE_ASSERT(mtop_ != nullptr, "Topology not initialized");
    int index = 0;
    mtop_->mols.nalloc_index = (mtop_->natoms + moleculeSize - 1) / moleculeSize + 1;
    srenew(mtop_->mols.index, mtop_->mols.nalloc_index);
    mtop_->mols.nr = 0;
    while (index < mtop_->natoms)
    {
        mtop_->mols.index[mtop_->mols.nr] = index;
        ++mtop_->mols.nr;
        index += moleculeSize;
    }
    mtop_->mols.index[mtop_->mols.nr] = mtop_->natoms;
}

void TopologyManager::initFrameIndices(const ArrayRef<const int> &index)
{
    GMX_RELEASE_ASSERT(frame_ != nullptr, "Frame not initialized");
    GMX_RELEASE_ASSERT(!frame_->bIndex, "Frame atom indices can only be set once");

    frame_->bIndex = TRUE;
    snew(frame_->index, index.size());
    std::copy(index.begin(), index.end(), frame_->index);

    frame_->natoms = index.size();
}

t_atoms &TopologyManager::atoms()
{
    GMX_RELEASE_ASSERT(mtop_ != nullptr, "Topology not initialized");
    GMX_RELEASE_ASSERT(mtop_->natoms == mtop_->moltype[0].atoms.nr,
                       "Test setup assumes all atoms in a single molecule type");
    return mtop_->moltype[0].atoms;
}

} // namespace test
} // namespace gmx
