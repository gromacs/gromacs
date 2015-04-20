/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

TopologyManager::TopologyManager()
    : top_(NULL), frame_(NULL)
{
}

TopologyManager::~TopologyManager()
{
    if (top_ != NULL)
    {
        free_t_atoms(&top_->atoms, TRUE);
        done_top(top_);
        sfree(top_);
    }

    if (frame_ != NULL)
    {
        sfree(frame_->x);
        sfree(frame_->v);
        sfree(frame_->f);
        sfree(frame_);
    }
}

void TopologyManager::requestFrame()
{
    GMX_RELEASE_ASSERT(top_ == NULL,
                       "Frame must be requested before initializing topology");
    if (frame_ == NULL)
    {
        snew(frame_, 1);
    }
}

void TopologyManager::requestVelocities()
{
    GMX_RELEASE_ASSERT(frame_ != NULL,
                       "Velocities requested before requesting a frame");
    frame_->bV = TRUE;
    if (frame_->natoms > 0)
    {
        snew(frame_->v, frame_->natoms);
    }
}

void TopologyManager::requestForces()
{
    GMX_RELEASE_ASSERT(frame_ != NULL,
                       "Forces requested before requesting a frame");
    frame_->bF = TRUE;
    if (frame_->natoms > 0)
    {
        snew(frame_->f, frame_->natoms);
    }
}

void TopologyManager::loadTopology(const char *filename)
{
    char    title[STRLEN];
    int     ePBC;
    rvec   *xtop = NULL;
    matrix  box;

    GMX_RELEASE_ASSERT(top_ == NULL, "Topology initialized more than once");
    snew(top_, 1);
    read_tps_conf(gmx::test::TestFileManager::getInputFilePath(filename).c_str(),
                  title, top_, &ePBC, frame_ != NULL ? &xtop : NULL,
                  NULL, box, FALSE);

    if (frame_ != NULL)
    {
        frame_->flags  = TRX_NEED_X;
        frame_->natoms = top_->atoms.nr;
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
    GMX_RELEASE_ASSERT(top_ == NULL, "Topology initialized more than once");
    snew(top_, 1);
    init_t_atoms(&top_->atoms, count, FALSE);
    for (int i = 0; i < count; ++i)
    {
        top_->atoms.atom[i].m = (i % 3 == 0 ? 2.0 : 1.0);
    }
    if (frame_ != NULL)
    {
        frame_->flags  = TRX_NEED_X;
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

void TopologyManager::initAtomTypes(int count, const char *const types[])
{
    GMX_RELEASE_ASSERT(top_ != NULL, "Topology not initialized");
    atomtypes_.reserve(count);
    for (int i = 0; i < count; ++i)
    {
        atomtypes_.push_back(gmx_strdup(types[i]));
    }
    snew(top_->atoms.atomtype, top_->atoms.nr);
    for (int i = 0, j = 0; i < top_->atoms.nr; ++i, ++j)
    {
        if (j == count)
        {
            j = 0;
        }
        top_->atoms.atomtype[i] = &atomtypes_[j];
    }
}

void TopologyManager::initUniformResidues(int residueSize)
{
    GMX_RELEASE_ASSERT(top_ != NULL, "Topology not initialized");
    int residueIndex = -1;
    for (int i = 0; i < top_->atoms.nr; ++i)
    {
        if (i % residueSize == 0)
        {
            ++residueIndex;
        }
        top_->atoms.atom[i].resind = residueIndex;
    }
}

void TopologyManager::initUniformMolecules(int moleculeSize)
{
    GMX_RELEASE_ASSERT(top_ != NULL, "Topology not initialized");
    int index = 0;
    top_->mols.nalloc_index = (top_->atoms.nr + moleculeSize - 1) / moleculeSize + 1;
    snew(top_->mols.index, top_->mols.nalloc_index);
    top_->mols.nr = 0;
    while (index < top_->atoms.nr)
    {
        top_->mols.index[top_->mols.nr] = index;
        ++top_->mols.nr;
        index += moleculeSize;
    }
    top_->mols.index[top_->mols.nr] = top_->atoms.nr;
}

} // namespace test
} // namespace gmx
