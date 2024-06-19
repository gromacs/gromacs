/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/topology/mtop_atomloops.h"

#include <cstddef>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

AtomIterator::AtomIterator(const gmx_mtop_t& mtop, int globalAtomNumber) :
    mtop_(&mtop),
    mblock_(0),
    atoms_(&mtop.moltype[mtop.molblock[0].type].atoms),
    currentMolecule_(0),
    highestResidueNumber_(mtop.maxResNumberNotRenumbered()),
    localAtomNumber_(0),
    globalAtomNumber_(globalAtomNumber)
{
    GMX_ASSERT(globalAtomNumber == 0 || globalAtomNumber == mtop.natoms,
               "Starting at other atoms not implemented yet");
}

AtomIterator& AtomIterator::operator++()
{
    localAtomNumber_++;
    globalAtomNumber_++;

    if (localAtomNumber_ >= atoms_->nr)
    {
        if (atoms_->nres <= mtop_->maxResiduesPerMoleculeToTriggerRenumber())
        {
            /* Single residue molecule, increase the count with one */
            highestResidueNumber_ += atoms_->nres;
        }
        currentMolecule_++;
        localAtomNumber_ = 0;
        if (currentMolecule_ >= mtop_->molblock[mblock_].nmol)
        {
            mblock_++;
            if (mblock_ >= mtop_->molblock.size())
            {
                return *this;
            }
            atoms_           = &mtop_->moltype[mtop_->molblock[mblock_].type].atoms;
            currentMolecule_ = 0;
        }
    }
    return *this;
}

bool AtomIterator::operator==(const AtomIterator& o) const
{
    return mtop_ == o.mtop_ && globalAtomNumber_ == o.globalAtomNumber_;
}

const t_atom& AtomProxy::atom() const
{
    return it_->atoms_->atom[it_->localAtomNumber_];
}

int AtomProxy::globalAtomNumber() const
{
    return it_->globalAtomNumber_;
}

const char* AtomProxy::atomName() const
{
    return *(it_->atoms_->atomname[it_->localAtomNumber_]);
}

const char* AtomProxy::residueName() const
{
    int residueIndexInMolecule = it_->atoms_->atom[it_->localAtomNumber_].resind;
    return *(it_->atoms_->resinfo[residueIndexInMolecule].name);
}

int AtomProxy::residueNumber() const
{
    int residueIndexInMolecule = it_->atoms_->atom[it_->localAtomNumber_].resind;
    if (it_->atoms_->nres <= it_->mtop_->maxResiduesPerMoleculeToTriggerRenumber())
    {
        return it_->highestResidueNumber_ + 1 + residueIndexInMolecule;
    }
    else
    {
        return it_->atoms_->resinfo[residueIndexInMolecule].nr;
    }
}

const gmx_moltype_t& AtomProxy::moleculeType() const
{
    return it_->mtop_->moltype[it_->mtop_->molblock[it_->mblock_].type];
}

int AtomProxy::atomNumberInMol() const
{
    return it_->localAtomNumber_;
}

struct gmx_mtop_atomloop_block
{
    const gmx_mtop_t* mtop;
    size_t            mblock;
    const t_atoms*    atoms;
    int               at_local;
};

gmx_mtop_atomloop_block_t gmx_mtop_atomloop_block_init(const gmx_mtop_t& mtop)
{
    struct gmx_mtop_atomloop_block* aloop = nullptr;

    snew(aloop, 1);

    aloop->mtop     = &mtop;
    aloop->mblock   = 0;
    aloop->atoms    = &mtop.moltype[mtop.molblock[aloop->mblock].type].atoms;
    aloop->at_local = -1;

    return aloop;
}

static void gmx_mtop_atomloop_block_destroy(gmx_mtop_atomloop_block_t aloop)
{
    sfree(aloop);
}

gmx_bool gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop, const t_atom** atom, int* nmol)
{
    if (aloop == nullptr)
    {
        gmx_incons("gmx_mtop_atomloop_all_next called without calling gmx_mtop_atomloop_all_init");
    }

    aloop->at_local++;

    if (aloop->at_local >= aloop->atoms->nr)
    {
        aloop->mblock++;
        if (aloop->mblock >= aloop->mtop->molblock.size())
        {
            gmx_mtop_atomloop_block_destroy(aloop);
            return FALSE;
        }
        aloop->atoms    = &aloop->mtop->moltype[aloop->mtop->molblock[aloop->mblock].type].atoms;
        aloop->at_local = 0;
    }

    *atom = &aloop->atoms->atom[aloop->at_local];
    *nmol = aloop->mtop->molblock[aloop->mblock].nmol;

    return TRUE;
}

IListIterator::IListIterator(const gmx_mtop_t& mtop, size_t molblock) :
    mtop_(&mtop), mblock_(molblock)
{
}

IListIterator& IListIterator::operator++()
{
    mblock_++;
    return *this;
}

bool IListIterator::operator==(const IListIterator& o) const
{
    return mtop_ == o.mtop_ && mblock_ == o.mblock_;
}

const InteractionLists& IListProxy::list() const
{
    // one past the end means we want to take the
    // intermolecular list instead.
    if (it_->mblock_ == it_->mtop_->molblock.size())
    {
        return *it_->mtop_->intermolecular_ilist;
    }
    else
    {
        return it_->mtop_->moltype[it_->mtop_->molblock[it_->mblock_].type].ilist;
    }
}

int IListProxy::nmol() const
{
    // one past the end means we want to take the
    // intermolecular list instead.
    if (it_->mblock_ == it_->mtop_->molblock.size())
    {
        return 1;
    }
    else
    {
        return it_->mtop_->molblock[it_->mblock_].nmol;
    }
}

IListRange::IListRange(const gmx_mtop_t& mtop) : begin_(mtop), end_(mtop, mtop.molblock.size())
{
    if (mtop.bIntermolecularInteractions)
    {
        end_ = IListIterator(mtop, mtop.molblock.size() + 1);
    }
}
