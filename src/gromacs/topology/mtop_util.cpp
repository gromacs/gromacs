/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "mtop_util.h"

#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static int gmx_mtop_maxresnr(const gmx_mtop_t* mtop, int maxres_renum)
{
    int maxresnr = 0;

    for (const gmx_moltype_t& moltype : mtop->moltype)
    {
        const t_atoms& atoms = moltype.atoms;
        if (atoms.nres > maxres_renum)
        {
            for (int r = 0; r < atoms.nres; r++)
            {
                if (atoms.resinfo[r].nr > maxresnr)
                {
                    maxresnr = atoms.resinfo[r].nr;
                }
            }
        }
    }

    return maxresnr;
}

static void buildMolblockIndices(gmx_mtop_t* mtop)
{
    mtop->moleculeBlockIndices.resize(mtop->molblock.size());

    int atomIndex          = 0;
    int residueIndex       = 0;
    int residueNumberStart = mtop->maxresnr + 1;
    int moleculeIndexStart = 0;
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
    {
        const gmx_molblock_t& molb         = mtop->molblock[mb];
        MoleculeBlockIndices& indices      = mtop->moleculeBlockIndices[mb];
        const int             numResPerMol = mtop->moltype[molb.type].atoms.nres;

        indices.numAtomsPerMolecule = mtop->moltype[molb.type].atoms.nr;
        indices.globalAtomStart     = atomIndex;
        indices.globalResidueStart  = residueIndex;
        atomIndex += molb.nmol * indices.numAtomsPerMolecule;
        residueIndex += molb.nmol * numResPerMol;
        indices.globalAtomEnd      = atomIndex;
        indices.residueNumberStart = residueNumberStart;
        if (numResPerMol <= mtop->maxres_renum)
        {
            residueNumberStart += molb.nmol * numResPerMol;
        }
        indices.moleculeIndexStart = moleculeIndexStart;
        moleculeIndexStart += molb.nmol;
    }
}

void gmx_mtop_finalize(gmx_mtop_t* mtop)
{
    char* env;

    if (mtop->molblock.size() == 1 && mtop->molblock[0].nmol == 1)
    {
        /* We have a single molecule only, no renumbering needed.
         * This case also covers an mtop converted from pdb/gro/... input,
         * so we retain the original residue numbering.
         */
        mtop->maxres_renum = 0;
    }
    else
    {
        /* We only renumber single residue molecules. Their intra-molecular
         * residue numbering is anyhow irrelevant.
         */
        mtop->maxres_renum = 1;
    }

    env = getenv("GMX_MAXRESRENUM");
    if (env != nullptr)
    {
        sscanf(env, "%d", &mtop->maxres_renum);
    }
    if (mtop->maxres_renum == -1)
    {
        /* -1 signals renumber residues in all molecules */
        mtop->maxres_renum = INT_MAX;
    }

    mtop->maxresnr = gmx_mtop_maxresnr(mtop, mtop->maxres_renum);

    buildMolblockIndices(mtop);
}

void gmx_mtop_count_atomtypes(const gmx_mtop_t* mtop, int state, int typecount[])
{
    for (int i = 0; i < mtop->ffparams.atnr; ++i)
    {
        typecount[i] = 0;
    }
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const t_atoms& atoms = mtop->moltype[molb.type].atoms;
        for (int i = 0; i < atoms.nr; ++i)
        {
            int tpi;
            if (state == 0)
            {
                tpi = atoms.atom[i].type;
            }
            else
            {
                tpi = atoms.atom[i].typeB;
            }
            typecount[tpi] += molb.nmol;
        }
    }
}

int gmx_mtop_num_molecules(const gmx_mtop_t& mtop)
{
    int numMolecules = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        numMolecules += molb.nmol;
    }
    return numMolecules;
}

int gmx_mtop_nres(const gmx_mtop_t* mtop)
{
    int nres = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        nres += molb.nmol * mtop->moltype[molb.type].atoms.nres;
    }
    return nres;
}

AtomIterator::AtomIterator(const gmx_mtop_t& mtop, int globalAtomNumber) :
    mtop_(&mtop),
    mblock_(0),
    atoms_(&mtop.moltype[mtop.molblock[0].type].atoms),
    currentMolecule_(0),
    highestResidueNumber_(mtop.maxresnr),
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
        if (atoms_->nres <= mtop_->maxresnr)
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

AtomIterator AtomIterator::operator++(int)
{
    AtomIterator temp = *this;
    ++(*this);
    return temp;
}

bool AtomIterator::operator==(const AtomIterator& o) const
{
    return mtop_ == o.mtop_ && globalAtomNumber_ == o.globalAtomNumber_;
}

bool AtomIterator::operator!=(const AtomIterator& o) const
{
    return !(*this == o);
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
    if (it_->atoms_->nres <= it_->mtop_->maxres_renum)
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

typedef struct gmx_mtop_atomloop_block
{
    const gmx_mtop_t* mtop;
    size_t            mblock;
    const t_atoms*    atoms;
    int               at_local;
} t_gmx_mtop_atomloop_block;

gmx_mtop_atomloop_block_t gmx_mtop_atomloop_block_init(const gmx_mtop_t* mtop)
{
    struct gmx_mtop_atomloop_block* aloop;

    snew(aloop, 1);

    aloop->mtop     = mtop;
    aloop->mblock   = 0;
    aloop->atoms    = &mtop->moltype[mtop->molblock[aloop->mblock].type].atoms;
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

typedef struct gmx_mtop_ilistloop
{
    const gmx_mtop_t* mtop;
    int               mblock;
} t_gmx_mtop_ilist;

gmx_mtop_ilistloop_t gmx_mtop_ilistloop_init(const gmx_mtop_t* mtop)
{
    struct gmx_mtop_ilistloop* iloop;

    snew(iloop, 1);

    iloop->mtop   = mtop;
    iloop->mblock = -1;

    return iloop;
}

gmx_mtop_ilistloop_t gmx_mtop_ilistloop_init(const gmx_mtop_t& mtop)
{
    return gmx_mtop_ilistloop_init(&mtop);
}

static void gmx_mtop_ilistloop_destroy(gmx_mtop_ilistloop_t iloop)
{
    sfree(iloop);
}

const InteractionLists* gmx_mtop_ilistloop_next(gmx_mtop_ilistloop_t iloop, int* nmol)
{
    if (iloop == nullptr)
    {
        gmx_incons("gmx_mtop_ilistloop_next called without calling gmx_mtop_ilistloop_init");
    }

    iloop->mblock++;
    if (iloop->mblock >= gmx::ssize(iloop->mtop->molblock))
    {
        if (iloop->mblock == gmx::ssize(iloop->mtop->molblock) && iloop->mtop->bIntermolecularInteractions)
        {
            *nmol = 1;
            return iloop->mtop->intermolecular_ilist.get();
        }

        gmx_mtop_ilistloop_destroy(iloop);
        return nullptr;
    }

    *nmol = iloop->mtop->molblock[iloop->mblock].nmol;

    return &iloop->mtop->moltype[iloop->mtop->molblock[iloop->mblock].type].ilist;
}
typedef struct gmx_mtop_ilistloop_all
{
    const gmx_mtop_t* mtop;
    size_t            mblock;
    int               mol;
    int               a_offset;
} t_gmx_mtop_ilist_all;

gmx_mtop_ilistloop_all_t gmx_mtop_ilistloop_all_init(const gmx_mtop_t* mtop)
{
    struct gmx_mtop_ilistloop_all* iloop;

    snew(iloop, 1);

    iloop->mtop     = mtop;
    iloop->mblock   = 0;
    iloop->mol      = -1;
    iloop->a_offset = 0;

    return iloop;
}

static void gmx_mtop_ilistloop_all_destroy(gmx_mtop_ilistloop_all_t iloop)
{
    sfree(iloop);
}

const InteractionLists* gmx_mtop_ilistloop_all_next(gmx_mtop_ilistloop_all_t iloop, int* atnr_offset)
{

    if (iloop == nullptr)
    {
        gmx_incons(
                "gmx_mtop_ilistloop_all_next called without calling gmx_mtop_ilistloop_all_init");
    }

    if (iloop->mol >= 0)
    {
        iloop->a_offset += iloop->mtop->moleculeBlockIndices[iloop->mblock].numAtomsPerMolecule;
    }

    iloop->mol++;

    /* Inter-molecular interactions, if present, are indexed with
     * iloop->mblock == iloop->mtop->nmolblock, thus we should separately
     * check for this value in this conditional.
     */
    if (iloop->mblock == iloop->mtop->molblock.size()
        || iloop->mol >= iloop->mtop->molblock[iloop->mblock].nmol)
    {
        iloop->mblock++;
        iloop->mol = 0;
        if (iloop->mblock >= iloop->mtop->molblock.size())
        {
            if (iloop->mblock == iloop->mtop->molblock.size() && iloop->mtop->bIntermolecularInteractions)
            {
                *atnr_offset = 0;
                return iloop->mtop->intermolecular_ilist.get();
            }

            gmx_mtop_ilistloop_all_destroy(iloop);
            return nullptr;
        }
    }

    *atnr_offset = iloop->a_offset;

    return &iloop->mtop->moltype[iloop->mtop->molblock[iloop->mblock].type].ilist;
}

int gmx_mtop_ftype_count(const gmx_mtop_t* mtop, int ftype)
{
    gmx_mtop_ilistloop_t iloop;
    int                  n, nmol;

    n = 0;

    iloop = gmx_mtop_ilistloop_init(mtop);
    while (const InteractionLists* il = gmx_mtop_ilistloop_next(iloop, &nmol))
    {
        n += nmol * (*il)[ftype].size() / (1 + NRAL(ftype));
    }

    if (mtop->bIntermolecularInteractions)
    {
        n += (*mtop->intermolecular_ilist)[ftype].size() / (1 + NRAL(ftype));
    }

    return n;
}

int gmx_mtop_ftype_count(const gmx_mtop_t& mtop, int ftype)
{
    return gmx_mtop_ftype_count(&mtop, ftype);
}

int gmx_mtop_interaction_count(const gmx_mtop_t& mtop, const int unsigned if_flags)
{
    int n = 0;

    gmx_mtop_ilistloop_t iloop = gmx_mtop_ilistloop_init(mtop);
    int                  nmol;
    while (const InteractionLists* il = gmx_mtop_ilistloop_next(iloop, &nmol))
    {
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & if_flags) == if_flags)
            {
                n += nmol * (*il)[ftype].size() / (1 + NRAL(ftype));
            }
        }
    }

    if (mtop.bIntermolecularInteractions)
    {
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & if_flags) == if_flags)
            {
                n += (*mtop.intermolecular_ilist)[ftype].size() / (1 + NRAL(ftype));
            }
        }
    }

    return n;
}

std::array<int, eptNR> gmx_mtop_particletype_count(const gmx_mtop_t& mtop)
{
    std::array<int, eptNR> count = { { 0 } };

    for (const auto& molblock : mtop.molblock)
    {
        const t_atoms& atoms = mtop.moltype[molblock.type].atoms;
        for (int a = 0; a < atoms.nr; a++)
        {
            count[atoms.atom[a].ptype] += molblock.nmol;
        }
    }

    return count;
}

static void atomcat(t_atoms* dest, const t_atoms* src, int copies, int maxres_renum, int* maxresnr)
{
    int i, j, l, size;
    int srcnr  = src->nr;
    int destnr = dest->nr;

    if (dest->nr == 0)
    {
        dest->haveMass    = src->haveMass;
        dest->haveType    = src->haveType;
        dest->haveCharge  = src->haveCharge;
        dest->haveBState  = src->haveBState;
        dest->havePdbInfo = src->havePdbInfo;
    }
    else
    {
        dest->haveMass    = dest->haveMass && src->haveMass;
        dest->haveType    = dest->haveType && src->haveType;
        dest->haveCharge  = dest->haveCharge && src->haveCharge;
        dest->haveBState  = dest->haveBState && src->haveBState;
        dest->havePdbInfo = dest->havePdbInfo && src->havePdbInfo;
    }

    if (srcnr)
    {
        size = destnr + copies * srcnr;
        srenew(dest->atom, size);
        srenew(dest->atomname, size);
        if (dest->haveType)
        {
            srenew(dest->atomtype, size);
            if (dest->haveBState)
            {
                srenew(dest->atomtypeB, size);
            }
        }
        if (dest->havePdbInfo)
        {
            srenew(dest->pdbinfo, size);
        }
    }
    if (src->nres)
    {
        size = dest->nres + copies * src->nres;
        srenew(dest->resinfo, size);
    }

    /* residue information */
    for (l = dest->nres, j = 0; (j < copies); j++, l += src->nres)
    {
        memcpy(reinterpret_cast<char*>(&(dest->resinfo[l])), reinterpret_cast<char*>(&(src->resinfo[0])),
               static_cast<size_t>(src->nres * sizeof(src->resinfo[0])));
    }

    for (l = destnr, j = 0; (j < copies); j++, l += srcnr)
    {
        memcpy(reinterpret_cast<char*>(&(dest->atom[l])), reinterpret_cast<char*>(&(src->atom[0])),
               static_cast<size_t>(srcnr * sizeof(src->atom[0])));
        memcpy(reinterpret_cast<char*>(&(dest->atomname[l])),
               reinterpret_cast<char*>(&(src->atomname[0])),
               static_cast<size_t>(srcnr * sizeof(src->atomname[0])));
        if (dest->haveType)
        {
            memcpy(reinterpret_cast<char*>(&(dest->atomtype[l])),
                   reinterpret_cast<char*>(&(src->atomtype[0])),
                   static_cast<size_t>(srcnr * sizeof(src->atomtype[0])));
            if (dest->haveBState)
            {
                memcpy(reinterpret_cast<char*>(&(dest->atomtypeB[l])),
                       reinterpret_cast<char*>(&(src->atomtypeB[0])),
                       static_cast<size_t>(srcnr * sizeof(src->atomtypeB[0])));
            }
        }
        if (dest->havePdbInfo)
        {
            memcpy(reinterpret_cast<char*>(&(dest->pdbinfo[l])),
                   reinterpret_cast<char*>(&(src->pdbinfo[0])),
                   static_cast<size_t>(srcnr * sizeof(src->pdbinfo[0])));
        }
    }

    /* Increment residue indices */
    for (l = destnr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < srcnr); i++, l++)
        {
            dest->atom[l].resind = dest->nres + j * src->nres + src->atom[i].resind;
        }
    }

    if (src->nres <= maxres_renum)
    {
        /* Single residue molecule, continue counting residues */
        for (j = 0; (j < copies); j++)
        {
            for (l = 0; l < src->nres; l++)
            {
                (*maxresnr)++;
                dest->resinfo[dest->nres + j * src->nres + l].nr = *maxresnr;
            }
        }
    }

    dest->nres += copies * src->nres;
    dest->nr += copies * src->nr;
}

t_atoms gmx_mtop_global_atoms(const gmx_mtop_t* mtop)
{
    t_atoms atoms;

    init_t_atoms(&atoms, 0, FALSE);

    int maxresnr = mtop->maxresnr;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        atomcat(&atoms, &mtop->moltype[molb.type].atoms, molb.nmol, mtop->maxres_renum, &maxresnr);
    }

    return atoms;
}

/*
 * The cat routines below are old code from src/kernel/topcat.c
 */

static void blockacat(t_blocka* dest, const t_blocka* src, int copies, int dnum, int snum)
{
    int i, j, l, size;
    int destnr  = dest->nr;
    int destnra = dest->nra;

    if (src->nr)
    {
        size = (dest->nr + copies * src->nr + 1);
        srenew(dest->index, size);
    }
    if (src->nra)
    {
        size = (dest->nra + copies * src->nra);
        srenew(dest->a, size);
    }

    for (l = destnr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < src->nr); i++)
        {
            dest->index[l++] = dest->nra + src->index[i];
        }
        dest->nra += src->nra;
    }
    for (l = destnra, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < src->nra); i++)
        {
            dest->a[l++] = dnum + src->a[i];
        }
        dnum += snum;
        dest->nr += src->nr;
    }
    dest->index[dest->nr] = dest->nra;
    dest->nalloc_index    = dest->nr;
    dest->nalloc_a        = dest->nra;
}

static void ilistcat(int ftype, t_ilist* dest, const InteractionList& src, int copies, int dnum, int snum)
{
    int nral, c, i, a;

    nral = NRAL(ftype);

    dest->nalloc = dest->nr + copies * src.size();
    srenew(dest->iatoms, dest->nalloc);

    for (c = 0; c < copies; c++)
    {
        for (i = 0; i < src.size();)
        {
            dest->iatoms[dest->nr++] = src.iatoms[i++];
            for (a = 0; a < nral; a++)
            {
                dest->iatoms[dest->nr++] = dnum + src.iatoms[i++];
            }
        }
        dnum += snum;
    }
}

static void set_posres_params(t_idef* idef, const gmx_molblock_t* molb, int i0, int a_offset)
{
    t_ilist*   il;
    int        i1, i, a_molb;
    t_iparams* ip;

    il                          = &idef->il[F_POSRES];
    i1                          = il->nr / 2;
    idef->iparams_posres_nalloc = i1;
    srenew(idef->iparams_posres, idef->iparams_posres_nalloc);
    for (i = i0; i < i1; i++)
    {
        ip = &idef->iparams_posres[i];
        /* Copy the force constants */
        *ip    = idef->iparams[il->iatoms[i * 2]];
        a_molb = il->iatoms[i * 2 + 1] - a_offset;
        if (molb->posres_xA.empty())
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        ip->posres.pos0A[XX] = molb->posres_xA[a_molb][XX];
        ip->posres.pos0A[YY] = molb->posres_xA[a_molb][YY];
        ip->posres.pos0A[ZZ] = molb->posres_xA[a_molb][ZZ];
        if (!molb->posres_xB.empty())
        {
            ip->posres.pos0B[XX] = molb->posres_xB[a_molb][XX];
            ip->posres.pos0B[YY] = molb->posres_xB[a_molb][YY];
            ip->posres.pos0B[ZZ] = molb->posres_xB[a_molb][ZZ];
        }
        else
        {
            ip->posres.pos0B[XX] = ip->posres.pos0A[XX];
            ip->posres.pos0B[YY] = ip->posres.pos0A[YY];
            ip->posres.pos0B[ZZ] = ip->posres.pos0A[ZZ];
        }
        /* Set the parameter index for idef->iparams_posre */
        il->iatoms[i * 2] = i;
    }
}

static void set_fbposres_params(t_idef* idef, const gmx_molblock_t* molb, int i0, int a_offset)
{
    t_ilist*   il;
    int        i1, i, a_molb;
    t_iparams* ip;

    il                            = &idef->il[F_FBPOSRES];
    i1                            = il->nr / 2;
    idef->iparams_fbposres_nalloc = i1;
    srenew(idef->iparams_fbposres, idef->iparams_fbposres_nalloc);
    for (i = i0; i < i1; i++)
    {
        ip = &idef->iparams_fbposres[i];
        /* Copy the force constants */
        *ip    = idef->iparams[il->iatoms[i * 2]];
        a_molb = il->iatoms[i * 2 + 1] - a_offset;
        if (molb->posres_xA.empty())
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        /* Take flat-bottom posres reference from normal position restraints */
        ip->fbposres.pos0[XX] = molb->posres_xA[a_molb][XX];
        ip->fbposres.pos0[YY] = molb->posres_xA[a_molb][YY];
        ip->fbposres.pos0[ZZ] = molb->posres_xA[a_molb][ZZ];
        /* Note: no B-type for flat-bottom posres */

        /* Set the parameter index for idef->iparams_posre */
        il->iatoms[i * 2] = i;
    }
}

/*! \brief Copy idef structure from mtop.
 *
 * Makes a deep copy of an idef data structure from a gmx_mtop_t.
 * Used to initialize legacy topology types.
 *
 * \param[in] mtop Reference to input mtop.
 * \param[in] idef Pointer to idef to populate.
 * \param[in] mergeConstr Decide if constraints will be merged.
 * \param[in] freeEnergyInteractionsAtEnd Decide if free energy stuff should
 *              be added at the end.
 */
static void copyIdefFromMtop(const gmx_mtop_t& mtop, t_idef* idef, bool freeEnergyInteractionsAtEnd, bool mergeConstr)
{
    const gmx_ffparams_t* ffp = &mtop.ffparams;

    idef->ntypes = ffp->numTypes();
    idef->atnr   = ffp->atnr;
    /* we can no longer copy the pointers to the mtop members,
     * because they will become invalid as soon as mtop gets free'd.
     * We also need to make sure to only operate on valid data!
     */

    if (!ffp->functype.empty())
    {
        snew(idef->functype, ffp->functype.size());
        std::copy(ffp->functype.data(), ffp->functype.data() + ffp->functype.size(), idef->functype);
    }
    else
    {
        idef->functype = nullptr;
    }
    if (!ffp->iparams.empty())
    {
        snew(idef->iparams, ffp->iparams.size());
        std::copy(ffp->iparams.data(), ffp->iparams.data() + ffp->iparams.size(), idef->iparams);
    }
    else
    {
        idef->iparams = nullptr;
    }
    idef->iparams_posres          = nullptr;
    idef->iparams_posres_nalloc   = 0;
    idef->iparams_fbposres        = nullptr;
    idef->iparams_fbposres_nalloc = 0;
    idef->fudgeQQ                 = ffp->fudgeQQ;
    idef->cmap_grid               = new gmx_cmap_t;
    *idef->cmap_grid              = ffp->cmap_grid;
    idef->ilsort                  = ilsortUNKNOWN;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        idef->il[ftype].nr     = 0;
        idef->il[ftype].nalloc = 0;
        idef->il[ftype].iatoms = nullptr;
    }

    int natoms = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        int srcnr  = molt.atoms.nr;
        int destnr = natoms;

        int nposre_old   = idef->il[F_POSRES].nr;
        int nfbposre_old = idef->il[F_FBPOSRES].nr;
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (mergeConstr && ftype == F_CONSTR && molt.ilist[F_CONSTRNC].size() > 0)
            {
                /* Merge all constrains into one ilist.
                 * This simplifies the constraint code.
                 */
                for (int mol = 0; mol < molb.nmol; mol++)
                {
                    ilistcat(ftype, &idef->il[F_CONSTR], molt.ilist[F_CONSTR], 1,
                             destnr + mol * srcnr, srcnr);
                    ilistcat(ftype, &idef->il[F_CONSTR], molt.ilist[F_CONSTRNC], 1,
                             destnr + mol * srcnr, srcnr);
                }
            }
            else if (!(mergeConstr && ftype == F_CONSTRNC))
            {
                ilistcat(ftype, &idef->il[ftype], molt.ilist[ftype], molb.nmol, destnr, srcnr);
            }
        }
        if (idef->il[F_POSRES].nr > nposre_old)
        {
            /* Executing this line line stops gmxdump -sys working
             * correctly. I'm not aware there's an elegant fix. */
            set_posres_params(idef, &molb, nposre_old / 2, natoms);
        }
        if (idef->il[F_FBPOSRES].nr > nfbposre_old)
        {
            set_fbposres_params(idef, &molb, nfbposre_old / 2, natoms);
        }

        natoms += molb.nmol * srcnr;
    }

    if (mtop.bIntermolecularInteractions)
    {
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            ilistcat(ftype, &idef->il[ftype], (*mtop.intermolecular_ilist)[ftype], 1, 0, mtop.natoms);
        }
    }

    if (freeEnergyInteractionsAtEnd && gmx_mtop_bondeds_free_energy(&mtop))
    {
        std::vector<real> qA(mtop.natoms);
        std::vector<real> qB(mtop.natoms);
        for (const AtomProxy atomP : AtomRange(mtop))
        {
            const t_atom& local = atomP.atom();
            int           index = atomP.globalAtomNumber();
            qA[index]           = local.q;
            qB[index]           = local.qB;
        }
        gmx_sort_ilist_fe(idef, qA.data(), qB.data());
    }
    else
    {
        idef->ilsort = ilsortNO_FE;
    }
}

/*! \brief Copy atomtypes from mtop
 *
 * Makes a deep copy of t_atomtypes from gmx_mtop_t.
 * Used to initialize legacy topology types.
 *
 * \param[in] mtop Reference to input mtop.
 * \param[in] atomtypes Pointer to atomtypes to populate.
 */
static void copyAtomtypesFromMtop(const gmx_mtop_t& mtop, t_atomtypes* atomtypes)
{
    atomtypes->nr = mtop.atomtypes.nr;
    if (mtop.atomtypes.atomnumber)
    {
        snew(atomtypes->atomnumber, mtop.atomtypes.nr);
        std::copy(mtop.atomtypes.atomnumber, mtop.atomtypes.atomnumber + mtop.atomtypes.nr,
                  atomtypes->atomnumber);
    }
    else
    {
        atomtypes->atomnumber = nullptr;
    }
}

/*! \brief Copy excls from mtop.
 *
 * Makes a deep copy of excls(t_blocka) from gmx_mtop_t.
 * Used to initialize legacy topology types.
 *
 * \param[in] mtop  Reference to input mtop.
 * \param[in] excls Pointer to final excls data structure.
 */
static void copyExclsFromMtop(const gmx_mtop_t& mtop, t_blocka* excls)
{
    init_blocka(excls);
    int natoms = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        int srcnr  = molt.atoms.nr;
        int destnr = natoms;

        blockacat(excls, &molt.excls, molb.nmol, destnr, srcnr);

        natoms += molb.nmol * srcnr;
    }
}

/*! \brief Updates inter-molecular exclusion lists
 *
 * This function updates inter-molecular exclusions to exclude all
 * non-bonded interactions between a given list of atoms
 *
 * \param[inout]    excls   existing exclusions in local topology
 * \param[in]       ids     list of global IDs of atoms
 */
static void addMimicExclusions(t_blocka* excls, const gmx::ArrayRef<const int> ids)
{
    t_blocka inter_excl{};
    init_blocka(&inter_excl);
    size_t n_q = ids.size();

    inter_excl.nr  = excls->nr;
    inter_excl.nra = n_q * n_q;

    size_t total_nra = n_q * n_q;

    snew(inter_excl.index, excls->nr + 1);
    snew(inter_excl.a, total_nra);

    for (int i = 0; i < excls->nr; ++i)
    {
        inter_excl.index[i] = 0;
    }

    /* Here we loop over the list of QM atom ids
     *  and create exclusions between all of them resulting in
     *  n_q * n_q sized exclusion list
     */
    int prev_index = 0;
    for (int k = 0; k < inter_excl.nr; ++k)
    {
        inter_excl.index[k] = prev_index;
        for (long i = 0; i < ids.ssize(); i++)
        {
            if (k != ids[i])
            {
                continue;
            }
            size_t index             = n_q * i;
            inter_excl.index[ids[i]] = index;
            prev_index               = index + n_q;
            for (size_t j = 0; j < n_q; ++j)
            {
                inter_excl.a[n_q * i + j] = ids[j];
            }
        }
    }
    inter_excl.index[ids[n_q - 1] + 1] = n_q * n_q;

    inter_excl.index[inter_excl.nr] = n_q * n_q;

    std::vector<gmx::ExclusionBlock> qmexcl2(excls->nr);
    gmx::blockaToExclusionBlocks(&inter_excl, qmexcl2);

    // Merge the created exclusion list with the existing one
    gmx::mergeExclusions(excls, qmexcl2);
}

static void gen_local_top(const gmx_mtop_t& mtop,
                          bool              freeEnergyInteractionsAtEnd,
                          bool              bMergeConstr,
                          gmx_localtop_t*   top)
{
    copyAtomtypesFromMtop(mtop, &top->atomtypes);
    copyIdefFromMtop(mtop, &top->idef, freeEnergyInteractionsAtEnd, bMergeConstr);
    copyExclsFromMtop(mtop, &top->excls);
    if (!mtop.intermolecularExclusionGroup.empty())
    {
        addMimicExclusions(&top->excls, mtop.intermolecularExclusionGroup);
    }
}

void gmx_mtop_generate_local_top(const gmx_mtop_t& mtop, gmx_localtop_t* top, bool freeEnergyInteractionsAtEnd)
{
    gen_local_top(mtop, freeEnergyInteractionsAtEnd, true, top);
}

/*! \brief Fills an array with molecule begin/end atom indices
 *
 * \param[in]  mtop   The global topology
 * \param[out] index  Array of size nr. of molecules + 1 to be filled with molecule begin/end indices
 */
static void fillMoleculeIndices(const gmx_mtop_t& mtop, gmx::ArrayRef<int> index)
{
    int globalAtomIndex   = 0;
    int globalMolIndex    = 0;
    index[globalMolIndex] = globalAtomIndex;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int numAtomsPerMolecule = mtop.moltype[molb.type].atoms.nr;
        for (int mol = 0; mol < molb.nmol; mol++)
        {
            globalAtomIndex += numAtomsPerMolecule;
            globalMolIndex += 1;
            index[globalMolIndex] = globalAtomIndex;
        }
    }
}

gmx::RangePartitioning gmx_mtop_molecules(const gmx_mtop_t& mtop)
{
    gmx::RangePartitioning mols;

    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int numAtomsPerMolecule = mtop.moltype[molb.type].atoms.nr;
        for (int mol = 0; mol < molb.nmol; mol++)
        {
            mols.appendBlock(numAtomsPerMolecule);
        }
    }

    return mols;
}

/*! \brief Creates and returns a deprecated t_block struct with molecule indices
 *
 * \param[in] mtop  The global topology
 */
static t_block gmx_mtop_molecules_t_block(const gmx_mtop_t& mtop)
{
    t_block mols;

    mols.nr           = gmx_mtop_num_molecules(mtop);
    mols.nalloc_index = mols.nr + 1;
    snew(mols.index, mols.nalloc_index);

    fillMoleculeIndices(mtop, gmx::arrayRefFromArray(mols.index, mols.nr + 1));

    return mols;
}

static void gen_t_topology(const gmx_mtop_t& mtop,
                           bool              freeEnergyInteractionsAtEnd,
                           bool              bMergeConstr,
                           t_topology*       top)
{
    copyAtomtypesFromMtop(mtop, &top->atomtypes);
    copyIdefFromMtop(mtop, &top->idef, freeEnergyInteractionsAtEnd, bMergeConstr);
    copyExclsFromMtop(mtop, &top->excls);

    top->name                        = mtop.name;
    top->atoms                       = gmx_mtop_global_atoms(&mtop);
    top->mols                        = gmx_mtop_molecules_t_block(mtop);
    top->bIntermolecularInteractions = mtop.bIntermolecularInteractions;
    top->symtab                      = mtop.symtab;
}

t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t* mtop, bool freeMTop)
{
    t_topology top;

    gen_t_topology(*mtop, false, false, &top);

    if (freeMTop)
    {
        // Clear pointers and counts, such that the pointers copied to top
        // keep pointing to valid data after destroying mtop.
        mtop->symtab.symbuf = nullptr;
        mtop->symtab.nr     = 0;
    }
    return top;
}

std::vector<int> get_atom_index(const gmx_mtop_t* mtop)
{

    std::vector<int> atom_index;
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           index = atomP.globalAtomNumber();
        if (local.ptype == eptAtom)
        {
            atom_index.push_back(index);
        }
    }
    return atom_index;
}

void convertAtomsToMtop(t_symtab* symtab, char** name, t_atoms* atoms, gmx_mtop_t* mtop)
{
    mtop->symtab = *symtab;

    mtop->name = name;

    mtop->moltype.clear();
    mtop->moltype.resize(1);
    mtop->moltype.back().atoms = *atoms;

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;

    mtop->bIntermolecularInteractions = FALSE;

    mtop->natoms = atoms->nr;

    mtop->haveMoleculeIndices = false;

    gmx_mtop_finalize(mtop);
}
