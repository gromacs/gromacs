/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

int gmx_mtop_maxresnr(const gmx_mtop_t *mtop, int maxres_renum)
{
    int            maxresnr, mt, r;
    const t_atoms *atoms;

    maxresnr = 0;

    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        atoms = &mtop->moltype[mt].atoms;
        if (atoms->nres > maxres_renum)
        {
            for (r = 0; r < atoms->nres; r++)
            {
                if (atoms->resinfo[r].nr > maxresnr)
                {
                    maxresnr = atoms->resinfo[r].nr;
                }
            }
        }
    }

    return maxresnr;
}

static void finalizeMolblocks(gmx_mtop_t *mtop)
{
    int atomIndex          = 0;
    int residueIndex       = 0;
    int residueNumberStart = mtop->maxresnr + 1;
    for (int mb = 0; mb < mtop->nmolblock; mb++)
    {
        gmx_molblock_t *molb          = &mtop->molblock[mb];
        int             numResPerMol  = mtop->moltype[molb->type].atoms.nres;
        molb->globalAtomStart         = atomIndex;
        molb->globalResidueStart      = residueIndex;
        atomIndex                    += molb->nmol*molb->natoms_mol;
        residueIndex                 += molb->nmol*numResPerMol;
        molb->globalAtomEnd           = atomIndex;
        molb->residueNumberStart      = residueNumberStart;
        if (numResPerMol <= mtop->maxres_renum)
        {
            residueNumberStart       += molb->nmol*numResPerMol;
        }
    }
}

void gmx_mtop_finalize(gmx_mtop_t *mtop)
{
    char *env;

    if (mtop->nmolblock == 1 && mtop->molblock[0].nmol == 1)
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

    finalizeMolblocks(mtop);
}

void gmx_mtop_count_atomtypes(const gmx_mtop_t *mtop, int state, int typecount[])
{
    int      i, mb, nmol, tpi;
    t_atoms *atoms;

    for (i = 0; i < mtop->ffparams.atnr; ++i)
    {
        typecount[i] = 0;
    }
    for (mb = 0; mb < mtop->nmolblock; ++mb)
    {
        nmol  = mtop->molblock[mb].nmol;
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        for (i = 0; i < atoms->nr; ++i)
        {
            if (state == 0)
            {
                tpi = atoms->atom[i].type;
            }
            else
            {
                tpi = atoms->atom[i].typeB;
            }
            typecount[tpi] += nmol;
        }
    }
}

int ncg_mtop(const gmx_mtop_t *mtop)
{
    int ncg;
    int mb;

    ncg = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        ncg +=
            mtop->molblock[mb].nmol*
            mtop->moltype[mtop->molblock[mb].type].cgs.nr;
    }

    return ncg;
}

int gmx_mtop_nres(const gmx_mtop_t *mtop)
{
    int nres = 0;
    for (int mb = 0; mb < mtop->nmolblock; ++mb)
    {
        nres +=
            mtop->molblock[mb].nmol*
            mtop->moltype[mtop->molblock[mb].type].atoms.nres;
    }
    return nres;
}

void gmx_mtop_remove_chargegroups(gmx_mtop_t *mtop)
{
    int      mt;
    t_block *cgs;
    int      i;

    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        cgs = &mtop->moltype[mt].cgs;
        if (cgs->nr < mtop->moltype[mt].atoms.nr)
        {
            cgs->nr = mtop->moltype[mt].atoms.nr;
            srenew(cgs->index, cgs->nr+1);
            for (i = 0; i < cgs->nr+1; i++)
            {
                cgs->index[i] = i;
            }
        }
    }
}

typedef struct gmx_mtop_atomloop_all
{
    const gmx_mtop_t *mtop;
    int               mblock;
    t_atoms          *atoms;
    int               mol;
    int               maxresnr;
    int               at_local;
    int               at_global;
} t_gmx_mtop_atomloop_all;

gmx_mtop_atomloop_all_t
gmx_mtop_atomloop_all_init(const gmx_mtop_t *mtop)
{
    struct gmx_mtop_atomloop_all *aloop;

    snew(aloop, 1);

    aloop->mtop         = mtop;
    aloop->mblock       = 0;
    aloop->atoms        =
        &mtop->moltype[mtop->molblock[aloop->mblock].type].atoms;
    aloop->mol          = 0;
    aloop->maxresnr     = mtop->maxresnr;
    aloop->at_local     = -1;
    aloop->at_global    = -1;

    return aloop;
}

static void gmx_mtop_atomloop_all_destroy(gmx_mtop_atomloop_all_t aloop)
{
    sfree(aloop);
}

gmx_bool gmx_mtop_atomloop_all_next(gmx_mtop_atomloop_all_t aloop,
                                    int *at_global, const t_atom **atom)
{
    if (aloop == nullptr)
    {
        gmx_incons("gmx_mtop_atomloop_all_next called without calling gmx_mtop_atomloop_all_init");
    }

    aloop->at_local++;
    aloop->at_global++;

    if (aloop->at_local >= aloop->atoms->nr)
    {
        if (aloop->atoms->nres <= aloop->mtop->maxres_renum)
        {
            /* Single residue molecule, increase the count with one */
            aloop->maxresnr += aloop->atoms->nres;
        }
        aloop->mol++;
        aloop->at_local = 0;
        if (aloop->mol >= aloop->mtop->molblock[aloop->mblock].nmol)
        {
            aloop->mblock++;
            if (aloop->mblock >= aloop->mtop->nmolblock)
            {
                gmx_mtop_atomloop_all_destroy(aloop);
                return FALSE;
            }
            aloop->atoms = &aloop->mtop->moltype[aloop->mtop->molblock[aloop->mblock].type].atoms;
            aloop->mol   = 0;
        }
    }

    *at_global = aloop->at_global;
    *atom      = &aloop->atoms->atom[aloop->at_local];

    return TRUE;
}

void gmx_mtop_atomloop_all_names(gmx_mtop_atomloop_all_t aloop,
                                 char **atomname, int *resnr, char **resname)
{
    int resind_mol;

    *atomname  = *(aloop->atoms->atomname[aloop->at_local]);
    resind_mol = aloop->atoms->atom[aloop->at_local].resind;
    *resnr     = aloop->atoms->resinfo[resind_mol].nr;
    if (aloop->atoms->nres <= aloop->mtop->maxres_renum)
    {
        *resnr = aloop->maxresnr + 1 + resind_mol;
    }
    *resname  = *(aloop->atoms->resinfo[resind_mol].name);
}

void gmx_mtop_atomloop_all_moltype(gmx_mtop_atomloop_all_t aloop,
                                   gmx_moltype_t **moltype, int *at_mol)
{
    *moltype = &aloop->mtop->moltype[aloop->mtop->molblock[aloop->mblock].type];
    *at_mol  = aloop->at_local;
}

typedef struct gmx_mtop_atomloop_block
{
    const gmx_mtop_t *mtop;
    int               mblock;
    t_atoms          *atoms;
    int               at_local;
} t_gmx_mtop_atomloop_block;

gmx_mtop_atomloop_block_t
gmx_mtop_atomloop_block_init(const gmx_mtop_t *mtop)
{
    struct gmx_mtop_atomloop_block *aloop;

    snew(aloop, 1);

    aloop->mtop      = mtop;
    aloop->mblock    = 0;
    aloop->atoms     = &mtop->moltype[mtop->molblock[aloop->mblock].type].atoms;
    aloop->at_local  = -1;

    return aloop;
}

static void gmx_mtop_atomloop_block_destroy(gmx_mtop_atomloop_block_t aloop)
{
    sfree(aloop);
}

gmx_bool gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop,
                                      const t_atom **atom, int *nmol)
{
    if (aloop == nullptr)
    {
        gmx_incons("gmx_mtop_atomloop_all_next called without calling gmx_mtop_atomloop_all_init");
    }

    aloop->at_local++;

    if (aloop->at_local >= aloop->atoms->nr)
    {
        aloop->mblock++;
        if (aloop->mblock >= aloop->mtop->nmolblock)
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
    const gmx_mtop_t *mtop;
    int               mblock;
} t_gmx_mtop_ilist;

gmx_mtop_ilistloop_t
gmx_mtop_ilistloop_init(const gmx_mtop_t *mtop)
{
    struct gmx_mtop_ilistloop *iloop;

    snew(iloop, 1);

    iloop->mtop      = mtop;
    iloop->mblock    = -1;

    return iloop;
}

static void gmx_mtop_ilistloop_destroy(gmx_mtop_ilistloop_t iloop)
{
    sfree(iloop);
}

gmx_bool gmx_mtop_ilistloop_next(gmx_mtop_ilistloop_t iloop,
                                 t_ilist **ilist_mol, int *nmol)
{
    if (iloop == nullptr)
    {
        gmx_incons("gmx_mtop_ilistloop_next called without calling gmx_mtop_ilistloop_init");
    }

    iloop->mblock++;
    if (iloop->mblock >= iloop->mtop->nmolblock)
    {
        if (iloop->mblock == iloop->mtop->nmolblock &&
            iloop->mtop->bIntermolecularInteractions)
        {
            *ilist_mol = iloop->mtop->intermolecular_ilist;
            *nmol      = 1;
            return TRUE;
        }

        gmx_mtop_ilistloop_destroy(iloop);
        return FALSE;
    }

    *ilist_mol =
        iloop->mtop->moltype[iloop->mtop->molblock[iloop->mblock].type].ilist;

    *nmol = iloop->mtop->molblock[iloop->mblock].nmol;

    return TRUE;
}
typedef struct gmx_mtop_ilistloop_all
{
    const gmx_mtop_t *mtop;
    int               mblock;
    int               mol;
    int               a_offset;
} t_gmx_mtop_ilist_all;

gmx_mtop_ilistloop_all_t
gmx_mtop_ilistloop_all_init(const gmx_mtop_t *mtop)
{
    struct gmx_mtop_ilistloop_all *iloop;

    snew(iloop, 1);

    iloop->mtop      = mtop;
    iloop->mblock    = 0;
    iloop->mol       = -1;
    iloop->a_offset  = 0;

    return iloop;
}

static void gmx_mtop_ilistloop_all_destroy(gmx_mtop_ilistloop_all_t iloop)
{
    sfree(iloop);
}

gmx_bool gmx_mtop_ilistloop_all_next(gmx_mtop_ilistloop_all_t iloop,
                                     t_ilist **ilist_mol, int *atnr_offset)
{

    if (iloop == nullptr)
    {
        gmx_incons("gmx_mtop_ilistloop_all_next called without calling gmx_mtop_ilistloop_all_init");
    }

    if (iloop->mol >= 0)
    {
        iloop->a_offset += iloop->mtop->molblock[iloop->mblock].natoms_mol;
    }

    iloop->mol++;

    /* Inter-molecular interactions, if present, are indexed with
     * iloop->mblock == iloop->mtop->nmolblock, thus we should separately
     * check for this value in this conditional.
     */
    if (iloop->mblock == iloop->mtop->nmolblock ||
        iloop->mol >= iloop->mtop->molblock[iloop->mblock].nmol)
    {
        iloop->mblock++;
        iloop->mol = 0;
        if (iloop->mblock >= iloop->mtop->nmolblock)
        {
            if (iloop->mblock == iloop->mtop->nmolblock &&
                iloop->mtop->bIntermolecularInteractions)
            {
                *ilist_mol   = iloop->mtop->intermolecular_ilist;
                *atnr_offset = 0;
                return TRUE;
            }

            gmx_mtop_ilistloop_all_destroy(iloop);
            return FALSE;
        }
    }

    *ilist_mol =
        iloop->mtop->moltype[iloop->mtop->molblock[iloop->mblock].type].ilist;

    *atnr_offset = iloop->a_offset;

    return TRUE;
}

int gmx_mtop_ftype_count(const gmx_mtop_t *mtop, int ftype)
{
    gmx_mtop_ilistloop_t iloop;
    t_ilist             *il;
    int                  n, nmol;

    n = 0;

    iloop = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        n += nmol*il[ftype].nr/(1+NRAL(ftype));
    }

    if (mtop->bIntermolecularInteractions)
    {
        n += mtop->intermolecular_ilist[ftype].nr/(1+NRAL(ftype));
    }

    return n;
}

t_block gmx_mtop_global_cgs(const gmx_mtop_t *mtop)
{
    t_block         cgs_gl, *cgs_mol;
    int             mb, mol, cg;
    gmx_molblock_t *molb;

    /* In most cases this is too much, but we realloc at the end */
    snew(cgs_gl.index, mtop->natoms+1);

    cgs_gl.nr       = 0;
    cgs_gl.index[0] = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb    = &mtop->molblock[mb];
        cgs_mol = &mtop->moltype[molb->type].cgs;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            for (cg = 0; cg < cgs_mol->nr; cg++)
            {
                cgs_gl.index[cgs_gl.nr+1] =
                    cgs_gl.index[cgs_gl.nr] +
                    cgs_mol->index[cg+1] - cgs_mol->index[cg];
                cgs_gl.nr++;
            }
        }
    }
    cgs_gl.nalloc_index = cgs_gl.nr + 1;
    srenew(cgs_gl.index, cgs_gl.nalloc_index);

    return cgs_gl;
}

static void atomcat(t_atoms *dest, t_atoms *src, int copies,
                    int maxres_renum, int *maxresnr)
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
        dest->haveMass    = dest->haveMass    && src->haveMass;
        dest->haveType    = dest->haveType    && src->haveType;
        dest->haveCharge  = dest->haveCharge  && src->haveCharge;
        dest->haveBState  = dest->haveBState  && src->haveBState;
        dest->havePdbInfo = dest->havePdbInfo && src->havePdbInfo;
    }

    if (srcnr)
    {
        size = destnr+copies*srcnr;
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
        size = dest->nres+copies*src->nres;
        srenew(dest->resinfo, size);
    }

    /* residue information */
    for (l = dest->nres, j = 0; (j < copies); j++, l += src->nres)
    {
        memcpy((char *) &(dest->resinfo[l]), (char *) &(src->resinfo[0]),
               (size_t)(src->nres*sizeof(src->resinfo[0])));
    }

    for (l = destnr, j = 0; (j < copies); j++, l += srcnr)
    {
        memcpy((char *) &(dest->atom[l]), (char *) &(src->atom[0]),
               (size_t)(srcnr*sizeof(src->atom[0])));
        memcpy((char *) &(dest->atomname[l]), (char *) &(src->atomname[0]),
               (size_t)(srcnr*sizeof(src->atomname[0])));
        if (dest->haveType)
        {
            memcpy((char *) &(dest->atomtype[l]), (char *) &(src->atomtype[0]),
                   (size_t)(srcnr*sizeof(src->atomtype[0])));
            if (dest->haveBState)
            {
                memcpy((char *) &(dest->atomtypeB[l]), (char *) &(src->atomtypeB[0]),
                       (size_t)(srcnr*sizeof(src->atomtypeB[0])));
            }
        }
        if (dest->havePdbInfo)
        {
            memcpy((char *) &(dest->pdbinfo[l]), (char *) &(src->pdbinfo[0]),
                   (size_t)(srcnr*sizeof(src->pdbinfo[0])));
        }
    }

    /* Increment residue indices */
    for (l = destnr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < srcnr); i++, l++)
        {
            dest->atom[l].resind = dest->nres+j*src->nres+src->atom[i].resind;
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
                dest->resinfo[dest->nres+j*src->nres+l].nr = *maxresnr;
            }
        }
    }

    dest->nres += copies*src->nres;
    dest->nr   += copies*src->nr;
}

t_atoms gmx_mtop_global_atoms(const gmx_mtop_t *mtop)
{
    t_atoms         atoms;
    int             maxresnr, mb;
    gmx_molblock_t *molb;

    init_t_atoms(&atoms, 0, FALSE);

    maxresnr = mtop->maxresnr;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        atomcat(&atoms, &mtop->moltype[molb->type].atoms, molb->nmol,
                mtop->maxres_renum, &maxresnr);
    }

    return atoms;
}

/*
 * The cat routines below are old code from src/kernel/topcat.c
 */

static void blockcat(t_block *dest, t_block *src, int copies)
{
    int i, j, l, nra, size;

    if (src->nr)
    {
        size = (dest->nr+copies*src->nr+1);
        srenew(dest->index, size);
    }

    nra = dest->index[dest->nr];
    for (l = dest->nr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < src->nr); i++)
        {
            dest->index[l++] = nra + src->index[i];
        }
        nra += src->index[src->nr];
    }
    dest->nr             += copies*src->nr;
    dest->index[dest->nr] = nra;
}

static void blockacat(t_blocka *dest, t_blocka *src, int copies,
                      int dnum, int snum)
{
    int i, j, l, size;
    int destnr  = dest->nr;
    int destnra = dest->nra;

    if (src->nr)
    {
        size = (dest->nr+copies*src->nr+1);
        srenew(dest->index, size);
    }
    if (src->nra)
    {
        size = (dest->nra+copies*src->nra);
        srenew(dest->a, size);
    }

    for (l = destnr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < src->nr); i++)
        {
            dest->index[l++] = dest->nra+src->index[i];
        }
        dest->nra += src->nra;
    }
    for (l = destnra, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < src->nra); i++)
        {
            dest->a[l++] = dnum+src->a[i];
        }
        dnum     += snum;
        dest->nr += src->nr;
    }
    dest->index[dest->nr] = dest->nra;
}

static void ilistcat(int ftype, t_ilist *dest, t_ilist *src, int copies,
                     int dnum, int snum)
{
    int nral, c, i, a;

    nral = NRAL(ftype);

    dest->nalloc = dest->nr + copies*src->nr;
    srenew(dest->iatoms, dest->nalloc);

    for (c = 0; c < copies; c++)
    {
        for (i = 0; i < src->nr; )
        {
            dest->iatoms[dest->nr++] = src->iatoms[i++];
            for (a = 0; a < nral; a++)
            {
                dest->iatoms[dest->nr++] = dnum + src->iatoms[i++];
            }
        }
        dnum += snum;
    }
}

#ifdef BUILD_WITH_FDA
static void pf_ilistcat(int ftype, t_ilist *dest, t_ilist *src, int copies,
                        int dnum, int snum, fda::FDASettings const& fda_settings)
{
    // Return if no bonded interaction is needed.
    if (!(fda_settings.type & (fda::InteractionType_BONDED + fda::InteractionType_NB14))) return;

    int nral, c, i, a, atomIdx;
    char needed;

    nral = NRAL(ftype);

    t_iatom *tmp;
    snew(tmp,copies*src->nr);
    int len = 0;

    int *g1atomsBeg = fda_settings.groups->a + fda_settings.groups->index[fda_settings.index_group1];
    int *g1atomsEnd = fda_settings.groups->a + fda_settings.groups->index[fda_settings.index_group1 + 1];
    int *g1atomsCur = NULL;
    int *g2atomsBeg = fda_settings.groups->a + fda_settings.groups->index[fda_settings.index_group2];
    int *g2atomsEnd = fda_settings.groups->a + fda_settings.groups->index[fda_settings.index_group2 + 1];
    int *g2atomsCur = NULL;

    for (c = 0; c < copies; c++)
    {
        for (i = 0; i < src->nr; )
        {
            needed = 0;
            for (a = 0; a < nral; a++)
            {
            	atomIdx = dnum + src->iatoms[i+a+1];
        		for (g1atomsCur = g1atomsBeg; g1atomsCur < g1atomsEnd; ++g1atomsCur) {
        			if (atomIdx == *g1atomsCur) needed = 1;
        		}
        		for (g2atomsCur = g2atomsBeg; g2atomsCur < g2atomsEnd; ++g2atomsCur) {
        			if (atomIdx == *g2atomsCur) needed = 1;
        		}
            }
            if (needed) {
                tmp[len++] = src->iatoms[i];
                for (a = 0; a < nral; a++) tmp[len++] = dnum + src->iatoms[i+a+1];

                #ifdef FDA_BONDEXCL_PRINT_DEBUG_ON
					fprintf(stderr, "=== DEBUG === bonded interaction %i", ftype);
					fprintf(stderr, " %i", src->iatoms[i]);
					for (a = 0; a < nral; a++) fprintf(stderr, " %i", dnum + src->iatoms[i + a + 1]);
					fprintf(stderr, " needed\n");
					fflush(stderr);
                #endif
            } else {
                #ifdef FDA_BONDEXCL_PRINT_DEBUG_ON
					fprintf(stderr, "=== DEBUG === bonded interaction %i", ftype);
					fprintf(stderr, " %i", src->iatoms[i]);
					for (a = 0; a < nral; a++) fprintf(stderr, " %i", dnum + src->iatoms[i + a + 1]);
					fprintf(stderr, " not needed\n");
					fflush(stderr);
                #endif
            }
            i += a + 1;
        }
        dnum += snum;
    }

    dest->nalloc = dest->nr + len;
    srenew(dest->iatoms, dest->nalloc);

    for (i = 0; i < len; i++) dest->iatoms[dest->nr++] = tmp[i];

    sfree(tmp);
}
#endif

static void set_posres_params(t_idef *idef, gmx_molblock_t *molb,
                              int i0, int a_offset)
{
    t_ilist   *il;
    int        i1, i, a_molb;
    t_iparams *ip;

    il = &idef->il[F_POSRES];
    i1 = il->nr/2;
    idef->iparams_posres_nalloc = i1;
    srenew(idef->iparams_posres, idef->iparams_posres_nalloc);
    for (i = i0; i < i1; i++)
    {
        ip = &idef->iparams_posres[i];
        /* Copy the force constants */
        *ip    = idef->iparams[il->iatoms[i*2]];
        a_molb = il->iatoms[i*2+1] - a_offset;
        if (molb->nposres_xA == 0)
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        ip->posres.pos0A[XX] = molb->posres_xA[a_molb][XX];
        ip->posres.pos0A[YY] = molb->posres_xA[a_molb][YY];
        ip->posres.pos0A[ZZ] = molb->posres_xA[a_molb][ZZ];
        if (molb->nposres_xB > 0)
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
        il->iatoms[i*2] = i;
    }
}

static void set_fbposres_params(t_idef *idef, gmx_molblock_t *molb,
                                int i0, int a_offset)
{
    t_ilist   *il;
    int        i1, i, a_molb;
    t_iparams *ip;

    il = &idef->il[F_FBPOSRES];
    i1 = il->nr/2;
    idef->iparams_fbposres_nalloc = i1;
    srenew(idef->iparams_fbposres, idef->iparams_fbposres_nalloc);
    for (i = i0; i < i1; i++)
    {
        ip = &idef->iparams_fbposres[i];
        /* Copy the force constants */
        *ip    = idef->iparams[il->iatoms[i*2]];
        a_molb = il->iatoms[i*2+1] - a_offset;
        if (molb->nposres_xA == 0)
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        /* Take flat-bottom posres reference from normal position restraints */
        ip->fbposres.pos0[XX] = molb->posres_xA[a_molb][XX];
        ip->fbposres.pos0[YY] = molb->posres_xA[a_molb][YY];
        ip->fbposres.pos0[ZZ] = molb->posres_xA[a_molb][ZZ];
        /* Note: no B-type for flat-bottom posres */

        /* Set the parameter index for idef->iparams_posre */
        il->iatoms[i*2] = i;
    }
}

static void gen_local_top(const gmx_mtop_t *mtop,
                          bool              freeEnergyInteractionsAtEnd,
                          bool              bMergeConstr,
                          gmx_localtop_t   *top
#ifdef BUILD_WITH_FDA
                          , fda::FDASettings *ptr_fda_settings
#endif
                         )
{
    int                     mb, srcnr, destnr, ftype, natoms, mol, nposre_old, nfbposre_old;
    gmx_molblock_t         *molb;
    gmx_moltype_t          *molt;
    const gmx_ffparams_t   *ffp;
    t_idef                 *idef;
    real                   *qA, *qB;
    gmx_mtop_atomloop_all_t aloop;
    int                     ag;

    top->atomtypes = mtop->atomtypes;

    ffp = &mtop->ffparams;

    idef                          = &top->idef;
    idef->ntypes                  = ffp->ntypes;
    idef->atnr                    = ffp->atnr;
    idef->functype                = ffp->functype;
    idef->iparams                 = ffp->iparams;
    idef->iparams_posres          = nullptr;
    idef->iparams_posres_nalloc   = 0;
    idef->iparams_fbposres        = nullptr;
    idef->iparams_fbposres_nalloc = 0;
    idef->fudgeQQ                 = ffp->fudgeQQ;
    idef->cmap_grid               = ffp->cmap_grid;
    idef->ilsort                  = ilsortUNKNOWN;

    init_block(&top->cgs);
    init_blocka(&top->excls);
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        idef->il[ftype].nr     = 0;
        idef->il[ftype].nalloc = 0;
        idef->il[ftype].iatoms = nullptr;
    }

    natoms = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];

        srcnr  = molt->atoms.nr;
        destnr = natoms;

        blockcat(&top->cgs, &molt->cgs, molb->nmol);

        blockacat(&top->excls, &molt->excls, molb->nmol, destnr, srcnr);

        nposre_old   = idef->il[F_POSRES].nr;
        nfbposre_old = idef->il[F_FBPOSRES].nr;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (bMergeConstr &&
                ftype == F_CONSTR && molt->ilist[F_CONSTRNC].nr > 0)
            {
                /* Merge all constrains into one ilist.
                 * This simplifies the constraint code.
                 */
                for (mol = 0; mol < molb->nmol; mol++)
                {
                    ilistcat(ftype, &idef->il[F_CONSTR], &molt->ilist[F_CONSTR],
                             1, destnr+mol*srcnr, srcnr);
                    ilistcat(ftype, &idef->il[F_CONSTR], &molt->ilist[F_CONSTRNC],
                             1, destnr+mol*srcnr, srcnr);
                }
            }
            else if (!(bMergeConstr && ftype == F_CONSTRNC))
            {
#ifdef BUILD_WITH_FDA
            	if (ptr_fda_settings and ptr_fda_settings->bonded_exclusion_on)
                    pf_ilistcat(ftype, &idef->il[ftype], &molt->ilist[ftype],
                             molb->nmol, destnr, srcnr, *ptr_fda_settings);
            	else
#endif
                    ilistcat(ftype, &idef->il[ftype], &molt->ilist[ftype],
                             molb->nmol, destnr, srcnr);
            }
        }
        if (idef->il[F_POSRES].nr > nposre_old)
        {
            /* Executing this line line stops gmxdump -sys working
             * correctly. I'm not aware there's an elegant fix. */
            set_posres_params(idef, molb, nposre_old/2, natoms);
        }
        if (idef->il[F_FBPOSRES].nr > nfbposre_old)
        {
            set_fbposres_params(idef, molb, nfbposre_old/2, natoms);
        }

        natoms += molb->nmol*srcnr;
    }

    if (mtop->bIntermolecularInteractions)
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            ilistcat(ftype, &idef->il[ftype], &mtop->intermolecular_ilist[ftype],
                     1, 0, mtop->natoms);
        }
    }

    if (freeEnergyInteractionsAtEnd && gmx_mtop_bondeds_free_energy(mtop))
    {
        snew(qA, mtop->natoms);
        snew(qB, mtop->natoms);
        aloop = gmx_mtop_atomloop_all_init(mtop);
        const t_atom *atom;
        while (gmx_mtop_atomloop_all_next(aloop, &ag, &atom))
        {
            qA[ag] = atom->q;
            qB[ag] = atom->qB;
        }
        gmx_sort_ilist_fe(&top->idef, qA, qB);
        sfree(qA);
        sfree(qB);
    }
    else
    {
        top->idef.ilsort = ilsortNO_FE;
    }
}

gmx_localtop_t *
gmx_mtop_generate_local_top(const gmx_mtop_t *mtop,
                            bool              freeEnergyInteractionsAtEnd
#ifdef BUILD_WITH_FDA
                            , fda::FDASettings *ptr_fda_settings
#endif
                           )
{
    gmx_localtop_t *top;

    snew(top, 1);

    gen_local_top(mtop, freeEnergyInteractionsAtEnd, true, top
#ifdef BUILD_WITH_FDA
                  , ptr_fda_settings
#endif
                 );

    return top;
}

t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t *mtop, bool freeMTop)
{
    int            mt, mb;
    gmx_localtop_t ltop;
    t_topology     top;

    gen_local_top(mtop, false, FALSE, &ltop

#ifdef BUILD_WITH_FDA
                  , nullptr
#endif
                 );
    ltop.idef.ilsort = ilsortUNKNOWN;

    top.name                        = mtop->name;
    top.idef                        = ltop.idef;
    top.atomtypes                   = ltop.atomtypes;
    top.cgs                         = ltop.cgs;
    top.excls                       = ltop.excls;
    top.atoms                       = gmx_mtop_global_atoms(mtop);
    top.mols                        = mtop->mols;
    top.bIntermolecularInteractions = mtop->bIntermolecularInteractions;
    top.symtab                      = mtop->symtab;

    if (freeMTop)
    {
        // Free pointers that have not been copied to top.
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            done_moltype(&mtop->moltype[mt]);
        }
        sfree(mtop->moltype);

        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            done_molblock(&mtop->molblock[mb]);
        }
        sfree(mtop->molblock);

        done_gmx_groups_t(&mtop->groups);
    }

    return top;
}

std::vector<size_t> get_atom_index(const gmx_mtop_t *mtop)
{

    std::vector<size_t>       atom_index;
    gmx_mtop_atomloop_block_t aloopb = gmx_mtop_atomloop_block_init(mtop);
    const t_atom             *atom;
    int                       nmol, j = 0;
    while (gmx_mtop_atomloop_block_next(aloopb, &atom, &nmol))
    {
        if (atom->ptype == eptAtom)
        {
            atom_index.push_back(j);
        }
        j++;
    }
    return atom_index;
}

void convertAtomsToMtop(t_symtab    *symtab,
                        char       **name,
                        t_atoms     *atoms,
                        gmx_mtop_t  *mtop)
{
    mtop->symtab                 = *symtab;

    mtop->name                   = name;

    mtop->nmoltype               = 1;
    // This snew clears all entries, we should replace it by an initializer
    snew(mtop->moltype, mtop->nmoltype);
    mtop->moltype[0].atoms       = *atoms;
    init_block(&mtop->moltype[0].cgs);
    init_blocka(&mtop->moltype[0].excls);

    mtop->nmolblock              = 1;
    // This snew clears all entries, we should replace it by an initializer
    snew(mtop->molblock, mtop->nmolblock);
    mtop->molblock[0].type       = 0;
    mtop->molblock[0].nmol       = 1;
    mtop->molblock[0].natoms_mol = atoms->nr;

    mtop->bIntermolecularInteractions = FALSE;

    mtop->natoms                 = atoms->nr;

    gmx_mtop_finalize(mtop);
}
