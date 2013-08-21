/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static int gmx_mtop_maxresnr(const gmx_mtop_t *mtop, int maxres_renum)
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

void gmx_mtop_finalize(gmx_mtop_t *mtop)
{
    char *env;

    mtop->maxres_renum = 1;

    env = getenv("GMX_MAXRESRENUM");
    if (env != NULL)
    {
        sscanf(env, "%d", &mtop->maxres_renum);
    }
    if (mtop->maxres_renum == -1)
    {
        /* -1 signals renumber residues in all molecules */
        mtop->maxres_renum = INT_MAX;
    }

    mtop->maxresnr = gmx_mtop_maxresnr(mtop, mtop->maxres_renum);
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


typedef struct
{
    int a_start;
    int a_end;
    int na_mol;
} mb_at_t;

typedef struct gmx_mtop_atomlookup
{
    const gmx_mtop_t *mtop;
    int               nmb;
    int               mb_start;
    mb_at_t          *mba;
} t_gmx_mtop_atomlookup;


gmx_mtop_atomlookup_t
gmx_mtop_atomlookup_init(const gmx_mtop_t *mtop)
{
    t_gmx_mtop_atomlookup *alook;
    int                    mb;
    int                    a_start, a_end, na, na_start = -1;

    snew(alook, 1);

    alook->mtop     = mtop;
    alook->nmb      = mtop->nmolblock;
    alook->mb_start = 0;
    snew(alook->mba, alook->nmb);

    a_start = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        na    = mtop->molblock[mb].nmol*mtop->molblock[mb].natoms_mol;
        a_end = a_start + na;

        alook->mba[mb].a_start = a_start;
        alook->mba[mb].a_end   = a_end;
        alook->mba[mb].na_mol  = mtop->molblock[mb].natoms_mol;

        /* We start the binary search with the largest block */
        if (mb == 0 || na > na_start)
        {
            alook->mb_start = mb;
            na_start        = na;
        }

        a_start = a_end;
    }

    return alook;
}

gmx_mtop_atomlookup_t
gmx_mtop_atomlookup_settle_init(const gmx_mtop_t *mtop)
{
    t_gmx_mtop_atomlookup *alook;
    int                    mb;
    int                    na, na_start = -1;

    alook = gmx_mtop_atomlookup_init(mtop);

    /* Check if the starting molblock has settle */
    if (mtop->moltype[mtop->molblock[alook->mb_start].type].ilist[F_SETTLE].nr  == 0)
    {
        /* Search the largest molblock with settle */
        alook->mb_start = -1;
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            if (mtop->moltype[mtop->molblock[mb].type].ilist[F_SETTLE].nr > 0)
            {
                na = alook->mba[mb].a_end - alook->mba[mb].a_start;
                if (alook->mb_start == -1 || na > na_start)
                {
                    alook->mb_start = mb;
                    na_start        = na;
                }
            }
        }

        if (alook->mb_start == -1)
        {
            gmx_incons("gmx_mtop_atomlookup_settle_init called without settles");
        }
    }

    return alook;
}

void
gmx_mtop_atomlookup_destroy(gmx_mtop_atomlookup_t alook)
{
    sfree(alook->mba);
    sfree(alook);
}

void gmx_mtop_atomnr_to_atom(const gmx_mtop_atomlookup_t alook,
                             int                         atnr_global,
                             t_atom                    **atom)
{
    int mb0, mb1, mb;
    int a_start, atnr_mol;

#ifdef DEBUG_MTOP
    if (atnr_global < 0 || atnr_global >= mtop->natoms)
    {
        gmx_fatal(FARGS, "gmx_mtop_atomnr_to_moltype was called with atnr_global=%d which is not in the atom range of this system (%d-%d)",
                  atnr_global, 0, mtop->natoms-1);
    }
#endif

    mb0 = -1;
    mb1 = alook->nmb;
    mb  = alook->mb_start;

    while (TRUE)
    {
        a_start = alook->mba[mb].a_start;
        if (atnr_global < a_start)
        {
            mb1 = mb;
        }
        else if (atnr_global >= alook->mba[mb].a_end)
        {
            mb0 = mb;
        }
        else
        {
            break;
        }
        mb = ((mb0 + mb1 + 1)>>1);
    }

    atnr_mol = (atnr_global - a_start) % alook->mba[mb].na_mol;

    *atom = &alook->mtop->moltype[alook->mtop->molblock[mb].type].atoms.atom[atnr_mol];
}

void gmx_mtop_atomnr_to_ilist(const gmx_mtop_atomlookup_t alook,
                              int atnr_global,
                              t_ilist **ilist_mol, int *atnr_offset)
{
    int mb0, mb1, mb;
    int a_start, atnr_local;

#ifdef DEBUG_MTOP
    if (atnr_global < 0 || atnr_global >= mtop->natoms)
    {
        gmx_fatal(FARGS, "gmx_mtop_atomnr_to_moltype was called with atnr_global=%d which is not in the atom range of this system (%d-%d)",
                  atnr_global, 0, mtop->natoms-1);
    }
#endif

    mb0 = -1;
    mb1 = alook->nmb;
    mb  = alook->mb_start;

    while (TRUE)
    {
        a_start = alook->mba[mb].a_start;
        if (atnr_global < a_start)
        {
            mb1 = mb;
        }
        else if (atnr_global >= alook->mba[mb].a_end)
        {
            mb0 = mb;
        }
        else
        {
            break;
        }
        mb = ((mb0 + mb1 + 1)>>1);
    }

    *ilist_mol = alook->mtop->moltype[alook->mtop->molblock[mb].type].ilist;

    atnr_local = (atnr_global - a_start) % alook->mba[mb].na_mol;

    *atnr_offset = atnr_global - atnr_local;
}

void gmx_mtop_atomnr_to_molblock_ind(const gmx_mtop_atomlookup_t alook,
                                     int atnr_global,
                                     int *molb, int *molnr, int *atnr_mol)
{
    int mb0, mb1, mb;
    int a_start;

#ifdef DEBUG_MTOP
    if (atnr_global < 0 || atnr_global >= mtop->natoms)
    {
        gmx_fatal(FARGS, "gmx_mtop_atomnr_to_moltype was called with atnr_global=%d which is not in the atom range of this system (%d-%d)",
                  atnr_global, 0, mtop->natoms-1);
    }
#endif

    mb0 = -1;
    mb1 = alook->nmb;
    mb  = alook->mb_start;

    while (TRUE)
    {
        a_start = alook->mba[mb].a_start;
        if (atnr_global < a_start)
        {
            mb1 = mb;
        }
        else if (atnr_global >= alook->mba[mb].a_end)
        {
            mb0 = mb;
        }
        else
        {
            break;
        }
        mb = ((mb0 + mb1 + 1)>>1);
    }

    *molb     = mb;
    *molnr    = (atnr_global - a_start) / alook->mba[mb].na_mol;
    *atnr_mol = atnr_global - a_start - (*molnr)*alook->mba[mb].na_mol;
}

void gmx_mtop_atominfo_global(const gmx_mtop_t *mtop, int atnr_global,
                              char **atomname, int *resnr, char **resname)
{
    int             mb, a_start, a_end, maxresnr, at_loc;
    gmx_molblock_t *molb;
    t_atoms        *atoms = NULL;

    if (atnr_global < 0 || atnr_global >= mtop->natoms)
    {
        gmx_fatal(FARGS, "gmx_mtop_atominfo_global was called with atnr_global=%d which is not in the atom range of this system (%d-%d)",
                  atnr_global, 0, mtop->natoms-1);
    }

    mb       = -1;
    a_end    = 0;
    maxresnr = mtop->maxresnr;
    do
    {
        if (mb >= 0)
        {
            /* cppcheck-suppress nullPointer #6330*/
            if (atoms->nres <= mtop->maxres_renum)
            {
                /* Single residue molecule, keep counting */
                /* cppcheck-suppress nullPointer #6330*/
                maxresnr += mtop->molblock[mb].nmol*atoms->nres;
            }
        }
        mb++;
        atoms   = &mtop->moltype[mtop->molblock[mb].type].atoms;
        a_start = a_end;
        a_end   = a_start + mtop->molblock[mb].nmol*atoms->nr;
    }
    while (atnr_global >= a_end);

    at_loc    = (atnr_global - a_start) % atoms->nr;
    *atomname = *(atoms->atomname[at_loc]);
    if (atoms->nres > mtop->maxres_renum)
    {
        *resnr = atoms->resinfo[atoms->atom[at_loc].resind].nr;
    }
    else
    {
        /* Single residue molecule, keep counting */
        *resnr = maxresnr + 1 + (atnr_global - a_start)/atoms->nr*atoms->nres + atoms->atom[at_loc].resind;
    }
    *resname  = *(atoms->resinfo[atoms->atom[at_loc].resind].name);
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
                                    int *at_global, t_atom **atom)
{
    if (aloop == NULL)
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
                                      t_atom **atom, int *nmol)
{
    if (aloop == NULL)
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
    if (iloop == NULL)
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
    gmx_molblock_t *molb;

    if (iloop == NULL)
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
    t_atoms        *atoms;

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

    if (srcnr)
    {
        size = destnr+copies*srcnr;
        srenew(dest->atom, size);
        srenew(dest->atomname, size);
        srenew(dest->atomtype, size);
        srenew(dest->atomtypeB, size);
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
        memcpy((char *) &(dest->atomname[l]), (char *) &(src->atomname[0]),
               (size_t)(srcnr*sizeof(src->atomname[0])));
        memcpy((char *) &(dest->atomtype[l]), (char *) &(src->atomtype[0]),
               (size_t)(srcnr*sizeof(src->atomtype[0])));
        memcpy((char *) &(dest->atomtypeB[l]), (char *) &(src->atomtypeB[0]),
               (size_t)(srcnr*sizeof(src->atomtypeB[0])));
        memcpy((char *) &(dest->atom[l]), (char *) &(src->atom[0]),
               (size_t)(srcnr*sizeof(src->atom[0])));
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

void gmx_mtop_make_atomic_charge_groups(gmx_mtop_t *mtop,
                                        gmx_bool    bKeepSingleMolCG)
{
    int      mb, cg;
    t_block *cgs_mol;

    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        cgs_mol = &mtop->moltype[mtop->molblock[mb].type].cgs;
        if (!(bKeepSingleMolCG && cgs_mol->nr == 1))
        {
            cgs_mol->nr           = mtop->molblock[mb].natoms_mol;
            cgs_mol->nalloc_index = cgs_mol->nr + 1;
            srenew(cgs_mol->index, cgs_mol->nalloc_index);
            for (cg = 0; cg < cgs_mol->nr+1; cg++)
            {
                cgs_mol->index[cg] = cg;
            }
        }
    }
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

static void gen_local_top(const gmx_mtop_t *mtop, const t_inputrec *ir,
                          gmx_bool bMergeConstr,
                          gmx_localtop_t *top)
{
    int                     mb, srcnr, destnr, ftype, ftype_dest, mt, natoms, mol, nposre_old, nfbposre_old;
    gmx_molblock_t         *molb;
    gmx_moltype_t          *molt;
    const gmx_ffparams_t   *ffp;
    t_idef                 *idef;
    real                   *qA, *qB;
    gmx_mtop_atomloop_all_t aloop;
    int                     ag;
    t_atom                 *atom;

    top->atomtypes = mtop->atomtypes;

    ffp = &mtop->ffparams;

    idef                          = &top->idef;
    idef->ntypes                  = ffp->ntypes;
    idef->atnr                    = ffp->atnr;
    idef->functype                = ffp->functype;
    idef->iparams                 = ffp->iparams;
    idef->iparams_posres          = NULL;
    idef->iparams_posres_nalloc   = 0;
    idef->iparams_fbposres        = NULL;
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
        idef->il[ftype].iatoms = NULL;
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

    if (ir == NULL)
    {
        top->idef.ilsort = ilsortUNKNOWN;
    }
    else
    {
        if (ir->efep != efepNO && gmx_mtop_bondeds_free_energy(mtop))
        {
            snew(qA, mtop->natoms);
            snew(qB, mtop->natoms);
            aloop = gmx_mtop_atomloop_all_init(mtop);
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
}

gmx_localtop_t *gmx_mtop_generate_local_top(const gmx_mtop_t *mtop,
                                            const t_inputrec *ir)
{
    gmx_localtop_t *top;

    snew(top, 1);

    gen_local_top(mtop, ir, TRUE, top);

    return top;
}

t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t *mtop)
{
    int            mt, mb;
    gmx_localtop_t ltop;
    t_topology     top;

    gen_local_top(mtop, NULL, FALSE, &ltop);

    top.name                        = mtop->name;
    top.idef                        = ltop.idef;
    top.atomtypes                   = ltop.atomtypes;
    top.cgs                         = ltop.cgs;
    top.excls                       = ltop.excls;
    top.atoms                       = gmx_mtop_global_atoms(mtop);
    top.mols                        = mtop->mols;
    top.bIntermolecularInteractions = mtop->bIntermolecularInteractions;
    top.symtab                      = mtop->symtab;

    /* We only need to free the moltype and molblock data,
     * all other pointers have been copied to top.
     *
     * Well, except for the group data, but we can't free those, because they
     * are used somewhere even after a call to this function.
     */
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

    return top;
}
