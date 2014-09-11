/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#include "gmxpre.h"

#include "genhydro.h"

#include <string.h>
#include <time.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/calch.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/ter_db.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static void copy_atom(t_atoms *atoms1, int a1, t_atoms *atoms2, int a2)
{
    atoms2->atom[a2] = atoms1->atom[a1];
    snew(atoms2->atomname[a2], 1);
    *atoms2->atomname[a2] = gmx_strdup(*atoms1->atomname[a1]);
}

static atom_id pdbasearch_atom(const char *name, int resind, t_atoms *pdba,
                               const char *searchtype, gmx_bool bAllowMissing)
{
    int  i;

    for (i = 0; (i < pdba->nr) && (pdba->atom[i].resind != resind); i++)
    {
        ;
    }

    return search_atom(name, i, pdba,
                       searchtype, bAllowMissing);
}

static void hacksearch_atom(int *ii, int *jj, char *name,
                            int nab[], t_hack *ab[],
                            int resind, t_atoms *pdba)
{
    int  i, j;

    *ii = -1;
    if (name[0] == '-')
    {
        name++;
        resind--;
    }
    for (i = 0; (i < pdba->nr) && (pdba->atom[i].resind != resind); i++)
    {
        ;
    }
    for (; (i < pdba->nr) && (pdba->atom[i].resind == resind) && (*ii < 0); i++)
    {
        for (j = 0; (j < nab[i]) && (*ii < 0); j++)
        {
            if (ab[i][j].nname && strcmp(name, ab[i][j].nname) == 0)
            {
                *ii = i;
                *jj = j;
            }
        }
    }

    return;
}

void dump_ab(FILE *out, int natom, int nab[], t_hack *ab[], gmx_bool bHeader)
{
    int i, j;

#define SS(s) (s) ? (s) : "-"
    /* dump ab */
    if (bHeader)
    {
        fprintf(out, "ADDBLOCK (t_hack) natom=%d\n"
                "%4s %2s %-4s %-4s %2s %-4s %-4s %-4s %-4s %1s %s\n",
                natom, "atom", "nr", "old", "new", "tp", "ai", "aj", "ak", "al", "a", "x");
    }
    for (i = 0; i < natom; i++)
    {
        for (j = 0; j < nab[i]; j++)
        {
            fprintf(out, "%4d %2d %-4s %-4s %2d %-4s %-4s %-4s %-4s %s %g %g %g\n",
                    i+1, ab[i][j].nr, SS(ab[i][j].oname), SS(ab[i][j].nname),
                    ab[i][j].tp,
                    SS(ab[i][j].AI), SS(ab[i][j].AJ),
                    SS(ab[i][j].AK), SS(ab[i][j].AL),
                    ab[i][j].atom ? "+" : "",
                    ab[i][j].newx[XX], ab[i][j].newx[YY], ab[i][j].newx[ZZ]);
        }
    }
#undef SS
}

static t_hackblock *get_hackblocks(t_atoms *pdba, int nah, t_hackblock ah[],
                                   int nterpairs,
                                   t_hackblock **ntdb, t_hackblock **ctdb,
                                   int *rN, int *rC)
{
    int          i, rnr;
    t_hackblock *hb, *ahptr;

    /* make space */
    snew(hb, pdba->nres);
    /* first the termini */
    for (i = 0; i < nterpairs; i++)
    {
        if (ntdb[i] != NULL)
        {
            copy_t_hackblock(ntdb[i], &hb[rN[i]]);
        }
        if (ctdb[i] != NULL)
        {
            merge_t_hackblock(ctdb[i], &hb[rC[i]]);
        }
    }
    /* then the whole hdb */
    for (rnr = 0; rnr < pdba->nres; rnr++)
    {
        ahptr = search_h_db(nah, ah, *pdba->resinfo[rnr].rtp);
        if (ahptr)
        {
            if (hb[rnr].name == NULL)
            {
                hb[rnr].name = gmx_strdup(ahptr->name);
            }
            merge_hacks(ahptr, &hb[rnr]);
        }
    }
    return hb;
}

static void expand_hackblocks_one(t_hackblock *hbr, char *atomname,
                                  int *nabi, t_hack **abi, gmx_bool bN, gmx_bool bC)
{
    int      j, k, l, d;
    gmx_bool bIgnore;

    /* we'll recursively add atoms to atoms */
    for (j = 0; j < hbr->nhack; j++)
    {
        /* first check if we're in the N- or C-terminus, then we should ignore
           all hacks involving atoms from resp. previous or next residue
           (i.e. which name begins with '-' (N) or '+' (C) */
        bIgnore = FALSE;
        if (bN) /* N-terminus: ignore '-' */
        {
            for (k = 0; k < 4 && hbr->hack[j].a[k] && !bIgnore; k++)
            {
                bIgnore = hbr->hack[j].a[k][0] == '-';
            }
        }
        if (bC) /* C-terminus: ignore '+' */
        {
            for (k = 0; k < 4 && hbr->hack[j].a[k] && !bIgnore; k++)
            {
                bIgnore = hbr->hack[j].a[k][0] == '+';
            }
        }
        /* must be either hdb entry (tp>0) or add from tdb (oname==NULL)
           and first control aton (AI) matches this atom or
           delete/replace from tdb (oname!=NULL) and oname matches this atom */
        if (debug)
        {
            fprintf(debug, " %s", hbr->hack[j].oname ? hbr->hack[j].oname : hbr->hack[j].AI);
        }

        if (!bIgnore &&
            ( ( ( hbr->hack[j].tp > 0 || hbr->hack[j].oname == NULL ) &&
                strcmp(atomname, hbr->hack[j].AI) == 0 ) ||
              ( hbr->hack[j].oname != NULL &&
                strcmp(atomname, hbr->hack[j].oname) == 0) ) )
        {
            /* now expand all hacks for this atom */
            if (debug)
            {
                fprintf(debug, " +%dh", hbr->hack[j].nr);
            }
            srenew(*abi, *nabi + hbr->hack[j].nr);
            for (k = 0; k < hbr->hack[j].nr; k++)
            {
                copy_t_hack(&hbr->hack[j], &(*abi)[*nabi + k]);
                (*abi)[*nabi + k].bXSet = FALSE;
                /* if we're adding (oname==NULL) and don't have a new name (nname)
                   yet, build it from atomname */
                if ( (*abi)[*nabi + k].nname == NULL)
                {
                    if ( (*abi)[*nabi + k].oname == NULL)
                    {
                        (*abi)[*nabi + k].nname    = gmx_strdup(atomname);
                        (*abi)[*nabi + k].nname[0] = 'H';
                    }
                }
                else
                {
                    if (gmx_debug_at)
                    {
                        fprintf(debug, "Hack '%s' %d, replacing nname '%s' with '%s' (old name '%s')\n",
                                atomname, j,
                                (*abi)[*nabi + k].nname, hbr->hack[j].nname,
                                (*abi)[*nabi + k].oname ? (*abi)[*nabi + k].oname : "");
                    }
                    sfree((*abi)[*nabi + k].nname);
                    (*abi)[*nabi + k].nname = gmx_strdup(hbr->hack[j].nname);
                }

                if (hbr->hack[j].tp == 10 && k == 2)
                {
                    /* This is a water virtual site, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    (*abi)[*nabi + k].nname[0] = 'M';
                }
                else if (hbr->hack[j].tp == 11 && k >= 2)
                {
                    /* This is a water lone pair, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    srenew((*abi)[*nabi + k].nname, 4);
                    (*abi)[*nabi + k].nname[0] = 'L';
                    (*abi)[*nabi + k].nname[1] = 'P';
                    (*abi)[*nabi + k].nname[2] = '1' + k - 2;
                    (*abi)[*nabi + k].nname[3] = '\0';
                }
                else if (hbr->hack[j].nr > 1)
                {
                    /* adding more than one atom, number them */
                    l = strlen((*abi)[*nabi + k].nname);
                    srenew((*abi)[*nabi + k].nname, l+2);
                    (*abi)[*nabi + k].nname[l]   = '1' + k;
                    (*abi)[*nabi + k].nname[l+1] = '\0';
                }
            }
            (*nabi) += hbr->hack[j].nr;

            /* add hacks to atoms we've just added */
            if (hbr->hack[j].tp > 0 || hbr->hack[j].oname == NULL)
            {
                for (k = 0; k < hbr->hack[j].nr; k++)
                {
                    expand_hackblocks_one(hbr, (*abi)[*nabi-hbr->hack[j].nr+k].nname,
                                          nabi, abi, bN, bC);
                }
            }
        }
    }
}

static void expand_hackblocks(t_atoms *pdba, t_hackblock hb[],
                              int nab[], t_hack *ab[],
                              int nterpairs, int *rN, int *rC)
{
    int      i, j;
    gmx_bool bN, bC;

    for (i = 0; i < pdba->nr; i++)
    {
        bN = FALSE;
        for (j = 0; j < nterpairs && !bN; j++)
        {
            bN = pdba->atom[i].resind == rN[j];
        }
        bC = FALSE;
        for (j = 0; j < nterpairs && !bC; j++)
        {
            bC = pdba->atom[i].resind == rC[j];
        }

        /* add hacks to this atom */
        expand_hackblocks_one(&hb[pdba->atom[i].resind], *pdba->atomname[i],
                              &nab[i], &ab[i], bN, bC);
    }
    if (debug)
    {
        fprintf(debug, "\n");
    }
}

static int check_atoms_present(t_atoms *pdba, int nab[], t_hack *ab[])
{
    int i, j, k, d, rnr, nadd;

    nadd = 0;
    for (i = 0; i < pdba->nr; i++)
    {
        rnr = pdba->atom[i].resind;
        for (j = 0; j < nab[i]; j++)
        {
            if (ab[i][j].oname == NULL)
            {
                /* we're adding */
                if (ab[i][j].nname == NULL)
                {
                    gmx_incons("ab[i][j].nname not allocated");
                }
                /* check if the atom is already present */
                k = pdbasearch_atom(ab[i][j].nname, rnr, pdba, "check", TRUE);
                if (k != -1)
                {
                    /* We found the added atom. */
                    ab[i][j].bAlreadyPresent = TRUE;
                    if (debug)
                    {
                        fprintf(debug, "Atom '%s' in residue '%s' %d is already present\n",
                                ab[i][j].nname,
                                *pdba->resinfo[rnr].name, pdba->resinfo[rnr].nr);
                    }
                }
                else
                {
                    ab[i][j].bAlreadyPresent = FALSE;
                    /* count how many atoms we'll add */
                    nadd++;
                }
            }
            else if (ab[i][j].nname == NULL)
            {
                /* we're deleting */
                nadd--;
            }
        }
    }

    return nadd;
}

static void calc_all_pos(t_atoms *pdba, rvec x[], int nab[], t_hack *ab[],
                         gmx_bool bCheckMissing)
{
    int      i, j, ii, jj, m, ia, d, rnr, l = 0;
#define MAXH 4
    rvec     xa[4];    /* control atoms for calc_h_pos */
    rvec     xh[MAXH]; /* hydrogen positions from calc_h_pos */
    gmx_bool bFoundAll;

    jj = 0;

    for (i = 0; i < pdba->nr; i++)
    {
        rnr   = pdba->atom[i].resind;
        for (j = 0; j < nab[i]; j += ab[i][j].nr)
        {
            /* check if we're adding: */
            if (ab[i][j].oname == NULL && ab[i][j].tp > 0)
            {
                bFoundAll = TRUE;
                for (m = 0; (m < ab[i][j].nctl && bFoundAll); m++)
                {
                    ia = pdbasearch_atom(ab[i][j].a[m], rnr, pdba,
                                         bCheckMissing ? "atom" : "check",
                                         !bCheckMissing);
                    if (ia < 0)
                    {
                        /* not found in original atoms, might still be in t_hack (ab) */
                        hacksearch_atom(&ii, &jj, ab[i][j].a[m], nab, ab, rnr, pdba);
                        if (ii >= 0)
                        {
                            copy_rvec(ab[ii][jj].newx, xa[m]);
                        }
                        else
                        {
                            bFoundAll = FALSE;
                            if (bCheckMissing)
                            {
                                gmx_fatal(FARGS, "Atom %s not found in residue %s %d"
                                          ", rtp entry %s"
                                          " while adding hydrogens",
                                          ab[i][j].a[m],
                                          *pdba->resinfo[rnr].name,
                                          pdba->resinfo[rnr].nr,
                                          *pdba->resinfo[rnr].rtp);
                            }
                        }
                    }
                    else
                    {
                        copy_rvec(x[ia], xa[m]);
                    }
                }
                if (bFoundAll)
                {
                    for (m = 0; (m < MAXH); m++)
                    {
                        for (d = 0; d < DIM; d++)
                        {
                            if (m < ab[i][j].nr)
                            {
                                xh[m][d] = 0;
                            }
                            else
                            {
                                xh[m][d] = NOTSET;
                            }
                        }
                    }
                    calc_h_pos(ab[i][j].tp, xa, xh, &l);
                    for (m = 0; m < ab[i][j].nr; m++)
                    {
                        copy_rvec(xh[m], ab[i][j+m].newx);
                        ab[i][j+m].bXSet = TRUE;
                    }
                }
            }
        }
    }
}

static void free_ab(int natoms, int *nab, t_hack **ab)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        free_t_hack(nab[i], &ab[i]);
    }
    sfree(nab);
    sfree(ab);
}

static int add_h_low(t_atoms **pdbaptr, rvec *xptr[],
                     int nah, t_hackblock ah[],
                     int nterpairs, t_hackblock **ntdb, t_hackblock **ctdb,
                     int *rN, int *rC, gmx_bool bCheckMissing,
                     int **nabptr, t_hack ***abptr,
                     gmx_bool bUpdate_pdba, gmx_bool bKeep_old_pdba)
{
    t_atoms        *newpdba = NULL, *pdba = NULL;
    int             nadd;
    int             i, newi, j, d, natoms, nalreadypresent;
    int            *nab = NULL;
    t_hack        **ab  = NULL;
    t_hackblock    *hb;
    rvec           *xn;
    gmx_bool        bKeep_ab;

    /* set flags for adding hydrogens (according to hdb) */
    pdba   = *pdbaptr;
    natoms = pdba->nr;

    if (nabptr && abptr)
    {
        /* the first time these will be pointers to NULL, but we will
           return in them the completed arrays, which we will get back
           the second time */
        nab      = *nabptr;
        ab       = *abptr;
        bKeep_ab = TRUE;
        if (debug)
        {
            fprintf(debug, "pointer to ab found\n");
        }
    }
    else
    {
        bKeep_ab = FALSE;
    }

    if (nab && ab)
    {
        /* WOW, everything was already figured out */
        bUpdate_pdba = FALSE;
        if (debug)
        {
            fprintf(debug, "pointer to non-null ab found\n");
        }
    }
    else
    {
        /* We'll have to do all the hard work */
        bUpdate_pdba = TRUE;
        /* first get all the hackblocks for each residue: */
        hb = get_hackblocks(pdba, nah, ah, nterpairs, ntdb, ctdb, rN, rC);
        if (debug)
        {
            dump_hb(debug, pdba->nres, hb);
        }

        /* expand the hackblocks to atom level */
        snew(nab, natoms);
        snew(ab, natoms);
        expand_hackblocks(pdba, hb, nab, ab, nterpairs, rN, rC);
        free_t_hackblock(pdba->nres, &hb);
    }

    if (debug)
    {
        fprintf(debug, "before calc_all_pos\n");
        dump_ab(debug, natoms, nab, ab, TRUE);
    }

    /* Now calc the positions */
    calc_all_pos(pdba, *xptr, nab, ab, bCheckMissing);

    if (debug)
    {
        fprintf(debug, "after calc_all_pos\n");
        dump_ab(debug, natoms, nab, ab, TRUE);
    }

    if (bUpdate_pdba)
    {
        /* we don't have to add atoms that are already present in pdba,
           so we will remove them from the ab (t_hack) */
        nadd = check_atoms_present(pdba, nab, ab);
        if (debug)
        {
            fprintf(debug, "removed add hacks that were already in pdba:\n");
            dump_ab(debug, natoms, nab, ab, TRUE);
            fprintf(debug, "will be adding %d atoms\n", nadd);
        }

        /* Copy old atoms, making space for new ones */
        snew(newpdba, 1);
        init_t_atoms(newpdba, natoms+nadd, FALSE);
        newpdba->nres    = pdba->nres;
        sfree(newpdba->resinfo);
        newpdba->resinfo = pdba->resinfo;
    }
    else
    {
        nadd = 0;
    }
    if (debug)
    {
        fprintf(debug, "snew xn for %d old + %d new atoms %d total)\n",
                natoms, nadd, natoms+nadd);
    }

    if (nadd == 0)
    {
        /* There is nothing to do: return now */
        if (!bKeep_ab)
        {
            free_ab(natoms, nab, ab);
        }

        return natoms;
    }

    snew(xn, natoms+nadd);
    newi = 0;
    for (i = 0; (i < natoms); i++)
    {
        /* check if this atom wasn't scheduled for deletion */
        if (nab[i] == 0 || (ab[i][0].nname != NULL) )
        {
            if (newi >= natoms+nadd)
            {
                /*gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                nadd += 10;
                srenew(xn, natoms+nadd);
                if (bUpdate_pdba)
                {
                    srenew(newpdba->atom, natoms+nadd);
                    srenew(newpdba->atomname, natoms+nadd);
                }
                debug_gmx();
            }
            if (debug)
            {
                fprintf(debug, "(%3d) %3d %4s %4s%3d %3d",
                        i+1, newi+1, *pdba->atomname[i],
                        *pdba->resinfo[pdba->atom[i].resind].name,
                        pdba->resinfo[pdba->atom[i].resind].nr, nab[i]);
            }
            if (bUpdate_pdba)
            {
                copy_atom(pdba, i, newpdba, newi);
            }
            copy_rvec((*xptr)[i], xn[newi]);
            /* process the hacks for this atom */
            nalreadypresent = 0;
            for (j = 0; j < nab[i]; j++)
            {
                if (ab[i][j].oname == NULL) /* add */
                {
                    newi++;
                    if (newi >= natoms+nadd)
                    {
                        /* gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                        nadd += 10;
                        srenew(xn, natoms+nadd);
                        if (bUpdate_pdba)
                        {
                            srenew(newpdba->atom, natoms+nadd);
                            srenew(newpdba->atomname, natoms+nadd);
                        }
                        debug_gmx();
                    }
                    if (bUpdate_pdba)
                    {
                        newpdba->atom[newi].resind = pdba->atom[i].resind;
                    }
                    if (debug)
                    {
                        fprintf(debug, " + %d", newi+1);
                    }
                }
                if (ab[i][j].nname != NULL &&
                    (ab[i][j].oname == NULL ||
                     strcmp(ab[i][j].oname, *newpdba->atomname[newi]) == 0))
                {
                    /* add or replace */
                    if (ab[i][j].oname == NULL && ab[i][j].bAlreadyPresent)
                    {
                        /* This atom is already present, copy it from the input. */
                        nalreadypresent++;
                        if (bUpdate_pdba)
                        {
                            copy_atom(pdba, i+nalreadypresent, newpdba, newi);
                        }
                        copy_rvec((*xptr)[i+nalreadypresent], xn[newi]);
                    }
                    else
                    {
                        if (bUpdate_pdba)
                        {
                            if (gmx_debug_at)
                            {
                                fprintf(debug, "Replacing %d '%s' with (old name '%s') %s\n",
                                        newi,
                                        (newpdba->atomname[newi] && *newpdba->atomname[newi]) ? *newpdba->atomname[newi] : "",
                                        ab[i][j].oname ? ab[i][j].oname : "",
                                        ab[i][j].nname);
                            }
                            snew(newpdba->atomname[newi], 1);
                            *newpdba->atomname[newi] = gmx_strdup(ab[i][j].nname);
                            if (ab[i][j].oname != NULL && ab[i][j].atom) /* replace */
                            {                                            /*          newpdba->atom[newi].m    = ab[i][j].atom->m; */
/*        newpdba->atom[newi].q    = ab[i][j].atom->q; */
/*        newpdba->atom[newi].type = ab[i][j].atom->type; */
                            }
                        }
                        if (ab[i][j].bXSet)
                        {
                            copy_rvec(ab[i][j].newx, xn[newi]);
                        }
                    }
                    if (bUpdate_pdba && debug)
                    {
                        fprintf(debug, " %s %g %g", *newpdba->atomname[newi],
                                newpdba->atom[newi].m, newpdba->atom[newi].q);
                    }
                }
            }
            newi++;
            i += nalreadypresent;
            if (debug)
            {
                fprintf(debug, "\n");
            }
        }
    }
    if (bUpdate_pdba)
    {
        newpdba->nr = newi;
    }

    if (bKeep_ab)
    {
        *nabptr = nab;
        *abptr  = ab;
    }
    else
    {
        /* Clean up */
        free_ab(natoms, nab, ab);
    }

    if (bUpdate_pdba)
    {
        if (!bKeep_old_pdba)
        {
            for (i = 0; i < natoms; i++)
            {
                /* Do not free the atomname string itself, it might be in symtab */
                /* sfree(*(pdba->atomname[i])); */
                /* sfree(pdba->atomname[i]); */
            }
            sfree(pdba->atomname);
            sfree(pdba->atom);
            sfree(pdba->pdbinfo);
            sfree(pdba);
        }
        *pdbaptr = newpdba;
    }
    else
    {
        nadd = newi-natoms;
    }

    sfree(*xptr);
    *xptr = xn;

    return newi;
}

void deprotonate(t_atoms *atoms, rvec *x)
{
    int  i, j;

    j = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        if ( (*atoms->atomname[i])[0] != 'H')
        {
            atoms->atomname[j] = atoms->atomname[i];
            atoms->atom[j]     = atoms->atom[i];
            copy_rvec(x[i], x[j]);
            j++;
        }
    }
    atoms->nr = j;
}

int add_h(t_atoms **pdbaptr, rvec *xptr[],
          int nah, t_hackblock ah[],
          int nterpairs, t_hackblock **ntdb, t_hackblock **ctdb,
          int *rN, int *rC, gmx_bool bAllowMissing,
          int **nabptr, t_hack ***abptr,
          gmx_bool bUpdate_pdba, gmx_bool bKeep_old_pdba)
{
    int nold, nnew, niter;

    /* Here we loop to be able to add atoms to added atoms.
     * We should not check for missing atoms here.
     */
    niter = 0;
    nnew  = 0;
    do
    {
        nold = nnew;
        nnew = add_h_low(pdbaptr, xptr, nah, ah, nterpairs, ntdb, ctdb, rN, rC, FALSE,
                         nabptr, abptr, bUpdate_pdba, bKeep_old_pdba);
        niter++;
        if (niter > 100)
        {
            gmx_fatal(FARGS, "More than 100 iterations of add_h. Maybe you are trying to replace an added atom (this is not supported)?");
        }
    }
    while (nnew > nold);

    if (!bAllowMissing)
    {
        /* Call add_h_low once more, now only for the missing atoms check */
        add_h_low(pdbaptr, xptr, nah, ah, nterpairs, ntdb, ctdb, rN, rC, TRUE,
                  nabptr, abptr, bUpdate_pdba, bKeep_old_pdba);
    }

    return nnew;
}

int protonate(t_atoms **atomsptr, rvec **xptr, t_protonate *protdata)
{
#define NTERPAIRS 1
    t_atoms    *atoms;
    gmx_bool    bUpdate_pdba, bKeep_old_pdba;
    int         nntdb, nctdb, nt, ct;
    int         nadd;

    atoms = NULL;
    if (!protdata->bInit)
    {
        if (debug)
        {
            fprintf(debug, "protonate: Initializing protdata\n");
        }

        /* set forcefield to use: */
        strcpy(protdata->FF, "oplsaa.ff");

        /* get the databases: */
        protdata->nah = read_h_db(protdata->FF, &protdata->ah);
        open_symtab(&protdata->tab);
        protdata->atype = read_atype(protdata->FF, &protdata->tab);
        nntdb           = read_ter_db(protdata->FF, 'n', &protdata->ntdb, protdata->atype);
        if (nntdb < 1)
        {
            gmx_fatal(FARGS, "no N-terminus database");
        }
        nctdb = read_ter_db(protdata->FF, 'c', &protdata->ctdb, protdata->atype);
        if (nctdb < 1)
        {
            gmx_fatal(FARGS, "no C-terminus database");
        }

        /* set terminus types: -NH3+ (different for Proline) and -COO- */
        atoms = *atomsptr;
        snew(protdata->sel_ntdb, NTERPAIRS);
        snew(protdata->sel_ctdb, NTERPAIRS);

        if (nntdb >= 4 && nctdb >= 2)
        {
            /* Yuk, yuk, hardcoded default termini selections !!! */
            if (strncmp(*atoms->resinfo[atoms->atom[atoms->nr-1].resind].name, "PRO", 3) == 0)
            {
                nt = 3;
            }
            else
            {
                nt = 1;
            }
            ct = 1;
        }
        else
        {
            nt = 0;
            ct = 0;
        }
        protdata->sel_ntdb[0] = &(protdata->ntdb[nt]);
        protdata->sel_ctdb[0] = &(protdata->ctdb[ct]);

        /* set terminal residue numbers: */
        snew(protdata->rN, NTERPAIRS);
        snew(protdata->rC, NTERPAIRS);
        protdata->rN[0] = 0;
        protdata->rC[0] = atoms->atom[atoms->nr-1].resind;

        /* keep unprotonated topology: */
        protdata->upatoms = atoms;
        /* we don't have these yet: */
        protdata->patoms = NULL;
        bUpdate_pdba     = TRUE;
        bKeep_old_pdba   = TRUE;

        /* clear hackblocks: */
        protdata->nab = NULL;
        protdata->ab  = NULL;

        /* set flag to show we're initialized: */
        protdata->bInit = TRUE;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "protonate: using available protdata\n");
        }
        /* add_h will need the unprotonated topology again: */
        atoms          = protdata->upatoms;
        bUpdate_pdba   = FALSE;
        bKeep_old_pdba = FALSE;
    }

    /* now protonate */
    nadd = add_h(&atoms, xptr, protdata->nah, protdata->ah,
                 NTERPAIRS, protdata->sel_ntdb, protdata->sel_ctdb,
                 protdata->rN, protdata->rC, TRUE,
                 &protdata->nab, &protdata->ab, bUpdate_pdba, bKeep_old_pdba);
    if (!protdata->patoms)
    {
        /* store protonated topology */
        protdata->patoms = atoms;
    }
    *atomsptr = protdata->patoms;
    if (debug)
    {
        fprintf(debug, "natoms: %d -> %d (nadd=%d)\n",
                protdata->upatoms->nr, protdata->patoms->nr, nadd);
    }
    return nadd;
#undef NTERPAIRS
}
