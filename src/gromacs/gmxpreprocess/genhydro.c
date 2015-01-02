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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <time.h>

#include "typedefs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "calch.h"
#include "genhydro.h"
#include "h_db.h"
#include "ter_db.h"
#include "resall.h"
#include "pgutil.h"
#include "network.h"

#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* new structure for lone pair/virtual site construction */
typedef struct
{
    int     nr;     /* number of vsites to add to a particular atom (dim: natoms) */
    char  **names;  /* names of the vsites/LP to add */
    rvec   *x;      /* constructed coordinates of the vsites/LP to add */  
} t_vsiteadd;

static void copy_atom(t_atoms *atoms1, int a1, t_atoms *atoms2, int a2)
{
    atoms2->atom[a2] = atoms1->atom[a1];
    snew(atoms2->atomname[a2], 1);
    *atoms2->atomname[a2] = strdup(*atoms1->atomname[a1]);
}

static void copy_atom_drude(t_atoms *atoms1, int a1, t_atoms *atoms2, int a2, const char *dname)
{
    atoms2->atom[a2].type = 'D';
    atoms2->atom[a2].ptype = eptShell;
    atoms2->atom[a2].resind = atoms1->atom[a1].resind;
    snew(atoms2->atomname[a2], 1);
    *atoms2->atomname[a2] = strdup(dname);
}

static void add_lonepair_pdba(t_atoms *atoms1, int a1, t_atoms *atoms2, int a2, const char *lpname)
{
    atoms2->atom[a2].type = 'V';
    atoms2->atom[a2].ptype = eptVSite;
    atoms2->atom[a2].resind = atoms1->atom[a1].resind;
    snew(atoms2->atomname[a2], 1);
    *atoms2->atomname[a2] = strdup(lpname);
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
                hb[rnr].name = strdup(ahptr->name);
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
           and first control atom (AI) matches this atom or
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
                        (*abi)[*nabi + k].nname    = strdup(atomname);
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
                    (*abi)[*nabi + k].nname = strdup(hbr->hack[j].nname);
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
                    calc_h_pos(ab[i][j].tp, ab[i][j].nr, xa, xh, &l);
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
            /* TODO: BROKEN */
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
                                fprintf(debug, "\nReplacing %d '%s' with (old name '%s') %s\n",
                                        newi,
                                        (newpdba->atomname[newi] && *newpdba->atomname[newi]) ? *newpdba->atomname[newi] : "",
                                        ab[i][j].oname ? ab[i][j].oname : "", ab[i][j].nname);
                            }
                            snew(newpdba->atomname[newi], 1);
                            *newpdba->atomname[newi] = strdup(ab[i][j].nname);
                            if (ab[i][j].oname != NULL && ab[i][j].atom) /* replace */
                            {                                            
                                /* newpdba->atom[newi].m    = ab[i][j].atom->m; */
                                /* newpdba->atom[newi].q    = ab[i][j].atom->q; */
                                /* newpdba->atom[newi].type = ab[i][j].atom->type; */
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

/* Adding Drudes is a lot like adding H, so I put it here */
void add_drudes(t_atoms **pdbaptr, rvec *xptr[])
{

    int         i, j, natoms, nadd;
    int         resnr;
    char       *dname;              /* Drude name, created from heavy atom */
    rvec       *xn;                 /* new atom positions */
    gmx_bool   *bHeavy;             /* if an atom is a heavy atom (non-H, non-Drude, non-LP) */
    gmx_bool   *bPresent;           /* flag to check if Drude is already present on its parent heavy atom in input pdba */
                                    /* NOTE that the indices in this array are for heavy atoms, not that index == Drude */
    t_atoms    *newpdba = NULL, *pdba = NULL;

    pdba = *pdbaptr;
    natoms = pdba->nr;

    /* ugly hack - names should be <5 char */
    snew(dname, 8);

    snew(xn, natoms);
    snew(newpdba, 1);
    snew(bHeavy, natoms);
    snew(bPresent, natoms);

    nadd = 0;
    /* loop over all atoms, identify heavy atoms for adding Drudes */
    for (i=0; i<natoms; i++)
    {
        bPresent[i] = FALSE;
        bHeavy[i] = FALSE;

        /* Step 1. If we have a heavy atom, look for its associated Drude */
        if (!is_hydrogen(*pdba->atomname[i]) && !is_drude(*pdba->atomname[i])
            && !is_lonepair(*pdba->atomname[i]) && !is_dummymass(*pdba->atomname[i]))
        {
            bHeavy[i] = TRUE;
            resnr = pdba->resinfo[pdba->atom[i].resind].nr;

            if (debug)
            {
                fprintf(debug, "ADD DRUDES: looking for Drude on atom %s, i = %d, res = %d\n", 
                        *pdba->atomname[i], i, resnr); 
            }

            /* build name of associated Drude */
            strcpy(dname, "D");
            strcat(dname, *pdba->atomname[i]);

            if (debug)
            {
                fprintf(debug, "ADD DRUDES: Searching for Drude on atom %i, dname = %s aname = %s\n", i, dname, *pdba->atomname[i]);
            }

            /* If H are added to a structure with Drudes, the Drude will no longer be located at
             * i+1 since the H follow consecutively after atom i.  So we need to loop over the residue
             * to find any instance of dname.  It will be pretty uncommon for an input structure to
             * have Drudes but not H, but this method is still much more robust than making assumptions */
            j = 0;
            while ((j<natoms) && (!bPresent[i]))
            {
                if ((strcmp(dname, *pdba->atomname[j])==0) && (pdba->resinfo[pdba->atom[j].resind].nr == resnr))
                {
                    if (debug)
                    {
                        fprintf(debug, "ADD DRUDES: %s (%s) found for atom %d, res %d\n",
                                        dname, *pdba->atomname[j], j, pdba->resinfo[pdba->atom[j].resind].nr);
                    }

                    bPresent[i] = TRUE; /* exit the loop */
                }
                else
                {
                    if (debug)
                    {
                        fprintf(debug, "ADD DRUDES: %s (%s) is not a match at atom %d, res %d\n",
                                        dname, *pdba->atomname[j], j, pdba->resinfo[pdba->atom[j].resind].nr);
                    }

                    /* we need to add a Drude to this heavy atom */
                    bPresent[i] = FALSE;
                    j++;    /* keep looking */
                }
            }

            if (!bPresent[i])
            {
                /* we have a Drude that needs to be added */
                nadd++;
            }
        }
    }

    /* set up structures for adding atoms */
    srenew(xn, natoms+nadd);
    init_t_atoms(newpdba, natoms+nadd, FALSE);
    newpdba->nres = pdba->nres;
    sfree(newpdba->resinfo);
    newpdba->resinfo = pdba->resinfo;

    if (nadd > 0)
    {
        if (debug)
        {
            fprintf(debug, "ADD DRUDES: need to add %d Drudes\n", nadd);
        }
        srenew(newpdba->atom, natoms+nadd);
        srenew(newpdba->atomname, natoms+nadd); 
    }

    /* Loop through pdba to either copy atoms to newpdba or add to newpdba */
    j = 0;
    for (i=0; i<natoms; i++)
    {
        /* If Drude is not present, add it */
        if (bHeavy[i])
        {
            if (debug)
            {
                fprintf(debug, "ADD DRUDES: heavy atom %s (%d) found in adding step\n", *pdba->atomname[i], i);
            }
            if (!bPresent[i]) /* the Drude on atom i is not already present */
            {
                /* keep the heavy atom */
                copy_atom(pdba, i, newpdba, i+j);
                /* add the Drude to newpdba */
                strcpy(dname, "D");
                strcat(dname, *pdba->atomname[i]);
                if (debug)
                {
                    fprintf(debug, "ADD DRUDES: adding Drude %s to newpdba at %d\n", dname, i+j+1);
                }
                copy_atom_drude(pdba, i, newpdba, i+j+1, dname);

                /* Copy coords of heavy atom to Drude */
                copy_rvec((*xptr)[i], xn[i+j]);   /* preserve heavy atom coords */
                copy_rvec((*xptr)[i], xn[i+j+1]); /* use heavy atom coords for Drude */
                j++;    /* if we make an addition, shift indices */
            }
            else    /* the Drude is already there, so just copy to new structures */
            {
                copy_atom(pdba, i, newpdba, i+j);
                copy_rvec((*xptr)[i], xn[i+j]);
            }
        }
        else    /* not a heavy atom, thus no Drude, either */
        {
            copy_atom(pdba, i, newpdba, i+j);
            copy_rvec((*xptr)[i], xn[i+j]);
        }
    }

    /* update coords */
    sfree(*xptr);
    *xptr = xn;

    /* update pdba */
    *pdbaptr = newpdba;

    sfree(dname);
    sfree(bHeavy);
    sfree(bPresent);
}

/* Build lone pair (virtual site) coordinates based on input constructing atoms */
static void build_lonepair_coords(int bt, int f, atom_id ai[], real *params, rvec *xptr[], rvec *xlp)
{

    real    a, b, c, d, a1, b1, c1, invdij;
    rvec    xi, xj, xk, xl;     /* coords of constructing atoms, for convenience */
    rvec    xij, xjk, xik, xil, xp, temp;
    rvec    ra, rb, rja, rjb, rm;

    switch (bt)
    {
        case ebtsVSITE2:
            copy_rvec((*xptr)[ai[1]], xi);
            copy_rvec((*xptr)[ai[2]], xj);
            a = params[0];
            b = 1.0 - a;
            (*xlp)[XX] = b*xi[XX] + a*xj[XX];
            (*xlp)[YY] = b*xi[YY] + a*xj[XX];
            (*xlp)[ZZ] = b*xi[ZZ] + a*xj[XX];
            break;
        case ebtsVSITE3:
            /* always three constructing atoms for type 3 vsites */
            copy_rvec((*xptr)[ai[1]], xi);
            copy_rvec((*xptr)[ai[2]], xj);
            copy_rvec((*xptr)[ai[3]], xk);
            /* always at least 2 constructing constants, only 3out supplies 3 */
            a = params[0];
            b = params[1];

            if (debug)
            {
                fprintf(debug, "BUILD LP: received a = %f b = %f for type 3 construction.\n", a, b);
            }
            switch (f)
            {
                case 1:     /* 3 */
                    c = 1.0 - a - b;
                    (*xlp)[XX] = c*xi[XX] + a*xj[XX] + b*xk[XX];
                    (*xlp)[YY] = c*xi[YY] + a*xj[YY] + b*xk[YY];
                    (*xlp)[ZZ] = c*xi[ZZ] + a*xj[ZZ] + b*xk[ZZ];
                    break;
                case 2:     /* 3fd */
                    rvec_sub(xj, xi, xij);
                    rvec_sub(xk, xj, xjk);
                    temp[XX] = xij[XX] + a*xjk[XX];
                    temp[YY] = xij[YY] + a*xjk[YY];
                    temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
                    c = b*gmx_invsqrt(iprod(temp, temp));
                    (*xlp)[XX] = xi[XX] + c*temp[XX];
                    (*xlp)[YY] = xi[YY] + c*temp[YY];
                    (*xlp)[ZZ] = xi[ZZ] + c*temp[ZZ];
                    break;
                case 3:     /* 3fad */
                    rvec_sub(xj, xi, xij);
                    rvec_sub(xk, xj, xjk);
                    /* special case for calculating a and b, so replace values from above (same as convparm.c) */
                    a = params[1] * cos(DEG2RAD * params[0]);
                    b = params[1] * sin(DEG2RAD * params[0]);
                    invdij = gmx_invsqrt(iprod(xij, xij));
                    c1     = invdij * invdij * iprod(xij, xjk);
                    xp[XX] = xjk[XX] - c1*xij[XX];
                    xp[YY] = xjk[YY] - c1*xij[YY];
                    xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
                    a1     = a*invdij;
                    b1     = b*gmx_invsqrt(iprod(xp, xp));
                    (*xlp)[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
                    (*xlp)[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
                    (*xlp)[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
                    break;
                case 4:     /* 3out */
                    c = params[2];
                    rvec_sub(xj, xi, xij);
                    rvec_sub(xk, xi, xik);
                    cprod(xij, xik, temp);
                    (*xlp)[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
                    (*xlp)[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
                    (*xlp)[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
                    break;
                default:
                    gmx_fatal(FARGS, "Unknown subtype of virtual_sites3 in build_lonepair_coords.\n");
            }
            break;
        case ebtsVSITE4:
            /* constructing atoms */
            copy_rvec((*xptr)[ai[1]], xi);
            copy_rvec((*xptr)[ai[2]], xj);
            copy_rvec((*xptr)[ai[3]], xk);
            copy_rvec((*xptr)[ai[4]], xl);

            /* constants */
            a = params[0];
            b = params[1];
            c = params[2];

            rvec_sub(xj, xi, xij);
            rvec_sub(xk, xi, xik);
            rvec_sub(xl, xi, xil);
            ra[XX] = a*xik[XX];
            ra[YY] = a*xik[YY];
            ra[ZZ] = a*xik[ZZ];
            rb[XX] = b*xil[XX];
            rb[YY] = b*xil[YY];
            rb[ZZ] = b*xil[ZZ];
            rvec_sub(ra, xij, rja);
            rvec_sub(rb, xij, rjb);
            cprod(rja, rjb, rm);
            d = c*gmx_invsqrt(norm2(rm));
            (*xlp)[XX] = xi[XX] + d*rm[XX];
            (*xlp)[YY] = xi[YY] + d*rm[YY];
            (*xlp)[ZZ] = xi[ZZ] + d*rm[ZZ];
            break;
        default:
            gmx_fatal(FARGS, "Something totally weird in build_lonepair_coords.\n");
    }

    /* print coords */
    if (debug)
    {
        fprintf(debug, "BUILD LP: lp coords = %f %f %f\n", (*xlp)[XX], (*xlp)[YY], (*xlp)[ZZ]);
    }
}

/* Adding lone pairs is a bit different from H and Drudes, but since the Drude function
 * is here, I will just keep everything in one place - jal */
void add_drude_lonepairs(t_atoms **pdbaptr, rvec *xptr[], t_restp rtp[], int nssbonds, t_ssbond *ssbonds)
{

    int         f, i, j, k, m, ss, index, natoms, nadd, bt, start;
    int         nat;                /* number of atoms in vsite interaction definition */
    int         r1, r2;             /* residues involved in disulfide */
    atom_id     cb1, sg1, sg2, cb2; /* atom indices of CB and SG in a disulfide */
    atom_id     lpsa, lpsb;         /* atom indices of LPSA and LPSB of SG */
    atom_id     ai[MAXATOMLIST];
    real       *params;             /* constructing constants for vsite */
    rvec       *xlp;                /* coordinates of added lone pair */
    rvec       *xn;                 /* new coordinate array */
    char       *rname;              /* residue name */
    char       *lpname;             /* name of virtual site (LP) */

    /* lone pair construction constant strings */
    char ss_lpa[] = "4   -0.135847248   -0.131015228   -2.467798394";
    char ss_lpb[] = "4   -0.135847248   -0.131015228    2.467798394";

    t_atoms    *newpdba = NULL, *pdba = NULL;
    t_rbonded  *bvsite;
    t_restp    *rtp_tmp;
    t_vsiteadd *lp;

    /* for parsing parameter strings, note that 3, 3fd, and 3fad have the same format */
    const char *vsite2fmt = "%d %f";            /* ftype a     */
    const char *vsite3fadfmt = "%d %f %f";      /* ftype a b   */
    const char *vsite3outfmt = "%d %f %f %f";   /* ftype a b c */
    const char *vsite4fdnfmt = "%d %f %f %f";   /* ftype a b c */

    snew(lpname, 8);    /* another ugly hack */
    snew(xlp, 1);
    snew(bvsite, 1);

    pdba = *pdbaptr;
    natoms = pdba->nr;

    snew(xn, natoms);
    snew(newpdba, 1);
    snew(lp, natoms);

    newpdba->nres = pdba->nres;
    sfree(newpdba->resinfo);
    newpdba->resinfo = pdba->resinfo;

    /* initialize lp struct */
    for (i=0; i<natoms; i++)
    {
        lp[i].nr = 0;
        snew((lp[i]).names, 1);
        snew((lp[i]).x, 1);
    }

    nat = 0;
    nadd = 0;
    start = 0;
    /* find all of the lone pairs specified in the .rtp entries */
    for (i=0; i < pdba->nres; i++)
    {
        rtp_tmp = &rtp[i];

        /* Search .rtp entry for this residue for any ebtsVSITE, increment nadd
         * for each vsite found if it is not present in pdba */
        for (bt = ebtsVSITE2; bt <= ebtsVSITE4; bt++)
        {
            /* set number of atoms in vsite type */
            switch (bt)
            {
                case ebtsVSITE2:
                    nat = 3;
                    break;
                case ebtsVSITE3:
                    nat = 4;
                    break;
                case ebtsVSITE4:
                    nat = 5;
                    break;
                default:
                    gmx_fatal(FARGS, "Unknown vsite type in add_drude_lonepairs.\n");
            }

            /* only do this if there are vsites in the rtp entry */
            if (rtp_tmp->rb[bt].nb > 0)
            {
                /* Loop over nb to find atoms and get their indices */
                for (k = rtp_tmp->rb[bt].nb - 1; k >= 0; k--)
                {
                    *bvsite = rtp_tmp->rb[bt].b[k];

                    /* determine vsite sub-type from param string in rtp entry */
                    /* only type 3 has multiple sub-types */
                    if (bt == ebtsVSITE3)
                    {
                        sscanf(bvsite->s, "%d", &f);
                        switch (f)
                        {
                            /* the first three types (3, 3fd, 3fad) have the same format */
                            case 1:
                            case 2:
                            case 3:
                                snew(params, 2);
                                sscanf(bvsite->s, vsite3fadfmt, &f, &(params[0]), &(params[1]));
                                break;
                            case 4:
                                snew(params, 3);
                                sscanf(bvsite->s, vsite3outfmt, &f, &(params[0]), &(params[1]), &(params[2]));
                                break;
                            default:
                                gmx_fatal(FARGS, "Unknown vsite sub-type for type 3: %d\n", f);
                        }
                    }
                    else if (bt == ebtsVSITE2)
                    {
                        /* one constructing parameter */
                        snew(params, 1);
                        sscanf(bvsite->s, vsite2fmt, &f, &(params[0]));
                    }
                    else if (bt == ebtsVSITE4)
                    {
                        /* three constructing parameters */
                        snew(params, 3);
                        sscanf(bvsite->s, vsite4fdnfmt, &f, &(params[0]), &(params[1]), &(params[2]));
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Unknown vsite type in add_drude_lonepairs.\n");
                    }

                    /* initialize */
                    for (m = 0; m < MAXATOMLIST; m++)
                    {
                        ai[m] = NO_ATID; 
                    }

                    /* find the LP and constructing atoms for each vsite in the residue */
                    for (m = 0; m < nat; m++)
                    {
                        ai[m] = search_atom(bvsite->a[m], start, pdba, "check", TRUE);
                    }

                    /* Loop back over to add the LP if it is missing */
                    for (m = 0; m < nat; m++)
                    {
                        /* At this point, the only missing atoms should be lone pairs, and the atom
                           to which it is "bonded" is the next index (m+1), so the LP needs to be inserted
                           in newpdba at the position of ai[m+1] */
                        /* NOTE: The array a[m] holds atom names, ai[m] holds their indices */
                        if (ai[m] == NO_ATID)
                        {
                            /* Here, we have to guard against atoms that have been removed/replaced when combining
                             * .rtp and .tdb entries.  For instance, for a backbone carbonyl LP at the C-terminus,
                             * O has been renamed (so O will be NO_ATID) and the LP have been built already, so we
                             * we don't need to do this, but we still need to build other LP that may be in the .rtp
                             * entry (e.g. side chain LP). */
                            if (is_lonepair(bvsite->a[m]) && ai[m+1] != NO_ATID)
                            {
                                if (debug)
                                {
                                    fprintf(debug, "ADD DRUDE LP: adding %s to %s\n", bvsite->a[m], bvsite->a[m+1]);
                                    fprintf(debug, "ADD DRUDE LP: attached to %d (%s) %d (%s) %d (%s)\n", ai[1], bvsite->a[m+1],
                                            ai[2], bvsite->a[m+2], ai[3], bvsite->a[m+3]);
                                    fprintf(debug, "ADD DRUDE LP: x[%d]: %f %f %f\n", ai[1], (*xptr)[ai[1]][XX], (*xptr)[ai[1]][YY], (*xptr)[ai[1]][ZZ]);
                                    fprintf(debug, "ADD DRUDE LP: x[%d]: %f %f %f\n", ai[2], (*xptr)[ai[2]][XX], (*xptr)[ai[2]][YY], (*xptr)[ai[2]][ZZ]);
                                    fprintf(debug, "ADD DRUDE LP: x[%d]: %f %f %f\n", ai[3], (*xptr)[ai[3]][XX], (*xptr)[ai[3]][YY], (*xptr)[ai[3]][ZZ]);
                                }

                                /* we have a lone pair to add, so keep track of total number of additions (nadd)
                                 * and the number of lone pairs added to a given atom (lp[ai[m+1]]->nr) */
                                nadd++;
                                lp[ai[m+1]].nr++;
                                srenew(lp[ai[m+1]].names, lp[ai[m+1]].nr);
                                srenew(lp[ai[m+1]].x, lp[ai[m+1]].nr);
                                index = lp[ai[m+1]].nr - 1;

                                /* construct LP coordinates and set name */
                                build_lonepair_coords(bt, f, ai, params, xptr, xlp);
                                lpname = strdup(bvsite->a[m]);

                                if (debug)
                                {
                                    fprintf(debug, "ADD DRUDE LP: adding lp %s to atom %d, lp[%d].nr = %d\n", lpname, ai[m+1], ai[m+1], lp[ai[m+1]].nr);
                                }

                                /* copy name and coords to lp struct */
                                copy_rvec((*xlp), lp[ai[m+1]].x[index]);
                                snew(lp[ai[m+1]].names[index], 8);
                                lp[ai[m+1]].names[index] = strdup(lpname);
                            }
                        }
                    }
                }
            }
        }

        while ((start < pdba->nr) && (pdba->atom[start].resind == i))
        {
            start++;
        }

    }   /* end of loop over nres */

    /* Disulfides need to be treated separately since the .rtp for CYS2
     * cannot specify the lone pair construction, as it depends on SG
     * from the linked CYS2 */
    if (nssbonds > 0)
    {
        for (m = 0; m < MAXATOMLIST; m++)
        {
            ai[m] = NO_ATID;
        }

        for (ss=0; ss<nssbonds; ss++)
        {
            /* get atom indices for each half of the disulfide */
            r1 = ssbonds[ss].res1;
            r2 = ssbonds[ss].res2;

            /* check to see if LPSA is already there */
            lpsa = search_res_atom("LPSA", r1, pdba, "check", TRUE);

            if (lpsa == NO_ATID)
            { 
                /* LPSA construction constants */
                sscanf(ss_lpa, "%d %f %f %f", &f, &(params[0]), &(params[1]), &(params[2]));

                /* get constructing indices */
                cb1 = search_res_atom("CB", r1, pdba, "special bond", FALSE);
                sg1 = search_res_atom("SG", r1, pdba, "special bond", FALSE);
                cb2 = search_res_atom("CB", r2, pdba, "special bond", FALSE);
                sg2 = search_res_atom("SG", r2, pdba, "special bond", FALSE);

                ai[0] = NO_ATID;
                ai[1] = sg1;
                ai[2] = cb1;
                ai[3] = sg2;

                /* build it */
                nadd++;
                lp[sg1].nr++;
                srenew(lp[sg1].names, lp[sg1].nr);
                srenew(lp[sg1].x, lp[sg1].nr);
                index = lp[sg1].nr - 1;
                build_lonepair_coords(ebtsVSITE3, f, ai, params, xptr, xlp);

                /* add it */
                lpname = "LPSA";
                copy_rvec((*xlp), lp[sg1].x[index]);        
                snew(lp[sg1].names[index], 8);              
                lp[sg1].names[index] = strdup(lpname); 

                /* Now, the other half of the disulfide */
                ai[1] = sg2;
                ai[2] = cb2;
                ai[3] = sg1;

                nadd++;
                lp[sg2].nr++;
                srenew(lp[sg2].names, lp[sg2].nr);
                srenew(lp[sg2].x, lp[sg2].nr);
                index = lp[sg2].nr - 1;
                build_lonepair_coords(ebtsVSITE3, f, ai, params, xptr, xlp);

                copy_rvec((*xlp), lp[sg2].x[index]);
                snew(lp[sg2].names[index], 8);
                lp[sg2].names[index] = strdup(lpname);
            }

            /* same for LPSB */
            lpsb = search_res_atom("LPSB", r1, pdba, "check", TRUE);

            if (lpsb == NO_ATID)
            {            
                /* LPSB construction constants */
                sscanf(ss_lpb, "%d %f %f %f", &f, &(params[0]), &(params[1]), &(params[2]));

                /* get constructing indices */
                cb1 = search_res_atom("CB", r1, pdba, "special bond", FALSE);
                sg1 = search_res_atom("SG", r1, pdba, "special bond", FALSE);
                cb2 = search_res_atom("CB", r2, pdba, "special bond", FALSE);
                sg2 = search_res_atom("SG", r2, pdba, "special bond", FALSE);

                ai[0] = NO_ATID;
                ai[1] = sg1;
                ai[2] = cb1;
                ai[3] = sg2;

                /* build it */
                nadd++;
                lp[sg1].nr++;
                srenew(lp[sg1].names, lp[sg1].nr);
                srenew(lp[sg1].x, lp[sg1].nr);
                index = lp[sg1].nr - 1;
                build_lonepair_coords(ebtsVSITE3, f, ai, params, xptr, xlp);

                /* add it */
                lpname = "LPSB";
                copy_rvec((*xlp), lp[sg1].x[index]);
                snew(lp[sg1].names[index], 8);
                lp[sg1].names[index] = strdup(lpname);

                /* Now, the other half of the disulfide */
                ai[1] = sg2;
                ai[2] = cb2;
                ai[3] = sg1;

                nadd++;
                lp[sg2].nr++;
                srenew(lp[sg2].names, lp[sg2].nr);
                srenew(lp[sg2].x, lp[sg2].nr);
                index = lp[sg2].nr - 1;
                build_lonepair_coords(ebtsVSITE3, f, ai, params, xptr, xlp);

                copy_rvec((*xlp), lp[sg2].x[index]);
                snew(lp[sg2].names[index], 8);
                lp[sg2].names[index] = strdup(lpname);
            }
        }
    }

    /* re-allocate memory to prepare for additions */
    init_t_atoms(newpdba, natoms+nadd, FALSE);
    newpdba->nres = pdba->nres;
    sfree(newpdba->resinfo);
    newpdba->resinfo = pdba->resinfo;

    if (nadd > 0)
    {
        if (debug)
        {
            fprintf(debug, "ADD DRUDE LP: need to add %d lone pairs.\n", nadd);
        }
        srenew(newpdba->atom, natoms+nadd);
        srenew(newpdba->atomname, natoms+nadd);
    }

    srenew(xn, natoms+nadd);

    newpdba->nr = natoms+nadd;

    j = 0;  /* number of additions made for each atom */
    k = 0;  /* total number of additions made to pdba */ 
    /* loop over natoms and add the lone pairs from the lp structure to those atoms that need them */
    for (i=0; i<natoms; i++)
    {
        if (lp[i].nr > 0)
        {
            /* copy the parent atom to which additions are being made */
            copy_atom(pdba, i, newpdba, i+k);
            copy_rvec((*xptr)[i], xn[i+k]);
            /* now, add the lone pairs */
            for (j=(lp[i].nr-1); j >= 0; j--)
            {
                /* increment counter of total additions made */
                k++;
                /* add lone pair to newpdba and add coords to xn */
                add_lonepair_pdba(pdba, i, newpdba, i+k, lp[i].names[j]);
                copy_rvec(lp[i].x[j], xn[i+k]);
            }
        }
        else
        {
            /* nothing to add to this atom, just copy */
            copy_atom(pdba, i, newpdba, i+k);
            copy_rvec((*xptr)[i], xn[i+k]);
        }
    }

    /* update coords */
    sfree(*xptr);
    *xptr = xn;

    /* update pdba */
    *pdbaptr = newpdba;

    /* status check */
    if (debug)
    {
        fprintf(debug, "ADD DRUDE LP: End of function check of pdbaptr\n");
        fprintf(debug, "ADD DRUDE LP: New no. of atoms: %d\n", newpdba->nr);
        for (i=0; i<(newpdba->nr); i++)
        {
            fprintf(debug, "ADD DRUDE LP: Atom %d: %s", i, *(newpdba->atomname[i]));
            fprintf(debug, " x: %f %f %f\n", (*xptr)[i][XX], (*xptr)[i][YY], (*xptr)[i][ZZ]);
        }
    }

    sfree(params);

}
