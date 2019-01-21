/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018, by the GROMACS development team, led by
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

#include <cstring>
#include <ctime>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/calch.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/ter_db.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static int pdbasearch_atom(const char *name, int resind, gmx::ArrayRef<const AtomInfo> pdba,
                           const char *searchtype, bool bAllowMissing)
{
    int  i;

    for (i = 0; (i < pdba.size()) && (pdba[i].resind_ != resind); i++)
    {
        ;
    }

    return search_atom(name, i, pdba,
                       searchtype, bAllowMissing);
}

static void hacksearch_atom(int *ii, int *jj, char *name,
                            gmx::ArrayRef<const int> nab, gmx::ArrayRef<std::vector<t_hack>> ab,
                            int resind, gmx::ArrayRef<const AtomInfo> pdba)
{
    int  i, j;

    *ii = -1;
    if (name[0] == '-')
    {
        name++;
        resind--;
    }
    for (i = 0; (i < pdba.size()) && (pdba[i].resind_ != resind); i++)
    {
        ;
    }
    for (; (i < pdba.size()) && (pdba[i].resind_ == resind) && (*ii < 0); i++)
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

}

static std::vector<t_hackblock> get_hackblocks(gmx::ArrayRef<const Residue> resinfo,
                                   int nah, gmx::ArrayRef<t_hackblock> ah,
                                   int nterpairs,
                                   gmx::ArrayRef<std::vector<t_hackblock>> ntdb,
                                   gmx::ArrayRef<std::vector<t_hackblock>> ctdb,
                                   const int *rN, const int *rC)
{
    int          i, rnr;
    std::vector<t_hackblock> hb, ahptr;

    /* make space */
    hb.resize(resinfo.size());
    /* first the termini */
    for (i = 0; i < nterpairs; i++)
    {
        if (!ntdb[i].empty())
        {
            copy_t_hackblock(ntdb[i], &hb[rN[i]]);
        }
        if (!ctdb[i].empty())
        {
            merge_t_hackblock(ctdb[i], &hb[rC[i]]);
        }
    }
    /* then the whole hdb */
    for (rnr = 0; rnr < resinfo.size(); rnr++)
    {
        t_hackblock *ahptr = search_h_db(nah, ah, *resinfo[rnr].rtp_);
        if (ahptr)
        {
            if (hb[rnr].name == nullptr)
            {
                hb[rnr].name = gmx_strdup(ahptr->name);
            }
            merge_hacks(ahptr, &hb[rnr]);
        }
    }
    return hb;
}

static void expand_hackblocks_one(t_hackblock *hbr, char *atomname,
                                  int *nabi, std::vector<t_hack> *abi, bool bN, bool bC)
{
    int      j, k, l;
    bool     bIgnore;

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

        if (!bIgnore &&
            ( ( ( hbr->hack[j].tp > 0 || hbr->hack[j].oname == nullptr ) &&
                strcmp(atomname, hbr->hack[j].ai()) == 0 ) ||
              ( hbr->hack[j].oname != nullptr &&
                strcmp(atomname, hbr->hack[j].oname) == 0) ) )
        {
            /* now expand all hacks for this atom */
            abi->resize(*nabi + hbr->hack[j].nr);
            for (k = 0; k < hbr->hack[j].nr; k++)
            {
                copy_t_hack(&hbr->hack[j], &(*abi)[*nabi + k]);
                (*abi)[*nabi + k].bXSet = FALSE;
                /* if we're adding (oname==NULL) and don't have a new name (nname)
                   yet, build it from atomname */
                if ( (*abi)[*nabi + k].nname == nullptr)
                {
                    if ( (*abi)[*nabi + k].oname == nullptr)
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
            if (hbr->hack[j].tp > 0 || hbr->hack[j].oname == nullptr)
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

static void expand_hackblocks(gmx::ArrayRef<const AtomInfo> pdba, t_hackblock hb[],
                              gmx::ArrayRef<int> nab, gmx::ArrayRef<std::vector<t_hack>> ab,
                              int nterpairs, const int *rN, const int *rC)
{
    int i = 0;
    for (const auto &a : pdba)
    {
        bool bN = FALSE;
        for (int j = 0; j < nterpairs && !bN; j++)
        {
            bN = a.resind_ == rN[j];
        }
        bool bC = FALSE;
        for (int j = 0; j < nterpairs && !bC; j++)
        {
            bC = a.resind_ == rC[j];
        }

        /* add hacks to this atom */
        expand_hackblocks_one(&hb[a.resind_], *a.atomname,
                              &nab[i], &ab[i], bN, bC);
        i++;
    }
}

static int check_atoms_present(gmx::ArrayRef<const AtomInfo> pdba,
                               gmx::ArrayRef<const int> nab,
                               gmx::ArrayRef<std::vector<t_hack>> ab)
{
    int i, j, k, rnr, nadd;

    nadd = 0;
    for (i = 0; i < pdba.size(); i++)
    {
        rnr = pdba[i].resind_;
        for (j = 0; j < nab[i]; j++)
        {
            if (ab[i][j].oname == nullptr)
            {
                /* we're adding */
                if (ab[i][j].nname == nullptr)
                {
                    gmx_incons("ab[i][j].nname not allocated");
                }
                /* check if the atom is already present */
                k = pdbasearch_atom(ab[i][j].nname, rnr, pdba, "check", TRUE);
                if (k != -1)
                {
                    /* We found the added atom. */
                    ab[i][j].bAlreadyPresent = TRUE;
                }
                else
                {
                    ab[i][j].bAlreadyPresent = FALSE;
                    /* count how many atoms we'll add */
                    nadd++;
                }
            }
            else if (ab[i][j].nname == nullptr)
            {
                /* we're deleting */
                nadd--;
            }
        }
    }

    return nadd;
}

static void calc_all_pos(gmx::ArrayRef<const AtomInfo> pdba,
                         gmx::ArrayRef<const Residue> resinfo,
                         rvec x[],
                         gmx::ArrayRef<const int> nab,
                         gmx::ArrayRef<std::vector<t_hack>> ab,
                         bool bCheckMissing)
{
    int      i, j, ii, jj, m, ia, d, rnr, l = 0;
#define MAXH 4
    rvec     xa[4];    /* control atoms for calc_h_pos */
    rvec     xh[MAXH]; /* hydrogen positions from calc_h_pos */
    bool     bFoundAll;

    jj = 0;

    for (i = 0; i < pdba.size(); i++)
    {
        rnr   = pdba[i].resind_;
        for (j = 0; j < nab[i]; j += ab[i][j].nr)
        {
            /* check if we're adding: */
            if (ab[i][j].oname == nullptr && ab[i][j].tp > 0)
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
                                          *resinfo[rnr].name_,
                                          resinfo[rnr].nr_,
                                          *resinfo[rnr].rtp_);
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

static void free_ab(int natoms, gmx::ArrayRef<const int> nab, gmx::ArrayRef<std::vector<t_hack>> ab)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        free_t_hack(nab[i], ab[i]);
    }
}

static int add_h_low(std::vector<AtomInfo> *pdba,
                     gmx::ArrayRef<const Residue> resRef,
                     rvec *xptr[],
                     int nah, t_hackblock ah[],
                     int nterpairs, t_hackblock **ntdb, t_hackblock **ctdb,
                     int *rN, int *rC, bool bCheckMissing,
                     std::vector<int> *nabptr, std::vector<std::vector<t_hack>> *abptr,
                     bool bUpdate_pdba, bool bKeep_old_pdba)
{
    int             nadd;
    int             i, newi, j, natoms, nalreadypresent;
    std::vector<int> nab;
    std::vector<std::vector<t_hack>> ab;
    t_hackblock    *hb;
    rvec           *xn;
    bool            bKeep_ab;

    /* set flags for adding hydrogens (according to hdb) */
    natoms = pdba->size();
    gmx::ArrayRef<AtomInfo> pdbaRef = gmx::arrayRefFromArray(pdba->data(), pdba->size());

    if (nabptr && abptr)
    {
        /* the first time these will be pointers to NULL, but we will
           return in them the completed arrays, which we will get back
           the second time */
        nab      = *nabptr;
        ab       = *abptr;
        bKeep_ab = TRUE;
    }
    else
    {
        bKeep_ab = FALSE;
    }

    if (!nab.empty() && !ab.empty())
    {
        /* WOW, everything was already figured out */
        bUpdate_pdba = FALSE;
    }
    else
    {
        /* We'll have to do all the hard work */
        bUpdate_pdba = TRUE;
        /* first get all the hackblocks for each residue: */
        hb = get_hackblocks(resRef, nah, ah, nterpairs, ntdb, ctdb, rN, rC);

        /* expand the hackblocks to atom level */
        nab.resize(natoms);
        ab.resize(natoms);
        expand_hackblocks(pdbaRef, hb, nab, ab, nterpairs, rN, rC);
        free_t_hackblock(resRef.size(), &hb);
    }

    /* Now calc the positions */
    calc_all_pos(pdbaRef, resRef, *xptr, nab, ab, bCheckMissing);

    std::vector<AtomInfo> newPdba;
    std::vector<Residue> newRes;
    if (bUpdate_pdba)
    {
        /* we don't have to add atoms that are already present in pdba,
           so we will remove them from the ab (t_hack) */
        nadd = check_atoms_present(pdbaRef, nab, ab);
    }
    else
    {
        nadd = 0;
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
        if (nab[i] == 0 || (ab[i][0].nname != nullptr) )
        {
            if (newi >= natoms+nadd)
            {
                /*gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                nadd += 10;
                srenew(xn, natoms+nadd);
            }
            if (bUpdate_pdba)
            {
                newPdba.push_back(pdbaRef[i]);
            }
            copy_rvec((*xptr)[i], xn[newi]);
            /* process the hacks for this atom */
            nalreadypresent = 0;
            for (j = 0; j < nab[i]; j++)
            {
                if (ab[i][j].oname == nullptr) /* add */
                {
                    newi++;
                    if (newi >= natoms+nadd)
                    {
                        /* gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                        nadd += 10;
                        srenew(xn, natoms+nadd);
                        if (bUpdate_pdba)
                        {
                        }
                    }
                    if (bUpdate_pdba)
                    {
                        newPdba.back().resind_ = pdbaRef[i].resind_;
                    }
                }
                if (ab[i][j].nname != nullptr &&
                    (ab[i][j].oname == nullptr ||
                     strcmp(ab[i][j].oname, *newPdba.back().atomname) == 0))
                {
                    /* add or replace */
                    if (ab[i][j].oname == nullptr && ab[i][j].bAlreadyPresent)
                    {
                        /* This atom is already present, copy it from the input. */
                        nalreadypresent++;
                        if (bUpdate_pdba)
                        {
                            newPdba.push_back(pdbaRef[i+nalreadypresent]);
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
                                        (newPdba.back().atomname && *newPdba.back().atomname) ? *newPdba.back().atomname : "",
                                        ab[i][j].oname ? ab[i][j].oname : "",
                                        ab[i][j].nname);
                            }
                            snew(newPdba.back().atomname, 1);
                            *newPdba.back().atomname = gmx_strdup(ab[i][j].nname);
                        }
                        if (ab[i][j].bXSet)
                        {
                            copy_rvec(ab[i][j].newx, xn[newi]);
                        }
                    }
                    if (bUpdate_pdba && debug)
                    {
                        fprintf(debug, " %s %g %g", *newPdba.back().atomname,
                                newPdba.back().m_, newPdba.back().q_);
                    }
                }
            }
            newi++;
            i += nalreadypresent;
        }
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
        }
        *pdba = newPdba;
    }

    sfree(*xptr);
    *xptr = xn;

    return newi;
}

int add_h(std::vector<AtomInfo> *pdba,
          gmx::ArrayRef<const Residue> resRef,
          rvec *xptr[],
          int nah, t_hackblock ah[],
          int nterpairs, t_hackblock **ntdb, t_hackblock **ctdb,
          int *rN, int *rC, bool bAllowMissing,
          std::vector<int> *nabptr, std::vector<std::vector<t_hack>> *abptr, 
          bool bUpdate_pdba, bool bKeep_old_pdba)
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
        nnew = add_h_low(pdba, resRef, xptr, nah, ah, nterpairs, ntdb, ctdb, rN, rC, FALSE,
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
        add_h_low(pdba, resRef, xptr, nah, ah, nterpairs, ntdb, ctdb, rN, rC, TRUE,
                  nabptr, abptr, bUpdate_pdba, bKeep_old_pdba);
    }

    return nnew;
}
