/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "hackblock.h"

static void copy_atom(t_atoms *atoms1, int a1, t_atoms *atoms2, int a2)
{
    atoms2->atom[a2] = atoms1->atom[a1];
    snew(atoms2->atomname[a2], 1);
    *atoms2->atomname[a2] = gmx_strdup(*atoms1->atomname[a1]);
}

static int pdbasearch_atom(const char *name, int resind, t_atoms *pdba,
                           const char *searchtype, bool bAllowMissing)
{
    int  i;

    for (i = 0; (i < pdba->nr) && (pdba->atom[i].resind != resind); i++)
    {
        ;
    }

    return search_atom(name, i, pdba,
                       searchtype, bAllowMissing);
}

static void hacksearch_atom(int *ii, int *jj, const char *name,
                            gmx::ArrayRef < std::vector < HackBlock>> ab,
                            int resind, t_atoms *pdba)
{
    int  i;

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
        int j = 0;
        for (const auto &h : ab[i])
        {
            if (h.nname.compare(name))
            {
                *ii = i;
                *jj = j;
            }
            j++;
        }
    }

}

static std::vector<AtomModificationBlock>
getAtomModificationBlocks(t_atoms *pdba,
                          gmx::ArrayRef<AtomModificationBlock> amb,
                          int nterpairs,
                          gmx::ArrayRef<AtomModificationBlock *> ntdb,
                          gmx::ArrayRef<AtomModificationBlock *> ctdb,
                          const int *rN, const int *rC)
{
    std::vector<AtomModificationBlock> modBlock(pdba->nres);
    /* make space */
    /* first the termini */
    for (int i = 0; i < nterpairs; i++)
    {
        if (ntdb[i] != nullptr)
        {
            copyModificationBlocks(*ntdb[i], &modBlock[rN[i]]);
        }
        if (ctdb[i] != nullptr)
        {
            mergeAtomAndBondModifications(*ctdb[i], &modBlock[rC[i]]);
        }
    }
    /* then the whole hdb */
    for (int rnr = 0; rnr < pdba->nres; rnr++)
    {
        auto ahptr = search_h_db(amb, *pdba->resinfo[rnr].rtp);
        if (ahptr != amb.end())
        {
            if (amb[rnr].name.empty())
            {
                amb[rnr].name = ahptr->name;
            }
            mergeAtomModifications(*ahptr, &modBlock[rnr]);
        }
    }
    return modBlock;
}

static void expand_hackblocks_one(const AtomModificationBlock &hbr, const char *atomname,
                                  std::vector<HackBlock> *abi, bool bN, bool bC)
{
    /* we'll recursively add atoms to atoms */
    for (auto it = hbr.hack.begin(); it != hbr.hack.end(); it++)
    {
        /* first check if we're in the N- or C-terminus, then we should ignore
           all hacks involving atoms from resp. previous or next residue
           (i.e. which name begins with '-' (N) or '+' (C) */
        bool bIgnore = false;
        if (bN) /* N-terminus: ignore '-' */
        {
            for (int k = 0; k < 4 && !it->a[k].empty() && !bIgnore; k++)
            {
                bIgnore = it->a[k][0] == '-';
            }
        }
        if (bC) /* C-terminus: ignore '+' */
        {
            for (int k = 0; k < 4 && !it->a[k].empty() && !bIgnore; k++)
            {
                bIgnore = it->a[k][0] == '+';
            }
        }
        /* must be either hdb entry (tp>0) or add from tdb (oname==NULL)
           and first control aton (AI) matches this atom or
           delete/replace from tdb (oname!=NULL) and oname matches this atom */

        if (!bIgnore &&
            ( ( ( it->tp > 0 || it->oname.empty() ) &&
                it->a[0].compare(atomname) ) ||
              it->oname.compare(atomname)))
        {
            /* now expand all hacks for this atom */
            for (int k = 0; k < it->nr; k++)
            {
                abi->push_back(*it);
                HackBlock *hack = &abi->back();
                hack->bXSet = false;
                /* if we're adding (oname==NULL) and don't have a new name (nname)
                   yet, build it from atomname */
                if (hack->nname.empty())
                {
                    if (hack->oname.empty())
                    {
                        hack->nname    = atomname;
                        hack->nname[0] = 'H';
                    }
                }
                else
                {
                    if (gmx_debug_at)
                    {
                        int pos = std::distance(hbr.hack.begin(), it);
                        fprintf(debug, "Hack '%s' %d, replacing nname '%s' with '%s' (old name '%s')\n",
                                atomname, pos,
                                hack->nname.c_str(), it->nname.c_str(),
                                hack->oname.empty() ? "" : hack->oname.c_str());
                    }
                    hack->nname = it->nname;
                }

                if (it->tp == 10 && k == 2)
                {
                    /* This is a water virtual site, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    hack->nname.assign("M");
                }
                else if (it->tp == 11 && k >= 2)
                {
                    /* This is a water lone pair, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    hack->nname.assign(gmx::formatString("LP%d", 1+k-2));
                }
                else if (it->nr > 1)
                {
                    /* adding more than one atom, number them */
                    hack->nname.append(gmx::formatString("%d", 1+k));
                }
            }

            /* add hacks to atoms we've just added */
            if (it->tp > 0 || it->oname.empty())
            {
                for (int k = 0; k < it->nr; k++)
                {
                    expand_hackblocks_one(hbr, abi->at(abi->size() - it->nr + k).nname.c_str(),
                                          abi, bN, bC);
                }
            }
        }
    }
}

static void expand_hackblocks(t_atoms *pdba, gmx::ArrayRef<const AtomModificationBlock> hb,
                              gmx::ArrayRef < std::vector < HackBlock>> ab,
                              int nterpairs, const int *rN, const int *rC)
{
    for (int i = 0; i < pdba->nr; i++)
    {
        bool bN = false;
        for (int j = 0; j < nterpairs && !bN; j++)
        {
            bN = pdba->atom[i].resind == rN[j];
        }
        bool bC = false;
        for (int j = 0; j < nterpairs && !bC; j++)
        {
            bC = pdba->atom[i].resind == rC[j];
        }

        /* add hacks to this atom */
        expand_hackblocks_one(hb[pdba->atom[i].resind], *pdba->atomname[i],
                              &ab[i], bN, bC);
    }
}

static int check_atoms_present(t_atoms *pdba, gmx::ArrayRef < std::vector < HackBlock>> ab)
{
    int nadd = 0;
    for (int i = 0; i < pdba->nr; i++)
    {
        int rnr = pdba->atom[i].resind;
        for (auto it = ab[i].begin(); it != ab[i].end(); it++)
        {
            switch (it->type())
            {
                case HackType::Add:
                {
                    /* we're adding */
                    /* check if the atom is already present */
                    int k = pdbasearch_atom(it->nname.c_str(), rnr, pdba, "check", TRUE);
                    if (k != -1)
                    {
                        /* We found the added atom. */
                        it->bAlreadyPresent = true;
                    }
                    else
                    {
                        it->bAlreadyPresent = false;
                        /* count how many atoms we'll add */
                        nadd++;
                    }
                    break;
                }
                case HackType::Delete:
                {
                    /* we're deleting */
                    nadd--;
                    break;
                }
                case HackType::Replace:
                {
                    break;
                }
                default:
                {
                    GMX_THROW(gmx::InternalError("Case not handled"));
                }
            }
        }
    }
    return nadd;
}

static void calc_all_pos(t_atoms *pdba, rvec x[], gmx::ArrayRef < std::vector < HackBlock>> ab,
                         bool bCheckMissing)
{
    int      ii, l = 0;
#define MAXH 4
    rvec     xa[4];    /* control atoms for calc_h_pos */
    rvec     xh[MAXH]; /* hydrogen positions from calc_h_pos */

    int      jj = 0;

    for (int i = 0; i < pdba->nr; i++)
    {
        int rnr   = pdba->atom[i].resind;
        for (auto jt = ab[i].begin(); jt != ab[i].end(); jt += jt->nr)
        {
            /* check if we're adding: */
            if (jt->type() == HackType::Add && jt->tp > 0)
            {
                bool bFoundAll = true;
                for (int m = 0; (m < jt->nctl && bFoundAll); m++)
                {
                    int ia = pdbasearch_atom(jt->a[m].c_str(), rnr, pdba,
                                             bCheckMissing ? "atom" : "check",
                                             !bCheckMissing);
                    if (ia < 0)
                    {
                        /* not found in original atoms, might still be in t_hack (ab) */
                        hacksearch_atom(&ii, &jj, jt->a[m].c_str(), ab, rnr, pdba);
                        if (ii >= 0)
                        {
                            copy_rvec(ab[ii][jj].newx, xa[m]);
                        }
                        else
                        {
                            bFoundAll = false;
                            if (bCheckMissing)
                            {
                                gmx_fatal(FARGS, "Atom %s not found in residue %s %d"
                                          ", rtp entry %s"
                                          " while adding hydrogens",
                                          jt->a[m].c_str(),
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
                    for (int m = 0; (m < MAXH); m++)
                    {
                        for (int d = 0; d < DIM; d++)
                        {
                            if (m < jt->nr)
                            {
                                xh[m][d] = 0;
                            }
                            else
                            {
                                xh[m][d] = NOTSET;
                            }
                        }
                    }
                    calc_h_pos(jt->tp, xa, xh, &l);
                    for (int m = 0; m < jt->nr; m++)
                    {
                        auto next = jt + m;
                        copy_rvec(xh[m], next->newx);
                        next->bXSet = true;
                    }
                }
            }
        }
    }
}

static int add_h_low(t_atoms **pdbaptr, rvec *xptr[],
                     gmx::ArrayRef<AtomModificationBlock> amb,
                     int nterpairs,
                     std::vector<AtomModificationBlock *> ntdb,
                     std::vector<AtomModificationBlock *> ctdb,
                     int *rN, int *rC, bool bCheckMissing,
                     bool bUpdate_pdba, bool bKeep_old_pdba)
{
    t_atoms        *newpdba = nullptr, *pdba = nullptr;
    int             nadd;
    int             newi, natoms, nalreadypresent;
    std::vector < std::vector < HackBlock>> ab;
    rvec           *xn;

    /* set flags for adding hydrogens (according to hdb) */
    pdba   = *pdbaptr;
    natoms = pdba->nr;

    {
        /* We'll have to do all the hard work */
        bUpdate_pdba = true;
        /* first get all the hackblocks for each residue: */
        std::vector<AtomModificationBlock> hb =
            getAtomModificationBlocks(pdba, amb, nterpairs, ntdb, ctdb, rN, rC);

        /* expand the hackblocks to atom level */
        ab.resize(natoms);
        expand_hackblocks(pdba, hb, ab, nterpairs, rN, rC);
        freeModificationBlock(hb);
    }

    /* Now calc the positions */
    calc_all_pos(pdba, *xptr, ab, bCheckMissing);

    if (bUpdate_pdba)
    {
        /* we don't have to add atoms that are already present in pdba,
           so we will remove them from the ab (t_hack) */
        nadd = check_atoms_present(pdba, ab);

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

    if (nadd == 0)
    {
        return natoms;
    }

    snew(xn, natoms+nadd);
    newi = 0;
    for (int i = 0; (i < natoms); i++)
    {
        /* check if this atom wasn't scheduled for deletion */
        if (ab[i].empty()  || (!ab[i][0].nname.empty()) )
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
            }
            if (bUpdate_pdba)
            {
                copy_atom(pdba, i, newpdba, newi);
            }
            copy_rvec((*xptr)[i], xn[newi]);
            /* process the hacks for this atom */
            nalreadypresent = 0;
            for (auto it = ab[i].begin(); it != ab[i].end(); it++)
            {
                if (it->type() == HackType::Add) /* add */
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
                    }
                    if (bUpdate_pdba)
                    {
                        newpdba->atom[newi].resind = pdba->atom[i].resind;
                    }
                }
                if (!it->nname.empty() &&
                    (it->oname.empty() ||
                     it->oname.compare(*newpdba->atomname[newi])))
                {
                    /* add or replace */
                    if (it->type() == HackType::Add && it->bAlreadyPresent)
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
                                        it->oname.empty() ? "" : it->oname.c_str(),
                                        it->nname.c_str());
                            }
                            snew(newpdba->atomname[newi], 1);
                            *newpdba->atomname[newi] = gmx_strdup(it->nname.c_str());
                        }
                        if (it->bXSet)
                        {
                            copy_rvec(it->newx, xn[newi]);
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
        }
    }
    if (bUpdate_pdba)
    {
        newpdba->nr = newi;
    }

    if (bUpdate_pdba)
    {
        if (!bKeep_old_pdba)
        {
            sfree(pdba->atomname);
            sfree(pdba->atom);
            sfree(pdba->pdbinfo);
            sfree(pdba);
        }
        *pdbaptr = newpdba;
    }

    sfree(*xptr);
    *xptr = xn;

    return newi;
}

int add_h(t_atoms **pdbaptr, rvec *xptr[],
          gmx::ArrayRef<AtomModificationBlock> amb,
          int nterpairs,
          const std::vector<AtomModificationBlock *> &ntdb,
          const std::vector<AtomModificationBlock *> &ctdb,
          int *rN, int *rC, bool bAllowMissing,
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
        nnew = add_h_low(pdbaptr, xptr, amb, nterpairs, ntdb, ctdb, rN, rC, FALSE,
                         bUpdate_pdba, bKeep_old_pdba);
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
        add_h_low(pdbaptr, xptr, amb, nterpairs, ntdb, ctdb, rN, rC, TRUE,
                  bUpdate_pdba, bKeep_old_pdba);
    }

    return nnew;
}
