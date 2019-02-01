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
                            gmx::ArrayRef < std::vector < ModificationInstruction>> modInstructions,
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
        for (const auto &h : modInstructions[i])
        {
            if (!h.nname.empty() && strcmp(name, h.nname.c_str()) == 0)
            {
                *ii = i;
                *jj = j;
            }
            j++;
        }
    }

}

static std::vector<SystemModificationInstructions>
getSystemModificationInstructionss(t_atoms *pdba,
                                   gmx::ArrayRef<SystemModificationInstructions> globalModifications,
                                   int nterpairs,
                                   gmx::ArrayRef<SystemModificationInstructions *> ntdb,
                                   gmx::ArrayRef<SystemModificationInstructions *> ctdb,
                                   const int *rN, const int *rC)
{
    std::vector<SystemModificationInstructions> modBlock(pdba->nres);
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
        auto ahptr = search_h_db(globalModifications, *pdba->resinfo[rnr].rtp);
        if (ahptr != globalModifications.end())
        {
            if (globalModifications[rnr].name.empty())
            {
                globalModifications[rnr].name = ahptr->name;
            }
            mergeAtomModifications(*ahptr, &modBlock[rnr]);
        }
    }
    return modBlock;
}

static void expand_hackblocks_one(const SystemModificationInstructions &modInstruction,
                                  const char *atomname,
                                  std::vector<ModificationInstruction> *globalModifications,
                                  bool bN, bool bC)
{
    /* we'll recursively add atoms to atoms */
    int pos = 0;
    for (auto &modifiedResidue : modInstruction.hack)
    {
        /* first check if we're in the N- or C-terminus, then we should ignore
           all hacks involving atoms from resp. previous or next residue
           (i.e. which name begins with '-' (N) or '+' (C) */
        bool bIgnore = false;
        if (bN) /* N-terminus: ignore '-' */
        {
            for (int k = 0; k < 4 && !modifiedResidue.a[k].empty() && !bIgnore; k++)
            {
                bIgnore = modifiedResidue.a[k][0] == '-';
            }
        }
        if (bC) /* C-terminus: ignore '+' */
        {
            for (int k = 0; k < 4 && !modifiedResidue.a[k].empty() && !bIgnore; k++)
            {
                bIgnore = modifiedResidue.a[k][0] == '+';
            }
        }
        /* must be either hdb entry (tp>0) or add from tdb (oname==NULL)
           and first control aton (AI) matches this atom or
           delete/replace from tdb (oname!=NULL) and oname matches this atom */

        if (!bIgnore &&
            ( ( ( modifiedResidue.tp > 0 || modifiedResidue.oname.empty() ) &&
                strcmp(atomname, modifiedResidue.ai()) == 0 ) ||
              ( !modifiedResidue.oname.empty() &&
                strcmp(atomname, modifiedResidue.oname.c_str()) == 0) ) )
        {
            /* now expand all hacks for this atom */
            for (int k = 0; k < modifiedResidue.nr; k++)
            {
                globalModifications->push_back(modifiedResidue);
                ModificationInstruction *hack = &globalModifications->back();
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
                        fprintf(debug, "Hack '%s' %d, replacing nname '%s' with '%s' (old name '%s')\n",
                                atomname, pos,
                                hack->nname.c_str(), modifiedResidue.nname.c_str(),
                                hack->oname.empty() ? "" : hack->oname.c_str());
                    }
                    hack->nname = modifiedResidue.nname;
                }

                if (modifiedResidue.tp == 10 && k == 2)
                {
                    /* This is a water virtual site, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    hack->nname.assign("M");
                }
                else if (modifiedResidue.tp == 11 && k >= 2)
                {
                    /* This is a water lone pair, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    hack->nname.assign(gmx::formatString("LP%d", 1+k-2));
                }
                else if (modifiedResidue.nr > 1)
                {
                    /* adding more than one atom, number them */
                    hack->nname.append(gmx::formatString("%d", 1+k));
                }
            }

            /* add hacks to atoms we've just added */
            if (modifiedResidue.tp > 0 || modifiedResidue.oname.empty())
            {
                for (int k = 0; k < modifiedResidue.nr; k++)
                {
                    expand_hackblocks_one(modInstruction,
                                          globalModifications->at(
                                                  globalModifications->size() -
                                                  modifiedResidue.nr +
                                                  k).nname.c_str(),
                                          globalModifications, bN, bC);
                }
            }
        }
        pos++;
    }
}

static void expand_hackblocks(t_atoms *pdba, gmx::ArrayRef<const SystemModificationInstructions> hb,
                              gmx::ArrayRef < std::vector < ModificationInstruction>> modInstructions,
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
                              &modInstructions[i], bN, bC);
    }
}

static int check_atoms_present(t_atoms *pdba, gmx::ArrayRef < std::vector < ModificationInstruction>> modInstructions)
{
    int nadd = 0;
    for (int i = 0; i < pdba->nr; i++)
    {
        int rnr = pdba->atom[i].resind;
        for (auto modification = modInstructions[i].begin();
             modification != modInstructions[i].end();
             modification++)
        {
            switch (modification->type())
            {
                case ModificationType::Add:
                {
                    /* we're adding */
                    /* check if the atom is already present */
                    int k = pdbasearch_atom(modification->nname.c_str(), rnr, pdba, "check", TRUE);
                    if (k != -1)
                    {
                        /* We found the added atom. */
                        modification->bAlreadyPresent = true;
                    }
                    else
                    {
                        modification->bAlreadyPresent = false;
                        /* count how many atoms we'll add */
                        nadd++;
                    }
                    break;
                }
                case ModificationType::Delete:
                {
                    /* we're deleting */
                    nadd--;
                    break;
                }
                case ModificationType::Replace:
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

static void calc_all_pos(t_atoms *pdba, rvec x[], gmx::ArrayRef < std::vector < ModificationInstruction>> modInstructions,
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
        for (auto modification = modInstructions[i].begin(); modification != modInstructions[i].end();
             modification += modification->nr)
        {
            /* check if we're adding: */
            if (modification->type() == ModificationType::Add && modification->tp > 0)
            {
                bool bFoundAll = true;
                for (int m = 0; (m < modification->nctl && bFoundAll); m++)
                {
                    int ia = pdbasearch_atom(modification->a[m].c_str(), rnr, pdba,
                                             bCheckMissing ? "atom" : "check",
                                             !bCheckMissing);
                    if (ia < 0)
                    {
                        /* not found in original atoms, might still be in
                         * the modification Instructions (modInstructions) */
                        hacksearch_atom(&ii, &jj, modification->a[m].c_str(), modInstructions, rnr, pdba);
                        if (ii >= 0)
                        {
                            copy_rvec(modInstructions[ii][jj].newx, xa[m]);
                        }
                        else
                        {
                            bFoundAll = false;
                            if (bCheckMissing)
                            {
                                gmx_fatal(FARGS, "Atom %s not found in residue %s %d"
                                          ", rtp entry %s"
                                          " while adding hydrogens",
                                          modification->a[m].c_str(),
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
                            if (m < modification->nr)
                            {
                                xh[m][d] = 0;
                            }
                            else
                            {
                                xh[m][d] = NOTSET;
                            }
                        }
                    }
                    calc_h_pos(modification->tp, xa, xh, &l);
                    for (int m = 0; m < modification->nr; m++)
                    {
                        auto next = modification + m;
                        copy_rvec(xh[m], next->newx);
                        next->bXSet = true;
                    }
                }
            }
        }
    }
}

static int add_h_low(t_atoms **pdbaptr, rvec *xptr[],
                     gmx::ArrayRef<SystemModificationInstructions> globalModifications,
                     int nterpairs,
                     std::vector<SystemModificationInstructions *> ntdb,
                     std::vector<SystemModificationInstructions *> ctdb,
                     int *rN, int *rC, bool bCheckMissing,
                     bool bUpdate_pdba, bool bKeep_old_pdba)
{
    t_atoms        *newpdba = nullptr, *pdba = nullptr;
    int             nadd;
    int             newi, natoms, nalreadypresent;
    std::vector < std::vector < ModificationInstruction>> modInstructions;
    rvec           *xn;

    /* set flags for adding hydrogens (according to hdb) */
    pdba   = *pdbaptr;
    natoms = pdba->nr;

    {
        /* We'll have to do all the hard work */
        bUpdate_pdba = true;
        /* first get all the hackblocks for each residue: */
        std::vector<SystemModificationInstructions> hb =
            getSystemModificationInstructionss(pdba, globalModifications, nterpairs, ntdb, ctdb, rN, rC);

        /* expand the hackblocks to atom level */
        modInstructions.resize(natoms);
        expand_hackblocks(pdba, hb, modInstructions, nterpairs, rN, rC);
        freeModificationBlock(hb);
    }

    /* Now calc the positions */
    calc_all_pos(pdba, *xptr, modInstructions, bCheckMissing);

    if (bUpdate_pdba)
    {
        /* we don't have to add atoms that are already present in pdba,
           so we will remove them from the modInstructions (t_hack) */
        nadd = check_atoms_present(pdba, modInstructions);

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
        if (modInstructions[i].empty()  || (!modInstructions[i][0].nname.empty()) )
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
            for (auto modification = modInstructions[i].begin();
                 modification != modInstructions[i].end();
                 modification++)
            {
                if (modification->type() == ModificationType::Add) /* add */
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
                if (!modification->nname.empty() &&
                    (modification->oname.empty() ||
                     strcmp(modification->oname.c_str(), *newpdba->atomname[newi]) == 0))
                {
                    /* add or replace */
                    if (modification->type() == ModificationType::Add && modification->bAlreadyPresent)
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
                                        modification->oname.empty() ? "" : modification->oname.c_str(),
                                        modification->nname.c_str());
                            }
                            snew(newpdba->atomname[newi], 1);
                            *newpdba->atomname[newi] = gmx_strdup(modification->nname.c_str());
                        }
                        if (modification->bXSet)
                        {
                            copy_rvec(modification->newx, xn[newi]);
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
          gmx::ArrayRef<SystemModificationInstructions> globalModifications,
          int nterpairs,
          const std::vector<SystemModificationInstructions *> &ntdb,
          const std::vector<SystemModificationInstructions *> &ctdb,
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
        nnew = add_h_low(pdbaptr, xptr, globalModifications, nterpairs, ntdb, ctdb, rN, rC, FALSE,
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
        add_h_low(pdbaptr, xptr, globalModifications, nterpairs, ntdb, ctdb, rN, rC, TRUE,
                  bUpdate_pdba, bKeep_old_pdba);
    }

    return nnew;
}
