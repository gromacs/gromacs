/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "genhydro.h"

#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/calch.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/ter_db.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"
#include "resall.h"

static void copy_atom(const t_atoms* atoms1, int a1, t_atoms* atoms2, int a2, t_symtab* symtab)
{
    atoms2->atom[a2]     = atoms1->atom[a1];
    atoms2->atomname[a2] = put_symtab(symtab, *atoms1->atomname[a1]);
}

static int pdbasearch_atom(const char*              name,
                           int                      resind,
                           const t_atoms*           pdba,
                           const char*              searchtype,
                           bool                     bAllowMissing,
                           gmx::ArrayRef<const int> cyclicBondsIndex)
{
    int i;

    for (i = 0; (i < pdba->nr) && (pdba->atom[i].resind != resind); i++) {}

    return search_atom(name, i, pdba, searchtype, bAllowMissing, cyclicBondsIndex);
}

/*! \brief Return the index of the first atom whose residue index
 * matches and which has a patch with the given name.
 *
 * \param[out] ii      Index of the first atom in the residue that matches,
 *                       -1 if no match occurs.
 * \param[out] jj      Index of the patch that matches,
 *                       unchanged if no match occurs.
 * \param[in]  name    Name of the desired patch to match
 * \param[in]  patches The patch database to search
 * \param[in]  resind  The residue index to match
 * \param[in]  pdba    The atoms to work with
 *
 * \todo The short-circuit logic will be simpler if this returned a
 * std::pair<int, int> as soon as the first double match is found.
 */
static void hacksearch_atom(int*                                            ii,
                            int*                                            jj,
                            const char*                                     name,
                            gmx::ArrayRef<const std::vector<MoleculePatch>> patches,
                            int                                             resind,
                            const t_atoms*                                  pdba)
{
    int i;

    *ii = -1;
    if (name[0] == '-')
    {
        name++;
        resind--;
    }
    for (i = 0; (i < pdba->nr) && (pdba->atom[i].resind != resind); i++) {}
    for (; (i < pdba->nr) && (pdba->atom[i].resind == resind) && (*ii < 0); i++)
    {
        int j = 0;
        for (const auto& patch : patches[i])
        {
            if (patch.nname == name)
            {
                *ii = i;
                *jj = j;
                if (*ii >= 0)
                {
                    break;
                }
            }
            j++;
        }
    }
}

static std::vector<MoleculePatchDatabase>
getMoleculePatchDatabases(const t_atoms*                              pdba,
                          gmx::ArrayRef<const MoleculePatchDatabase>  globalPatches,
                          int                                         nterpairs,
                          gmx::ArrayRef<MoleculePatchDatabase* const> ntdb,
                          gmx::ArrayRef<MoleculePatchDatabase* const> ctdb,
                          gmx::ArrayRef<const int>                    rN,
                          gmx::ArrayRef<const int>                    rC)
{
    std::vector<MoleculePatchDatabase> modBlock(pdba->nres);
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
        auto ahptr = search_h_db(globalPatches, *pdba->resinfo[rnr].rtp);
        if (ahptr != globalPatches.end())
        {
            if (modBlock[rnr].name.empty())
            {
                modBlock[rnr].name = ahptr->name;
            }
            mergeAtomModifications(*ahptr, &modBlock[rnr]);
        }
    }
    return modBlock;
}

static void expand_hackblocks_one(const MoleculePatchDatabase& newPatch,
                                  const std::string localAtomName, //NOLINT(performance-unnecessary-value-param)
                                  std::vector<MoleculePatch>* globalPatches,
                                  bool                        bN,
                                  bool                        bC)
{
    /* we'll recursively add atoms to atoms */
    int pos = 0;
    for (const auto& singlePatch : newPatch.hack)
    {
        /* first check if we're in the N- or C-terminus, then we should ignore
           all hacks involving atoms from resp. previous or next residue
           (i.e. which name begins with '-' (N) or '+' (C) */
        bool bIgnore = false;
        if (bN) /* N-terminus: ignore '-' */
        {
            for (int k = 0; k < 4 && !singlePatch.a[k].empty() && !bIgnore; k++)
            {
                bIgnore = singlePatch.a[k][0] == '-';
            }
        }
        if (bC) /* C-terminus: ignore '+' */
        {
            for (int k = 0; k < 4 && !singlePatch.a[k].empty() && !bIgnore; k++)
            {
                bIgnore = singlePatch.a[k][0] == '+';
            }
        }
        /* must be either hdb entry (tp>0) or add from tdb (oname==NULL)
           and first control aton (AI) matches this atom or
           delete/replace from tdb (oname!=NULL) and oname matches this atom */

        if (!bIgnore
            && (((singlePatch.tp > 0 || singlePatch.oname.empty()) && singlePatch.a[0] == localAtomName)
                || (singlePatch.oname == localAtomName)))
        {
            /* now expand all hacks for this atom */
            for (int k = 0; k < singlePatch.nr; k++)
            {
                globalPatches->push_back(singlePatch);
                MoleculePatch* patch = &globalPatches->back();
                patch->bXSet         = false;
                /* if we're adding (oname==NULL) and don't have a new name (nname)
                   yet, build it from localAtomName */
                if (patch->nname.empty())
                {
                    if (patch->oname.empty())
                    {
                        patch->nname    = localAtomName;
                        patch->nname[0] = 'H';
                    }
                }
                else
                {
                    if (gmx_debug_at)
                    {
                        fprintf(debug,
                                "Hack '%s' %d, replacing nname '%s' with '%s' (old name '%s')\n",
                                localAtomName.c_str(),
                                pos,
                                patch->nname.c_str(),
                                singlePatch.nname.c_str(),
                                patch->oname.empty() ? "" : patch->oname.c_str());
                    }
                    patch->nname = singlePatch.nname;
                }

                if (singlePatch.tp == 10 && k == 2)
                {
                    /* This is a water virtual site, not a hydrogen */
                    /* Ugly hardcoded name hack to replace 'H' with 'M' */
                    GMX_RELEASE_ASSERT(
                            !patch->nname.empty() && patch->nname[0] == 'H',
                            "Water virtual site should be named starting with H at this point");
                    patch->nname[0] = 'M';
                }
                else if (singlePatch.tp == 11 && k >= 2)
                {
                    /* This is a water lone pair, not a hydrogen */
                    /* Ugly hardcoded name hack */
                    patch->nname.assign(gmx::formatString("LP%d", 1 + k - 2));
                }
                else if (singlePatch.nr > 1)
                {
                    /* adding more than one atom, number them */
                    patch->nname.append(gmx::formatString("%d", 1 + k));
                }
            }

            /* add hacks to atoms we've just added */
            if (singlePatch.tp > 0 || singlePatch.oname.empty())
            {
                for (int k = 0; k < singlePatch.nr; k++)
                {
                    expand_hackblocks_one(
                            newPatch,
                            globalPatches->at(globalPatches->size() - singlePatch.nr + k).nname,
                            globalPatches,
                            bN,
                            bC);
                }
            }
        }
        pos++;
    }
}

static void expand_hackblocks(const t_atoms*                             pdba,
                              gmx::ArrayRef<const MoleculePatchDatabase> hb,
                              gmx::ArrayRef<std::vector<MoleculePatch>>  patches,
                              int                                        nterpairs,
                              gmx::ArrayRef<const int>                   rN,
                              gmx::ArrayRef<const int>                   rC)
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
        expand_hackblocks_one(hb[pdba->atom[i].resind], *pdba->atomname[i], &patches[i], bN, bC);
    }
}

static int check_atoms_present(const t_atoms*                            pdba,
                               gmx::ArrayRef<std::vector<MoleculePatch>> patches,
                               gmx::ArrayRef<const int>                  cyclicBondsIndex)
{
    int nadd = 0;
    for (int i = 0; i < pdba->nr; i++)
    {
        int rnr = pdba->atom[i].resind;
        for (auto patch = patches[i].begin(); patch != patches[i].end(); patch++)
        {
            switch (patch->type())
            {
                case MoleculePatchType::Add:
                {
                    /* we're adding */
                    /* check if the atom is already present */
                    int k = pdbasearch_atom(patch->nname.c_str(), rnr, pdba, "check", TRUE, cyclicBondsIndex);
                    if (k != -1)
                    {
                        /* We found the added atom. */
                        patch->bAlreadyPresent = true;
                    }
                    else
                    {
                        patch->bAlreadyPresent = false;
                        /* count how many atoms we'll add */
                        nadd++;
                    }
                    break;
                }
                case MoleculePatchType::Delete:
                {
                    /* we're deleting */
                    nadd--;
                    break;
                }
                case MoleculePatchType::Replace:
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

static void calc_all_pos(const t_atoms*                            pdba,
                         gmx::ArrayRef<const gmx::RVec>            x,
                         gmx::ArrayRef<std::vector<MoleculePatch>> patches,
                         bool                                      bCheckMissing,
                         gmx::ArrayRef<const int>                  cyclicBondsIndex)
{
    int ii, l = 0;
#define MAXH 4
    rvec xa[4];    /* control atoms for calc_h_pos */
    rvec xh[MAXH]; /* hydrogen positions from calc_h_pos */

    int jj = 0;

    for (int i = 0; i < pdba->nr; i++)
    {
        int rnr = pdba->atom[i].resind;
        for (auto patch = patches[i].begin(); patch != patches[i].end(); patch += patch->nr)
        {
            GMX_RELEASE_ASSERT(patch < patches[i].end(),
                               "The number of patches in the last patch can not exceed the total "
                               "number of patches");
            /* check if we're adding: */
            if (patch->type() == MoleculePatchType::Add && patch->tp > 0)
            {
                bool bFoundAll = true;
                for (int m = 0; (m < patch->nctl && bFoundAll); m++)
                {
                    int ia = pdbasearch_atom(patch->a[m].c_str(),
                                             rnr,
                                             pdba,
                                             bCheckMissing ? "atom" : "check",
                                             !bCheckMissing,
                                             cyclicBondsIndex);
                    if (ia < 0)
                    {
                        /* not found in original atoms, might still be in
                         * the patch Instructions (patches) */
                        hacksearch_atom(&ii, &jj, patch->a[m].c_str(), patches, rnr, pdba);
                        if (ii >= 0)
                        {
                            copy_rvec(patches[ii][jj].newx, xa[m]);
                        }
                        else
                        {
                            bFoundAll = false;
                            if (bCheckMissing)
                            {
                                gmx_fatal(FARGS,
                                          "Atom %s not found in residue %s %d"
                                          ", rtp entry %s"
                                          " while adding hydrogens",
                                          patch->a[m].c_str(),
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
                            if (m < patch->nr)
                            {
                                xh[m][d] = 0;
                            }
                            else
                            {
                                xh[m][d] = NOTSET;
                            }
                        }
                    }
                    calc_h_pos(patch->tp, xa, xh, &l);
                    for (int m = 0; m < patch->nr; m++)
                    {
                        auto next = patch + m;
                        copy_rvec(xh[m], next->newx);
                        next->bXSet = true;
                    }
                }
            }
        }
    }
}

static int add_h_low(t_atoms**                                   initialAtoms,
                     t_atoms**                                   modifiedAtoms,
                     std::vector<gmx::RVec>*                     xptr,
                     gmx::ArrayRef<const MoleculePatchDatabase>  globalPatches,
                     t_symtab*                                   symtab,
                     const int                                   nterpairs,
                     gmx::ArrayRef<MoleculePatchDatabase* const> ntdb,
                     gmx::ArrayRef<MoleculePatchDatabase* const> ctdb,
                     gmx::ArrayRef<const int>                    rN,
                     gmx::ArrayRef<const int>                    rC,
                     const bool                                  bCheckMissing,
                     gmx::ArrayRef<const int>                    cyclicBondsIndex)
{
    int                                     nadd;
    int                                     newi, natoms, nalreadypresent;
    std::vector<std::vector<MoleculePatch>> patches;
    std::vector<gmx::RVec>                  xn;

    t_atoms* pdba = *initialAtoms;

    /* set flags for adding hydrogens (according to hdb) */
    natoms = pdba->nr;

    {
        /* We'll have to do all the hard work */
        /* first get all the hackblocks for each residue: */
        std::vector<MoleculePatchDatabase> hb =
                getMoleculePatchDatabases(pdba, globalPatches, nterpairs, ntdb, ctdb, rN, rC);

        /* expand the hackblocks to atom level */
        patches.resize(natoms);
        expand_hackblocks(pdba, hb, patches, nterpairs, rN, rC);
    }

    /* Now calc the positions */
    calc_all_pos(pdba, *xptr, patches, bCheckMissing, cyclicBondsIndex);

    /* we don't have to add atoms that are already present in initialAtoms,
       so we will remove them from the patches (MoleculePatch) */
    nadd = check_atoms_present(pdba, patches, cyclicBondsIndex);

    /* Copy old atoms, making space for new ones */
    if (nadd > 0)
    {
        srenew(*modifiedAtoms, 1);
        init_t_atoms(*modifiedAtoms, natoms + nadd, FALSE);
        (*modifiedAtoms)->nres = pdba->nres;
        srenew((*modifiedAtoms)->resinfo, pdba->nres);
        std::copy(pdba->resinfo, pdba->resinfo + pdba->nres, (*modifiedAtoms)->resinfo);
    }
    if (nadd == 0)
    {
        return natoms;
    }

    xn.resize(natoms + nadd);
    newi = 0;
    for (int i = 0; (i < natoms); i++)
    {
        /* check if this atom wasn't scheduled for deletion */
        if (patches[i].empty() || (!patches[i][0].nname.empty()))
        {
            if (newi >= natoms + nadd)
            {
                /*gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                nadd += 10;
                xn.resize(natoms + nadd);
                srenew((*modifiedAtoms)->atom, natoms + nadd);
                srenew((*modifiedAtoms)->atomname, natoms + nadd);
            }
            copy_atom(pdba, i, (*modifiedAtoms), newi, symtab);
            copy_rvec((*xptr)[i], xn[newi]);
            /* process the hacks for this atom */
            nalreadypresent = 0;
            for (auto patch = patches[i].begin(); patch != patches[i].end(); patch++)
            {
                if (patch->type() == MoleculePatchType::Add) /* add */
                {
                    newi++;
                    if (newi >= natoms + nadd)
                    {
                        /* gmx_fatal(FARGS,"Not enough space for adding atoms");*/
                        nadd += 10;
                        xn.resize(natoms + nadd);
                        srenew((*modifiedAtoms)->atom, natoms + nadd);
                        srenew((*modifiedAtoms)->atomname, natoms + nadd);
                    }
                    (*modifiedAtoms)->atom[newi].resind = pdba->atom[i].resind;
                }
                if (!patch->nname.empty()
                    && (patch->oname.empty() || patch->oname == *(*modifiedAtoms)->atomname[newi]))
                {
                    /* add or replace */
                    if (patch->type() == MoleculePatchType::Add && patch->bAlreadyPresent)
                    {
                        /* This atom is already present, copy it from the input. */
                        nalreadypresent++;
                        copy_atom(pdba, i + nalreadypresent, (*modifiedAtoms), newi, symtab);
                        copy_rvec((*xptr)[i + nalreadypresent], xn[newi]);
                    }
                    else
                    {
                        if (gmx_debug_at)
                        {
                            fprintf(debug,
                                    "Replacing %d '%s' with (old name '%s') %s\n",
                                    newi,
                                    ((*modifiedAtoms)->atomname[newi] && *(*modifiedAtoms)->atomname[newi])
                                            ? *(*modifiedAtoms)->atomname[newi]
                                            : "",
                                    patch->oname.empty() ? "" : patch->oname.c_str(),
                                    patch->nname.c_str());
                        }
                        (*modifiedAtoms)->atomname[newi] = put_symtab(symtab, patch->nname.c_str());
                        if (patch->bXSet)
                        {
                            copy_rvec(patch->newx, xn[newi]);
                        }
                    }
                    if (debug)
                    {
                        fprintf(debug,
                                " %s %g %g",
                                *(*modifiedAtoms)->atomname[newi],
                                (*modifiedAtoms)->atom[newi].m,
                                (*modifiedAtoms)->atom[newi].q);
                    }
                }
            }
            newi++;
            i += nalreadypresent;
        }
    }
    (*modifiedAtoms)->nr = newi;

    done_atom(pdba);
    *initialAtoms = *modifiedAtoms;

    *xptr = xn;

    return newi;
}

int add_h(t_atoms**                                   initialAtoms,
          t_atoms**                                   localAtoms,
          std::vector<gmx::RVec>*                     xptr,
          gmx::ArrayRef<const MoleculePatchDatabase>  globalPatches,
          t_symtab*                                   symtab,
          const int                                   nterpairs,
          gmx::ArrayRef<MoleculePatchDatabase* const> ntdb,
          gmx::ArrayRef<MoleculePatchDatabase* const> ctdb,
          gmx::ArrayRef<const int>                    rN,
          gmx::ArrayRef<const int>                    rC,
          const bool                                  bAllowMissing,
          gmx::ArrayRef<const int>                    cyclicBondsIndex)
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
        nnew = add_h_low(
                initialAtoms, localAtoms, xptr, globalPatches, symtab, nterpairs, ntdb, ctdb, rN, rC, FALSE, cyclicBondsIndex);
        niter++;
        if (niter > 100)
        {
            gmx_fatal(FARGS,
                      "More than 100 iterations of add_h. Maybe you are trying to replace an added "
                      "atom (this is not supported)?");
        }
    } while (nnew > nold);

    if (!bAllowMissing)
    {
        /* Call add_h_low once more, now only for the missing atoms check */
        add_h_low(initialAtoms, localAtoms, xptr, globalPatches, symtab, nterpairs, ntdb, ctdb, rN, rC, TRUE, cyclicBondsIndex);
    }

    return nnew;
}
