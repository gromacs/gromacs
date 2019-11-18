/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "gen_ad.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "hackblock.h"
#include "resall.h"

#define DIHEDRAL_WAS_SET_IN_RTP 0
static bool was_dihedral_set_in_rtp(const InteractionOfType& dih)
{
    // This is a bad way to check this, but I don't know how to make this better now.
    gmx::ArrayRef<const real> forceParam = dih.forceParam();
    return forceParam[MAXFORCEPARAM - 1] == DIHEDRAL_WAS_SET_IN_RTP;
}

typedef bool (*peq)(const InteractionOfType& p1, const InteractionOfType& p2);

static bool acomp(const InteractionOfType& a1, const InteractionOfType& a2)
{
    int ac;

    if ((ac = (a1.aj() - a2.aj())) != 0)
    {
        return ac < 0;
    }
    else if ((ac = (a1.ai() - a2.ai())) != 0)
    {
        return ac < 0;
    }
    else
    {
        return (a1.ak() < a2.ak());
    }
}

static bool pcomp(const InteractionOfType& a1, const InteractionOfType& a2)
{
    int pc;

    if ((pc = (a1.ai() - a2.ai())) != 0)
    {
        return pc < 0;
    }
    else
    {
        return (a1.aj() < a2.aj());
    }
}

static bool dcomp(const InteractionOfType& d1, const InteractionOfType& d2)
{
    int dc;

    /* First sort by J & K (the two central) atoms */
    if ((dc = (d1.aj() - d2.aj())) != 0)
    {
        return dc < 0;
    }
    else if ((dc = (d1.ak() - d2.ak())) != 0)
    {
        return dc < 0;
    }
    /* Then make sure to put rtp dihedrals before generated ones */
    else if (was_dihedral_set_in_rtp(d1) && !was_dihedral_set_in_rtp(d2))
    {
        return true;
    }
    else if (!was_dihedral_set_in_rtp(d1) && was_dihedral_set_in_rtp(d2))
    {
        return false;
    }
    /* Then sort by I and J (two outer) atoms */
    else if ((dc = (d1.ai() - d2.ai())) != 0)
    {
        return dc < 0;
    }
    else if ((dc = (d1.al() - d2.al())) != 0)
    {
        return dc < 0;
    }
    else
    {
        // AMBER force fields with type 9 dihedrals can reach here, where we sort on
        // the contents of the string that names the macro for the parameters.
        return std::lexicographical_compare(
                d1.interactionTypeName().begin(), d1.interactionTypeName().end(),
                d2.interactionTypeName().begin(), d2.interactionTypeName().end());
    }
}


static bool is_dihedral_on_same_bond(const InteractionOfType& p1, const InteractionOfType& p2)
{
    return ((p1.aj() == p2.aj()) && (p1.ak() == p2.ak()))
           || ((p1.aj() == p2.ak()) && (p1.ak() == p2.aj()));
}


static bool preq(const InteractionOfType& p1, const InteractionOfType& p2)
{
    return (p1.ai() == p2.ai()) && (p1.aj() == p2.aj());
}

static void rm2par(std::vector<InteractionOfType>* p, peq eq)
{
    if (p->empty())
    {
        return;
    }

    for (auto param = p->begin() + 1; param != p->end();)
    {
        auto prev = param - 1;
        if (eq(*param, *prev))
        {
            param = p->erase(param);
        }
        else
        {
            ++param;
        }
    }
}

static void cppar(gmx::ArrayRef<const InteractionOfType> types, gmx::ArrayRef<InteractionsOfType> plist, int ftype)
{
    /* Keep old stuff */
    for (const auto& type : types)
    {
        plist[ftype].interactionTypes.push_back(type);
    }
}

static bool idcomp(const InteractionOfType& a, const InteractionOfType& b)
{
    int d;

    if ((d = (a.ai() - b.ai())) != 0)
    {
        return d < 0;
    }
    else if ((d = (a.al() - b.al())) != 0)
    {
        return d < 0;
    }
    else if ((d = (a.aj() - b.aj())) != 0)
    {
        return d < 0;
    }
    else
    {
        return (a.ak() < b.ak());
    }
}

static void sort_id(gmx::ArrayRef<InteractionOfType> ps)
{
    if (ps.size() > 1)
    {
        for (auto& parm : ps)
        {
            parm.sortAtomIds();
        }
        std::sort(ps.begin(), ps.end(), idcomp);
    }
}

static int n_hydro(gmx::ArrayRef<const int> a, char*** atomname)
{
    int nh = 0;

    for (auto atom = a.begin(); atom < a.end(); atom += 3)
    {
        const char* aname = *atomname[*atom];
        char        c0    = toupper(aname[0]);
        if (c0 == 'H')
        {
            nh++;
        }
        else if ((static_cast<int>(strlen(aname)) > 1) && (c0 >= '0') && (c0 <= '9'))
        {
            char c1 = toupper(aname[1]);
            if (c1 == 'H')
            {
                nh++;
            }
        }
    }
    return nh;
}

/* Clean up the dihedrals (both generated and read from the .rtp
 * file). */
static std::vector<InteractionOfType> clean_dih(gmx::ArrayRef<const InteractionOfType> dih,
                                                gmx::ArrayRef<const InteractionOfType> improper,
                                                t_atoms*                               atoms,
                                                bool bKeepAllGeneratedDihedrals,
                                                bool bRemoveDihedralIfWithImproper)
{
    /* Construct the list of the indices of the dihedrals
     * (i.e. generated or read) that might be kept. */
    std::vector<std::pair<InteractionOfType, int>> newDihedrals;
    if (bKeepAllGeneratedDihedrals)
    {
        fprintf(stderr, "Keeping all generated dihedrals\n");
        int i = 0;
        for (const auto& dihedral : dih)
        {
            newDihedrals.emplace_back(std::pair<InteractionOfType, int>(dihedral, i++));
        }
    }
    else
    {
        /* Check if generated dihedral i should be removed. The
         * dihedrals have been sorted by dcomp() above, so all those
         * on the same two central atoms are together, with those from
         * the .rtp file preceding those that were automatically
         * generated. We remove the latter if the former exist. */
        int i = 0;
        for (auto dihedral = dih.begin(); dihedral != dih.end(); dihedral++)
        {
            /* Keep the dihedrals that were defined in the .rtp file,
             * and the dihedrals that were generated and different
             * from the last one (whether it was generated or not). */
            if (was_dihedral_set_in_rtp(*dihedral) || dihedral == dih.begin()
                || !is_dihedral_on_same_bond(*dihedral, *(dihedral - 1)))
            {
                newDihedrals.emplace_back(std::pair<InteractionOfType, int>(*dihedral, i++));
            }
        }
    }
    int k = 0;
    for (auto dihedral = newDihedrals.begin(); dihedral != newDihedrals.end();)
    {
        bool bWasSetInRTP = was_dihedral_set_in_rtp(dihedral->first);
        bool bKeep        = true;
        if (!bWasSetInRTP && bRemoveDihedralIfWithImproper)
        {
            /* Remove the dihedral if there is an improper on the same
             * bond. */
            for (auto imp = improper.begin(); imp != improper.end() && bKeep; ++imp)
            {
                bKeep = !is_dihedral_on_same_bond(dihedral->first, *imp);
            }
        }

        if (bKeep)
        {
            /* If we don't want all dihedrals, we want to select the
             * ones with the fewest hydrogens. Note that any generated
             * dihedrals on the same bond as an .rtp dihedral may have
             * been already pruned above in the construction of
             * index[]. However, their parameters are still present,
             * and l is looping over this dihedral and all of its
             * pruned siblings. */
            int bestl = dihedral->second;
            if (!bKeepAllGeneratedDihedrals && !bWasSetInRTP)
            {
                /* Minimum number of hydrogens for i and l atoms */
                int minh = 2;
                int next = dihedral->second + 1;
                for (int l = dihedral->second;
                     (l < next && is_dihedral_on_same_bond(dihedral->first, dih[l])); l++)
                {
                    int nh = n_hydro(dih[l].atoms(), atoms->atomname);
                    if (nh < minh)
                    {
                        minh  = nh;
                        bestl = l;
                    }
                    if (0 == minh)
                    {
                        break;
                    }
                }
                dihedral->first = dih[bestl];
            }
            if (k == bestl)
            {
                ++dihedral;
            }
            k++;
        }
        else
        {
            dihedral = newDihedrals.erase(dihedral);
        }
    }
    std::vector<InteractionOfType> finalDihedrals;
    finalDihedrals.reserve(newDihedrals.size());
    for (const auto& param : newDihedrals)
    {
        finalDihedrals.emplace_back(param.first);
    }
    return finalDihedrals;
}

static std::vector<InteractionOfType> get_impropers(t_atoms*                             atoms,
                                                    gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
                                                    bool bAllowMissing)
{
    std::vector<InteractionOfType> improper;

    /* Add all the impropers from the residue database to the list. */
    int start = 0;
    if (!globalPatches.empty())
    {
        for (int i = 0; (i < atoms->nres); i++)
        {
            BondedInteractionList* impropers = &globalPatches[i].rb[ebtsIDIHS];
            for (const auto& bondeds : impropers->b)
            {
                bool             bStop = false;
                std::vector<int> ai;
                for (int k = 0; (k < 4) && !bStop; k++)
                {
                    int entry = search_atom(bondeds.a[k].c_str(), start, atoms, "improper", bAllowMissing);

                    if (entry != -1)
                    {
                        ai.emplace_back(entry);
                    }
                    else
                    {
                        bStop = true;
                    }
                }
                if (!bStop)
                {
                    /* Not broken out */
                    improper.emplace_back(InteractionOfType(ai, {}, bondeds.s));
                }
            }
            while ((start < atoms->nr) && (atoms->atom[start].resind == i))
            {
                start++;
            }
        }
    }

    return improper;
}

static int nb_dist(t_nextnb* nnb, int ai, int aj)
{
    int  nre, nrx, NRE;
    int* nrexcl;
    int* a;

    if (ai == aj)
    {
        return 0;
    }

    NRE    = -1;
    nrexcl = nnb->nrexcl[ai];
    for (nre = 1; (nre < nnb->nrex); nre++)
    {
        a = nnb->a[ai][nre];
        for (nrx = 0; (nrx < nrexcl[nre]); nrx++)
        {
            if ((aj == a[nrx]) && (NRE == -1))
            {
                NRE = nre;
            }
        }
    }
    return NRE;
}

static bool is_hydro(t_atoms* atoms, int ai)
{
    return ((*(atoms->atomname[ai]))[0] == 'H');
}

static void get_atomnames_min(int                        n,
                              gmx::ArrayRef<std::string> anm,
                              int                        resind,
                              t_atoms*                   atoms,
                              gmx::ArrayRef<const int>   a)
{
    /* Assume ascending residue numbering */
    for (int m = 0; m < n; m++)
    {
        if (atoms->atom[a[m]].resind < resind)
        {
            anm[m] = "-";
        }
        else if (atoms->atom[a[m]].resind > resind)
        {
            anm[m] = "+";
        }
        else
        {
            anm[m] = "";
        }
        anm[m].append(*(atoms->atomname[a[m]]));
    }
}

static void gen_excls(t_atoms* atoms, t_excls* excls, gmx::ArrayRef<MoleculePatchDatabase> globalPatches, bool bAllowMissing)
{
    int astart = 0;
    for (int a = 0; a < atoms->nr; a++)
    {
        int r = atoms->atom[a].resind;
        if (a == atoms->nr - 1 || atoms->atom[a + 1].resind != r)
        {
            BondedInteractionList* hbexcl = &globalPatches[r].rb[ebtsEXCLS];

            for (const auto& bondeds : hbexcl->b)
            {
                const char* anm = bondeds.a[0].c_str();
                int         i1  = search_atom(anm, astart, atoms, "exclusion", bAllowMissing);
                anm             = bondeds.a[1].c_str();
                int i2          = search_atom(anm, astart, atoms, "exclusion", bAllowMissing);
                if (i1 != -1 && i2 != -1)
                {
                    if (i1 > i2)
                    {
                        int itmp = i1;
                        i1       = i2;
                        i2       = itmp;
                    }
                    srenew(excls[i1].e, excls[i1].nr + 1);
                    excls[i1].e[excls[i1].nr] = i2;
                    excls[i1].nr++;
                }
            }

            astart = a + 1;
        }
    }

    for (int a = 0; a < atoms->nr; a++)
    {
        if (excls[a].nr > 1)
        {
            std::sort(excls[a].e, excls[a].e + excls[a].nr);
        }
    }
}

static void remove_excl(t_excls* excls, int remove)
{
    int i;

    for (i = remove + 1; i < excls->nr; i++)
    {
        excls->e[i - 1] = excls->e[i];
    }

    excls->nr--;
}

static void clean_excls(t_nextnb* nnb, int nrexcl, t_excls excls[])
{
    int      i, j, j1, k, k1, l, l1, e;
    t_excls* excl;

    if (nrexcl >= 1)
    {
        /* extract all i-j-k-l neighbours from nnb struct */
        for (i = 0; (i < nnb->nr); i++)
        {
            /* For all particles */
            excl = &excls[i];

            for (j = 0; (j < nnb->nrexcl[i][1]); j++)
            {
                /* For all first neighbours */
                j1 = nnb->a[i][1][j];

                for (e = 0; e < excl->nr; e++)
                {
                    if (excl->e[e] == j1)
                    {
                        remove_excl(excl, e);
                    }
                }

                if (nrexcl >= 2)
                {
                    for (k = 0; (k < nnb->nrexcl[j1][1]); k++)
                    {
                        /* For all first neighbours of j1 */
                        k1 = nnb->a[j1][1][k];

                        for (e = 0; e < excl->nr; e++)
                        {
                            if (excl->e[e] == k1)
                            {
                                remove_excl(excl, e);
                            }
                        }

                        if (nrexcl >= 3)
                        {
                            for (l = 0; (l < nnb->nrexcl[k1][1]); l++)
                            {
                                /* For all first neighbours of k1 */
                                l1 = nnb->a[k1][1][l];

                                for (e = 0; e < excl->nr; e++)
                                {
                                    if (excl->e[e] == l1)
                                    {
                                        remove_excl(excl, e);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/* Generate pairs, angles and dihedrals from .rtp settings */
void gen_pad(t_atoms*                               atoms,
             gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
             gmx::ArrayRef<InteractionsOfType>      plist,
             t_excls                                excls[],
             gmx::ArrayRef<MoleculePatchDatabase>   globalPatches,
             bool                                   bAllowMissing)
{
    t_nextnb nnb;
    init_nnb(&nnb, atoms->nr, 4);
    gen_nnb(&nnb, plist);
    print_nnb(&nnb, "NNB");

    /* These are the angles, dihedrals and pairs that we generate
     * from the bonds. The ones that are already there from the rtp file
     * will be retained.
     */
    std::vector<InteractionOfType> ang;
    std::vector<InteractionOfType> dih;
    std::vector<InteractionOfType> pai;
    std::vector<InteractionOfType> improper;

    std::array<std::string, 4> anm;

    if (!globalPatches.empty())
    {
        gen_excls(atoms, excls, globalPatches, bAllowMissing);
        /* mark all entries as not matched yet */
        for (int i = 0; i < atoms->nres; i++)
        {
            for (int j = 0; j < ebtsNR; j++)
            {
                for (auto& bondeds : globalPatches[i].rb[j].b)
                {
                    bondeds.match = false;
                }
            }
        }
    }

    /* Extract all i-j-k-l neighbours from nnb struct to generate all
     * angles and dihedrals. */
    for (int i = 0; (i < nnb.nr); i++)
    {
        /* For all particles */
        for (int j = 0; (j < nnb.nrexcl[i][1]); j++)
        {
            /* For all first neighbours */
            int j1 = nnb.a[i][1][j];
            for (int k = 0; (k < nnb.nrexcl[j1][1]); k++)
            {
                /* For all first neighbours of j1 */
                int k1 = nnb.a[j1][1][k];
                if (k1 != i)
                {
                    /* Generate every angle only once */
                    if (i < k1)
                    {
                        std::vector<int> atomNumbers = { i, j1, k1 };
                        std::string      name;
                        if (!globalPatches.empty())
                        {
                            int minres = atoms->atom[i].resind;
                            int maxres = minres;
                            for (int m = 1; m < 3; m++)
                            {
                                minres = std::min(minres, atoms->atom[atomNumbers[m]].resind);
                                maxres = std::max(maxres, atoms->atom[atomNumbers[m]].resind);
                            }
                            int res = 2 * minres - maxres;
                            do
                            {
                                res += maxres - minres;
                                get_atomnames_min(3, anm, res, atoms, atomNumbers);
                                BondedInteractionList* hbang = &globalPatches[res].rb[ebtsANGLES];
                                for (auto& bondeds : hbang->b)
                                {
                                    if (anm[1] == bondeds.aj())
                                    {
                                        bool bFound = false;
                                        for (int m = 0; m < 3; m += 2)
                                        {
                                            bFound = (bFound
                                                      || ((anm[m] == bondeds.ai())
                                                          && (anm[2 - m] == bondeds.ak())));
                                        }
                                        if (bFound)
                                        {
                                            name = bondeds.s;
                                            /* Mark that we found a match for this entry */
                                            bondeds.match = true;
                                        }
                                    }
                                }
                            } while (res < maxres);
                        }
                        ang.push_back(InteractionOfType(atomNumbers, {}, name));
                    }
                    /* Generate every dihedral, 1-4 exclusion and 1-4 interaction
                       only once */
                    if (j1 < k1)
                    {
                        for (int l = 0; (l < nnb.nrexcl[k1][1]); l++)
                        {
                            /* For all first neighbours of k1 */
                            int l1 = nnb.a[k1][1][l];
                            if ((l1 != i) && (l1 != j1))
                            {
                                std::vector<int> atomNumbers = { i, j1, k1, l1 };
                                std::string      name;
                                int              nFound = 0;
                                if (!globalPatches.empty())
                                {
                                    int minres = atoms->atom[i].resind;
                                    int maxres = minres;
                                    for (int m = 1; m < 4; m++)
                                    {
                                        minres = std::min(minres, atoms->atom[atomNumbers[m]].resind);
                                        maxres = std::max(maxres, atoms->atom[atomNumbers[m]].resind);
                                    }
                                    int res = 2 * minres - maxres;
                                    do
                                    {
                                        res += maxres - minres;
                                        get_atomnames_min(4, anm, res, atoms, atomNumbers);
                                        BondedInteractionList* hbdih = &globalPatches[res].rb[ebtsPDIHS];
                                        for (auto& bondeds : hbdih->b)
                                        {
                                            bool bFound = false;
                                            for (int m = 0; m < 2; m++)
                                            {
                                                bFound = (bFound
                                                          || ((anm[3 * m] == bondeds.ai())
                                                              && (anm[1 + m] == bondeds.aj())
                                                              && (anm[2 - m] == bondeds.ak())
                                                              && (anm[3 - 3 * m] == bondeds.al())));
                                            }
                                            if (bFound)
                                            {
                                                name = bondeds.s;
                                                /* Mark that we found a match for this entry */
                                                bondeds.match = true;

                                                /* Set the last parameter to be able to see
                                                   if the dihedral was in the rtp list.
                                                 */
                                                nFound++;
                                                dih.push_back(InteractionOfType(atomNumbers, {}, name));
                                                dih.back().setForceParameter(
                                                        MAXFORCEPARAM - 1, DIHEDRAL_WAS_SET_IN_RTP);
                                            }
                                        }
                                    } while (res < maxres);
                                }
                                if (nFound == 0)
                                {
                                    std::vector<int> atoms = { i, j1, k1, l1 };
                                    dih.push_back(InteractionOfType(atoms, {}, ""));
                                }

                                int nbd = nb_dist(&nnb, i, l1);
                                if (nbd == 3)
                                {
                                    int  i1    = std::min(i, l1);
                                    int  i2    = std::max(i, l1);
                                    bool bExcl = false;
                                    for (int m = 0; m < excls[i1].nr; m++)
                                    {
                                        bExcl = bExcl || excls[i1].e[m] == i2;
                                    }
                                    if (!bExcl)
                                    {
                                        if (rtpFFDB[0].bGenerateHH14Interactions
                                            || !(is_hydro(atoms, i1) && is_hydro(atoms, i2)))
                                        {
                                            std::vector<int> atoms = { i1, i2 };
                                            pai.push_back(InteractionOfType(atoms, {}, ""));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (!globalPatches.empty())
    {
        /* The above approach is great in that we double-check that e.g. an angle
         * really corresponds to three atoms connected by bonds, but this is not
         * generally true. Go through the angle and dihedral hackblocks to add
         * entries that we have not yet marked as matched when going through bonds.
         */
        for (int i = 0; i < atoms->nres; i++)
        {
            /* Add remaining angles from hackblock */
            BondedInteractionList* hbang = &globalPatches[i].rb[ebtsANGLES];
            for (auto& bondeds : hbang->b)
            {
                if (bondeds.match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                std::vector<int> atomNumbers;
                bool             bFound = true;
                for (int k = 0; k < 3 && bFound; k++)
                {
                    const char* p   = bondeds.a[k].c_str();
                    int         res = i;
                    if (p[0] == '-')
                    {
                        p++;
                        res--;
                    }
                    else if (p[0] == '+')
                    {
                        p++;
                        res++;
                    }
                    atomNumbers.emplace_back(search_res_atom(p, res, atoms, "angle", TRUE));
                    bFound = (atomNumbers.back() != -1);
                }

                if (bFound)
                {
                    bondeds.match = true;
                    /* Incrementing nang means we save this angle */
                    ang.push_back(InteractionOfType(atomNumbers, {}, bondeds.s));
                }
            }

            /* Add remaining dihedrals from hackblock */
            BondedInteractionList* hbdih = &globalPatches[i].rb[ebtsPDIHS];
            for (auto& bondeds : hbdih->b)
            {
                if (bondeds.match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                std::vector<int> atomNumbers;
                bool             bFound = true;
                for (int k = 0; k < 4 && bFound; k++)
                {
                    const char* p   = bondeds.a[k].c_str();
                    int         res = i;
                    if (p[0] == '-')
                    {
                        p++;
                        res--;
                    }
                    else if (p[0] == '+')
                    {
                        p++;
                        res++;
                    }
                    atomNumbers.emplace_back(search_res_atom(p, res, atoms, "dihedral", TRUE));
                    bFound = (atomNumbers.back() != -1);
                }

                if (bFound)
                {
                    bondeds.match = true;
                    /* Incrementing ndih means we save this dihedral */
                    dih.push_back(InteractionOfType(atomNumbers, {}, bondeds.s));
                }
            }
        }
    }

    /* Sort angles with respect to j-i-k (middle atom first) */
    if (ang.size() > 1)
    {
        std::sort(ang.begin(), ang.end(), acomp);
    }

    /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
    if (dih.size() > 1)
    {
        std::sort(dih.begin(), dih.end(), dcomp);
    }

    /* Sort the pairs */
    if (pai.size() > 1)
    {
        std::sort(pai.begin(), pai.end(), pcomp);
    }
    if (!pai.empty())
    {
        /* Remove doubles, could occur in 6-rings, such as phenyls,
           maybe one does not want this when fudgeQQ < 1.
         */
        fprintf(stderr, "Before cleaning: %zu pairs\n", pai.size());
        rm2par(&pai, preq);
    }

    /* Get the impropers from the database */
    improper = get_impropers(atoms, globalPatches, bAllowMissing);

    /* Sort the impropers */
    sort_id(improper);

    if (!dih.empty())
    {
        fprintf(stderr, "Before cleaning: %zu dihedrals\n", dih.size());
        dih = clean_dih(dih, improper, atoms, rtpFFDB[0].bKeepAllGeneratedDihedrals,
                        rtpFFDB[0].bRemoveDihedralIfWithImproper);
    }

    /* Now we have unique lists of angles and dihedrals
     * Copy them into the destination struct
     */
    cppar(ang, plist, F_ANGLES);
    cppar(dih, plist, F_PDIHS);
    cppar(improper, plist, F_IDIHS);
    cppar(pai, plist, F_LJ14);

    /* Remove all exclusions which are within nrexcl */
    clean_excls(&nnb, rtpFFDB[0].nrexcl, excls);
    done_nnb(&nnb);
}
