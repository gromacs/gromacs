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

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define DIHEDRAL_WAS_SET_IN_RTP 0
static bool was_dihedral_set_in_rtp(const t_param &dih)
{
    return dih.c[MAXFORCEPARAM-1] == DIHEDRAL_WAS_SET_IN_RTP;
}

typedef bool (*peq)(const t_param &p1, const t_param &p2);

static int acomp(const t_param &a1, const t_param &a2)
{
    int            ac;

    if ((ac = (a1.aj()-a2.aj())) != 0)
    {
        return ac;
    }
    else if ((ac = (a1.ai()-a2.ai())) != 0)
    {
        return ac;
    }
    else
    {
        return (a1.ak()-a2.ak());
    }
}

static int pcomp(const t_param &a1, const t_param &a2)
{
    int            pc;

    if ((pc = (a1.ai()-a2.ai())) != 0)
    {
        return pc;
    }
    else
    {
        return (a1.aj()-a2.aj());
    }
}

static int dcomp(const t_param &d1, const t_param &d2)
{
    int            dc;

    /* First sort by J & K (the two central) atoms */
    if ((dc = (d1.aj()-d2.aj())) != 0)
    {
        return dc;
    }
    else if ((dc = (d1.ak()-d2.ak())) != 0)
    {
        return dc;
    }
    /* Then make sure to put rtp dihedrals before generated ones */
    else if (was_dihedral_set_in_rtp(d1) &&
             !was_dihedral_set_in_rtp(d2))
    {
        return -1;
    }
    else if (!was_dihedral_set_in_rtp(d1) &&
             was_dihedral_set_in_rtp(d2))
    {
        return 1;
    }
    /* Then sort by I and J (two outer) atoms */
    else if ((dc = (d1.ai()-d2.ai())) != 0)
    {
        return dc;
    }
    else if ((dc = (d1.al()-d2.al())) != 0)
    {
        return dc;
    }
    else
    {
        // AMBER force fields with type 9 dihedrals can reach here, where we sort on
        // the contents of the string that names the macro for the parameters.
        return d1.s == d2.s;
    }
}


static bool is_dihedral_on_same_bond(const t_param &p1, const t_param &p2)
{
    return ((p1.aj() == p2.aj()) && (p1.ak() == p2.ak())) ||
           ((p1.aj() == p2.ak()) && (p1.ak() == p2.aj()));
}


static bool preq(const t_param &p1, const t_param &p2)
{
    return (p1.ai() == p2.ai()) && (p1.aj() == p2.aj());
}

static void rm2par(std::vector<t_param> *p, peq eq)
{
    if (p->empty())
    {
        return;
    }

    for (auto it = p->begin(); it != p->end(); )
    {
        auto next = it;
        if (eq(*it, *(++next)))
        {
            it = p->erase(it);
        }
        else
        {
            ++it;
        }
    }
}

static void cppar(gmx::ArrayRef<t_param> p, gmx::ArrayRef<t_params> plist, int ftype)
{
    t_params *ps   = &plist[ftype];
    int       nral = NRAL(ftype);
    int       nrfp = NRFP(ftype);

    /* Keep old stuff */
    for (int i = 0; (i < p.size()); i++)
    {
        ps->param.push_back(t_param());
        for (int j = 0; (j < nral); j++)
        {
            ps->param.back().a[j] = p[i].a[j];
        }
        for (int j = 0; (j < nrfp); j++)
        {
            ps->param.back().c[j] = p[i].c[j];
        }
        ps->param.back().s = p[i].s;
    }
}

static void set_p(t_param *p, const int ai[4], const real *c, const char *s)
{
    for (int j = 0; (j < 4); j++)
    {
        p->a[j] = ai[j];
    }
    for (int j = 0; (j < MAXFORCEPARAM); j++)
    {
        if (c)
        {
            p->c[j] = c[j];
        }
        else
        {
            p->c[j] = NOTSET;
        }
    }

    set_p_string(p, s);
}

static int idcomp(const t_param &a, const t_param &b)
{
    int            d;

    if ((d = (a.a[0]-b.a[0])) != 0)
    {
        return d;
    }
    else if ((d = (a.a[3]-b.a[3])) != 0)
    {
        return d;
    }
    else if ((d = (a.a[1]-b.a[1])) != 0)
    {
        return d;
    }
    else
    {
        return (a.a[2]-b.a[2]);
    }
}

static void sort_id(gmx::ArrayRef<t_param> ps)
{
    /* First swap order of atoms around if necessary */
    for (auto it = ps.begin(); it != ps.end(); it++)
    {
        int tmp;
        if (it->a[3] < it->a[0])
        {
            tmp = it->a[3]; it->a[3] = it->a[0]; it->a[0] = tmp;
            tmp = it->a[2]; it->a[2] = it->a[1]; it->a[1] = tmp;
        }
    }
    std::sort(ps.begin(), ps.end(), idcomp);
}

static int n_hydro(gmx::ArrayRef<const int> a, char ***atomname)
{
    int  i, nh = 0;
    char c0, c1, *aname;

    for (i = 0; (i < 4); i += 3)
    {
        aname = *atomname[a[i]];
        c0    = toupper(aname[0]);
        if (c0 == 'H')
        {
            nh++;
        }
        else if ((static_cast<int>(strlen(aname)) > 1) && (c0 >= '0') && (c0 <= '9'))
        {
            c1 = toupper(aname[1]);
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
static void clean_dih(std::vector<t_param> *dih,
                      gmx::ArrayRef<t_param> improper,
                      const t_atoms &atoms, bool bKeepAllGeneratedDihedrals,
                      bool bRemoveDihedralIfWithImproper)
{
    /* Construct the list of the indices of the dihedrals
     * (i.e. generated or read) that might be kept. */
    std::vector<int> index;
    if (bKeepAllGeneratedDihedrals)
    {
        fprintf(stderr, "Keeping all generated dihedrals\n");
        for (int i = 0; i < gmx::index(dih->size()); i++)
        {
            index.push_back(i);
        }
    }
    else
    {
        /* Check if generated dihedral i should be removed. The
         * dihedrals have been sorted by dcomp() above, so all those
         * on the same two central atoms are together, with those from
         * the .rtp file preceding those that were automatically
         * generated. We remove the latter if the former exist. */
        for (int i = 0; i < gmx::index(dih->size()); i++)
        {
            /* Keep the dihedrals that were defined in the .rtp file,
             * and the dihedrals that were generated and different
             * from the last one (whether it was generated or not). */
            if (was_dihedral_set_in_rtp(dih->at(i)) ||
                0 == i ||
                !is_dihedral_on_same_bond(dih->at(i), dih->at(i-1)))
            {
                index.push_back(i);
            }
        }
    }

    std::vector<t_param> newDih;
    int                  k = 0;
    for (int i = 0; i < gmx::index(index.size()); i++)
    {
        bool bWasSetInRTP = was_dihedral_set_in_rtp(dih->at(index[i]));
        bool bKeep        = TRUE;
        if (!bWasSetInRTP && bRemoveDihedralIfWithImproper)
        {
            /* Remove the dihedral if there is an improper on the same
             * bond. */
            for (int j = 0; j < gmx::index(improper.size()) && bKeep; j++)
            {
                bKeep = !is_dihedral_on_same_bond(dih->at(index[i]), improper[j]);
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
            int bestl = index[i];
            if (!bKeepAllGeneratedDihedrals && !bWasSetInRTP)
            {
                /* Minimum number of hydrogens for i and l atoms */
                int minh = 2;
                for (int l = index[i];
                     (l < index[i+1] &&
                      is_dihedral_on_same_bond(dih->at(index[i]), dih->at(l)));
                     l++)
                {
                    int nh = n_hydro(dih->at(l).a, atoms.atomname);
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
            }
            if (k != bestl)
            {
                newDih.push_back(dih->at(bestl));
            }
            k++;
        }
    }
    *dih = newDih;
}

static void get_impropers(const t_atoms             &atoms,
                          gmx::ArrayRef<t_hackblock> hb,
                          std::vector<t_param>      *improper,
                          bool                       bAllowMissing)
{
    int           start;
    int           ai[MAXATOMLIST];
    bool          bStop;

    /* Add all the impropers from the residue database to the list. */
    start     = 0;
    if (!hb.empty())
    {
        for (int i = 0; (i < atoms.nres); i++)
        {
            t_rbondeds *impropers = &hb[i].rb[ebtsIDIHS];
            for (int j = 0; (j < impropers->nb()); j++)
            {
                bStop = FALSE;
                for (int k = 0; (k < 4) && !bStop; k++)
                {
                    ai[k] = search_atom(impropers->b[j].a[k].c_str(), start,
                                        atoms,
                                        "improper", bAllowMissing);
                    if (ai[k] == -1)
                    {
                        bStop = TRUE;
                    }
                }
                if (!bStop)
                {
                    improper->push_back(t_param());
                    /* Not broken out */
                    set_p(&improper->back(), ai, nullptr, impropers->b[j].s.c_str());
                }
            }
            while ((start < atoms.nr) && (atoms.atom[start].resind == i))
            {
                start++;
            }
        }
    }
}

static int nb_dist(t_nextnb *nnb, int ai, int aj)
{
    int  nre, nrx, NRE;
    int *nrexcl;
    int *a;

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

static bool is_hydro(const t_atoms &atoms, int ai)
{
    return ((*(atoms.atomname[ai]))[0] == 'H');
}

static void get_atomnames_min(int n, gmx::ArrayRef<std::string> anm,
                              int resind, const t_atoms &atoms, gmx::ArrayRef<const int> a)
{
    /* Assume ascending residue numbering */
    for (int m = 0; m < n; m++)
    {
        if (atoms.atom[a[m]].resind < resind)
        {
            anm[m] = "-";
        }
        else if (atoms.atom[a[m]].resind > resind)
        {
            anm[m] = "+";
        }
        else
        {
            anm[m] = "";
        }
        anm[m].append(*(atoms.atomname[a[m]]));
    }
}

static void gen_excls(const t_atoms &atoms, t_excls *excls, gmx::ArrayRef<t_hackblock> hb,
                      bool bAllowMissing)
{
    int         itmp;
    t_rbondeds *hbexcl;

    int         astart = 0;
    for (int a = 0; a < atoms.nr; a++)
    {
        int r = atoms.atom[a].resind;
        if (a == atoms.nr-1 || atoms.atom[a+1].resind != r)
        {
            hbexcl = &hb[r].rb[ebtsEXCLS];

            for (int e = 0; e < hbexcl->nb(); e++)
            {
                std::string anm = hbexcl->b[e].a[0].c_str();
                int         i1  = search_atom(anm.c_str(), astart, atoms,
                                              "exclusion", bAllowMissing);
                anm = hbexcl->b[e].a[1];
                int i2  = search_atom(anm.c_str(), astart, atoms,
                                      "exclusion", bAllowMissing);
                if (i1 != -1 && i2 != -1)
                {
                    if (i1 > i2)
                    {
                        itmp = i1;
                        i1   = i2;
                        i2   = itmp;
                    }
                    srenew(excls[i1].e, excls[i1].nr+1);
                    excls[i1].e[excls[i1].nr] = i2;
                    excls[i1].nr++;
                }
            }

            astart = a+1;
        }
    }

    for (int a = 0; a < atoms.nr; a++)
    {
        if (excls[a].nr > 1)
        {
            std::sort(excls[a].e, excls[a].e+excls[a].nr);
        }
    }
}

static void remove_excl(t_excls *excls, int remove)
{
    int i;

    for (i = remove+1; i < excls->nr; i++)
    {
        excls->e[i-1] = excls->e[i];
    }

    excls->nr--;
}

void clean_excls(t_nextnb *nnb, int nrexcl, t_excls excls[])
{
    int      i, j, j1, k, k1, l, l1, e;
    t_excls *excl;

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

void generate_excls(t_nextnb *nnb, int nrexcl, t_excls excls[])
{
    int      i, j, n, N;
    t_excls *excl;

    for (N = 1; (N < std::min(nrexcl, nnb->nrex)); N++)
    {
        /* extract all i-j-k-l neighbours from nnb struct */
        for (i = 0; (i < nnb->nr); i++)
        {
            /* For all particles */
            excl      = &excls[i];
            n         = excl->nr;
            excl->nr += nnb->nrexcl[i][N];
            srenew(excl->e, excl->nr);
            for (j = 0; (j < nnb->nrexcl[i][N]); j++)
            {
                /* For all first neighbours */
                if (nnb->a[i][N][j] != i)
                {
                    excl->e[n++] = nnb->a[i][N][j];
                }
            }
        }
    }
}

/* Generate pairs, angles and dihedrals from .rtp settings */
void gen_pad(t_nextnb *nnb, const t_atoms &atoms, gmx::ArrayRef<t_restp> rtp,
             gmx::ArrayRef<t_params> plist, t_excls excls[], gmx::ArrayRef<t_hackblock> hb,
             bool bAllowMissing)
{
    t_rbondeds *hbang, *hbdih;
    int         res, minres, maxres;
    int         i1, i2;
    int         nbd;
    int         nFound;
    bool        bFound, bExcl;

    /* These are the angles, dihedrals and pairs that we generate
     * from the bonds. The ones that are already there from the rtp file
     * will be retained.
     */
    std::vector<t_param>     ang;
    std::vector<t_param>     dih;
    std::vector<t_param>     pai;
    std::vector<t_param>     improper;

    std::vector<std::string> anm(4);

    if (!hb.empty())
    {
        gen_excls(atoms, excls, hb, bAllowMissing);
        /* mark all entries as not matched yet */
        for (int i = 0; i < atoms.nres; i++)
        {
            for (int j = 0; j < ebtsNR; j++)
            {
                for (int k = 0; k < hb[i].rb[j].nb(); k++)
                {
                    hb[i].rb[j].b[k].match = FALSE;
                }
            }
        }
    }

    /* Extract all i-j-k-l neighbours from nnb struct to generate all
     * angles and dihedrals. */
    for (int i = 0; (i < nnb->nr); i++)
    {
        /* For all particles */
        for (int j = 0; (j < nnb->nrexcl[i][1]); j++)
        {
            /* For all first neighbours */
            int j1 = nnb->a[i][1][j];
            for (int k = 0; (k < nnb->nrexcl[j1][1]); k++)
            {
                /* For all first neighbours of j1 */
                int k1 = nnb->a[j1][1][k];
                if (k1 != i)
                {
                    /* Generate every angle only once */
                    if (i < k1)
                    {
                        ang.push_back(t_param());
                        ang.back().ai() = i;
                        ang.back().aj() = j1;
                        ang.back().ak() = k1;
                        set_p_string(&(ang.back()), "");
                        if (!hb.empty())
                        {
                            minres = atoms.atom[ang.back().a[0]].resind;
                            maxres = minres;
                            for (int m = 1; m < 3; m++)
                            {
                                minres = std::min(minres, atoms.atom[ang.back().a[m]].resind);
                                maxres = std::max(maxres, atoms.atom[ang.back().a[m]].resind);
                            }
                            res = 2*minres-maxres;
                            do
                            {
                                res += maxres-minres;
                                get_atomnames_min(3, anm, res, atoms, ang.back().a);
                                hbang = &hb[res].rb[ebtsANGLES];
                                for (int l = 0; (l < hbang->nb()); l++)
                                {
                                    if (anm[1] == hbang->b[l].aj())
                                    {
                                        bFound = FALSE;
                                        for (int m = 0; m < 3; m += 2)
                                        {
                                            bFound = (bFound ||
                                                      ((anm[m] == hbang->b[l].ai()) &&
                                                       (anm[2-m] == hbang->b[l].ak())));
                                        }
                                        if (bFound)
                                        {
                                            set_p_string(&(ang.back()), hbang->b[l].s.c_str());
                                            /* Mark that we found a match for this entry */
                                            hbang->b[l].match = TRUE;
                                        }
                                    }
                                }
                            }
                            while (res < maxres);
                        }
                    }
                    /* Generate every dihedral, 1-4 exclusion and 1-4 interaction
                       only once */
                    if (j1 < k1)
                    {
                        for (int l = 0; (l < nnb->nrexcl[k1][1]); l++)
                        {
                            /* For all first neighbours of k1 */
                            int l1 = nnb->a[k1][1][l];
                            if ((l1 != i) && (l1 != j1))
                            {
                                dih.push_back(t_param());
                                dih.back().ai() = i;
                                dih.back().aj() = j1;
                                dih.back().ak() = k1;
                                dih.back().al() = l1;
                                set_p_string(&(dih.back()), "");
                                nFound = 0;
                                if (!hb.empty())
                                {
                                    minres = atoms.atom[dih.back().a[0]].resind;
                                    maxres = minres;
                                    for (int m = 1; m < 4; m++)
                                    {
                                        minres = std::min(minres, atoms.atom[dih.back().a[m]].resind);
                                        maxres = std::max(maxres, atoms.atom[dih.back().a[m]].resind);
                                    }
                                    res = 2*minres-maxres;
                                    do
                                    {
                                        res += maxres-minres;
                                        get_atomnames_min(4, anm, res, atoms, dih.back().a);
                                        hbdih = &hb[res].rb[ebtsPDIHS];
                                        for (int n = 0; (n < hbdih->nb()); n++)
                                        {
                                            bFound = FALSE;
                                            for (int m = 0; m < 2; m++)
                                            {
                                                bFound = (bFound ||
                                                          ((anm[3*m] == hbdih->b[n].ai()) &&
                                                           (anm[1+m] == hbdih->b[n].aj()) &&
                                                           (anm[2-m] == hbdih->b[n].ak()) &&
                                                           (anm[3-3*m] == hbdih->b[n].al())));
                                            }
                                            if (bFound)
                                            {
                                                set_p_string(&dih.back(), hbdih->b[n].s);
                                                /* Mark that we found a match for this entry */
                                                hbdih->b[n].match = TRUE;

                                                /* Set the last parameter to be able to see
                                                   if the dihedral was in the rtp list.
                                                 */
                                                dih.back().c[MAXFORCEPARAM-1] = DIHEDRAL_WAS_SET_IN_RTP;
                                                nFound++;
                                                /* Set the next direct in case the rtp contains
                                                   multiple entries for this dihedral.
                                                 */
                                                dih.push_back(t_param());
                                                dih.back().ai() = i;
                                                dih.back().aj() = j1;
                                                dih.back().ak() = k1;
                                                dih.back().al() = l1;
                                            }
                                        }
                                    }
                                    while (res < maxres);
                                }
                                if (nFound == 0)
                                {
                                    dih.push_back(t_param());
                                    dih.back().ai() = i;
                                    dih.back().aj() = j1;
                                    dih.back().ak() = k1;
                                    dih.back().al() = l1;
                                    set_p_string(&(dih.back()), "");
                                }

                                nbd = nb_dist(nnb, i, l1);
                                if (nbd == 3)
                                {
                                    i1    = std::min(i, l1);
                                    i2    = std::max(i, l1);
                                    bExcl = FALSE;
                                    for (int m = 0; m < excls[i1].nr; m++)
                                    {
                                        bExcl = bExcl || excls[i1].e[m] == i2;
                                    }
                                    if (!bExcl)
                                    {
                                        if (rtp[0].bGenerateHH14Interactions ||
                                            !(is_hydro(atoms, i1) && is_hydro(atoms, i2)))
                                        {
                                            pai.push_back(t_param());
                                            pai.back().ai() = i1;
                                            pai.back().aj() = i2;
                                            set_p_string(&(pai.back()), "");
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

    if (!hb.empty())
    {
        /* The above approach is great in that we double-check that e.g. an angle
         * really corresponds to three atoms connected by bonds, but this is not
         * generally true. Go through the angle and dihedral hackblocks to add
         * entries that we have not yet marked as matched when going through bonds.
         */
        for (int i = 0; i < atoms.nres; i++)
        {
            /* Add remaining angles from hackblock */
            hbang = &hb[i].rb[ebtsANGLES];
            for (int j = 0; j < hbang->nb(); j++)
            {
                if (hbang->b[j].match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                t_param tmp;
                bFound = TRUE;
                for (int k = 0; k < 3 && bFound; k++)
                {
                    std::string p   = hbang->b[j].a[k];
                    res = i;
                    if (p[0] == '-')
                    {
                        p.erase(p.begin());
                        res--;
                    }
                    else if (p[0] == '+')
                    {
                        p.erase(p.begin());
                        res++;
                    }
                    tmp.a[k]       = search_res_atom(p, res, atoms, "angle", TRUE);
                    bFound         = (tmp.a[k] != -1);
                }

                if (bFound)
                {
                    set_p_string(&(tmp), hbang->b[j].s);
                    hbang->b[j].match = TRUE;
                    ang.push_back(tmp);
                }
            }

            /* Add remaining dihedrals from hackblock */
            hbdih = &hb[i].rb[ebtsPDIHS];
            for (int j = 0; j < hbdih->nb(); j++)
            {
                if (hbdih->b[j].match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                bFound = TRUE;
                t_param tmp;
                for (int k = 0; k < 4 && bFound; k++)
                {
                    std::string p   = hbdih->b[j].a[k];
                    res = i;
                    if (p[0] == '-')
                    {
                        p.erase(p.begin());
                        res--;
                    }
                    else if (p[0] == '+')
                    {
                        p.erase(p.begin());
                        res++;
                    }
                    tmp.a[k]       = search_res_atom(p, res, atoms, "dihedral", TRUE);
                    bFound         = (tmp.a[k] != -1);
                }

                if (bFound)
                {
                    set_p_string(&(tmp), hbdih->b[j].s);
                    hbdih->b[j].match = TRUE;
                    /* Incrementing ndih means we save this dihedral */
                    dih.push_back(tmp);
                }
            }
        }
    }

    /* Sort angles with respect to j-i-k (middle atom first) */
    std::sort(ang.begin(), ang.end(), acomp);

    /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
    std::sort(dih.begin(), dih.end(), dcomp);

    /* Sort the pairs */
    std::sort(pai.begin(), pai.end(), pcomp);
    if (!pai.empty())
    {
        /* Remove doubles, could occur in 6-rings, such as phenyls,
           maybe one does not want this when fudgeQQ < 1.
         */
        int size = pai.size();
        fprintf(stderr, "Before cleaning: %d pairs\n", size);
        rm2par(&pai, preq);
    }

    /* Get the impropers from the database */
    get_impropers(atoms, hb, &improper, bAllowMissing);

    /* Sort the impropers */
    sort_id(improper);

    if (dih.size() > 0)
    {
        int size = dih.size();
        fprintf(stderr, "Before cleaning: %d dihedrals\n", size);
        clean_dih(&dih, improper, atoms,
                  rtp[0].bKeepAllGeneratedDihedrals,
                  rtp[0].bRemoveDihedralIfWithImproper);
    }

    /* Now we have unique lists of angles and dihedrals
     * Copy them into the destination struct
     */
    cppar(ang, plist, F_ANGLES);
    cppar(dih, plist, F_PDIHS);
    cppar(improper, plist, F_IDIHS);
    cppar(pai, plist, F_LJ14);

    /* Remove all exclusions which are within nrexcl */
    clean_excls(nnb, rtp[0].nrexcl, excls);
}
