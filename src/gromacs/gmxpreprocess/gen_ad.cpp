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
static bool was_dihedral_set_in_rtp(const t_param *dih)
{
    return dih->c[MAXFORCEPARAM-1] == DIHEDRAL_WAS_SET_IN_RTP;
}

typedef bool (*peq)(t_param *p1, t_param *p2);

static int acomp(const void *a1, const void *a2)
{
    const t_param *p1, *p2;
    int            ac;

    p1 = static_cast<const t_param *>(a1);
    p2 = static_cast<const t_param *>(a2);
    if ((ac = (p1->aj()-p2->aj())) != 0)
    {
        return ac;
    }
    else if ((ac = (p1->ai()-p2->ai())) != 0)
    {
        return ac;
    }
    else
    {
        return (p1->ak()-p2->ak());
    }
}

static int pcomp(const void *a1, const void *a2)
{
    const t_param *p1, *p2;
    int            pc;

    p1 = static_cast<const t_param *>(a1);
    p2 = static_cast<const t_param *>(a2);
    if ((pc = (p1->ai()-p2->ai())) != 0)
    {
        return pc;
    }
    else
    {
        return (p1->aj()-p2->aj());
    }
}

static int dcomp(const void *d1, const void *d2)
{
    const t_param *p1, *p2;
    int            dc;

    p1 = static_cast<const t_param *>(d1);
    p2 = static_cast<const t_param *>(d2);
    /* First sort by J & K (the two central) atoms */
    if ((dc = (p1->aj()-p2->aj())) != 0)
    {
        return dc;
    }
    else if ((dc = (p1->ak()-p2->ak())) != 0)
    {
        return dc;
    }
    /* Then make sure to put rtp dihedrals before generated ones */
    else if (was_dihedral_set_in_rtp(p1) &&
             !was_dihedral_set_in_rtp(p2))
    {
        return -1;
    }
    else if (!was_dihedral_set_in_rtp(p1) &&
             was_dihedral_set_in_rtp(p2))
    {
        return 1;
    }
    /* Then sort by I and J (two outer) atoms */
    else if ((dc = (p1->ai()-p2->ai())) != 0)
    {
        return dc;
    }
    else if ((dc = (p1->al()-p2->al())) != 0)
    {
        return dc;
    }
    else
    {
        // AMBER force fields with type 9 dihedrals can reach here, where we sort on
        // the contents of the string that names the macro for the parameters.
        return strcmp(p1->s, p2->s);
    }
}


static bool is_dihedral_on_same_bond(t_param *p1, t_param *p2)
{
    return ((p1->aj() == p2->aj()) && (p1->ak() == p2->ak())) ||
           ((p1->aj() == p2->ak()) && (p1->ak() == p2->aj()));
}


static bool preq(t_param *p1, t_param *p2)
{
    return (p1->ai() == p2->ai()) && (p1->aj() == p2->aj());
}

static void rm2par(t_param p[], int *np, peq eq)
{
    int *index, nind;
    int  i, j;

    if ((*np) == 0)
    {
        return;
    }

    snew(index, *np);
    nind          = 0;
    index[nind++] = 0;
    for (i = 1; (i < (*np)); i++)
    {
        if (!eq(&p[i], &p[i-1]))
        {
            index[nind++] = i;
        }
    }
    /* Index now holds pointers to all the non-equal params,
     * this only works when p is sorted of course
     */
    for (i = 0; (i < nind); i++)
    {
        for (j = 0; (j < MAXATOMLIST); j++)
        {
            p[i].a[j] = p[index[i]].a[j];
        }
        for (j = 0; (j < MAXFORCEPARAM); j++)
        {
            p[i].c[j] = p[index[i]].c[j];
        }
        if (p[index[i]].a[0] == p[index[i]].a[1])
        {
            strcpy(p[i].s, "");
        }
        else if (index[i] > i)
        {
            /* Copy the string only if it comes from somewhere else
             * otherwise we will end up copying a random (newly freed) pointer.
             * Since the index is sorted we only have to test for index[i] > i.
             */
            strcpy(p[i].s, p[index[i]].s);
        }
    }
    (*np) = nind;

    sfree(index);
}

static void cppar(t_param p[], int np, t_params plist[], int ftype)
{
    int       i, j, nral, nrfp;
    t_params *ps;

    ps   = &plist[ftype];
    nral = NRAL(ftype);
    nrfp = NRFP(ftype);

    /* Keep old stuff */
    pr_alloc(np, ps);
    for (i = 0; (i < np); i++)
    {
        for (j = 0; (j < nral); j++)
        {
            ps->param[ps->nr].a[j] = p[i].a[j];
        }
        for (j = 0; (j < nrfp); j++)
        {
            ps->param[ps->nr].c[j] = p[i].c[j];
        }
        for (j = 0; (j < MAXSLEN); j++)
        {
            ps->param[ps->nr].s[j] = p[i].s[j];
        }
        ps->nr++;
    }
}

static void cpparam(t_param *dest, t_param *src)
{
    int j;

    for (j = 0; (j < MAXATOMLIST); j++)
    {
        dest->a[j] = src->a[j];
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        dest->c[j] = src->c[j];
    }
    for (j = 0; (j < MAXSLEN); j++)
    {
        dest->s[j] = src->s[j];
    }
}

static void set_p(t_param *p, const int ai[4], const real *c, const char *s)
{
    int j;

    for (j = 0; (j < 4); j++)
    {
        p->a[j] = ai[j];
    }
    for (j = 0; (j < MAXFORCEPARAM); j++)
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

static int idcomp(const void *a, const void *b)
{
    const t_param *pa, *pb;
    int            d;

    pa = static_cast<const t_param *>(a);
    pb = static_cast<const t_param *>(b);
    if ((d = (pa->a[0]-pb->a[0])) != 0)
    {
        return d;
    }
    else if ((d = (pa->a[3]-pb->a[3])) != 0)
    {
        return d;
    }
    else if ((d = (pa->a[1]-pb->a[1])) != 0)
    {
        return d;
    }
    else
    {
        return (pa->a[2]-pb->a[2]);
    }
}

static void sort_id(int nr, t_param ps[])
{
    int i, tmp;

    /* First swap order of atoms around if necessary */
    for (i = 0; (i < nr); i++)
    {
        if (ps[i].a[3] < ps[i].a[0])
        {
            tmp = ps[i].a[3]; ps[i].a[3] = ps[i].a[0]; ps[i].a[0] = tmp;
            tmp = ps[i].a[2]; ps[i].a[2] = ps[i].a[1]; ps[i].a[1] = tmp;
        }
    }
    /* Now sort it */
    if (nr > 1)
    {
        qsort(ps, nr, static_cast<size_t>(sizeof(ps[0])), idcomp);
    }
}

static int n_hydro(const int a[], char ***atomname)
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
static void clean_dih(t_param *dih, int *ndih, t_param improper[], int nimproper,
                      t_atoms *atoms, bool bKeepAllGeneratedDihedrals,
                      bool bRemoveDihedralIfWithImproper)
{
    int   i, j, k, l;
    int  *index, nind;

    /* Construct the list of the indices of the dihedrals
     * (i.e. generated or read) that might be kept. */
    snew(index, *ndih+1);
    if (bKeepAllGeneratedDihedrals)
    {
        fprintf(stderr, "Keeping all generated dihedrals\n");
        nind = *ndih;
        for (i = 0; i < nind; i++)
        {
            index[i] = i;
        }
        index[nind] = *ndih;
    }
    else
    {
        nind = 0;
        /* Check if generated dihedral i should be removed. The
         * dihedrals have been sorted by dcomp() above, so all those
         * on the same two central atoms are together, with those from
         * the .rtp file preceding those that were automatically
         * generated. We remove the latter if the former exist. */
        for (i = 0; i < *ndih; i++)
        {
            /* Keep the dihedrals that were defined in the .rtp file,
             * and the dihedrals that were generated and different
             * from the last one (whether it was generated or not). */
            if (was_dihedral_set_in_rtp(&dih[i]) ||
                0 == i ||
                !is_dihedral_on_same_bond(&dih[i], &dih[i-1]))
            {
                index[nind++] = i;
            }
        }
        index[nind] = *ndih;
    }

    k = 0;
    for (i = 0; i < nind; i++)
    {
        bool bWasSetInRTP = was_dihedral_set_in_rtp(&dih[index[i]]);
        bool bKeep        = TRUE;
        if (!bWasSetInRTP && bRemoveDihedralIfWithImproper)
        {
            /* Remove the dihedral if there is an improper on the same
             * bond. */
            for (j = 0; j < nimproper && bKeep; j++)
            {
                bKeep = !is_dihedral_on_same_bond(&dih[index[i]], &improper[j]);
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
                for (l = index[i];
                     (l < index[i+1] &&
                      is_dihedral_on_same_bond(&dih[index[i]], &dih[l]));
                     l++)
                {
                    int nh = n_hydro(dih[l].a, atoms->atomname);
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
                cpparam(&dih[k], &dih[bestl]);
            }
            k++;
        }
    }

    for (i = k; i < *ndih; i++)
    {
        strcpy(dih[i].s, "");
    }
    *ndih = k;

    sfree(index);
}

static int get_impropers(t_atoms *atoms, gmx::ArrayRef<MoleculePatchDatabase> globalPatches, t_param **improper,
                         bool bAllowMissing)
{
    int           nimproper, start, ninc, nalloc;
    int           ai[MAXATOMLIST];

    ninc   = 500;
    nalloc = ninc;
    snew(*improper, nalloc);

    /* Add all the impropers from the residue database to the list. */
    nimproper = 0;
    start     = 0;
    if (!globalPatches.empty())
    {
        for (int i = 0; (i < atoms->nres); i++)
        {
            BondedInteractionList *impropers = &globalPatches[i].rb[ebtsIDIHS];
            for (const auto &bondeds : impropers->b)
            {
                bool bStop = false;
                for (int k = 0; (k < 4) && !bStop; k++)
                {
                    ai[k] = search_atom(bondeds.a[k].c_str(), start,
                                        atoms,
                                        "improper", bAllowMissing);
                    if (ai[k] == -1)
                    {
                        bStop = true;
                    }
                }
                if (!bStop)
                {
                    if (nimproper == nalloc)
                    {
                        nalloc += ninc;
                        srenew(*improper, nalloc);
                    }
                    /* Not broken out */
                    set_p(&((*improper)[nimproper]), ai, nullptr, bondeds.s.c_str());
                    nimproper++;
                }
            }
            while ((start < atoms->nr) && (atoms->atom[start].resind == i))
            {
                start++;
            }
        }
    }

    return nimproper;
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

static bool is_hydro(t_atoms *atoms, int ai)
{
    return ((*(atoms->atomname[ai]))[0] == 'H');
}

static void get_atomnames_min(int n, char **anm,
                              int resind, t_atoms *atoms, const int *a)
{
    int m;

    /* Assume ascending residue numbering */
    for (m = 0; m < n; m++)
    {
        if (atoms->atom[a[m]].resind < resind)
        {
            strcpy(anm[m], "-");
        }
        else if (atoms->atom[a[m]].resind > resind)
        {
            strcpy(anm[m], "+");
        }
        else
        {
            strcpy(anm[m], "");
        }
        strcat(anm[m], *(atoms->atomname[a[m]]));
    }
}

static void gen_excls(t_atoms                             *atoms,
                      t_excls                             *excls,
                      gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
                      bool                                 bAllowMissing)
{
    int astart = 0;
    for (int a = 0; a < atoms->nr; a++)
    {
        int r = atoms->atom[a].resind;
        if (a == atoms->nr-1 || atoms->atom[a+1].resind != r)
        {
            BondedInteractionList *hbexcl = &globalPatches[r].rb[ebtsEXCLS];

            for (const auto &bondeds : hbexcl->b)
            {
                const char *anm = bondeds.a[0].c_str();
                int         i1  = search_atom(anm, astart, atoms,
                                              "exclusion", bAllowMissing);
                anm = bondeds.a[1].c_str();
                int i2  = search_atom(anm, astart, atoms,
                                      "exclusion", bAllowMissing);
                if (i1 != -1 && i2 != -1)
                {
                    if (i1 > i2)
                    {
                        int itmp = i1;
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

    for (int a = 0; a < atoms->nr; a++)
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
void gen_pad(t_nextnb *nnb, t_atoms *atoms, gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
             t_params plist[], t_excls excls[], gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
             bool bAllowMissing)
{
    t_param    *ang, *dih, *pai, *improper;
    char      **anm;
    int         ninc, maxang, maxdih, maxpai;
    int         nang, ndih, npai, nimproper;

    /* These are the angles, dihedrals and pairs that we generate
     * from the bonds. The ones that are already there from the rtp file
     * will be retained.
     */
    nang   = 0;
    npai   = 0;
    ndih   = 0;
    ninc   = 500;
    maxang = maxdih = maxpai = ninc;
    snew(ang, maxang);
    snew(dih, maxdih);
    snew(pai, maxpai);

    snew(anm, 4);
    for (int i = 0; i < 4; i++)
    {
        snew(anm[i], 12);
    }

    if (!globalPatches.empty())
    {
        gen_excls(atoms, excls, globalPatches, bAllowMissing);
        /* mark all entries as not matched yet */
        for (int i = 0; i < atoms->nres; i++)
        {
            for (int j = 0; j < ebtsNR; j++)
            {
                for (auto &bondeds : globalPatches[i].rb[j].b)
                {
                    bondeds.match = false;
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
                        if (nang == maxang)
                        {
                            maxang += ninc;
                            srenew(ang, maxang);
                        }
                        ang[nang].ai() = i;
                        ang[nang].aj() = j1;
                        ang[nang].ak() = k1;
                        ang[nang].c0() = NOTSET;
                        ang[nang].c1() = NOTSET;
                        set_p_string(&(ang[nang]), "");
                        if (!globalPatches.empty())
                        {
                            int minres = atoms->atom[ang[nang].a[0]].resind;
                            int maxres = minres;
                            for (int m = 1; m < 3; m++)
                            {
                                minres = std::min(minres, atoms->atom[ang[nang].a[m]].resind);
                                maxres = std::max(maxres, atoms->atom[ang[nang].a[m]].resind);
                            }
                            int res = 2*minres-maxres;
                            do
                            {
                                res += maxres-minres;
                                get_atomnames_min(3, anm, res, atoms, ang[nang].a);
                                BondedInteractionList *hbang = &globalPatches[res].rb[ebtsANGLES];
                                for (auto &bondeds : hbang->b)
                                {
                                    if (anm[1] == bondeds.aj())
                                    {
                                        bool bFound = false;
                                        for (int m = 0; m < 3; m += 2)
                                        {
                                            bFound = (bFound ||
                                                      ((anm[m] == bondeds.ai()) &&
                                                       (anm[2-m] == bondeds.ak())));
                                        }
                                        if (bFound)
                                        {
                                            set_p_string(&(ang[nang]), bondeds.s.c_str());
                                            /* Mark that we found a match for this entry */
                                            bondeds.match = true;
                                        }
                                    }
                                }
                            }
                            while (res < maxres);
                        }
                        nang++;
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
                                if (ndih == maxdih)
                                {
                                    maxdih += ninc;
                                    srenew(dih, maxdih);
                                }
                                dih[ndih].ai() = i;
                                dih[ndih].aj() = j1;
                                dih[ndih].ak() = k1;
                                dih[ndih].al() = l1;
                                for (int m = 0; m < MAXFORCEPARAM; m++)
                                {
                                    dih[ndih].c[m] = NOTSET;
                                }
                                set_p_string(&(dih[ndih]), "");
                                int nFound = 0;
                                if (!globalPatches.empty())
                                {
                                    int minres = atoms->atom[dih[ndih].a[0]].resind;
                                    int maxres = minres;
                                    for (int m = 1; m < 4; m++)
                                    {
                                        minres = std::min(minres, atoms->atom[dih[ndih].a[m]].resind);
                                        maxres = std::max(maxres, atoms->atom[dih[ndih].a[m]].resind);
                                    }
                                    int res = 2*minres-maxres;
                                    do
                                    {
                                        res += maxres-minres;
                                        get_atomnames_min(4, anm, res, atoms, dih[ndih].a);
                                        BondedInteractionList *hbdih = &globalPatches[res].rb[ebtsPDIHS];
                                        for (auto &bondeds : hbdih->b)
                                        {
                                            bool bFound = false;
                                            for (int m = 0; m < 2; m++)
                                            {
                                                bFound = (bFound ||
                                                          ((anm[3*m] == bondeds.ai()) &&
                                                           (anm[1+m] == bondeds.aj()) &&
                                                           (anm[2-m] == bondeds.ak()) &&
                                                           (anm[3-3*m] == bondeds.al())));
                                            }
                                            if (bFound)
                                            {
                                                set_p_string(&dih[ndih], bondeds.s.c_str());
                                                /* Mark that we found a match for this entry */
                                                bondeds.match = true;

                                                /* Set the last parameter to be able to see
                                                   if the dihedral was in the rtp list.
                                                 */
                                                dih[ndih].c[MAXFORCEPARAM-1] = DIHEDRAL_WAS_SET_IN_RTP;
                                                nFound++;
                                                ndih++;
                                                /* Set the next direct in case the rtp contains
                                                   multiple entries for this dihedral.
                                                 */
                                                if (ndih == maxdih)
                                                {
                                                    maxdih += ninc;
                                                    srenew(dih, maxdih);
                                                }
                                                dih[ndih].ai() = i;
                                                dih[ndih].aj() = j1;
                                                dih[ndih].ak() = k1;
                                                dih[ndih].al() = l1;
                                                for (int m = 0; m < MAXFORCEPARAM; m++)
                                                {
                                                    dih[ndih].c[m] = NOTSET;
                                                }
                                            }
                                        }
                                    }
                                    while (res < maxres);
                                }
                                if (nFound == 0)
                                {
                                    if (ndih == maxdih)
                                    {
                                        maxdih += ninc;
                                        srenew(dih, maxdih);
                                    }
                                    dih[ndih].ai() = i;
                                    dih[ndih].aj() = j1;
                                    dih[ndih].ak() = k1;
                                    dih[ndih].al() = l1;
                                    for (int m = 0; m < MAXFORCEPARAM; m++)
                                    {
                                        dih[ndih].c[m] = NOTSET;
                                    }
                                    set_p_string(&(dih[ndih]), "");
                                    ndih++;
                                }

                                int nbd = nb_dist(nnb, i, l1);
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
                                        if (rtpFFDB[0].bGenerateHH14Interactions ||
                                            !(is_hydro(atoms, i1) && is_hydro(atoms, i2)))
                                        {
                                            if (npai == maxpai)
                                            {
                                                maxpai += ninc;
                                                srenew(pai, maxpai);
                                            }
                                            pai[npai].ai() = i1;
                                            pai[npai].aj() = i2;
                                            pai[npai].c0() = NOTSET;
                                            pai[npai].c1() = NOTSET;
                                            set_p_string(&(pai[npai]), "");
                                            npai++;
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
            BondedInteractionList *hbang = &globalPatches[i].rb[ebtsANGLES];
            for (auto &bondeds : hbang->b)
            {
                if (bondeds.match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                if (nang == maxang)
                {
                    maxang += ninc;
                    srenew(ang, maxang);
                }
                bool bFound = true;
                for (int k = 0; k < 3 && bFound; k++)
                {
                    const char *p   = bondeds.a[k].c_str();
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
                    ang[nang].a[k] = search_res_atom(p, res, atoms, "angle", TRUE);
                    bFound         = (ang[nang].a[k] != -1);
                }
                ang[nang].c0() = NOTSET;
                ang[nang].c1() = NOTSET;

                if (bFound)
                {
                    set_p_string(&(ang[nang]), bondeds.s.c_str());
                    bondeds.match = true;
                    /* Incrementing nang means we save this angle */
                    nang++;
                }
            }

            /* Add remaining dihedrals from hackblock */
            BondedInteractionList *hbdih = &globalPatches[i].rb[ebtsPDIHS];
            for (auto &bondeds : hbdih->b)
            {
                if (bondeds.match)
                {
                    /* We already used this entry, continue to the next */
                    continue;
                }
                /* Hm - entry not used, let's see if we can find all atoms */
                if (ndih == maxdih)
                {
                    maxdih += ninc;
                    srenew(dih, maxdih);
                }
                bool bFound = true;
                for (int k = 0; k < 4 && bFound; k++)
                {
                    const char *p   = bondeds.a[k].c_str();
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
                    dih[ndih].a[k] = search_res_atom(p, res, atoms, "dihedral", TRUE);
                    bFound         = (dih[ndih].a[k] != -1);
                }
                for (int m = 0; m < MAXFORCEPARAM; m++)
                {
                    dih[ndih].c[m] = NOTSET;
                }

                if (bFound)
                {
                    set_p_string(&(dih[ndih]), bondeds.s.c_str());
                    bondeds.match = true;
                    /* Incrementing ndih means we save this dihedral */
                    ndih++;
                }
            }
        }
    }

    /* Sort angles with respect to j-i-k (middle atom first) */
    if (nang > 1)
    {
        qsort(ang, nang, static_cast<size_t>(sizeof(ang[0])), acomp);
    }

    /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
    if (ndih > 1)
    {
        qsort(dih, ndih, static_cast<size_t>(sizeof(dih[0])), dcomp);
    }

    /* Sort the pairs */
    if (npai > 1)
    {
        qsort(pai, npai, static_cast<size_t>(sizeof(pai[0])), pcomp);
    }
    if (npai > 0)
    {
        /* Remove doubles, could occur in 6-rings, such as phenyls,
           maybe one does not want this when fudgeQQ < 1.
         */
        fprintf(stderr, "Before cleaning: %d pairs\n", npai);
        rm2par(pai, &npai, preq);
    }

    /* Get the impropers from the database */
    nimproper = get_impropers(atoms, globalPatches, &improper, bAllowMissing);

    /* Sort the impropers */
    sort_id(nimproper, improper);

    if (ndih > 0)
    {
        fprintf(stderr, "Before cleaning: %d dihedrals\n", ndih);
        clean_dih(dih, &ndih, improper, nimproper, atoms,
                  rtpFFDB[0].bKeepAllGeneratedDihedrals,
                  rtpFFDB[0].bRemoveDihedralIfWithImproper);
    }

    /* Now we have unique lists of angles and dihedrals
     * Copy them into the destination struct
     */
    cppar(ang, nang, plist, F_ANGLES);
    cppar(dih, ndih, plist, F_PDIHS);
    cppar(improper, nimproper, plist, F_IDIHS);
    cppar(pai, npai, plist, F_LJ14);

    /* Remove all exclusions which are within nrexcl */
    clean_excls(nnb, rtpFFDB[0].nrexcl, excls);

    sfree(ang);
    sfree(dih);
    sfree(improper);
    sfree(pai);
}
