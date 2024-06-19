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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "hizzie.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <filesystem>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int in_strings(char* key, int nstr, const char** str)
{
    int j;

    for (j = 0; (j < nstr); j++)
    {
        if (strcmp(str[j], key) == 0)
        {
            return j;
        }
    }

    return -1;
}

static bool hbond(rvec x[], int i, int j, real distance)
{
    real tol = distance * distance;
    rvec tmp;

    rvec_sub(x[i], x[j], tmp);

    return (iprod(tmp, tmp) < tol);
}

static void chk_allhb(t_atoms* pdba, rvec x[], t_blocka* hb, const bool donor[], const bool accept[], real dist)
{
    int i, j, k, ii, natom;

    natom = pdba->nr;
    snew(hb->index, natom + 1);
    snew(hb->a, 6 * natom);
    hb->nr  = natom;
    hb->nra = 6 * natom;

    k = ii          = 0;
    hb->index[ii++] = 0;
    for (i = 0; (i < natom); i++)
    {
        if (donor[i])
        {
            for (j = i + 1; (j < natom); j++)
            {
                if ((accept[j]) && (hbond(x, i, j, dist)))
                {
                    hb->a[k++] = j;
                }
            }
        }
        else if (accept[i])
        {
            for (j = i + 1; (j < natom); j++)
            {
                if ((donor[j]) && (hbond(x, i, j, dist)))
                {
                    hb->a[k++] = j;
                }
            }
        }
        hb->index[ii++] = k;
    }
    hb->nra = k;
}

static bool chk_hbonds(int i, t_atoms* pdba, rvec x[], const bool ad[], bool hbond[], rvec xh, real angle, real dist)
{
    bool bHB;
    int  j, aj, ri, natom;
    real d2, dist2, a;
    rvec nh, oh;

    natom = pdba->nr;
    bHB   = FALSE;
    ri    = pdba->atom[i].resind;
    dist2 = gmx::square(dist);
    for (j = 0; (j < natom); j++)
    {
        /* Check whether the other atom is a donor/acceptor and not i */
        if ((ad[j]) && (j != i))
        {
            /* Check whether the other atom is on the same ring as well */
            if ((pdba->atom[j].resind != ri)
                || ((strcmp(*pdba->atomname[j], "ND1") != 0) && (strcmp(*pdba->atomname[j], "NE2") != 0)))
            {
                aj = j;
                d2 = distance2(x[i], x[j]);
                rvec_sub(x[i], xh, nh);
                rvec_sub(x[aj], xh, oh);
                a = gmx::c_rad2Deg * std::acos(cos_angle(nh, oh));
                if ((d2 < dist2) && (a > angle))
                {
                    hbond[i] = TRUE;
                    bHB      = TRUE;
                }
            }
        }
    }
    return bHB;
}

static void calc_ringh(rvec xattach, rvec xb, rvec xc, rvec xh)
{
    rvec tab, tac;
    real n;

    /* Add a proton on a ring to atom attach at distance 0.1 nm */
    rvec_sub(xattach, xb, tab);
    rvec_sub(xattach, xc, tac);
    rvec_add(tab, tac, xh);
    n = 0.1 / norm(xh);
    svmul(n, xh, xh);
    rvec_inc(xh, xattach);
}

void set_histp(t_atoms* pdba, rvec* x, t_symtab* symtab, real angle, real dist)
{
    static const char* prot_acc[] = { "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OW" };
#define NPA asize(prot_acc)
    static const char* prot_don[] = { "N",  "NH1", "NH2", "NE", "ND1", "ND2", "NE2",
                                      "NZ", "OG",  "OG1", "OH", "NE1", "OW" };
#define NPD asize(prot_don)

    bool *          donor, *acceptor;
    bool*           hbond;
    bool            bHDd, bHEd;
    rvec            xh1, xh2;
    int             natom;
    int             i, j, nd, na, hisind;
    HistidineStates type = HistidineStates::Count;
    int             nd1, ne2, cg, cd2, ce1;
    t_blocka*       hb;
    char*           atomnm;

    natom = pdba->nr;

    i = 0;
    while (i < natom && gmx_strcasecmp(*pdba->resinfo[pdba->atom[i].resind].name, "HIS") != 0)
    {
        i++;
    }
    if (natom == i)
    {
        return;
    }

    /* A histidine residue exists that requires automated assignment, so
     * doing the analysis of donors and acceptors is worthwhile. */
    fprintf(stderr,
            "Analysing hydrogen-bonding network for automated assignment of histidine\n"
            " protonation.");

    snew(donor, natom);
    snew(acceptor, natom);
    snew(hbond, natom);
    snew(hb, 1);

    nd = na = 0;
    for (j = 0; (j < natom); j++)
    {
        if (in_strings(*pdba->atomname[j], NPA, prot_acc) != -1)
        {
            acceptor[j] = TRUE;
            na++;
        }
        if (in_strings(*pdba->atomname[j], NPD, prot_don) != -1)
        {
            donor[j] = TRUE;
            nd++;
        }
    }
    fprintf(stderr, " %d donors and %d acceptors were found.\n", nd, na);
    chk_allhb(pdba, x, hb, donor, acceptor, dist);
    fprintf(stderr, "There are %d hydrogen bonds\n", hb->nra);

    /* Now do the HIS stuff */
    hisind = -1;
    while (i < natom)
    {
        if (gmx_strcasecmp(*pdba->resinfo[pdba->atom[i].resind].name, "HIS") != 0)
        {
            i++;
        }
        else
        {
            if (pdba->atom[i].resind != hisind)
            {
                hisind = pdba->atom[i].resind;

                /* Find the  atoms in the ring */
                nd1 = ne2 = cg = cd2 = ce1 = -1;
                while (i < natom && pdba->atom[i].resind == hisind)
                {
                    atomnm = *pdba->atomname[i];
                    if (strcmp(atomnm, "CD2") == 0)
                    {
                        cd2 = i;
                    }
                    else if (strcmp(atomnm, "CG") == 0)
                    {
                        cg = i;
                    }
                    else if (strcmp(atomnm, "CE1") == 0)
                    {
                        ce1 = i;
                    }
                    else if (strcmp(atomnm, "ND1") == 0)
                    {
                        nd1 = i;
                    }
                    else if (strcmp(atomnm, "NE2") == 0)
                    {
                        ne2 = i;
                    }

                    i++;
                }

                if (!((cg == -1) || (cd2 == -1) || (ce1 == -1) || (nd1 == -1) || (ne2 == -1)))
                {
                    calc_ringh(x[nd1], x[cg], x[ce1], xh1);
                    calc_ringh(x[ne2], x[ce1], x[cd2], xh2);

                    bHDd = chk_hbonds(nd1, pdba, x, acceptor, hbond, xh1, angle, dist);
                    chk_hbonds(nd1, pdba, x, donor, hbond, xh1, angle, dist);
                    bHEd = chk_hbonds(ne2, pdba, x, acceptor, hbond, xh2, angle, dist);
                    chk_hbonds(ne2, pdba, x, donor, hbond, xh2, angle, dist);

                    if (bHDd)
                    {
                        if (bHEd)
                        {
                            type = HistidineStates::H;
                        }
                        else
                        {
                            type = HistidineStates::A;
                        }
                    }
                    else
                    {
                        type = HistidineStates::B;
                    }
                    fprintf(stderr,
                            "Will use %s for residue %d\n",
                            enumValueToString(type),
                            pdba->resinfo[hisind].nr);
                }
                else
                {
                    gmx_fatal(FARGS, "Incomplete ring in HIS%d", pdba->resinfo[hisind].nr);
                }

                pdba->resinfo[hisind].rtp = put_symtab(symtab, enumValueToString(type));
            }
        }
    }
    done_blocka(hb);
    sfree(hb);
    sfree(donor);
    sfree(acceptor);
    sfree(hbond);
}
