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

#include "hxprops.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <filesystem>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

real ellipticity(int nres, t_bb bb[])
{
    typedef struct
    {
        real phi, psi, w;
    } t_ppwstr;
    // Avoid warnings about narrowing conversions from double to real
#ifdef _MSC_VER
#    pragma warning(disable : 4838)
#endif
    static const t_ppwstr ppw[] = { { -67, -44, 0.31 },     { -66, -41, 0.31 }, { -59, -44, 0.44 },
                                    { -57, -47, 0.56 },     { -53, -52, 0.78 }, { -48, -57, 1.00 },
                                    { -70.5, -35.8, 0.15 }, { -57, -79, 0.23 }, { -38, -78, 1.20 },
                                    { -60, -30, 0.24 },     { -54, -28, 0.46 }, { -44, -33, 0.68 } };
#ifdef _MSC_VER
#    pragma warning(default : 4838)
#endif
#define NPPW asize(ppw)

    int  i, j;
    real ell, pp2, phi, psi;

    ell = 0;
    for (i = 0; (i < nres); i++)
    {
        phi = bb[i].phi;
        psi = bb[i].psi;
        for (j = 0; (j < NPPW); j++)
        {
            pp2 = gmx::square(phi - ppw[j].phi) + gmx::square(psi - ppw[j].psi);
            if (pp2 < 64)
            {
                bb[i].nhx++;
                ell += ppw[j].w;
                break;
            }
        }
    }
    return ell;
}

real ahx_len(int gnx, const int index[], rvec x[])
/* Assume we have a list of Calpha atoms only! */
{
    rvec dx;

    rvec_sub(x[index[0]], x[index[gnx - 1]], dx);

    return norm(dx);
}

real radius(FILE* fp, int nca, const int ca_index[], rvec x[])
/* Assume we have all the backbone */
{
    real dl2, dlt;
    int  i, ai;

    dlt = 0;
    for (i = 0; (i < nca); i++)
    {
        ai  = ca_index[i];
        dl2 = gmx::square(x[ai][XX]) + gmx::square(x[ai][YY]);

        if (fp)
        {
            fprintf(fp, "  %10g", dl2);
        }

        dlt += dl2;
    }
    if (fp)
    {
        fprintf(fp, "\n");
    }

    return std::sqrt(dlt / nca);
}

static real rot(rvec x1, const rvec x2)
{
    real phi1, dphi, cp, sp;
    real xx, yy;

    phi1 = std::atan2(x1[YY], x1[XX]);
    cp   = std::cos(phi1);
    sp   = std::sin(phi1);
    xx   = cp * x2[XX] + sp * x2[YY];
    yy   = -sp * x2[XX] + cp * x2[YY];

    dphi = gmx::c_rad2Deg * std::atan2(yy, xx);

    return dphi;
}

real twist(int nca, const int caindex[], rvec x[])
{
    real pt, dphi;
    int  i, a0, a1;

    pt = 0;
    a0 = caindex[0];
    for (i = 1; (i < nca); i++)
    {
        a1 = caindex[i];

        dphi = rot(x[a0], x[a1]);
        if (dphi < -90)
        {
            dphi += 360;
        }
        pt += dphi;
        a0 = a1;
    }

    return (pt / (nca - 1));
}

real ca_phi(int gnx, const int index[], rvec x[])
/* Assume we have a list of Calpha atoms only! */
{
    real phi, phitot;
    int  i, ai, aj, ak, al, t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;

    if (gnx <= 4)
    {
        return 0;
    }

    phitot = 0;
    for (i = 0; (i < gnx - 4); i++)
    {
        ai  = index[i + 0];
        aj  = index[i + 1];
        ak  = index[i + 2];
        al  = index[i + 3];
        phi = gmx::c_rad2Deg
              * dih_angle(x[ai], x[aj], x[ak], x[al], nullptr, r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);
        phitot += phi;
    }

    return (phitot / (gnx - 4.0));
}

real dip(int nbb, int const bbind[], const rvec x[], const t_atom atom[])
{
    int  i, m, ai;
    rvec dipje;
    real q;

    clear_rvec(dipje);
    for (i = 0; (i < nbb); i++)
    {
        ai = bbind[i];
        q  = atom[ai].q;
        for (m = 0; (m < DIM); m++)
        {
            dipje[m] += x[ai][m] * q;
        }
    }
    return norm(dipje);
}

real rise(int gnx, const int index[], rvec x[])
/* Assume we have a list of Calpha atoms only! */
{
    real z, z0, ztot;
    int  i, ai;

    ai   = index[0];
    z0   = x[ai][ZZ];
    ztot = 0;
    for (i = 1; (i < gnx); i++)
    {
        ai = index[i];
        z  = x[ai][ZZ];
        ztot += (z - z0);
        z0 = z;
    }

    return (ztot / (gnx - 1.0));
}

void av_hblen(FILE* fp3, FILE* fp3a, FILE* fp4, FILE* fp4a, FILE* fp5, FILE* fp5a, real t, int nres, t_bb bb[])
{
    int  i, n3 = 0, n4 = 0, n5 = 0;
    real d3 = 0, d4 = 0, d5 = 0;

    for (i = 0; (i < nres - 3); i++)
    {
        if (bb[i].bHelix)
        {
            fprintf(fp3a, "%10g", bb[i].d3);
            n3++;
            d3 += bb[i].d3;
            if (i < nres - 4)
            {
                fprintf(fp4a, "%10g", bb[i].d4);
                n4++;
                d4 += bb[i].d4;
            }
            if (i < nres - 5)
            {
                fprintf(fp5a, "%10g", bb[i].d5);
                n5++;
                d5 += bb[i].d5;
            }
        }
    }
    fprintf(fp3, "%10g  %10g\n", t, d3 / n3);
    fprintf(fp4, "%10g  %10g\n", t, d4 / n4);
    fprintf(fp5, "%10g  %10g\n", t, d5 / n5);
    fprintf(fp3a, "\n");
    fprintf(fp4a, "\n");
    fprintf(fp5a, "\n");
}

void av_phipsi(FILE* fphi, FILE* fpsi, FILE* fphi2, FILE* fpsi2, real t, int nres, t_bb bb[])
{
    int  i, n = 0;
    real phi = 0, psi = 0;

    fprintf(fphi2, "%10g", t);
    fprintf(fpsi2, "%10g", t);
    for (i = 0; (i < nres); i++)
    {
        if (bb[i].bHelix)
        {
            phi += bb[i].phi;
            psi += bb[i].psi;
            fprintf(fphi2, "  %10g", bb[i].phi);
            fprintf(fpsi2, "  %10g", bb[i].psi);
            n++;
        }
    }
    fprintf(fphi, "%10g  %10g\n", t, (phi / n));
    fprintf(fpsi, "%10g  %10g\n", t, (psi / n));
    fprintf(fphi2, "\n");
    fprintf(fpsi2, "\n");
}

static void set_ahcity(int nbb, t_bb bb[])
{
    real pp2;
    int  n;

    for (n = 0; (n < nbb); n++)
    {
        pp2 = gmx::square(bb[n].phi - PHI_AHX) + gmx::square(bb[n].psi - PSI_AHX);

        bb[n].bHelix = FALSE;
        if (pp2 < 2500)
        {
            if ((bb[n].d4 < 0.36) || ((n > 0) && bb[n - 1].bHelix))
            {
                bb[n].bHelix = TRUE;
            }
        }
    }
}

t_bb* mkbbind(const char* fn,
              int*        nres,
              int*        nbb,
              int         res0,
              int*        nall,
              int**       index,
              char***     atomname,
              t_atom      atom[],
              t_resinfo*  resinfo)
{
    static const char* bb_nm[] = { "N", "H", "CA", "C", "O", "HN" };
#define NBB asize(bb_nm)
    t_bb* bb;
    char* grpname;
    int   gnx, r0, r1;

    fprintf(stderr, "Please select a group containing the entire backbone\n");
    rd_index(fn, 1, &gnx, index, &grpname);
    *nall = gnx;
    fprintf(stderr, "Checking group %s\n", grpname);
    r0 = r1 = atom[(*index)[0]].resind;
    for (int i = 1; (i < gnx); i++)
    {
        r0 = std::min(r0, atom[(*index)[i]].resind);
        r1 = std::max(r1, atom[(*index)[i]].resind);
    }
    int rnr = r1 - r0 + 1;
    fprintf(stderr, "There are %d residues\n", rnr);
    snew(bb, rnr);
    for (int i = 0; (i < rnr); i++)
    {
        bb[i].N = bb[i].H = bb[i].CA = bb[i].C = bb[i].O = -1;
        bb[i].resno                                      = res0 + i;
    }

    for (int i = 0; (i < gnx); i++)
    {
        int ai = (*index)[i];
        // Create an index into the residue index for the topology.
        int resindex = atom[ai].resind;
        // Create an index into the residues present in the selected
        // index group.
        int bbindex = resindex - r0;
        if (std::strcmp(*(resinfo[resindex].name), "PRO") == 0)
        {
            // For PRO in a peptide, there is no H bound to backbone
            // N, so use CD instead.
            if (std::strcmp(*(atomname[ai]), "CD") == 0)
            {
                bb[bbindex].H = ai;
            }
        }
        int k = 0;
        for (; (k < NBB); k++)
        {
            if (std::strcmp(bb_nm[k], *(atomname[ai])) == 0)
            {
                break;
            }
        }
        switch (k)
        {
            case 0: bb[bbindex].N = ai; break;
            case 1:
            case 5:
                /* No attempt to address the case where some weird input has both H and HN atoms in the group */
                bb[bbindex].H = ai;
                break;
            case 2: bb[bbindex].CA = ai; break;
            case 3: bb[bbindex].C = ai; break;
            case 4: bb[bbindex].O = ai; break;
            default: break;
        }
    }

    int i0 = 0;
    for (; (i0 < rnr); i0++)
    {
        if ((bb[i0].N != -1) && (bb[i0].H != -1) && (bb[i0].CA != -1) && (bb[i0].C != -1)
            && (bb[i0].O != -1))
        {
            break;
        }
    }
    int i1 = rnr - 1;
    for (; (i1 >= 0); i1--)
    {
        if ((bb[i1].N != -1) && (bb[i1].H != -1) && (bb[i1].CA != -1) && (bb[i1].C != -1)
            && (bb[i1].O != -1))
        {
            break;
        }
    }
    if (i0 == 0)
    {
        i0++;
    }
    if (i1 == rnr - 1)
    {
        i1--;
    }

    for (int i = i0; (i < i1); i++)
    {
        bb[i].Cprev = bb[i - 1].C;
        bb[i].Nnext = bb[i + 1].N;
    }
    rnr = std::max(0, i1 - i0 + 1);
    fprintf(stderr, "There are %d complete backbone residues (from %d to %d)\n", rnr, bb[i0].resno, bb[i1].resno);
    if (rnr == 0)
    {
        gmx_fatal(FARGS, "Zero complete backbone residues were found, cannot proceed");
    }
    for (int i = 0; (i < rnr); i++, i0++)
    {
        bb[i] = bb[i0];
    }

    /* Set the labels */
    for (int i = 0; (i < rnr); i++)
    {
        int resindex = atom[bb[i].CA].resind;
        sprintf(bb[i].label, "%s%d", *(resinfo[resindex].name), resinfo[resindex].nr);
    }

    *nres = rnr;
    *nbb  = rnr * asize(bb_nm);

    return bb;
}

real pprms(FILE* fp, int nbb, t_bb bb[])
{
    int  i, n;
    real rms, rmst, rms2;

    rmst = rms2 = 0;
    for (i = n = 0; (i < nbb); i++)
    {
        if (bb[i].bHelix)
        {
            rms = std::sqrt(bb[i].pprms2);
            rmst += rms;
            rms2 += bb[i].pprms2;
            fprintf(fp, "%10g  ", rms);
            n++;
        }
    }
    fprintf(fp, "\n");
    rms = std::sqrt(rms2 / n - gmx::square(rmst / n));

    return rms;
}

void calc_hxprops(int nres, t_bb bb[], const rvec x[])
{
    int  i, ao, an, t1, t2, t3;
    rvec dx, r_ij, r_kj, r_kl, m, n;

    for (i = 0; (i < nres); i++)
    {
        ao       = bb[i].O;
        bb[i].d4 = bb[i].d3 = bb[i].d5 = 0;
        if (i < nres - 3)
        {
            an = bb[i + 3].N;
            rvec_sub(x[ao], x[an], dx);
            bb[i].d3 = norm(dx);
        }
        if (i < nres - 4)
        {
            an = bb[i + 4].N;
            rvec_sub(x[ao], x[an], dx);
            bb[i].d4 = norm(dx);
        }
        if (i < nres - 5)
        {
            an = bb[i + 5].N;
            rvec_sub(x[ao], x[an], dx);
            bb[i].d5 = norm(dx);
        }

        bb[i].phi = gmx::c_rad2Deg
                    * dih_angle(x[bb[i].Cprev],
                                x[bb[i].N],
                                x[bb[i].CA],
                                x[bb[i].C],
                                nullptr,
                                r_ij,
                                r_kj,
                                r_kl,
                                m,
                                n,
                                &t1,
                                &t2,
                                &t3);
        bb[i].psi = gmx::c_rad2Deg
                    * dih_angle(x[bb[i].N],
                                x[bb[i].CA],
                                x[bb[i].C],
                                x[bb[i].Nnext],
                                nullptr,
                                r_ij,
                                r_kj,
                                r_kl,
                                m,
                                n,
                                &t1,
                                &t2,
                                &t3);
        bb[i].pprms2 = gmx::square(bb[i].phi - PHI_AHX) + gmx::square(bb[i].psi - PSI_AHX);

        bb[i].jcaha += 1.4 * std::sin((bb[i].psi + 138.0) * gmx::c_deg2Rad)
                       - 4.1 * std::cos(2.0 * gmx::c_deg2Rad * (bb[i].psi + 138.0))
                       + 2.0 * std::cos(2.0 * gmx::c_deg2Rad * (bb[i].phi + 30.0));
    }
}

static void check_ahx(int nres, t_bb bb[], int* hstart, int* hend)
{
    int h0, h1, h0sav, h1sav;

    set_ahcity(nres, bb);
    h0 = h0sav = h1sav = 0;
    do
    {
        for (; (!bb[h0].bHelix) && (h0 < nres - 4); h0++) {}
        for (h1 = h0; bb[h1 + 1].bHelix && (h1 < nres - 1); h1++) {}
        if (h1 > h0)
        {
            /*fprintf(stderr,"Helix from %d to %d\n",h0,h1);*/
            if (h1 - h0 > h1sav - h0sav)
            {
                h0sav = h0;
                h1sav = h1;
            }
        }
        h0 = h1 + 1;
    } while (h1 < nres - 1);
    *hstart = h0sav;
    *hend   = h1sav;
}

void do_start_end(int nres, t_bb bb[], int* nbb, int bbindex[], int* nca, int caindex[], gmx_bool bRange, int rStart, int rEnd)
{
    int i, j, hstart = 0, hend = 0;

    if (bRange)
    {
        for (i = 0; (i < nres); i++)
        {
            if ((bb[i].resno >= rStart) && (bb[i].resno <= rEnd))
            {
                bb[i].bHelix = TRUE;
            }
            if (bb[i].resno == rStart)
            {
                hstart = i;
            }
            if (bb[i].resno == rEnd)
            {
                hend = i;
            }
        }
    }
    else
    {
        /* Find start and end of longest helix fragment */
        check_ahx(nres, bb, &hstart, &hend);
    }
    fprintf(stderr, "helix from: %d through %d\n", bb[hstart].resno, bb[hend].resno);

    for (j = 0, i = hstart; (i <= hend); i++)
    {
        bbindex[j++]        = bb[i].N;
        bbindex[j++]        = bb[i].H;
        bbindex[j++]        = bb[i].CA;
        bbindex[j++]        = bb[i].C;
        bbindex[j++]        = bb[i].O;
        caindex[i - hstart] = bb[i].CA;
    }
    *nbb = j;
    *nca = (hend - hstart + 1);
}

void pr_bb(FILE* fp, int nres, t_bb bb[])
{
    int i;

    fprintf(fp, "\n");
    fprintf(fp,
            "%3s %3s %3s %3s %3s %7s %7s %7s %7s %7s %3s\n",
            "AA",
            "N",
            "Ca",
            "C",
            "O",
            "Phi",
            "Psi",
            "D3",
            "D4",
            "D5",
            "Hx?");
    for (i = 0; (i < nres); i++)
    {
        fprintf(fp,
                "%3d %3d %3d %3d %3d %7.2f %7.2f %7.3f %7.3f %7.3f %3s\n",
                bb[i].resno,
                bb[i].N,
                bb[i].CA,
                bb[i].C,
                bb[i].O,
                bb[i].phi,
                bb[i].psi,
                bb[i].d3,
                bb[i].d4,
                bb[i].d5,
                bb[i].bHelix ? "Yes" : "No");
    }
    fprintf(fp, "\n");
}
