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

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

static void calc_com_pbc(int nrefat, t_topology* top, rvec x[], t_pbc* pbc, const int index[], rvec xref, gmx_bool bPBC)
{
    const real tol = 1e-4;
    gmx_bool   bChanged;
    int        m, j, ai, iter;
    real       mass, mtot;
    rvec       dx, xtest;

    /* First simple calculation */
    clear_rvec(xref);
    mtot = 0;
    for (m = 0; (m < nrefat); m++)
    {
        ai   = index[m];
        mass = top->atoms.atom[ai].m;
        for (j = 0; (j < DIM); j++)
        {
            xref[j] += mass * x[ai][j];
        }
        mtot += mass;
    }
    svmul(1 / mtot, xref, xref);
    /* Now check if any atom is more than half the box from the COM */
    if (bPBC)
    {
        iter = 0;
        do
        {
            bChanged = FALSE;
            for (m = 0; (m < nrefat); m++)
            {
                ai   = index[m];
                mass = top->atoms.atom[ai].m / mtot;
                pbc_dx(pbc, x[ai], xref, dx);
                rvec_add(xref, dx, xtest);
                for (j = 0; (j < DIM); j++)
                {
                    if (std::abs(xtest[j] - x[ai][j]) > tol)
                    {
                        /* Here we have used the wrong image for contributing to the COM */
                        xref[j] += mass * (xtest[j] - x[ai][j]);
                        x[ai][j] = xtest[j];
                        bChanged = TRUE;
                    }
                }
            }
            if (bChanged)
            {
                printf("COM: %8.3f  %8.3f  %8.3f  iter = %d\n", xref[XX], xref[YY], xref[ZZ], iter);
            }
            iter++;
        } while (bChanged);
    }
}

int gmx_sorient(int argc, char* argv[])
{
    t_topology   top;
    PbcType      pbcType = PbcType::Unset;
    t_trxstatus* status;
    int          natoms;
    real         t;
    rvec *       xtop, *x;
    matrix       box;

    FILE*       fp;
    int         i, p, sa0, sa1, sa2, n, ntot, nf, m, *hist1, *hist2, *histn, nbin1, nbin2, nrbin;
    real *      histi1, *histi2, invbw, invrbw;
    double      sum1, sum2;
    int *       isize, nrefgrp, nrefat;
    int**       index;
    char**      grpname;
    real        inp, outp, nav, normfac, rmin2, rmax2, rcut, rcut2, r2, r;
    real        c1, c2;
    char        str[STRLEN];
    gmx_bool    bTPS;
    rvec        xref, dx, dxh1, dxh2, outer;
    gmx_rmpbc_t gpbc = nullptr;
    t_pbc       pbc;
    std::array<std::string, 2> legr = { "<cos(\\8q\\4\\s1\\N)>", "<3cos\\S2\\N(\\8q\\4\\s2\\N)-1>" };
    std::array<std::string, 2> legc = { "cos(\\8q\\4\\s1\\N)", "3cos\\S2\\N(\\8q\\4\\s2\\N)-1" };

    const char* desc[] = {
        "[THISMODULE] analyzes solvent orientation around solutes.",
        "It calculates two angles between the vector from one or more",
        "reference positions to the first atom of each solvent molecule:",
        "",
        " * [GRK]theta[grk][SUB]1[sub]: the angle with the vector from the first atom of "
        "the solvent",
        "   molecule to the midpoint between atoms 2 and 3.",
        " * [GRK]theta[grk][SUB]2[sub]: the angle with the normal of the solvent plane, "
        "defined by the",
        "   same three atoms, or, when the option [TT]-v23[tt] is set, ",
        "   the angle with the vector between atoms 2 and 3.",
        "",
        "The reference can be a set of atoms or",
        "the center of mass of a set of atoms. The group of solvent atoms should",
        "consist of 3 atoms per solvent molecule.",
        "Only solvent molecules between [TT]-rmin[tt] and [TT]-rmax[tt] are",
        "considered for [TT]-o[tt] and [TT]-no[tt] each frame.[PAR]",
        "[TT]-o[tt]: distribution of [MATH][COS][GRK]theta[grk][SUB]1[sub][cos][math] for "
        "rmin<=r<=rmax.[PAR]",
        "[TT]-no[tt]: distribution of [MATH][COS][GRK]theta[grk][SUB]2[sub][cos][math] for "
        "rmin<=r<=rmax.[PAR]",
        "[TT]-ro[tt]: [MATH][CHEVRON][COS][GRK]theta[grk][SUB]1[sub][cos][chevron][math] "
        "and [MATH][CHEVRON]3[COS]^2[GRK]theta[grk][SUB]2[sub][cos]-1[chevron][math] as a "
        "function of the",
        "distance.[PAR]",
        "[TT]-co[tt]: the sum over all solvent molecules within distance r",
        "of [MATH][COS][GRK]theta[grk][SUB]1[sub][cos][math] and "
        "[MATH]3[COS]^2([GRK]theta[grk][SUB]2[sub])-1[cos][math] as a function of r.[PAR]",
        "[TT]-rc[tt]: the distribution of the solvent molecules as a function of r"
    };

    gmx_output_env_t* oenv;
    static gmx_bool   bCom = FALSE, bVec23 = FALSE, bPBC = FALSE;
    static real       rmin = 0.0, rmax = 0.5, binwidth = 0.02, rbinw = 0.02;
    t_pargs           pa[] = {
        { "-com", FALSE, etBOOL, { &bCom }, "Use the center of mass as the reference position" },
        { "-v23", FALSE, etBOOL, { &bVec23 }, "Use the vector between atoms 2 and 3" },
        { "-rmin", FALSE, etREAL, { &rmin }, "Minimum distance (nm)" },
        { "-rmax", FALSE, etREAL, { &rmax }, "Maximum distance (nm)" },
        { "-cbin", FALSE, etREAL, { &binwidth }, "Binwidth for the cosine" },
        { "-rbin", FALSE, etREAL, { &rbinw }, "Binwidth for r (nm)" },
        { "-pbc",
          FALSE,
          etBOOL,
          { &bPBC },
          "Check PBC for the center of mass calculation. Only necessary when your reference group "
          "consists of several molecules." }
    };

    t_filenm fnm[] = { { efTRX, nullptr, nullptr, ffREAD },  { efTPS, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffOPTRD }, { efXVG, nullptr, "sori", ffWRITE },
                       { efXVG, "-no", "snor", ffWRITE },    { efXVG, "-ro", "sord", ffWRITE },
                       { efXVG, "-co", "scum", ffWRITE },    { efXVG, "-rc", "scount", ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    bTPS = (opt2bSet("-s", NFILE, fnm) || !opt2bSet("-n", NFILE, fnm) || bCom);
    if (bTPS)
    {
        read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, box, bCom);
    }

    /* get index groups */
    printf("Select a group of reference particles and a solvent group:\n");
    snew(grpname, 2);
    snew(index, 2);
    snew(isize, 2);
    if (bTPS)
    {
        get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 2, isize, index, grpname);
    }
    else
    {
        get_index(nullptr, ftp2fn(efNDX, NFILE, fnm), 2, isize, index, grpname);
    }

    if (bCom)
    {
        nrefgrp = 1;
        nrefat  = isize[0];
    }
    else
    {
        nrefgrp = isize[0];
        nrefat  = 1;
    }

    if (isize[1] % 3)
    {
        gmx_fatal(FARGS, "The number of solvent atoms (%d) is not a multiple of 3", isize[1]);
    }

    /* initialize reading trajectory:                         */
    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    rmin2 = gmx::square(rmin);
    rmax2 = gmx::square(rmax);
    rcut  = 0.99 * std::sqrt(max_cutoff2(guessPbcType(box), box));
    if (rcut == 0)
    {
        rcut = 10 * rmax;
    }
    rcut2 = gmx::square(rcut);

    invbw = 1 / binwidth;
    nbin1 = 1 + gmx::roundToInt(2 * invbw);
    nbin2 = 1 + gmx::roundToInt(invbw);

    invrbw = 1 / rbinw;

    snew(hist1, nbin1);
    snew(hist2, nbin2);
    nrbin = 1 + static_cast<int>(rcut / rbinw);
    if (nrbin == 0)
    {
        nrbin = 1;
    }
    snew(histi1, nrbin);
    snew(histi2, nrbin);
    snew(histn, nrbin);

    ntot = 0;
    nf   = 0;
    sum1 = 0;
    sum2 = 0;

    if (bTPS)
    {
        /* make molecules whole again */
        gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    }
    /* start analysis of trajectory */
    do
    {
        if (bTPS)
        {
            /* make molecules whole again */
            gmx_rmpbc_apply(gpbc, natoms, box, x);
        }

        set_pbc(&pbc, pbcType, box);
        n   = 0;
        inp = 0;
        for (p = 0; (p < nrefgrp); p++)
        {
            if (bCom)
            {
                calc_com_pbc(nrefat, &top, x, &pbc, index[0], xref, bPBC);
            }
            else
            {
                copy_rvec(x[index[0][p]], xref);
            }

            for (m = 0; m < isize[1]; m += 3)
            {
                sa0 = index[1][m];
                sa1 = index[1][m + 1];
                sa2 = index[1][m + 2];
                range_check(sa0, 0, natoms);
                range_check(sa1, 0, natoms);
                range_check(sa2, 0, natoms);
                pbc_dx(&pbc, x[sa0], xref, dx);
                r2 = norm2(dx);
                if (r2 < rcut2)
                {
                    r = std::sqrt(r2);
                    if (!bVec23)
                    {
                        /* Determine the normal to the plain */
                        rvec_sub(x[sa1], x[sa0], dxh1);
                        rvec_sub(x[sa2], x[sa0], dxh2);
                        rvec_inc(dxh1, dxh2);
                        svmul(1 / r, dx, dx);
                        unitv(dxh1, dxh1);
                        inp = iprod(dx, dxh1);
                        cprod(dxh1, dxh2, outer);
                        unitv(outer, outer);
                        outp = iprod(dx, outer);
                    }
                    else
                    {
                        /* Use the vector between the 2nd and 3rd atom */
                        rvec_sub(x[sa2], x[sa1], dxh2);
                        unitv(dxh2, dxh2);
                        outp = iprod(dx, dxh2) / r;
                    }
                    {
                        int ii = static_cast<int>(invrbw * r);
                        range_check(ii, 0, nrbin);
                        histi1[ii] += inp;
                        histi2[ii] += 3 * gmx::square(outp) - 1;
                        histn[ii]++;
                    }
                    if ((r2 >= rmin2) && (r2 < rmax2))
                    {
                        int ii1 = static_cast<int>(invbw * (inp + 1));
                        int ii2 = static_cast<int>(invbw * std::abs(outp));

                        range_check(ii1, 0, nbin1);
                        range_check(ii2, 0, nbin2);
                        hist1[ii1]++;
                        hist2[ii2]++;
                        sum1 += inp;
                        sum2 += outp;
                        n++;
                    }
                }
            }
        }
        ntot += n;
        nf++;

    } while (read_next_x(oenv, status, &t, x, box));

    /* clean up */
    sfree(x);
    close_trx(status);
    gmx_rmpbc_done(gpbc);

    /* Add the bin for the exact maximum to the previous bin */
    hist1[nbin1 - 1] += hist1[nbin1];
    hist2[nbin2 - 1] += hist2[nbin2];

    nav     = static_cast<real>(ntot) / (nrefgrp * nf);
    normfac = invbw / ntot;

    fprintf(stderr, "Average nr of molecules between %g and %g nm: %.1f\n", rmin, rmax, nav);
    if (ntot > 0)
    {
        sum1 /= ntot;
        sum2 /= ntot;
        fprintf(stderr, "Average cos(theta1)     between %g and %g nm: %6.3f\n", rmin, rmax, sum1);
        fprintf(stderr, "Average 3cos2(theta2)-1 between %g and %g nm: %6.3f\n", rmin, rmax, sum2);
    }

    sprintf(str, "Solvent orientation between %g and %g nm", rmin, rmax);
    fp = xvgropen(opt2fn("-o", NFILE, fnm), str, "cos(\\8q\\4\\s1\\N)", "", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"average shell size %.1f molecules\"\n", nav);
    }
    for (i = 0; i < nbin1; i++)
    {
        fprintf(fp, "%g %g\n", (i + 0.5) * binwidth - 1, 2 * normfac * hist1[i]);
    }
    xvgrclose(fp);

    sprintf(str, "Solvent normal orientation between %g and %g nm", rmin, rmax);
    fp = xvgropen(opt2fn("-no", NFILE, fnm), str, "cos(\\8q\\4\\s2\\N)", "", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"average shell size %.1f molecules\"\n", nav);
    }
    for (i = 0; i < nbin2; i++)
    {
        fprintf(fp, "%g %g\n", (i + 0.5) * binwidth, normfac * hist2[i]);
    }
    xvgrclose(fp);


    sprintf(str, "Solvent orientation");
    fp = xvgropen(opt2fn("-ro", NFILE, fnm), str, "r (nm)", "", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"as a function of distance\"\n");
    }
    xvgrLegend(fp, legr, oenv);
    for (i = 0; i < nrbin; i++)
    {
        fprintf(fp,
                "%g %g %g\n",
                (i + 0.5) * rbinw,
                histn[i] ? histi1[i] / histn[i] : 0,
                histn[i] ? histi2[i] / histn[i] : 0);
    }
    xvgrclose(fp);

    sprintf(str, "Cumulative solvent orientation");
    fp = xvgropen(opt2fn("-co", NFILE, fnm), str, "r (nm)", "", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"as a function of distance\"\n");
    }
    xvgrLegend(fp, legc, oenv);
    normfac = 1.0 / (nrefgrp * nf);
    c1      = 0;
    c2      = 0;
    fprintf(fp, "%g %g %g\n", 0.0, c1, c2);
    for (i = 0; i < nrbin; i++)
    {
        c1 += histi1[i] * normfac;
        c2 += histi2[i] * normfac;
        fprintf(fp, "%g %g %g\n", (i + 1) * rbinw, c1, c2);
    }
    xvgrclose(fp);

    sprintf(str, "Solvent distribution");
    fp = xvgropen(opt2fn("-rc", NFILE, fnm), str, "r (nm)", "molecules/nm", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"as a function of distance\"\n");
    }
    normfac = 1.0 / (rbinw * nf);
    for (i = 0; i < nrbin; i++)
    {
        fprintf(fp, "%g %g\n", (i + 0.5) * rbinw, histn[i] * normfac);
    }
    xvgrclose(fp);

    do_view(oenv, opt2fn("-o", NFILE, fnm), nullptr);
    do_view(oenv, opt2fn("-no", NFILE, fnm), nullptr);
    do_view(oenv, opt2fn("-ro", NFILE, fnm), "-nxy");
    do_view(oenv, opt2fn("-co", NFILE, fnm), "-nxy");

    return 0;
}
