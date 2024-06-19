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
#include <cstring>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;


int gmx_vanhove(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes the Van Hove correlation function.",
        "The Van Hove G(r,t) is the probability that a particle that is at r[SUB]0[sub]",
        "at time zero can be found at position r[SUB]0[sub]+r at time t.",
        "[THISMODULE] determines G not for a vector r, but for the length of r.",
        "Thus it gives the probability that a particle moves a distance of r",
        "in time t.",
        "Jumps across the periodic boundaries are removed.",
        "Corrections are made for scaling due to isotropic",
        "or anisotropic pressure coupling.",
        "[PAR]",
        "With option [TT]-om[tt] the whole matrix can be written as a function",
        "of t and r or as a function of [SQRT]t[sqrt] and r (option [TT]-sqrt[tt]).",
        "[PAR]",
        "With option [TT]-or[tt] the Van Hove function is plotted for one",
        "or more values of t. Option [TT]-nr[tt] sets the number of times,",
        "option [TT]-fr[tt] the number spacing between the times.",
        "The binwidth is set with option [TT]-rbin[tt]. The number of bins",
        "is determined automatically.",
        "[PAR]",
        "With option [TT]-ot[tt] the integral up to a certain distance",
        "(option [TT]-rt[tt]) is plotted as a function of time.",
        "[PAR]",
        "For all frames that are read the coordinates of the selected particles",
        "are stored in memory. Therefore the program may use a lot of memory.",
        "For options [TT]-om[tt] and [TT]-ot[tt] the program may be slow.",
        "This is because the calculation scales as the number of frames times",
        "[TT]-fm[tt] or [TT]-ft[tt].",
        "Note that with the [TT]-dt[tt] option the memory usage and calculation",
        "time can be reduced."
    };
    static int  fmmax = 0, ftmax = 0, nlev = 81, nr = 1, fshift = 0;
    static real sbin = 0, rmax = 2, rbin = 0.01, mmax = 0, rint = 0;
    t_pargs     pa[] = {
        { "-sqrt",
          FALSE,
          etREAL,
          { &sbin },
          "Use [SQRT]t[sqrt] on the matrix axis which binspacing # in [SQRT]ps[sqrt]" },
        { "-fm", FALSE, etINT, { &fmmax }, "Number of frames in the matrix, 0 is plot all" },
        { "-rmax", FALSE, etREAL, { &rmax }, "Maximum r in the matrix (nm)" },
        { "-rbin", FALSE, etREAL, { &rbin }, "Binwidth in the matrix and for [TT]-or[tt] (nm)" },
        { "-mmax",
          FALSE,
          etREAL,
          { &mmax },
          "Maximum density in the matrix, 0 is calculate (1/nm)" },
        { "-nlevels", FALSE, etINT, { &nlev }, "Number of levels in the matrix" },
        { "-nr", FALSE, etINT, { &nr }, "Number of curves for the [TT]-or[tt] output" },
        { "-fr", FALSE, etINT, { &fshift }, "Frame spacing for the [TT]-or[tt] output" },
        { "-rt", FALSE, etREAL, { &rint }, "Integration limit for the [TT]-ot[tt] output (nm)" },
        { "-ft",
          FALSE,
          etINT,
          { &ftmax },
          "Number of frames in the [TT]-ot[tt] output, 0 is plot all" }
    };
#define NPA asize(pa)

    t_filenm fnm[] = {
        { efTRX, nullptr, nullptr, ffREAD },    { efTPS, nullptr, nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffOPTRD },   { efXPM, "-om", "vanhove", ffOPTWR },
        { efXVG, "-or", "vanhove_r", ffOPTWR }, { efXVG, "-ot", "vanhove_t", ffOPTWR }
    };
#define NFILE asize(fnm)

    gmx_output_env_t*        oenv;
    const char *             matfile, *otfile, *orfile;
    t_topology               top;
    PbcType                  pbcType;
    matrix                   boxtop, box, *sbox, avbox, corr;
    rvec *                   xtop, *x, **sx;
    int                      isize, nalloc, nallocn;
    t_trxstatus*             status;
    int*                     index;
    char*                    grpname;
    int                      nfr, f, ff, i, m, mat_nx = 0, nbin = 0, bin, mbin, fbin;
    real *                   time, t, invbin = 0, rmax2 = 0, rint2 = 0, d2;
    real                     invsbin = 0, matmax, normfac, dt, *tickx, *ticky;
    std::vector<std::string> legend;
    real**                   mat = nullptr;
    int * pt = nullptr, **pr = nullptr, *mcount = nullptr, *tcount = nullptr, *rcount = nullptr;
    FILE* fp;
    t_rgb rlo = { 1, 1, 1 }, rhi = { 0, 0, 0 };

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    matfile = opt2fn_null("-om", NFILE, fnm);
    if (opt2parg_bSet("-fr", NPA, pa))
    {
        orfile = opt2fn("-or", NFILE, fnm);
    }
    else
    {
        orfile = opt2fn_null("-or", NFILE, fnm);
    }
    if (opt2parg_bSet("-rt", NPA, pa))
    {
        otfile = opt2fn("-ot", NFILE, fnm);
    }
    else
    {
        otfile = opt2fn_null("-ot", NFILE, fnm);
    }

    if (!matfile && !otfile && !orfile)
    {
        fprintf(stderr, "For output set one (or more) of the output file options\n");
        exit(0);
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, boxtop, FALSE);
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);

    nalloc = 0;
    time   = nullptr;
    sbox   = nullptr;
    sx     = nullptr;
    clear_mat(avbox);

    read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    nfr = 0;
    do
    {
        if (nfr >= nalloc)
        {
            nalloc += 100;
            srenew(time, nalloc);
            srenew(sbox, nalloc);
            srenew(sx, nalloc);
        }
        GMX_RELEASE_ASSERT(time != nullptr, "Memory allocation failure; time array is NULL");
        GMX_RELEASE_ASSERT(sbox != nullptr, "Memory allocation failure; sbox array is NULL");

        time[nfr] = t;
        copy_mat(box, sbox[nfr]);
        /* This assumes that the off-diagonal box elements
         * are not affected by jumps across the periodic boundaries.
         */
        m_add(avbox, box, avbox);
        snew(sx[nfr], isize);
        for (i = 0; i < isize; i++)
        {
            copy_rvec(x[index[i]], sx[nfr][i]);
        }

        nfr++;
    } while (read_next_x(oenv, status, &t, x, box));

    /* clean up */
    sfree(x);
    close_trx(status);

    fprintf(stderr, "Read %d frames\n", nfr);

    dt = (time[nfr - 1] - time[0]) / (nfr - 1);
    /* Some ugly rounding to get nice nice times in the output */
    dt = std::round(10000.0 * dt) / 10000.0;

    invbin = 1.0 / rbin;

    if (matfile)
    {
        if (fmmax <= 0 || fmmax >= nfr)
        {
            fmmax = nfr - 1;
        }
        snew(mcount, fmmax);
        nbin = gmx::roundToInt(rmax * invbin);
        if (sbin == 0)
        {
            mat_nx = fmmax + 1;
        }
        else
        {
            invsbin = 1.0 / sbin;
            mat_nx  = static_cast<int>(std::sqrt(fmmax * dt) * invsbin + 1);
        }
        snew(mat, mat_nx);
        for (f = 0; f < mat_nx; f++)
        {
            snew(mat[f], nbin);
        }
        rmax2 = gmx::square(nbin * rbin);
        /* Initialize time zero */
        mat[0][0] = nfr * isize;
        mcount[0] += nfr;
    }
    else
    {
        fmmax = 0;
    }

    if (orfile)
    {
        snew(pr, nr);
        nalloc = 0;
        snew(rcount, nr);
    }

    if (otfile)
    {
        if (ftmax <= 0)
        {
            ftmax = nfr - 1;
        }
        snew(tcount, ftmax);
        snew(pt, nfr);
        rint2 = rint * rint;
        /* Initialize time zero */
        pt[0] = nfr * isize;
        tcount[0] += nfr;
    }
    else
    {
        ftmax = 0;
    }

    msmul(avbox, 1.0 / nfr, avbox);
    for (f = 0; f < nfr; f++)
    {
        if (f % 100 == 0)
        {
            fprintf(stderr, "\rProcessing frame %d", f);
            fflush(stderr);
        }
        if (pbcType != PbcType::No)
        {
            /* Scale all the configuration to the average box */
            gmx::invertBoxMatrix(sbox[f], corr);
            mmul_ur0(avbox, corr, corr);
            for (i = 0; i < isize; i++)
            {
                mvmul_ur0(corr, sx[f][i], sx[f][i]);
                if (f > 0)
                {
                    /* Correct for periodic jumps */
                    for (m = DIM - 1; m >= 0; m--)
                    {
                        while (sx[f][i][m] - sx[f - 1][i][m] > 0.5 * avbox[m][m])
                        {
                            rvec_dec(sx[f][i], avbox[m]);
                        }
                        while (sx[f][i][m] - sx[f - 1][i][m] <= -0.5 * avbox[m][m])
                        {
                            rvec_inc(sx[f][i], avbox[m]);
                        }
                    }
                }
            }
        }
        for (ff = 0; ff < f; ff++)
        {
            fbin = f - ff;
            if (fbin <= fmmax || fbin <= ftmax)
            {
                if (sbin == 0)
                {
                    mbin = fbin;
                }
                else
                {
                    mbin = gmx::roundToInt(std::sqrt(fbin * dt) * invsbin);
                }
                for (i = 0; i < isize; i++)
                {
                    d2 = distance2(sx[f][i], sx[ff][i]);
                    if (mbin < mat_nx && d2 < rmax2)
                    {
                        bin = gmx::roundToInt(std::sqrt(d2) * invbin);
                        if (bin < nbin)
                        {
                            mat[mbin][bin] += 1;
                        }
                    }
                    if (fbin <= ftmax && d2 <= rint2)
                    {
                        pt[fbin]++;
                    }
                }
                if (matfile)
                {
                    mcount[mbin]++;
                }
                if (otfile)
                {
                    tcount[fbin]++;
                }
            }
        }
        if (orfile)
        {
            for (fbin = 0; fbin < nr; fbin++)
            {
                ff = f - (fbin + 1) * fshift;
                if (ff >= 0)
                {
                    for (i = 0; i < isize; i++)
                    {
                        d2  = distance2(sx[f][i], sx[ff][i]);
                        bin = gmx::roundToInt(std::sqrt(d2) * invbin);
                        if (bin >= nalloc)
                        {
                            nallocn = 10 * (bin / 10) + 11;
                            for (m = 0; m < nr; m++)
                            {
                                srenew(pr[m], nallocn);
                                for (i = nalloc; i < nallocn; i++)
                                {
                                    pr[m][i] = 0;
                                }
                            }
                            nalloc = nallocn;
                        }
                        pr[fbin][bin]++;
                    }
                    rcount[fbin]++;
                }
            }
        }
    }
    fprintf(stderr, "\n");

    if (matfile)
    {
        matmax = 0;
        for (f = 0; f < mat_nx; f++)
        {
            normfac = 1.0 / (mcount[f] * isize * rbin);
            for (i = 0; i < nbin; i++)
            {
                mat[f][i] *= normfac;
                if (mat[f][i] > matmax && (f != 0 || i != 0))
                {
                    matmax = mat[f][i];
                }
            }
        }
        fprintf(stdout, "Value at (0,0): %.3f, maximum of the rest %.3f\n", mat[0][0], matmax);
        if (mmax > 0)
        {
            matmax = mmax;
        }
        snew(tickx, mat_nx);
        for (f = 0; f < mat_nx; f++)
        {
            if (sbin == 0)
            {
                tickx[f] = f * dt;
            }
            else
            {
                tickx[f] = f * sbin;
            }
        }
        snew(ticky, nbin + 1);
        for (i = 0; i <= nbin; i++)
        {
            ticky[i] = i * rbin;
        }
        fp = gmx_ffopen(matfile, "w");
        write_xpm(fp,
                  MAT_SPATIAL_Y,
                  "Van Hove function",
                  "G (1/nm)",
                  sbin == 0 ? "time (ps)" : "sqrt(time) (ps^1/2)",
                  "r (nm)",
                  mat_nx,
                  nbin,
                  tickx,
                  ticky,
                  mat,
                  0,
                  matmax,
                  rlo,
                  rhi,
                  &nlev);
        gmx_ffclose(fp);
    }

    if (orfile)
    {
        fp = xvgropen(orfile, "Van Hove function", "r (nm)", "G (nm\\S-1\\N)", oenv);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fp, "@ subtitle \"for particles in group %s\"\n", grpname);
        }
        for (fbin = 0; fbin < nr; fbin++)
        {
            legend.emplace_back(gmx::formatString("%g ps", (fbin + 1) * fshift * dt));
        }
        xvgrLegend(fp, legend, oenv);
        for (i = 0; i < nalloc; i++)
        {
            fprintf(fp, "%g", i * rbin);
            for (fbin = 0; fbin < nr; fbin++)
            {
                fprintf(fp,
                        " %g",
                        static_cast<real>(pr[fbin][i]
                                          / (rcount[fbin] * isize * rbin * (i == 0 ? 0.5 : 1.0))));
            }
            fprintf(fp, "\n");
        }
        xvgrclose(fp);
    }

    if (otfile)
    {
        auto buf = gmx::formatString("Probability of moving less than %g nm", rint);
        fp       = xvgropen(otfile, buf.c_str(), "t (ps)", "", oenv);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fp, "@ subtitle \"for particles in group %s\"\n", grpname);
        }
        for (f = 0; f <= ftmax; f++)
        {
            fprintf(fp, "%g %g\n", f * dt, static_cast<real>(pt[f]) / (tcount[f] * isize));
        }
        xvgrclose(fp);
    }

    do_view(oenv, matfile, nullptr);
    do_view(oenv, orfile, nullptr);
    do_view(oenv, otfile, nullptr);

    return 0;
}
