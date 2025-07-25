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

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/angle_correction.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

struct gmx_output_env_t;

static void dump_dih_trr(int nframes, int nangles, real** dih, const char* fn, real* time)
{
    int              i, j, k, l, m, na;
    struct t_fileio* fio;
    rvec*            x;
    matrix           box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };

    na = (nangles * 2);
    if ((na % 3) != 0)
    {
        na = 1 + na / 3;
    }
    else
    {
        na = na / 3;
    }
    printf("There are %d dihedrals. Will fill %d atom positions with cos/sin\n", nangles, na);
    snew(x, na);
    fio = gmx_trr_open(fn, "w");
    for (i = 0; (i < nframes); i++)
    {
        k = l = 0;
        for (j = 0; (j < nangles); j++)
        {
            for (m = 0; (m < 2); m++)
            {
                // This is just because the compler and static-analyzer cannot
                // know that dih[j][i] is always valid. Since it occurs in the innermost
                // loop over angles and will only trigger on coding errors, we
                // only enable it for debug builds.
                GMX_ASSERT(dih != nullptr && dih[j] != nullptr, "Incorrect dihedral array data");
                x[k][l] = (m == 0) ? std::cos(dih[j][i]) : std::sin(dih[j][i]);
                l++;
                if (l == DIM)
                {
                    l = 0;
                    k++;
                }
            }
        }
        gmx_trr_write_frame(fio, i, time[i], 0, box, na, x, nullptr, nullptr);
    }
    gmx_trr_close(fio);
    sfree(x);
}

int gmx_g_angle(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] computes the angle distribution for a number of angles",
        "or dihedrals.[PAR]",
        "With option [TT]-ov[tt], you can plot the average angle of",
        "a group of angles as a function of time. With the [TT]-all[tt] option,",
        "the first graph is the average and the rest are the individual angles.[PAR]",
        "With the [TT]-of[tt] option, [THISMODULE] also calculates the fraction of trans",
        "dihedrals (only for dihedrals) as function of time, but this is",
        "probably only fun for a select few.[PAR]",
        "With option [TT]-oc[tt], a dihedral correlation function is calculated.[PAR]",
        "It should be noted that the index file must contain",
        "atom triplets for angles or atom quadruplets for dihedrals.",
        "If this is not the case, the program will crash.[PAR]",
        "With option [TT]-or[tt], a trajectory file is dumped containing cos and",
        "sin of selected dihedral angles, which subsequently can be used as",
        "input for a principal components analysis using [gmx-covar].[PAR]",
        "Option [TT]-ot[tt] plots when transitions occur between",
        "dihedral rotamers of multiplicity 3 and [TT]-oh[tt]",
        "records a histogram of the times between such transitions,",
        "assuming the input trajectory frames are equally spaced in time."
    };
    static const char* opt[] = { nullptr, "angle", "dihedral", "improper", "ryckaert-bellemans",
                                 nullptr };
    static gmx_bool    bALL = FALSE, bChandler = FALSE, bAverCorr = FALSE, bPBC = TRUE;
    static real        binwidth = 1;
    t_pargs            pa[]     = {
        { "-type", FALSE, etENUM, { opt }, "Type of angle to analyse" },
        { "-all",
                         FALSE,
                         etBOOL,
                         { &bALL },
                         "Plot all angles separately in the averages file, in the order of appearance in the "
                                        "index file." },
        { "-binwidth",
                         FALSE,
                         etREAL,
                         { &binwidth },
                         "binwidth (degrees) for calculating the distribution" },
        { "-periodic", FALSE, etBOOL, { &bPBC }, "Print dihedral angles modulo 360 degrees" },
        { "-chandler",
                         FALSE,
                         etBOOL,
                         { &bChandler },
                         "Use Chandler correlation function (N[trans] = 1, N[gauche] = 0) rather than cosine "
                                        "correlation function. Trans is defined as phi < -60 or phi > 60." },
        { "-avercorr",
                         FALSE,
                         etBOOL,
                         { &bAverCorr },
                         "Average the correlation functions for the individual angles/dihedrals" }
    };
    static const char* bugs[] = {
        "Counting transitions only works for dihedrals with multiplicity 3"
    };

    FILE*         out;
    real          dt;
    int           isize;
    int*          index;
    char*         grpname;
    real          maxang, S2, norm_fac, maxstat;
    unsigned long mode;
    int           nframes, maxangstat, mult, *angstat;
    int           i, j, nangles, first, last;
    gmx_bool      bAver, bRb, bPeriodic, bFrac, /* calculate fraction too?  */
            bTrans,                             /* worry about transtions too? */
            bCorr;                              /* correlation function ? */
    double   tfrac = 0;
    char     title[256];
    real**   dih = nullptr; /* mega array with all dih. angles at all times*/
    real *   time, *trans_frac, *aver_angle;
    t_filenm fnm[] = { { efTRX, "-f", nullptr, ffREAD },     { efNDX, nullptr, "angle", ffREAD },
                       { efXVG, "-od", "angdist", ffWRITE }, { efXVG, "-ov", "angaver", ffOPTWR },
                       { efXVG, "-of", "dihfrac", ffOPTWR }, { efXVG, "-ot", "dihtrans", ffOPTWR },
                       { efXVG, "-oh", "trhisto", ffOPTWR }, { efXVG, "-oc", "dihcorr", ffOPTWR },
                       { efTRR, "-or", nullptr, ffOPTWR } };
#define NFILE asize(fnm)
    int               npargs;
    t_pargs*          ppa;
    gmx_output_env_t* oenv;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, npargs, ppa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    mult   = 4;
    maxang = 360.0;
    bRb    = FALSE;

    GMX_RELEASE_ASSERT(opt[0] != nullptr,
                       "Internal option inconsistency; opt[0]==NULL after processing");

    switch (opt[0][0])
    {
        case 'a':
            mult   = 3;
            maxang = 180.0;
            break;
        case 'd': // Intended fall through
        case 'i': break;
        case 'r': bRb = TRUE; break;
    }

    if (opt2bSet("-or", NFILE, fnm))
    {
        if (mult != 4)
        {
            gmx_fatal(FARGS, "Can not combine angles with trr dump");
        }
        else
        {
            please_cite(stdout, "Mu2005a");
        }
    }

    /* Calculate bin size */
    maxangstat = gmx::roundToInt(maxang / binwidth);
    binwidth   = maxang / maxangstat;

    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
    nangles = isize / mult;
    if ((isize % mult) != 0)
    {
        gmx_fatal(FARGS,
                  "number of index elements not multiple of %d, "
                  "these can not be %s\n",
                  mult,
                  (mult == 3) ? "angle triplets" : "dihedral quadruplets");
    }


    /* Check whether specific analysis has to be performed */
    bCorr  = opt2bSet("-oc", NFILE, fnm);
    bAver  = opt2bSet("-ov", NFILE, fnm);
    bTrans = opt2bSet("-ot", NFILE, fnm);
    bFrac  = opt2bSet("-of", NFILE, fnm);
    if (bTrans && opt[0][0] != 'd')
    {
        fprintf(stderr, "Option -ot should only accompany -type dihedral. Disabling -ot.\n");
        bTrans = FALSE;
    }

    if (bChandler && !bCorr)
    {
        bCorr = TRUE;
    }

    if (bFrac && !bRb)
    {
        fprintf(stderr,
                "Warning:"
                " calculating fractions as defined in this program\n"
                "makes sense for Ryckaert Bellemans dihs. only. Ignoring -of\n\n");
        bFrac = FALSE;
    }

    if ((bTrans || bFrac || bCorr) && mult == 3)
    {
        gmx_fatal(FARGS,
                  "Can only do transition, fraction or correlation\n"
                  "on dihedrals. Select -d\n");
    }

    /*
     * We need to know the nr of frames so we can allocate memory for an array
     * with all dihedral angles at all timesteps. Works for me.
     */
    if (bTrans || bCorr || bALL || opt2bSet("-or", NFILE, fnm))
    {
        snew(dih, nangles);
    }

    snew(angstat, maxangstat);

    read_ang_dih(ftp2fn(efTRX, NFILE, fnm),
                 (mult == 3),
                 bALL || bCorr || bTrans || opt2bSet("-or", NFILE, fnm),
                 bRb,
                 bPBC,
                 maxangstat,
                 angstat,
                 &nframes,
                 &time,
                 isize,
                 index,
                 &trans_frac,
                 &aver_angle,
                 dih,
                 oenv);

    dt = (time[nframes - 1] - time[0]) / (nframes - 1);

    if (bAver)
    {
        sprintf(title, "Average Angle: %s", grpname);
        out = xvgropen(opt2fn("-ov", NFILE, fnm), title, "Time (ps)", "Angle (degrees)", oenv);
        for (i = 0; (i < nframes); i++)
        {
            fprintf(out, "%10.5f  %8.3f", time[i], aver_angle[i] * gmx::c_rad2Deg);
            if (bALL)
            {
                for (j = 0; (j < nangles); j++)
                {
                    if (bPBC)
                    {
                        real dd = dih[j][i];
                        fprintf(out, "  %8.3f", std::atan2(std::sin(dd), std::cos(dd)) * gmx::c_rad2Deg);
                    }
                    else
                    {
                        fprintf(out, "  %8.3f", dih[j][i] * gmx::c_rad2Deg);
                    }
                }
            }
            fprintf(out, "\n");
        }
        xvgrclose(out);
    }
    if (opt2bSet("-or", NFILE, fnm))
    {
        dump_dih_trr(nframes, nangles, dih, opt2fn("-or", NFILE, fnm), time);
    }

    if (bFrac)
    {
        sprintf(title, "Trans fraction: %s", grpname);
        out   = xvgropen(opt2fn("-of", NFILE, fnm), title, "Time (ps)", "Fraction", oenv);
        tfrac = 0.0;
        for (i = 0; (i < nframes); i++)
        {
            fprintf(out, "%10.5f  %10.3f\n", time[i], trans_frac[i]);
            tfrac += trans_frac[i];
        }
        xvgrclose(out);

        tfrac /= nframes;
        fprintf(stderr, "Average trans fraction: %g\n", tfrac);
    }
    sfree(trans_frac);

    if (bTrans)
    {
        ana_dih_trans(
                opt2fn("-ot", NFILE, fnm), opt2fn("-oh", NFILE, fnm), dih, nframes, nangles, grpname, time, bRb, oenv);
    }

    if (bCorr)
    {
        /* Autocorrelation function */
        if (nframes < 2)
        {
            fprintf(stderr, "Not enough frames for correlation function\n");
        }
        else
        {

            if (bChandler)
            {
                real     dval, sixty = gmx::c_deg2Rad * 60;
                gmx_bool bTest;

                for (i = 0; (i < nangles); i++)
                {
                    for (j = 0; (j < nframes); j++)
                    {
                        dval = dih[i][j];
                        if (bRb)
                        {
                            bTest = (dval > -sixty) && (dval < sixty);
                        }
                        else
                        {
                            bTest = (dval < -sixty) || (dval > sixty);
                        }
                        if (bTest)
                        {
                            dih[i][j] = dval - tfrac;
                        }
                        else
                        {
                            dih[i][j] = -tfrac;
                        }
                    }
                }
            }
            if (bChandler)
            {
                mode = eacNormal;
            }
            else
            {
                mode = eacCos;
            }
            do_autocorr(opt2fn("-oc", NFILE, fnm),
                        oenv,
                        "Dihedral Autocorrelation Function",
                        nframes,
                        nangles,
                        dih,
                        dt,
                        mode,
                        bAverCorr);
        }
    }


    /* Determine the non-zero part of the distribution */
    for (first = 0; (first < maxangstat - 1) && (angstat[first + 1] == 0); first++) {}
    for (last = maxangstat - 1; (last > 0) && (angstat[last - 1] == 0); last--) {}

    double aver = 0;
    printf("Found points in the range from %d to %d (max %d)\n", first, last, maxangstat);
    if (bTrans || bCorr || bALL || opt2bSet("-or", NFILE, fnm))
    { /* It's better to re-calculate Std. Dev per sample */
        real b_aver = aver_angle[0];
        real b      = dih[0][0];
        real delta;
        for (int i = 0; (i < nframes); i++)
        {
            delta = correctRadianAngleRange(aver_angle[i] - b_aver);
            b_aver += delta;
            aver += b_aver;
            for (int j = 0; (j < nangles); j++)
            {
                delta = correctRadianAngleRange(dih[j][i] - b);
                b += delta;
            }
        }
    }
    else
    { /* Incorrect  for Std. Dev. */
        real delta, b_aver = aver_angle[0];
        for (i = 0; (i < nframes); i++)
        {
            delta = correctRadianAngleRange(aver_angle[i] - b_aver);
            b_aver += delta;
            aver += b_aver;
        }
    }
    aver /= nframes;
    double aversig = correctRadianAngleRange(aver);
    aversig *= gmx::c_rad2Deg;
    aver *= gmx::c_rad2Deg;
    printf(" < angle >  = %g\n", aversig);

    if (mult == 3)
    {
        sprintf(title, "Angle Distribution: %s", grpname);
    }
    else
    {
        sprintf(title, "Dihedral Distribution: %s", grpname);

        calc_distribution_props(maxangstat, angstat, -180.0, 0, nullptr, &S2);
        fprintf(stderr, "Order parameter S^2 = %g\n", S2);
    }

    bPeriodic = (mult == 4) && (first == 0) && (last == maxangstat - 1);

    out = xvgropen(opt2fn("-od", NFILE, fnm), title, "Degrees", "", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@    subtitle \"average angle: %g\\So\\N\"\n", aver);
    }
    norm_fac = 1.0 / (nangles * nframes * binwidth);
    if (bPeriodic)
    {
        maxstat = 0;
        for (i = first; (i <= last); i++)
        {
            maxstat = std::max(maxstat, angstat[i] * norm_fac);
        }
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(out, "@with g0\n");
            fprintf(out, "@    world xmin -180\n");
            fprintf(out, "@    world xmax  180\n");
            fprintf(out, "@    world ymin 0\n");
            fprintf(out, "@    world ymax %g\n", maxstat * 1.05);
            fprintf(out, "@    xaxis  tick major 60\n");
            fprintf(out, "@    xaxis  tick minor 30\n");
            fprintf(out, "@    yaxis  tick major 0.005\n");
            fprintf(out, "@    yaxis  tick minor 0.0025\n");
        }
    }
    for (i = first; (i <= last); i++)
    {
        fprintf(out, "%10g  %10f\n", i * binwidth + 180.0 - maxang, angstat[i] * norm_fac);
    }
    if (bPeriodic)
    {
        /* print first bin again as last one */
        fprintf(out, "%10g  %10f\n", 180.0, angstat[0] * norm_fac);
    }

    xvgrclose(out);

    do_view(oenv, opt2fn("-od", NFILE, fnm), "-nxy");
    if (bAver)
    {
        do_view(oenv, opt2fn("-ov", NFILE, fnm), "-nxy");
    }

    return 0;
}
