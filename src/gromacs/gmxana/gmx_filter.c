/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"

int gmx_filter(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] performs frequency filtering on a trajectory.",
        "The filter shape is cos([GRK]pi[grk] t/A) + 1 from -A to +A, where A is given",
        "by the option [TT]-nf[tt] times the time step in the input trajectory.",
        "This filter reduces fluctuations with period A by 85%, with period",
        "2*A by 50% and with period 3*A by 17% for low-pass filtering.",
        "Both a low-pass and high-pass filtered trajectory can be written.[PAR]",

        "Option [TT]-ol[tt] writes a low-pass filtered trajectory.",
        "A frame is written every [TT]-nf[tt] input frames.",
        "This ratio of filter length and output interval ensures a good",
        "suppression of aliasing of high-frequency motion, which is useful for",
        "making smooth movies. Also averages of properties which are linear",
        "in the coordinates are preserved, since all input frames are weighted",
        "equally in the output.",
        "When all frames are needed, use the [TT]-all[tt] option.[PAR]",

        "Option [TT]-oh[tt] writes a high-pass filtered trajectory.",
        "The high-pass filtered coordinates are added to the coordinates",
        "from the structure file. When using high-pass filtering use [TT]-fit[tt]",
        "or make sure you use a trajectory that has been fitted on",
        "the coordinates in the structure file."
    };

    static int      nf      = 10;
    static gmx_bool bNoJump = TRUE, bFit = FALSE, bLowAll = FALSE;
    t_pargs         pa[]    = {
        { "-nf", FALSE, etINT, {&nf},
          "Sets the filter length as well as the output interval for low-pass filtering" },
        { "-all", FALSE, etBOOL, {&bLowAll},
          "Write all low-pass filtered frames" },
        { "-nojump", FALSE, etBOOL, {&bNoJump},
          "Remove jumps of atoms across the box" },
        { "-fit", FALSE, etBOOL, {&bFit},
          "Fit all frames to a reference structure" }
    };
    const char     *topfile, *lowfile, *highfile;
    gmx_bool        bTop = FALSE;
    t_topology      top;
    int             ePBC = -1;
    rvec           *xtop;
    matrix          topbox, *box, boxf;
    char            title[256], *grpname;
    int             isize;
    atom_id        *index;
    real           *w_rls = NULL;
    t_trxstatus    *in;
    t_trxstatus    *outl, *outh;
    int             nffr, i, fr, nat, j, d, m;
    atom_id        *ind;
    real            flen, *filt, sum, *t;
    rvec            xcmtop, xcm, **x, *ptr, *xf, *xn, *xp, hbox;
    output_env_t    oenv;
    gmx_rmpbc_t     gpbc = NULL;

#define NLEG asize(leg)
    t_filenm fnm[] = {
        { efTRX, "-f", NULL, ffREAD  },
        { efTPS, NULL, NULL, ffOPTRD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efTRO, "-ol", "lowpass",  ffOPTWR },
        { efTRO, "-oh", "highpass", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    highfile = opt2fn_null("-oh", NFILE, fnm);
    if (highfile)
    {
        topfile = ftp2fn(efTPS, NFILE, fnm);
        lowfile = opt2fn_null("-ol", NFILE, fnm);
    }
    else
    {
        topfile = ftp2fn_null(efTPS, NFILE, fnm);
        lowfile = opt2fn("-ol", NFILE, fnm);
    }
    if (topfile)
    {
        bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC,
                             &xtop, NULL, topbox, TRUE);
        if (bTop)
        {
            gpbc = gmx_rmpbc_init(&top.idef, ePBC, top.atoms.nr);
            gmx_rmpbc(gpbc, top.atoms.nr, topbox, xtop);
        }
    }

    clear_rvec(xcmtop);
    if (bFit)
    {
        fprintf(stderr, "Select group for least squares fit\n");
        get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
        /* Set the weight */
        snew(w_rls, top.atoms.nr);
        for (i = 0; i < isize; i++)
        {
            w_rls[index[i]] = top.atoms.atom[index[i]].m;
        }
        calc_xcm(xtop, isize, index, top.atoms.atom, xcmtop, FALSE);
        for (j = 0; j < top.atoms.nr; j++)
        {
            rvec_dec(xtop[j], xcmtop);
        }
    }

    /* The actual filter length flen can actually be any real number */
    flen = 2*nf;
    /* nffr is the number of frames that we filter over */
    nffr = 2*nf - 1;
    snew(filt, nffr);
    sum = 0;
    for (i = 0; i < nffr; i++)
    {
        filt[i] = cos(2*M_PI*(i - nf + 1)/(real)flen) + 1;
        sum    += filt[i];
    }
    fprintf(stdout, "filter weights:");
    for (i = 0; i < nffr; i++)
    {
        filt[i] /= sum;
        fprintf(stdout, " %5.3f", filt[i]);
    }
    fprintf(stdout, "\n");

    snew(t, nffr);
    snew(x, nffr);
    snew(box, nffr);

    nat = read_first_x(oenv, &in, opt2fn("-f", NFILE, fnm),
                       &(t[nffr - 1]), &(x[nffr - 1]), box[nffr - 1]);
    snew(ind, nat);
    for (i = 0; i < nat; i++)
    {
        ind[i] = i;
    }
    /* x[nffr - 1] was already allocated by read_first_x */
    for (i = 0; i < nffr-1; i++)
    {
        snew(x[i], nat);
    }
    snew(xf, nat);
    if (lowfile)
    {
        outl = open_trx(lowfile, "w");
    }
    else
    {
        outl = 0;
    }
    if (highfile)
    {
        outh = open_trx(highfile, "w");
    }
    else
    {
        outh = 0;
    }

    fr = 0;
    do
    {
        xn = x[nffr - 1];
        if (bNoJump && fr > 0)
        {
            xp = x[nffr - 2];
            for (j = 0; j < nat; j++)
            {
                for (d = 0; d < DIM; d++)
                {
                    hbox[d] = 0.5*box[nffr - 1][d][d];
                }
            }
            for (i = 0; i < nat; i++)
            {
                for (m = DIM-1; m >= 0; m--)
                {
                    if (hbox[m] > 0)
                    {
                        while (xn[i][m] - xp[i][m] <= -hbox[m])
                        {
                            for (d = 0; d <= m; d++)
                            {
                                xn[i][d] += box[nffr - 1][m][d];
                            }
                        }
                        while (xn[i][m] - xp[i][m] > hbox[m])
                        {
                            for (d = 0; d <= m; d++)
                            {
                                xn[i][d] -= box[nffr - 1][m][d];
                            }
                        }
                    }
                }
            }
        }
        if (bTop)
        {
            gmx_rmpbc(gpbc, nat, box[nffr - 1], xn);
        }
        if (bFit)
        {
            calc_xcm(xn, isize, index, top.atoms.atom, xcm, FALSE);
            for (j = 0; j < nat; j++)
            {
                rvec_dec(xn[j], xcm);
            }
            do_fit(nat, w_rls, xtop, xn);
            for (j = 0; j < nat; j++)
            {
                rvec_inc(xn[j], xcmtop);
            }
        }
        if (fr >= nffr && (outh || bLowAll || fr % nf == nf - 1))
        {
            /* Lowpass filtering */
            for (j = 0; j < nat; j++)
            {
                clear_rvec(xf[j]);
            }
            clear_mat(boxf);
            for (i = 0; i < nffr; i++)
            {
                for (j = 0; j < nat; j++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        xf[j][d] += filt[i]*x[i][j][d];
                    }
                }
                for (j = 0; j < DIM; j++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        boxf[j][d] += filt[i]*box[i][j][d];
                    }
                }
            }
            if (outl && (bLowAll || fr % nf == nf - 1))
            {
                write_trx(outl, nat, ind, topfile ? &(top.atoms) : NULL,
                          0, t[nf - 1], bFit ? topbox : boxf, xf, NULL, NULL);
            }
            if (outh)
            {
                /* Highpass filtering */
                for (j = 0; j < nat; j++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        xf[j][d] = xtop[j][d] + x[nf - 1][j][d] - xf[j][d];
                    }
                }
                if (bFit)
                {
                    for (j = 0; j < nat; j++)
                    {
                        rvec_inc(xf[j], xcmtop);
                    }
                }
                for (j = 0; j < DIM; j++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        boxf[j][d] = topbox[j][d] + box[nf - 1][j][d] - boxf[j][d];
                    }
                }
                write_trx(outh, nat, ind, topfile ? &(top.atoms) : NULL,
                          0, t[nf - 1], bFit ? topbox : boxf, xf, NULL, NULL);
            }
        }
        /* Cycle all the pointer and the box by one */
        ptr = x[0];
        for (i = 0; i < nffr-1; i++)
        {
            t[i] = t[i+1];
            x[i] = x[i+1];
            copy_mat(box[i+1], box[i]);
        }
        x[nffr - 1] = ptr;
        fr++;
    }
    while (read_next_x(oenv, in, &(t[nffr - 1]), x[nffr - 1], box[nffr - 1]));

    if (bTop)
    {
        gmx_rmpbc_done(gpbc);
    }

    if (outh)
    {
        close_trx(outh);
    }
    if (outl)
    {
        close_trx(outl);
    }
    close_trx(in);

    return 0;
}
