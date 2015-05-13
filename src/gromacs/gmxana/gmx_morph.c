/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static real dointerp(int n, rvec x1[], rvec x2[], rvec xx[],
                     int I, int N, real first, real last)
{
    int    i, j;
    double fac, fac0, fac1;

    fac  = first + (I*(last-first))/(N-1);
    fac0 = 1-fac;
    fac1 = fac;
    for (i = 0; (i < n); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            xx[i][j] = fac0*x1[i][j] + fac1*x2[i][j];
        }
    }

    return fac;
}

int gmx_morph(int argc, char *argv[])
{
    const char      *desc[] = {
        "[THISMODULE] does a linear interpolation of conformations in order to",
        "create intermediates. Of course these are completely unphysical, but",
        "that you may try to justify yourself. Output is in the form of a ",
        "generic trajectory. The number of intermediates can be controlled with",
        "the [TT]-ninterm[tt] flag. The first and last flag correspond to the way of",
        "interpolating: 0 corresponds to input structure 1 while",
        "1 corresponds to input structure 2.",
        "If you specify [TT]-first[tt] < 0 or [TT]-last[tt] > 1 extrapolation will be",
        "on the path from input structure x[SUB]1[sub] to x[SUB]2[sub]. In general, the coordinates",
        "of the intermediate x(i) out of N total intermediates correspond to:[PAR]",
        "x(i) = x[SUB]1[sub] + (first+(i/(N-1))*(last-first))*(x[SUB]2[sub]-x[SUB]1[sub])[PAR]",
        "Finally the RMSD with respect to both input structures can be computed",
        "if explicitly selected ([TT]-or[tt] option). In that case, an index file may be",
        "read to select the group from which the RMS is computed."
    };
    t_filenm         fnm[] = {
        { efSTX, "-f1", "conf1",  ffREAD },
        { efSTX, "-f2", "conf2",  ffREAD },
        { efTRX, "-o",  "interm", ffWRITE },
        { efXVG, "-or", "rms-interm", ffOPTWR },
        { efNDX, "-n",  "index",  ffOPTRD }
    };
#define NFILE asize(fnm)
    static  int      ninterm = 11;
    static  real     first   = 0.0;
    static  real     last    = 1.0;
    static  gmx_bool bFit    = TRUE;
    t_pargs          pa []   = {
        { "-ninterm", FALSE, etINT,  {&ninterm},
          "Number of intermediates" },
        { "-first",   FALSE, etREAL, {&first},
          "Corresponds to first generated structure (0 is input x[SUB]1[sub], see above)" },
        { "-last",    FALSE, etREAL, {&last},
          "Corresponds to last generated structure (1 is input x[SUB]2[sub], see above)" },
        { "-fit",     FALSE, etBOOL, {&bFit},
          "Do a least squares fit of the second to the first structure before interpolating" }
    };
    const char      *leg[] = { "Ref = 1\\Sst\\N conf", "Ref = 2\\Snd\\N conf" };
    FILE            *fp    = NULL;
    int              i, isize, is_lsq, nat1, nat2;
    t_trxstatus     *status;
    atom_id         *index, *index_lsq, *index_all, *dummy;
    t_atoms          atoms;
    rvec            *x1, *x2, *xx, *v;
    matrix           box;
    real             rms1, rms2, fac, *mass;
    char             title[STRLEN], *grpname;
    gmx_bool         bRMS;
    output_env_t     oenv;

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }
    get_stx_coordnum (opt2fn("-f1", NFILE, fnm), &nat1);
    get_stx_coordnum (opt2fn("-f2", NFILE, fnm), &nat2);
    if (nat1 != nat2)
    {
        gmx_fatal(FARGS, "Number of atoms in first structure is %d, in second %d",
                  nat1, nat2);
    }

    init_t_atoms(&atoms, nat1, TRUE);
    snew(x1, nat1);
    snew(x2, nat1);
    snew(xx, nat1);
    snew(v, nat1);

    read_stx_conf(opt2fn("-f1", NFILE, fnm), title, &atoms, x1, v, NULL, box);
    read_stx_conf(opt2fn("-f2", NFILE, fnm), title, &atoms, x2, v, NULL, box);

    snew(mass, nat1);
    snew(index_all, nat1);
    for (i = 0; (i < nat1); i++)
    {
        mass[i]      = 1;
        index_all[i] = i;
    }
    if (bFit)
    {
        printf("Select group for LSQ superposition:\n");
        get_index(&atoms, opt2fn_null("-n", NFILE, fnm), 1, &is_lsq, &index_lsq,
                  &grpname);
        reset_x(is_lsq, index_lsq, nat1, index_all, x1, mass);
        reset_x(is_lsq, index_lsq, nat1, index_all, x2, mass);
        do_fit(nat1, mass, x1, x2);
    }

    bRMS = opt2bSet("-or", NFILE, fnm);
    if (bRMS)
    {
        fp = xvgropen(opt2fn("-or", NFILE, fnm), "RMSD", "Conf", "(nm)", oenv);
        xvgr_legend(fp, asize(leg), leg, oenv);
        printf("Select group for RMSD calculation:\n");
        get_index(&atoms, opt2fn_null("-n", NFILE, fnm), 1, &isize, &index, &grpname);
        printf("You selected group %s, containing %d atoms\n", grpname, isize);
        rms1 = rmsdev_ind(isize, index, mass, x1, x2);
        fprintf(stderr, "RMSD between input conformations is %g nm\n", rms1);
    }

    snew(dummy, nat1);
    for (i = 0; (i < nat1); i++)
    {
        dummy[i] = i;
    }
    status = open_trx(ftp2fn(efTRX, NFILE, fnm), "w");

    for (i = 0; (i < ninterm); i++)
    {
        fac = dointerp(nat1, x1, x2, xx, i, ninterm, first, last);
        write_trx(status, nat1, dummy, &atoms, i, fac, box, xx, NULL, NULL);
        if (bRMS)
        {
            rms1 = rmsdev_ind(isize, index, mass, x1, xx);
            rms2 = rmsdev_ind(isize, index, mass, x2, xx);
            fprintf(fp, "%10g  %10g  %10g\n", fac, rms1, rms2);
        }
    }

    close_trx(status);

    if (bRMS)
    {
        xvgrclose(fp);
        do_view(oenv, opt2fn("-or", NFILE, fnm), "-nxy");
    }

    return 0;
}
