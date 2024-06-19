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
#include <cstring>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;

int gmx_rotacf(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] calculates the rotational correlation function",
        "for molecules. Atom triplets (i,j,k) must be given in the index",
        "file, defining two vectors ij and jk. The rotational ACF",
        "is calculated as the autocorrelation function of the vector",
        "n = ij x jk, i.e. the cross product of the two vectors.",
        "Since three atoms span a plane, the order of the three atoms",
        "does not matter. Optionally, by invoking the [TT]-d[tt] switch, you can",
        "calculate the rotational correlation function for linear molecules",
        "by specifying atom pairs (i,j) in the index file.",
        "[PAR]",
        "EXAMPLES[PAR]",
        "[TT]gmx rotacf -P 1 -nparm 2 -fft -n index -o rotacf-x-P1",
        "-fa expfit-x-P1 -beginfit 2.5 -endfit 20.0[tt][PAR]",
        "This will calculate the rotational correlation function using a first",
        "order Legendre polynomial of the angle of a vector defined by the index",
        "file. The correlation function will be fitted from 2.5 ps until 20.0 ps",
        "to a two-parameter exponential."
    };
    static gmx_bool bVec = FALSE, bAver = TRUE;

    t_pargs pa[] = {
        { "-d",
          FALSE,
          etBOOL,
          { &bVec },
          "Use index doublets (vectors) for correlation function instead of triplets (planes)" },
        { "-aver", FALSE, etBOOL, { &bAver }, "Average over molecules" }
    };

    t_trxstatus*  status;
    int           isize;
    int*          index;
    char*         grpname;
    rvec *        x, *x_s;
    matrix        box;
    real**        c1;
    rvec          xij, xjk, n;
    int           i, m, teller, n_alloc, natoms, nvec, ai, aj, ak;
    unsigned long mode;
    real          t, t0, t1, dt;
    gmx_rmpbc_t   gpbc = nullptr;
    t_topology*   top;
    PbcType       pbcType;
    t_filenm      fnm[] = { { efTRX, "-f", nullptr, ffREAD },
                       { efTPR, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffREAD },
                       { efXVG, "-o", "rotacf", ffWRITE } };
#define NFILE asize(fnm)
    int      npargs;
    t_pargs* ppa;

    gmx_output_env_t* oenv;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);

    if (bVec)
    {
        nvec = isize / 2;
    }
    else
    {
        nvec = isize / 3;
    }

    if (((isize % 3) != 0) && !bVec)
    {
        gmx_fatal(FARGS,
                  "number of index elements not multiple of 3, "
                  "these can not be atom triplets\n");
    }
    if (((isize % 2) != 0) && bVec)
    {
        gmx_fatal(FARGS,
                  "number of index elements not multiple of 2, "
                  "these can not be atom doublets\n");
    }

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType);

    snew(c1, nvec);
    for (i = 0; (i < nvec); i++)
    {
        c1[i] = nullptr;
    }
    n_alloc = 0;

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    snew(x_s, natoms);

    gpbc = gmx_rmpbc_init(&(top->idef), pbcType, natoms);

    /* Start the loop over frames */
    t0     = t;
    teller = 0;
    do
    {
        if (teller >= n_alloc)
        {
            n_alloc += 100;
            for (i = 0; (i < nvec); i++)
            {
                srenew(c1[i], DIM * n_alloc);
            }
        }
        t1 = t;

        /* Remove periodicity */
        gmx_rmpbc_copy(gpbc, natoms, box, x, x_s);

        /* Compute crossproducts for all vectors, if triplets.
         * else, just get the vectors in case of doublets.
         */
        if (!bVec)
        {
            for (i = 0; (i < nvec); i++)
            {
                ai = index[3 * i];
                aj = index[3 * i + 1];
                ak = index[3 * i + 2];
                rvec_sub(x_s[ai], x_s[aj], xij);
                rvec_sub(x_s[aj], x_s[ak], xjk);
                cprod(xij, xjk, n);
                for (m = 0; (m < DIM); m++)
                {
                    c1[i][DIM * teller + m] = n[m];
                }
            }
        }
        else
        {
            for (i = 0; (i < nvec); i++)
            {
                ai = index[2 * i];
                aj = index[2 * i + 1];
                rvec_sub(x_s[ai], x_s[aj], n);
                for (m = 0; (m < DIM); m++)
                {
                    c1[i][DIM * teller + m] = n[m];
                }
            }
        }
        /* Increment loop counter */
        teller++;
    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);
    fprintf(stderr, "\nDone with trajectory\n");

    gmx_rmpbc_done(gpbc);


    /* Autocorrelation function */
    if (teller < 2)
    {
        fprintf(stderr, "Not enough frames for correlation function\n");
    }
    else
    {
        dt = (t1 - t0) / (static_cast<real>(teller - 1));

        mode = eacVector;

        do_autocorr(
                ftp2fn(efXVG, NFILE, fnm), oenv, "Rotational Correlation Function", teller, nvec, c1, dt, mode, bAver);
    }

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), nullptr);

    return 0;
}
