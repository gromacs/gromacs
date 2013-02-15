/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "mshift.h"
#include "xvgr.h"
#include "princ.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "tpxio.h"
#include "gstat.h"
#include "gmx_ana.h"


void
calc_principal_axes(t_topology *   top,
                    rvec *         x,
                    atom_id *      index,
                    int            n,
                    matrix         axes,
                    rvec           inertia)
{
    rvec   xcm;

    sub_xcm(x, n, index, top->atoms.atom, xcm, FALSE);
    principal_comp(n, index, top->atoms.atom, x, axes, inertia);
}

int gmx_principal(int argc, char *argv[])
{
    const char     *desc[] = {
        "[TT]g_principal[tt] calculates the three principal axes of inertia for a group",
        "of atoms.",
    };
    static gmx_bool foo = FALSE;

    t_pargs         pa[] = {
        { "-foo",      FALSE, etBOOL, {&foo}, "Dummy option to avoid empty array" }
    };
    t_trxstatus    *status;
    t_topology      top;
    int             ePBC;
    real            t;
    rvec      *     x;

    int             natoms;
    char           *grpname, title[256];
    int             i, j, m, gnx, nam, mol;
    atom_id        *index;
    rvec            a1, a2, a3, moi;
    FILE      *     axis1;
    FILE      *     axis2;
    FILE      *     axis3;
    FILE      *     fmoi;
    matrix          axes, box;
    output_env_t    oenv;
    gmx_rmpbc_t     gpbc = NULL;


    t_filenm fnm[] = {
        { efTRX, "-f",   NULL,       ffREAD },
        { efTPS, NULL,   NULL,       ffREAD },
        { efNDX, NULL,   NULL,       ffOPTRD },
        { efDAT, "-a1",  "axis1",    ffWRITE },
        { efDAT, "-a2",  "axis2",    ffWRITE },
        { efDAT, "-a3",  "axis3",    ffWRITE },
        { efDAT, "-om",  "moi",      ffWRITE }
    };
#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv,
                      PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    axis1 = ffopen(opt2fn("-a1", NFILE, fnm), "w");
    axis2 = ffopen(opt2fn("-a2", NFILE, fnm), "w");
    axis3 = ffopen(opt2fn("-a3", NFILE, fnm), "w");
    fmoi  = ffopen(opt2fn("-om", NFILE, fnm), "w");

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, NULL, NULL, box, TRUE);

    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);

    do
    {
        gmx_rmpbc(gpbc, natoms, box, x);

        calc_principal_axes(&top, x, index, gnx, axes, moi);

        fprintf(axis1, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[XX][XX], axes[YY][XX], axes[ZZ][XX]);
        fprintf(axis2, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[XX][YY], axes[YY][YY], axes[ZZ][YY]);
        fprintf(axis3, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[XX][ZZ], axes[YY][ZZ], axes[ZZ][ZZ]);
        fprintf(fmoi,  "%15.10f     %15.10f  %15.10f  %15.10f\n", t, moi[XX], moi[YY], moi[ZZ]);
    }
    while (read_next_x(oenv, status, &t, natoms, x, box));

    gmx_rmpbc_done(gpbc);


    close_trj(status);
    ffclose(axis1);
    ffclose(axis2);
    ffclose(axis3);
    ffclose(fmoi);

    thanx(stderr);

    return 0;
}
