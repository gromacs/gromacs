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

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/nrama.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


static void plot_rama(FILE *out, t_xrama *xr)
{
    int  i;
    real phi, psi;

    for (i = 0; (i < xr->npp); i++)
    {
        phi = xr->dih[xr->pp[i].iphi].ang*RAD2DEG;
        psi = xr->dih[xr->pp[i].ipsi].ang*RAD2DEG;
        fprintf(out, "%g  %g  %s\n", phi, psi, xr->pp[i].label);
    }
}

int gmx_rama(int argc, char *argv[])
{
    const char  *desc[] = {
        "[THISMODULE] selects the [GRK]phi[grk]/[GRK]psi[grk] dihedral combinations from your topology file",
        "and computes these as a function of time.",
        "Using simple Unix tools such as [IT]grep[it] you can select out",
        "specific residues."
    };

    FILE        *out;
    t_xrama     *xr;
    int          j;
    output_env_t oenv;
    t_filenm     fnm[] = {
        { efTRX, "-f", NULL,  ffREAD },
        { efTPR, NULL, NULL,  ffREAD },
        { efXVG, NULL, "rama", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, 0, NULL, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }


    snew(xr, 1);
    init_rama(oenv, ftp2fn(efTRX, NFILE, fnm), ftp2fn(efTPR, NFILE, fnm), xr, 3);

    out = xvgropen(ftp2fn(efXVG, NFILE, fnm), "Ramachandran Plot", "Phi", "Psi", oenv);
    xvgr_line_props(out, 0, elNone, ecFrank, oenv);
    xvgr_view(out, 0.2, 0.2, 0.8, 0.8, oenv);
    xvgr_world(out, -180, -180, 180, 180, oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@    xaxis  tick on\n@    xaxis  tick major 60\n@    xaxis  tick minor 30\n");
        fprintf(out, "@    yaxis  tick on\n@    yaxis  tick major 60\n@    yaxis  tick minor 30\n");
        fprintf(out, "@ s0 symbol 2\n@ s0 symbol size 0.4\n@ s0 symbol fill 1\n");
    }
    j = 0;
    do
    {
        plot_rama(out, xr);
        j++;
    }
    while (new_data(xr));
    fprintf(stderr, "\n");
    xvgrclose(out);

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), NULL);

    return 0;
}
