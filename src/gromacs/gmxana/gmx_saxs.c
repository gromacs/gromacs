/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/sfactor.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/smalloc.h"

int gmx_saxs(int argc, char *argv[])
{
    const char  *desc[] = {
        "[THISMODULE] calculates SAXS structure factors for given index",
        "groups based on Cromer's method.",
        "Both topology and trajectory files are required."
    };

    static real  start_q = 0.0, end_q = 60.0, energy = 12.0;
    static int   ngroups = 1;

    t_pargs      pa[] = {
        { "-ng",       FALSE, etINT, {&ngroups},
          "Number of groups to compute SAXS" },
        {"-startq", FALSE, etREAL, {&start_q},
         "Starting q (1/nm) "},
        {"-endq", FALSE, etREAL, {&end_q},
         "Ending q (1/nm)"},
        {"-energy", FALSE, etREAL, {&energy},
         "Energy of the incoming X-ray (keV) "}
    };
#define NPA asize(pa)
    const char  *fnTPS, *fnTRX, *fnNDX, *fnDAT = NULL;
    output_env_t oenv;

    t_filenm     fnm[] = {
        { efTRX, "-f",  NULL,      ffREAD },
        { efTPS, NULL,  NULL,      ffREAD },
        { efNDX, NULL,  NULL,      ffOPTRD },
        { efDAT, "-d",  "sfactor", ffOPTRD },
        { efXVG, "-sq", "sq",      ffWRITE },
    };
#define NFILE asize(fnm)
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    fnTPS = ftp2fn(efTPS, NFILE, fnm);
    fnTRX = ftp2fn(efTRX, NFILE, fnm);
    fnDAT = ftp2fn(efDAT, NFILE, fnm);
    fnNDX = ftp2fn_null(efNDX, NFILE, fnm);

    do_scattering_intensity(fnTPS, fnNDX, opt2fn("-sq", NFILE, fnm),
                            fnTRX, fnDAT,
                            start_q, end_q, energy, ngroups, oenv);

    please_cite(stdout, "Cromer1968a");

    return 0;
}
