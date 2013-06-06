/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "gstat.h"
#include "matio.h"
#include "gmx_ana.h"
#include "names.h"
#include "sfactor.h"

int gmx_saxs(int argc, char *argv[])
{
    const char  *desc[] = {
        "g_saxs calculates SAXS structure factors for given index groups based on Cromer's method.",
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
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv);

    fnTPS = ftp2fn(efTPS, NFILE, fnm);
    fnTRX = ftp2fn(efTRX, NFILE, fnm);
    fnDAT = ftp2fn(efDAT, NFILE, fnm);
    fnNDX = ftp2fn_null(efNDX, NFILE, fnm);

    do_scattering_intensity(fnTPS, fnNDX, opt2fn("-sq", NFILE, fnm),
                            fnTRX, fnDAT,
                            start_q, end_q, energy, ngroups, oenv);

    please_cite(stdout, "Cromer1968a");

    thanx(stderr);

    return 0;
}
