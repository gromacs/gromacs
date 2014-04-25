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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/gmxana/gstat.h"
#include "gmx_viscosity.h"
#include "viscosity.h"
#include "handler.h"

int gmx_viscosity(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] computes shear- and bulk viscosity using a",
        "number of different algorithms."
    };

    t_filenm           fnm[] = {
        { efEDR, "-f",    NULL,       ffRDMULT },
        { efTPX, "-s",    NULL,       ffOPTRD },
        { efXVG, "-vis",  "visco",    ffWRITE },
        { efXVG, "-evis", "evisco",   ffOPTWR },
        { efXVG, "-eivis", "eivisco", ffOPTWR }
    };
#define NFILE asize(fnm)
    int                npargs;
    t_pargs           *ppa;
    output_env_t       oenv;

    npargs = 0;
    ppa    = add_acf_pargs(&npargs, NULL);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_BE_NICE,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    gmx::EnergyHandler eh;
    {
        char **fnms;
        int    nfile = opt2fns(&fnms, "-f", NFILE, fnm);
        if (0 == nfile)
        {
            fprintf(stderr, "No input files, nothing to do!\n");
            return 0;
        }
        for (int i = 0; (i < nfile); i++)
        {
            eh.addEnergyFile(fnms[i]);
        }
    }
    const char     *fn;
    // Viscosity analysis
    gmx::Viscosity *vis = new(gmx::Viscosity);

    vis->helper()->setOutputEnvironment(oenv);
    fn = opt2fn_null("-vis", NFILE, fnm);
    if (NULL != fn)
    {
        vis->helper()->setOutputFile(fn);
    }
    fn = opt2fn_null("-evis", NFILE, fnm);
    if (NULL != fn)
    {
        vis->setFnEinstein(fn);
    }
    fn = opt2fn_null("-eivis", NFILE, fnm);
    if (NULL != fn)
    {
        vis->setFnEinsteinIntegral(fn);
    }
    eh.addAnalysisTool(vis);

    if (!eh.readFiles())
    {
        fprintf(stderr, "Too bad. Something went wrong reading %s.\n",
                opt2fn("-f", NFILE, fnm));
    }

    return 0;
}
