/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "statutil.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "gmx_ana.h"
#include "string2.h"

int gmx_membed(int argc, char *argv[])
{
    const char *desc[] = {
        "[TT]g_membed[tt] embeds a membrane protein into an equilibrated lipid bilayer at the position",
        "and orientation specified by the user.[PAR]",
        "SHORT MANUAL[BR]------------[BR]",
        "The user should merge the structure files of the protein and membrane (+solvent), creating a",
        "single structure file with the protein overlapping the membrane at the desired position and",
        "orientation. The box size is taken from the membrane structure file. The corresponding topology",
        "files should also be merged. Consecutively, create a [TT].tpr[tt] file (input for [TT]g_membed[tt]) from these files,"
        "with the following options included in the [TT].mdp[tt] file.[BR]",
        " - [TT]integrator      = md[tt][BR]",
        " - [TT]energygrps      = Protein[tt] (or other group that you want to insert)[BR]",
        " - [TT]freezegrps      = Protein[tt][BR]",
        " - [TT]freezedim       = Y Y Y[tt][BR]",
        " - [TT]energygrp_excl  = Protein Protein[tt][BR]",
        "The output is a structure file containing the protein embedded in the membrane. If a topology",
        "file is provided, the number of lipid and ",
        "solvent molecules will be updated to match the new structure file.[BR]",
        "For a more extensive manual see Wolf et al, J Comp Chem 31 (2010) 2169-2174, Appendix.[PAR]",
        "SHORT METHOD DESCRIPTION[BR]",
        "------------------------[BR]",
        "1. The protein is resized around its center of mass by a factor [TT]-xy[tt] in the xy-plane",
        "(the membrane plane) and a factor [TT]-z[tt] in the [IT]z[it]-direction (if the size of the",
        "protein in the z-direction is the same or smaller than the width of the membrane, a",
        "[TT]-z[tt] value larger than 1 can prevent that the protein will be enveloped by the lipids).[BR]",
        "2. All lipid and solvent molecules overlapping with the resized protein are removed. All",
        "intraprotein interactions are turned off to prevent numerical issues for small values of [TT]-xy[tt]",
        " or [TT]-z[tt][BR]",
        "3. One md step is performed.[BR]",
        "4. The resize factor ([TT]-xy[tt] or [TT]-z[tt]) is incremented by a small amount ((1-xy)/nxy or (1-z)/nz) and the",
        "protein is resized again around its center of mass. The resize factor for the xy-plane",
        "is incremented first. The resize factor for the z-direction is not changed until the [TT]-xy[tt] factor",
        "is 1 (thus after [TT]-nxy[tt] iterations).[BR]",
        "5. Repeat step 3 and 4 until the protein reaches its original size ([TT]-nxy[tt] + [TT]-nz[tt] iterations).[BR]",
        "For a more extensive method description see Wolf et al, J Comp Chem, 31 (2010) 2169-2174.[PAR]",
        "NOTE[BR]----[BR]",
        " - Protein can be any molecule you want to insert in the membrane.[BR]",
        " - It is recommended to perform a short equilibration run after the embedding",
        "(see Wolf et al, J Comp Chem 31 (2010) 2169-2174), to re-equilibrate the membrane. Clearly",
        "protein equilibration might require longer.[PAR]"
    };
    t_filenm    fnm[] = {
        { efTPX, "-f",      "into_mem", ffREAD },
        { efNDX, "-n",      "index",    ffOPTRD },
        { efTOP, "-p",      "topol",    ffOPTRW },
        { efTRN, "-o",      NULL,       ffWRITE },
        { efXTC, "-x",      NULL,       ffOPTWR },
        { efSTO, "-c",      "membedded",  ffWRITE },
        { efEDR, "-e",      "ener",     ffWRITE },
        { efDAT, "-dat",    "membed",   ffWRITE }
    };
#define NFILE asize(fnm)

    /* Command line options ! */
    real         xy_fac           = 0.5;
    real         xy_max           = 1.0;
    real         z_fac            = 1.0;
    real         z_max            = 1.0;
    int          it_xy            = 1000;
    int          it_z             = 0;
    real         probe_rad        = 0.22;
    int          low_up_rm        = 0;
    int          maxwarn          = 0;
    int          pieces           = 1;
    gmx_bool     bALLOW_ASYMMETRY = FALSE;
    gmx_bool     bStart           = FALSE;
    int          nstepout         = 100;
    gmx_bool     bVerbose         = FALSE;
    char        *mdrun_path       = NULL;

    t_pargs      pa[] = {
        { "-xyinit",   FALSE, etREAL,  {&xy_fac},
          "Resize factor for the protein in the xy dimension before starting embedding" },
        { "-xyend",   FALSE, etREAL,  {&xy_max},
          "Final resize factor in the xy dimension" },
        { "-zinit",    FALSE, etREAL,  {&z_fac},
          "Resize factor for the protein in the z dimension before starting embedding" },
        { "-zend",    FALSE, etREAL,  {&z_max},
          "Final resize faction in the z dimension" },
        { "-nxy",     FALSE,  etINT,  {&it_xy},
          "Number of iteration for the xy dimension" },
        { "-nz",      FALSE,  etINT,  {&it_z},
          "Number of iterations for the z dimension" },
        { "-rad",     FALSE, etREAL,  {&probe_rad},
          "Probe radius to check for overlap between the group to embed and the membrane"},
        { "-pieces",  FALSE,  etINT,  {&pieces},
          "Perform piecewise resize. Select parts of the group to insert and resize these with respect to their own geometrical center." },
        { "-asymmetry", FALSE, etBOOL, {&bALLOW_ASYMMETRY},
          "Allow asymmetric insertion, i.e. the number of lipids removed from the upper and lower leaflet will not be checked." },
        { "-ndiff",  FALSE, etINT, {&low_up_rm},
          "Number of lipids that will additionally be removed from the lower (negative number) or upper (positive number) membrane leaflet." },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Maximum number of warning allowed" },
        { "-start",   FALSE, etBOOL, {&bStart},
          "Call mdrun with membed options" },
        { "-stepout", FALSE, etINT, {&nstepout},
          "HIDDENFrequency of writing the remaining runtime" },
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-mdrun_path", FALSE, etSTR, {&mdrun_path},
          "Path to the mdrun executable compiled with this g_membed version" }
    };

    FILE        *data_out;
    output_env_t oenv;
    char         buf[256], buf2[64];
    gmx_bool     bSucces;

    parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                      asize(desc), desc, 0, NULL, &oenv);

    data_out = ffopen(opt2fn("-dat", NFILE, fnm), "w");

    fprintf(data_out, "nxy = %d\nnz = %d\nxyinit = %f\nxyend = %f\nzinit = %f\nzend = %f\n"
            "rad = %f\npieces = %d\nasymmetry = %s\nndiff = %d\nmaxwarn = %d\n",
            it_xy, it_z, xy_fac, xy_max, z_fac, z_max, probe_rad, pieces,
            bALLOW_ASYMMETRY ? "yes" : "no", low_up_rm, maxwarn);

    fclose(data_out);

    sprintf(buf, "%s -s %s -membed %s -o %s -c %s -e %s -nt 1 -cpt -1",
            (mdrun_path == NULL) ? "mdrun" : mdrun_path,
            opt2fn("-f", NFILE, fnm), opt2fn("-dat", NFILE, fnm), opt2fn("-o", NFILE, fnm),
            opt2fn("-c", NFILE, fnm), opt2fn("-e", NFILE, fnm));

    if (opt2bSet("-n", NFILE, fnm))
    {
        sprintf(buf2, " -mn %s", opt2fn("-n", NFILE, fnm));
        strcat(buf, buf2);
    }

    if (opt2bSet("-x", NFILE, fnm))
    {
        sprintf(buf2, " -x %s", opt2fn("-x", NFILE, fnm));
        strcat(buf, buf2);
    }

    if (opt2bSet("-p", NFILE, fnm))
    {
        sprintf(buf2, " -mp %s", opt2fn("-p", NFILE, fnm));
        strcat(buf, buf2);
    }

    if (bVerbose)
    {
        sprintf(buf2, " -v -stepout %d", nstepout);
        strcat(buf, buf2);
    }

    if (bStart)
    {
        printf("Start run with:\n%s\n", buf);
        bSucces = system(buf);
    }
    else
    {
        printf("You can membed your protein now by:\n%s\n", buf);
    }

    fprintf(stderr, "Please cite:\nWolf et al, J Comp Chem 31 (2010) 2169-2174.\n");

    return 0;
}
