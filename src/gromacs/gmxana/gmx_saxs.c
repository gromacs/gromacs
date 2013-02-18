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
    const char             *desc[] = {
        "g_saxs calculates SAXS structure factors for given index groups based on Cromer's method.",
        "Both topology and trajectory files are required."
    };

    static real             start_q = 0.0, end_q = 60.0, energy = 12.0;
    static int              ngroups = 1;

    t_pargs                 pa[] = {
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
    const char             *fnTPS, *fnTRX, *fnNDX, *fnDAT = NULL;
    int                     i, *isize, flags = TRX_READ_X, **index_atp;
    t_trxstatus            *status;
    char                  **grpname, title[STRLEN];
    atom_id               **index;
    t_topology              top;
    int                     ePBC;
    t_trxframe              fr;
    reduced_atom_t        **red;
    structure_factor_t     *sf;
    rvec                   *xtop;
    real                  **sf_table;
    int                     nsftable;
    matrix                  box;
    double                  r_tmp;

    gmx_structurefactors_t *gmx_sf;
    real                   *a, *b, c;
    int                     success;

    output_env_t            oenv;

    t_filenm                fnm[] = {
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

    snew(a, 4);
    snew(b, 4);


    gmx_sf = gmx_structurefactors_init(fnDAT);

    success = gmx_structurefactors_get_sf(gmx_sf, 0, a, b, &c);

    snew (sf, 1);
    sf->energy = energy;

    /* Read the topology informations */
    read_tps_conf (fnTPS, title, &top, &ePBC, &xtop, NULL, box, TRUE);
    sfree (xtop);

    /* groups stuff... */
    snew (isize, ngroups);
    snew (index, ngroups);
    snew (grpname, ngroups);

    fprintf (stderr, "\nSelect %d group%s\n", ngroups,
             ngroups == 1 ? "" : "s");
    if (fnTPS)
    {
        get_index (&top.atoms, fnNDX, ngroups, isize, index, grpname);
    }
    else
    {
        rd_index (fnNDX, ngroups, isize, index, grpname);
    }

    /* The first time we read data is a little special */
    read_first_frame (oenv, &status, fnTRX, &fr, flags);

    sf->total_n_atoms = fr.natoms;

    snew (red, ngroups);
    snew (index_atp, ngroups);

    r_tmp = max (box[XX][XX], box[YY][YY]);
    r_tmp = (double) max (box[ZZ][ZZ], r_tmp);

    sf->ref_k = (2.0 * M_PI) / (r_tmp);
    /* ref_k will be the reference momentum unit */
    sf->n_angles = (int) (end_q / sf->ref_k + 0.5);

    snew (sf->F, ngroups);
    for (i = 0; i < ngroups; i++)
    {
        snew (sf->F[i], sf->n_angles);
    }
    for (i = 0; i < ngroups; i++)
    {
        snew (red[i], isize[i]);
        rearrange_atoms (red[i], &fr, index[i], isize[i], &top, TRUE, gmx_sf);
        index_atp[i] = create_indexed_atom_type (red[i], isize[i]);
    }

    sf_table = compute_scattering_factor_table (gmx_sf, (structure_factor_t *)sf, &nsftable);


    /* This is the main loop over frames */

    do
    {
        sf->nSteps++;
        for (i = 0; i < ngroups; i++)
        {
            rearrange_atoms (red[i], &fr, index[i], isize[i], &top, FALSE, gmx_sf);

            compute_structure_factor ((structure_factor_t *)sf, box, red[i], isize[i],
                                      start_q, end_q, i, sf_table);
        }
    }

    while (read_next_frame (oenv, status, &fr));

    save_data ((structure_factor_t *)sf, opt2fn_null("-sq", NFILE, fnm), ngroups, start_q, end_q, oenv);


    sfree(a);
    sfree(b);

    done_gmx_structurefactors(gmx_sf);

    please_cite(stdout, "Cromer1968a");

    thanx(stderr);

    return 0;
}
