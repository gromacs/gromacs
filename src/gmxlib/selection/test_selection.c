/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
 */
/*! \internal \file
 * \brief Testing/debugging tool for the selection engine.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <copyrite.h>
#include <filenm.h>
#include <macros.h>
#include <smalloc.h>
#include <statutil.h>

#include <trajana.h>

typedef struct
{
    bool                     bFrameTree;
    gmx_ana_selcollection_t *sc;
} t_dumpdata;

int
dump_frame(t_topology * top, t_trxframe * fr, t_pbc * pbc,
           int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_dumpdata         *d = (t_dumpdata *)data;
    int                 g, i;

    fprintf(stderr, "\n");
    if (d->bFrameTree)
    {
        gmx_ana_selcollection_print_tree(stderr, d->sc, TRUE);
    }
    for (g = 0; g < nr; ++g)
    {
        gmx_ana_index_dump(sel[g]->g, g);
        fprintf(stderr, "  Positions (%d pcs):\n", sel[g]->p.nr);
        for (i = 0; i < sel[g]->p.nr; ++i)
        {
            fprintf(stderr, "    (%.2f,%.2f,%.2f) r=%d, m=%d, b=%d-%d\n",
                    sel[g]->p.x[i][XX], sel[g]->p.x[i][YY], sel[g]->p.x[i][ZZ],
                    sel[g]->p.m.refid[i], sel[g]->p.m.mapid[i],
                    sel[g]->p.m.mapb.index[i]+1,
                    sel[g]->p.m.mapb.index[i+1]);
        }
    }
    fprintf(stderr, "\n");
    return 0;
}

static void
print_selections(int nr, gmx_ana_selection_t **sel)
{
    int                 g, i;

    fprintf(stderr, "\nSelections:\n");
    for (g = 0; g < nr; ++g)
    {
        fprintf(stderr, "  ");
        gmx_ana_selection_print_info(sel[g]);
        fprintf(stderr, "    ");
        gmx_ana_index_dump(sel[g]->g, g);

        fprintf(stderr, "    Block (size=%d):", sel[g]->p.m.mapb.nr);
        if (!sel[g]->p.m.mapb.index)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            for (i = 0; i <= sel[g]->p.m.mapb.nr; ++i)
                fprintf(stderr, " %d", sel[g]->p.m.mapb.index[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "    RefId:");
        if (!sel[g]->p.m.refid)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            for (i = 0; i < sel[g]->p.m.nr; ++i)
                fprintf(stderr, " %d", sel[g]->p.m.refid[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "    MapId:");
        if (!sel[g]->p.m.mapid)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            for (i = 0; i < sel[g]->p.m.nr; ++i)
                fprintf(stderr, " %d", sel[g]->p.m.mapid[i]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

int
gmx_test_selection(int argc, char *argv[])
{
    const char         *desc[] = {
        "This is a test program for selections.",
    };

    bool                bMaskOnly  = FALSE;
    bool                bFrameTree = FALSE;
    int                 nref       = 0;
    t_pargs             pa[] = {
        {"-mask",   FALSE, etBOOL, {&bMaskOnly},
         "Test position mask functionality"},
        {"-frtree", FALSE, etBOOL, {&bFrameTree},
         "Print the whole evaluation tree for each frame"},
        {"-nref",   FALSE, etINT,  {&nref},
         "Number of reference selections to ask for"},
    };

    t_filenm            fnm[] = {
        {efDAT, "-o", "debug", ffOPTWR},
    };

    gmx_ana_traj_t       *trj;
    t_dumpdata            d;
    int                   ngrps;
    gmx_ana_selection_t **sel;
    output_env_t          oenv;

#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);

    gmx_ana_traj_create(&trj, ANA_DEBUG_SELECTION | ANA_USER_SELINIT | ANA_USE_FULLGRPS);
    gmx_ana_get_selcollection(trj, &d.sc);
    gmx_ana_set_nanagrps(trj, -1);
    parse_trjana_args(trj, &argc, argv, 0,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    if (bMaskOnly)
    {
        gmx_ana_add_flags(trj, ANA_USE_POSMASK);
        gmx_ana_selcollection_set_outpostype(d.sc, NULL, TRUE);
    }
    gmx_ana_set_nrefgrps(trj, nref);
    gmx_ana_init_selections(trj);
    gmx_ana_get_ngrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    print_selections(ngrps, sel);

    d.bFrameTree = bFrameTree;
    gmx_ana_do(trj, 0, &dump_frame, &d);

    print_selections(ngrps, sel);

    gmx_ana_traj_free(trj);

    return 0;
}

int
main(int argc, char *argv[])
{
    gmx_test_selection(argc, argv);
    return 0;
}
