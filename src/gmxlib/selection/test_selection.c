/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
    gmx_bool                     bFrameTree;
    int                          nmaxind;
    gmx_ana_selcollection_t     *sc;
} t_dumpdata;

int
dump_frame(t_topology * top, t_trxframe * fr, t_pbc * pbc,
           int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_dumpdata         *d = (t_dumpdata *)data;
    int                 g, i, n;

    fprintf(stderr, "\n");
    if (d->bFrameTree)
    {
        gmx_ana_selcollection_print_tree(stderr, d->sc, TRUE);
    }
    for (g = 0; g < nr; ++g)
    {
        gmx_ana_index_dump(sel[g]->g, g, d->nmaxind);
        fprintf(stderr, "  Positions (%d pcs):\n", sel[g]->p.nr);
        n = sel[g]->p.nr;
        if (d->nmaxind >= 0 && n > d->nmaxind)
        {
            n = d->nmaxind;
        }
        for (i = 0; i < n; ++i)
        {
            fprintf(stderr, "    (%.2f,%.2f,%.2f) r=%d, m=%d, b=%d-%d\n",
                    sel[g]->p.x[i][XX], sel[g]->p.x[i][YY], sel[g]->p.x[i][ZZ],
                    sel[g]->p.m.refid[i], sel[g]->p.m.mapid[i],
                    sel[g]->p.m.mapb.index[i]+1,
                    sel[g]->p.m.mapb.index[i+1]);
        }
        if (n < sel[g]->p.nr)
        {
            fprintf(stderr, "    ...\n");
        }
    }
    fprintf(stderr, "\n");
    return 0;
}

static void
print_selections(int nr, gmx_ana_selection_t **sel, int nmaxind)
{
    int                 g, i, n;

    fprintf(stderr, "\nSelections:\n");
    for (g = 0; g < nr; ++g)
    {
        fprintf(stderr, "  ");
        gmx_ana_selection_print_info(sel[g]);
        fprintf(stderr, "    ");
        gmx_ana_index_dump(sel[g]->g, g, nmaxind);

        fprintf(stderr, "    Block (size=%d):", sel[g]->p.m.mapb.nr);
        if (!sel[g]->p.m.mapb.index)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            n = sel[g]->p.m.mapb.nr;
            if (nmaxind >= 0 && n > nmaxind)
            {
                n = nmaxind;
            }
            for (i = 0; i <= n; ++i)
            {
                fprintf(stderr, " %d", sel[g]->p.m.mapb.index[i]);
            }
            if (n < sel[g]->p.m.mapb.nr)
            {
                fprintf(stderr, " ...");
            }
        }
        fprintf(stderr, "\n");

        n = sel[g]->p.m.nr;
        if (nmaxind >= 0 && n > nmaxind)
        {
            n = nmaxind;
        }
        fprintf(stderr, "    RefId:");
        if (!sel[g]->p.m.refid)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            for (i = 0; i < n; ++i)
            {
                fprintf(stderr, " %d", sel[g]->p.m.refid[i]);
            }
            if (n < sel[g]->p.m.nr)
            {
                fprintf(stderr, " ...");
            }
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "    MapId:");
        if (!sel[g]->p.m.mapid)
        {
            fprintf(stderr, " (null)");
        }
        else
        {
            for (i = 0; i < n; ++i)
            {
                fprintf(stderr, " %d", sel[g]->p.m.mapid[i]);
            }
            if (n < sel[g]->p.m.nr)
            {
                fprintf(stderr, " ...");
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

int
gmx_test_selection(int argc, char *argv[])
{
    const char             *desc[] = {
        "This is a test program for selections.",
    };

    gmx_bool                bMaskOnly     = FALSE;
    gmx_bool                bFrameTree    = FALSE;
    gmx_bool                bDebugCompile = FALSE;
    int                     nref          = 0;
    int                     nmaxind       = 20;
    t_pargs                 pa[]          = {
        {"-mask",   FALSE, etBOOL, {&bMaskOnly},
         "Test position mask functionality"},
        {"-compdebug", FALSE, etBOOL, {&bDebugCompile},
         "Print intermediate trees during compilation"},
        {"-frtree", FALSE, etBOOL, {&bFrameTree},
         "Print the whole evaluation tree for each frame"},
        {"-nref",   FALSE, etINT,  {&nref},
         "Number of reference selections to ask for"},
        {"-pmax",   FALSE, etINT,  {&nmaxind},
         "Maximum number of indices to print in lists (-1 = print all)"},
    };

    t_filenm                fnm[] = {
        {efDAT, "-o", "debug", ffOPTWR},
    };

    gmx_ana_traj_t         *trj;
    t_dumpdata              d;
    int                     ngrps;
    gmx_ana_selection_t   **sel;
    output_env_t            oenv;

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
    gmx_ana_selcollection_set_compile_debug(d.sc, bDebugCompile);
    gmx_ana_set_nrefgrps(trj, nref);
    gmx_ana_init_selections(trj);
    gmx_ana_get_ngrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    d.bFrameTree = bFrameTree;
    d.nmaxind    = nmaxind;

    print_selections(ngrps, sel, d.nmaxind);

    gmx_ana_do(trj, 0, &dump_frame, &d);

    print_selections(ngrps, sel, d.nmaxind);

    gmx_ana_traj_free(trj);
    done_filenms(NFILE, fnm);

    return 0;
}

int
main(int argc, char *argv[])
{
    gmx_test_selection(argc, argv);
    return 0;
}
