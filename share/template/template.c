/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
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
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>

#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>

/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    gmx_ana_selection_t *refsel;
    FILE                *fp;
    real                *ave;
    real                *n;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata     *d = (t_analysisdata *)data;
    int                 g, i;
    real                frave;
    int                 rc;

    /* Here, you can do whatever analysis your program requires for a frame. */
    if (d->fp)
    {
        fprintf(d->fp, "%10.3f", fr->time);
    }
    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }
    for (g = 0; g < nr; ++g)
    {
        frave = 0;
        for (i = 0; i < sel[g]->p.nr; ++i)
        {
            frave += gmx_ana_nbsearch_pos_mindist(d->nb, &sel[g]->p, i);
        }
        d->ave[g] += frave;
        d->n[g]   += sel[g]->p.nr;
        if (d->fp)
        {
            frave /= sel[g]->p.nr;
            fprintf(d->fp, " %.3f", frave);
        }
    }
    if (d->fp)
    {
        fprintf(d->fp, "\n");
    }
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */
int
gmx_template(int argc, char *argv[])
{
    const char         *desc[] = {
        "This is a template for writing your own analysis tools for",
        "Gromacs. The advantage of using Gromacs for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by Gromacs. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from the a reference group to one or more",
        "analysis groups.",
    };

    /* Command-line arguments */
    real                cutoff = 0;
    gmx_bool                bArg   = FALSE;
    t_pargs             pa[] = {
        {"-cutoff", FALSE, etREAL, {&cutoff},
         "Cutoff for distance calculation (0 = no cutoff)"},
        {"-arg2",   FALSE, etBOOL, {&bArg},
         "Example argument 2"},
    };
    /* The second argument is for demonstration purposes only */

    /* Output files */
    t_filenm            fnm[] = {
        {efXVG, "-o", "avedist", ffOPTWR},
    };
#define NFILE asize(fnm)

    gmx_ana_traj_t       *trj;
    output_env_t          oenv;
    t_analysisdata        d;
    int                   ngrps;
    gmx_ana_selection_t **sel;
    int                   g;
    int                   rc;

    CopyRight(stderr, argv[0]);
    /* Here, we can use flags to specify requirements for the selections and/or
     * other features of the library. */
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj, -1);
    /* If required, other functions can also be used to configure the library
     * before calling parse_trjana_args(). */
    parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);

    /* You can now do any initialization you wish, using the information
     * from the trj structure.
     * In particular, you should store any command-line parameter values that
     * the analysis part requires into d for them to be accessible. */

    /* First, we get some selection information from the structure */
    gmx_ana_get_refsel(trj, 0, &d.refsel);
    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    /* First, we initialize the neighborhood search for the first index
     * group. */
    rc = gmx_ana_nbsearch_create(&d.nb, cutoff, d.refsel->p.nr);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    /* We then allocate memory in d to store intermediate data
     * and initialize the counters to zero */
    snew(d.ave, ngrps);
    snew(d.n,   ngrps);

    /* We also open the output file if the user provided it */
    d.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Average distance",
                        "Time [ps]", "Distance [nm]", oenv);
        xvgr_selections(d.fp, trj);
    }

    /* Now, we do the actual analysis */
    gmx_ana_do(trj, 0, &analyze_frame, &d);

    /* Now, the analysis has been done for all frames, and you can access the
     * results in d. Here, you should post-process your data and write out any
     * averaged properties. */

    /* For the template, we close the output file if one was opened */
    if (d.fp)
    {
        ffclose(d.fp);
    }

    /* We also calculate the average distance for each group */
    for (g = 0; g < ngrps; ++g)
    {
        d.ave[g] /= d.n[g];
        fprintf(stderr, "Average distance for '%s': %.3f nm\n",
                sel[g]->name, d.ave[g]);
    }

    /* Here, we could free some memory, but this usually not necessary as we
     * are going to quit anyways. */

    thanx(stderr);
    return 0;
}

/*! \brief
 * The main function.
 *
 * In Gromacs, most analysis programs are implemented such that the \p main
 * function is only a wrapper for a \p gmx_something function that does all
 * the work, and that convention is also followed here. */
int
main(int argc, char *argv[])
{
    gmx_template(argc, argv);
    return 0;
}
