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
#include "gmxpre.h"

#include <stdlib.h>
#include <string.h>

#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


static char pulldim[STRLEN];

static void string2dvec(const char buf[], dvec nums)
{
    double dum;

    if (sscanf(buf, "%lf%lf%lf%lf", &nums[0], &nums[1], &nums[2], &dum) != 3)
    {
        gmx_fatal(FARGS, "Expected three numbers at input line %s", buf);
    }
}

static void init_pull_group(t_pull_group *pg,
                            const char   *wbuf)
{
    double d;
    int    n, m;

    pg->nweight = 0;
    while (sscanf(wbuf, "%lf %n", &d, &n) == 1)
    {
        if (pg->nweight % 100 == 0)
        {
            srenew(pg->weight, pg->nweight+100);
        }
        pg->weight[pg->nweight++] = d;
        wbuf += n;
    }
}

static void init_pull_coord(t_pull_coord *pcrd, int eGeom,
                            const char *origin_buf, const char *vec_buf)
{
    int    m;
    dvec   origin, vec;

    string2dvec(origin_buf, origin);
    if (pcrd->group[0] != 0 && dnorm(origin) > 0)
    {
        gmx_fatal(FARGS, "The pull origin can only be set with an absolute reference");
    }

    if (eGeom == epullgDIST)
    {
        clear_dvec(vec);
    }
    else
    {
        string2dvec(vec_buf, vec);
        if (eGeom == epullgDIR || eGeom == epullgCYL)
        {
            /* Normalize the direction vector */
            dsvmul(1/dnorm(vec), vec, vec);
        }
    }
    for (m = 0; m < DIM; m++)
    {
        pcrd->origin[m] = origin[m];
        pcrd->vec[m]    = vec[m];
    }
}

char **read_pullparams(int *ninp_p, t_inpfile **inp_p,
                       t_pull *pull, gmx_bool *bStart,
                       warninp_t wi)
{
    int           ninp, nerror = 0, i, nchar, nscan, m, idum;
    t_inpfile    *inp;
    const char   *tmp;
    char        **grpbuf;
    char          dummy[STRLEN], buf[STRLEN], groups[STRLEN], init[STRLEN];
    const char   *init_def1 = "0.0", *init_def3 = "0.0 0.0 0.0";
    char          wbuf[STRLEN], origin_buf[STRLEN], vec_buf[STRLEN];

    t_pull_group *pgrp;
    t_pull_coord *pcrd;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read pull parameters */
    CTYPE("Pull geometry: distance, direction, direction-periodic or cylinder");
    EETYPE("pull-geometry",   pull->eGeom, epullg_names);
    CTYPE("Select components for the pull vector. default: Y Y Y");
    STYPE("pull-dim",         pulldim, "Y Y Y");
    CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
    RTYPE("pull-r1",          pull->cyl_r1, 1.0);
    CTYPE("Switch from r1 to r0 in case of dynamic reaction force");
    RTYPE("pull-r0",          pull->cyl_r0, 1.5);
    RTYPE("pull-constr-tol",  pull->constr_tol, 1E-6);
    EETYPE("pull-start",      *bStart, yesno_names);
    EETYPE("pull-print-reference", pull->bPrintRef, yesno_names);
    ITYPE("pull-nstxout",     pull->nstxout, 10);
    ITYPE("pull-nstfout",     pull->nstfout,  1);
    CTYPE("Number of pull groups");
    ITYPE("pull-ngroups",     pull->ngroup, 1);
    CTYPE("Number of pull coordinates");
    ITYPE("pull-ncoords",     pull->ncoord, 1);

    if (pull->cyl_r1 > pull->cyl_r0)
    {
        warning_error(wi, "pull-r1 > pull_r0");
    }

    if (pull->ngroup < 1)
    {
        gmx_fatal(FARGS, "pull-ngroups should be >= 1");
    }
    /* We always add an absolute reference group (index 0), even if not used */
    pull->ngroup += 1;

    if (pull->ncoord < 1)
    {
        gmx_fatal(FARGS, "pull-ncoords should be >= 1");
    }

    snew(pull->group, pull->ngroup);

    snew(pull->coord, pull->ncoord);

    /* pull group options */
    CTYPE("Group name, weight (default all 1), vector, init, rate (nm/ps), kJ/(mol*nm^2)");

    /* Read the pull groups */
    snew(grpbuf, pull->ngroup);
    /* Group 0 is the absolute reference, we don't read anything for 0 */
    for (i = 1; i < pull->ngroup; i++)
    {
        pgrp = &pull->group[i];
        snew(grpbuf[i], STRLEN);
        sprintf(buf, "pull-group%d-name", i);
        STYPE(buf,              grpbuf[i], "");
        sprintf(buf, "pull-group%d-weights", i);
        STYPE(buf,              wbuf, "");
        sprintf(buf, "pull-group%d-pbcatom", i);
        ITYPE(buf,              pgrp->pbcatom, 0);

        /* Initialize the pull group */
        init_pull_group(pgrp, wbuf);
    }

    /* Read the pull coordinates */
    for (i = 1; i < pull->ncoord + 1; i++)
    {
        pcrd = &pull->coord[i-1];
        sprintf(buf, "pull-coord%d-groups", i);
        STYPE(buf,              groups, "");
        nscan = sscanf(groups, "%d %d %d", &pcrd->group[0], &pcrd->group[1], &idum);
        if (nscan != 2)
        {
            fprintf(stderr, "ERROR: %s should have %d components\n", buf, 2);
            nerror++;
        }
        sprintf(buf, "pull-coord%d-origin", i);
        STYPE(buf,              origin_buf, "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-vec", i);
        STYPE(buf,              vec_buf,    "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-init", i);
        RTYPE(buf,              pcrd->init, 0.0);
        sprintf(buf, "pull-coord%d-rate", i);
        RTYPE(buf,              pcrd->rate, 0.0);
        sprintf(buf, "pull-coord%d-k", i);
        RTYPE(buf,              pcrd->k, 0.0);
        sprintf(buf, "pull-coord%d-kB", i);
        RTYPE(buf,              pcrd->kB, pcrd->k);

        /* Initialize the pull coordinate */
        init_pull_coord(pcrd, pull->eGeom, origin_buf, vec_buf);
    }

    *ninp_p   = ninp;
    *inp_p    = inp;

    return grpbuf;
}

void make_pull_groups(t_pull *pull,
                      char **pgnames,
                      const t_blocka *grps, char **gnames)
{
    int           g, ig = -1, i;
    t_pull_group *pgrp;

    /* Absolute reference group (might not be used) is special */
    pgrp          = &pull->group[0];
    pgrp->nat     = 0;
    pgrp->pbcatom = -1;

    for (g = 1; g < pull->ngroup; g++)
    {
        pgrp = &pull->group[g];

        if (strcmp(pgnames[g], "") == 0)
        {
            gmx_fatal(FARGS, "Group pull_group%d required by grompp was undefined.", g);
        }

        ig        = search_string(pgnames[g], grps->nr, gnames);
        pgrp->nat = grps->index[ig+1] - grps->index[ig];

        fprintf(stderr, "Pull group %d '%s' has %d atoms\n",
                g, pgnames[g], pgrp->nat);

        if (pgrp->nat == 0)
        {
            gmx_fatal(FARGS, "Pull group %d '%s' is empty", g, pgnames[g]);
        }

        snew(pgrp->ind, pgrp->nat);
        for (i = 0; i < pgrp->nat; i++)
        {
            pgrp->ind[i] = grps->a[grps->index[ig]+i];
        }

        if (pull->eGeom == epullgCYL && g == 1 && pgrp->nweight > 0)
        {
            gmx_fatal(FARGS, "Weights are not supported for the reference group with cylinder pulling");
        }
        if (pgrp->nweight > 0 && pgrp->nweight != pgrp->nat)
        {
            gmx_fatal(FARGS, "Number of weights (%d) for pull group %d '%s' does not match the number of atoms (%d)",
                      pgrp->nweight, g, pgnames[g], pgrp->nat);
        }

        if (pgrp->nat == 1)
        {
            /* No pbc is required for this group */
            pgrp->pbcatom = -1;
        }
        else
        {
            if (pgrp->pbcatom > 0)
            {
                pgrp->pbcatom -= 1;
            }
            else if (pgrp->pbcatom == 0)
            {
                pgrp->pbcatom = pgrp->ind[(pgrp->nat-1)/2];
            }
            else
            {
                /* Use cosine weighting */
                pgrp->pbcatom = -1;
            }
        }
    }
}

void make_pull_coords(t_pull *pull)
{
    int           ndim, d, nchar, c;
    char         *ptr, pulldim1[STRLEN];
    t_pull_coord *pcrd;

    ptr  = pulldim;
    ndim = 0;
    for (d = 0; d < DIM; d++)
    {
        if (sscanf(ptr, "%s%n", pulldim1, &nchar) != 1)
        {
            gmx_fatal(FARGS, "Less than 3 pull dimensions given in pull_dim: '%s'",
                      pulldim);
        }

        if (gmx_strncasecmp(pulldim1, "N", 1) == 0)
        {
            pull->dim[d] = 0;
        }
        else if (gmx_strncasecmp(pulldim1, "Y", 1) == 0)
        {
            pull->dim[d] = 1;
            ndim++;
        }
        else
        {
            gmx_fatal(FARGS, "Please use Y(ES) or N(O) for pull_dim only (not %s)",
                      pulldim1);
        }
        ptr += nchar;
    }
    if (ndim == 0)
    {
        gmx_fatal(FARGS, "All entries in pull_dim are N");
    }

    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd = &pull->coord[c];

        if (pcrd->group[0] < 0 || pcrd->group[0] >= pull->ngroup ||
            pcrd->group[1] < 0 || pcrd->group[1] >= pull->ngroup)
        {
            gmx_fatal(FARGS, "Pull group index in pull-coord%d-groups out of range, should be between %d and %d", c+1, 0, pull->ngroup+1);
        }

        if (pcrd->group[0] == pcrd->group[1])
        {
            gmx_fatal(FARGS, "Identical pull group indices in pull-coord%d-groups", c+1);
        }

        if (pull->eGeom == epullgCYL && pcrd->group[0] != 1)
        {
            gmx_fatal(FARGS, "With pull geometry %s, the first pull group should always be 1", EPULLGEOM(pull->eGeom));
        }

        if (pull->eGeom != epullgDIST)
        {
            for (d = 0; d < DIM; d++)
            {
                if (pcrd->vec[d] != 0 && pull->dim[d] == 0)
                {
                    gmx_fatal(FARGS, "ERROR: pull-group%d-vec has non-zero %c-component while pull_dim for the %c-dimension is N\n", c+1, 'x'+d, 'x'+d);
                }
            }
        }

        if ((pull->eGeom == epullgDIR || pull->eGeom == epullgCYL) &&
            norm2(pcrd->vec) == 0)
        {
            gmx_fatal(FARGS, "pull-group%d-vec can not be zero with geometry %s",
                      c+1, EPULLGEOM(pull->eGeom));
        }
    }
}

void set_pull_init(t_inputrec *ir, gmx_mtop_t *mtop, rvec *x, matrix box, real lambda,
                   const output_env_t oenv, gmx_bool bStart)
{
    t_mdatoms    *md;
    t_pull       *pull;
    t_pull_coord *pcrd;
    t_pull_group *pgrp0, *pgrp1;
    t_pbc         pbc;
    int           c, m;
    double        t_start, tinvrate;
    real          init;
    dvec          dr;
    double        dev;

    init_pull(NULL, ir, 0, NULL, mtop, NULL, oenv, lambda, FALSE, 0);
    md = init_mdatoms(NULL, mtop, ir->efep);
    atoms2md(mtop, ir, 0, NULL, mtop->natoms, md);
    if (ir->efep)
    {
        update_mdatoms(md, lambda);
    }
    pull = ir->pull;

    set_pbc(&pbc, ir->ePBC, box);

    t_start = ir->init_t + ir->init_step*ir->delta_t;

    pull_calc_coms(NULL, pull, md, &pbc, t_start, x, NULL);

    fprintf(stderr, "Pull group  natoms  pbc atom  distance at start     reference at t=0\n");
    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd  = &pull->coord[c];

        pgrp0 = &pull->group[pcrd->group[0]];
        pgrp1 = &pull->group[pcrd->group[1]];
        fprintf(stderr, "%8d  %8d  %8d\n",
                pcrd->group[0], pgrp0->nat, pgrp0->pbcatom+1);
        fprintf(stderr, "%8d  %8d  %8d ",
                pcrd->group[1], pgrp1->nat, pgrp1->pbcatom+1);

        init       = pcrd->init;
        pcrd->init = 0;

        if (pcrd->rate == 0)
        {
            tinvrate = 0;
        }
        else
        {
            tinvrate = t_start/pcrd->rate;
        }
        get_pull_coord_distance(pull, c, &pbc, 0, dr, &dev);
        fprintf(stderr, " %6.3f ", dev);

        if (bStart)
        {
            pcrd->init = init + dev - tinvrate;
        }
        else
        {
            pcrd->init = init;
        }
        fprintf(stderr, " %6.3f\n", pcrd->init);
    }
}
