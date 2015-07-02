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

static void process_pull_dim(char *dim_buf, ivec dim)
{
    int           ndim, d, nchar, c;
    char         *ptr, pulldim1[STRLEN];
    t_pull_coord *pcrd;

    ptr  = dim_buf;
    ndim = 0;
    for (d = 0; d < DIM; d++)
    {
        if (sscanf(ptr, "%s%n", pulldim1, &nchar) != 1)
        {
            gmx_fatal(FARGS, "Less than 3 pull dimensions given in pull_dim: '%s'",
                      dim_buf);
        }

        if (gmx_strncasecmp(pulldim1, "N", 1) == 0)
        {
            dim[d] = 0;
        }
        else if (gmx_strncasecmp(pulldim1, "Y", 1) == 0)
        {
            dim[d] = 1;
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
        gmx_fatal(FARGS, "All entries in pull dim are N");
    }
}

static void init_pull_coord(t_pull_coord *pcrd,
                            char *dim_buf,
                            const char *origin_buf, const char *vec_buf,
                            warninp_t wi)
{
    int    m;
    dvec   origin, vec;
    char   buf[STRLEN];

    if (pcrd->eType == epullCONSTRAINT && (pcrd->eGeom == epullgCYL ||
                                           pcrd->eGeom == epullgDIRRELATIVE))
    {
        gmx_fatal(FARGS, "Pulling of type %s can not be combined with geometry %s. Consider using pull type %s.",
                  epull_names[pcrd->eType],
                  epullg_names[pcrd->eGeom],
                  epull_names[epullUMBRELLA]);
    }

    process_pull_dim(dim_buf, pcrd->dim);

    string2dvec(origin_buf, origin);
    if (pcrd->group[0] != 0 && dnorm(origin) > 0)
    {
        gmx_fatal(FARGS, "The pull origin can only be set with an absolute reference");
    }

    /* Check and set the pull vector */
    clear_dvec(vec);
    if (pcrd->eGeom == epullgDIST)
    {
        if (pcrd->init < 0)
        {
            sprintf(buf, "The initial pull distance is negative with geometry %s, while a distance can not be negative. Use geometry %s instead.",
                    EPULLGEOM(pcrd->eGeom), EPULLGEOM(epullgDIR));
            warning_error(wi, buf);
        }
        /* TODO: With a positive init but a negative rate things could still
         * go wrong, but it might be fine if you don't pull too far.
         * We should give a warning or note when there is only one pull dim
         * active, since that is usually the problematic case when you should
         * be using direction. We will do this later, since an already planned
         * generalization of the pull code makes pull dim available here.
         */
    }
    else if (pcrd->eGeom != epullgDIRRELATIVE)
    {
        string2dvec(vec_buf, vec);
        if (dnorm2(vec) == 0)
        {
            gmx_fatal(FARGS, "With pull geometry %s the pull vector can not be 0,0,0",
                      epullg_names[pcrd->eGeom]);
        }
        if (pcrd->eGeom == epullgDIR || pcrd->eGeom == epullgCYL)
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
                       pull_params_t *pull,
                       warninp_t wi)
{
    int           ninp, i, nchar, nscan, m, idum;
    t_inpfile    *inp;
    const char   *tmp;
    char        **grpbuf;
    char          dummy[STRLEN], buf[STRLEN], groups[STRLEN], dim_buf[STRLEN];
    char          init[STRLEN];
    const char   *init_def1 = "0.0", *init_def3 = "0.0 0.0 0.0";
    char          wbuf[STRLEN], origin_buf[STRLEN], vec_buf[STRLEN];

    t_pull_group *pgrp;
    t_pull_coord *pcrd;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read pull parameters */
    CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
    RTYPE("pull-cylinder-r",  pull->cylinder_r, 1.5);
    RTYPE("pull-constr-tol",  pull->constr_tol, 1E-6);
    EETYPE("pull-print-com1", pull->bPrintCOM1, yesno_names);
    EETYPE("pull-print-com2", pull->bPrintCOM2, yesno_names);
    EETYPE("pull-print-ref-value", pull->bPrintRefValue, yesno_names);
    EETYPE("pull-print-components", pull->bPrintComp, yesno_names);
    ITYPE("pull-nstxout",     pull->nstxout, 50);
    ITYPE("pull-nstfout",     pull->nstfout, 50);
    CTYPE("Number of pull groups");
    ITYPE("pull-ngroups",     pull->ngroup, 1);
    CTYPE("Number of pull coordinates");
    ITYPE("pull-ncoords",     pull->ncoord, 1);

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
        int ngroup;

        pcrd = &pull->coord[i-1];
        sprintf(buf, "pull-coord%d-type", i);
        EETYPE(buf,             pcrd->eType, epull_names);
        sprintf(buf, "pull-coord%d-geometry", i);
        EETYPE(buf,             pcrd->eGeom, epullg_names);
        sprintf(buf, "pull-coord%d-groups", i);
        STYPE(buf,              groups, "");

        nscan  = sscanf(groups, "%d %d %d %d %d", &pcrd->group[0], &pcrd->group[1],  &pcrd->group[2], &pcrd->group[3], &idum);
        ngroup = (pcrd->eGeom == epullgDIRRELATIVE) ? 4 : 2;
        if (nscan != ngroup)
        {
            sprintf(wbuf, "%s should contain %d pull group indices with geometry %s",
                    buf, ngroup, epullg_names[pcrd->eGeom]);
            set_warning_line(wi, NULL, -1);
            warning_error(wi, wbuf);
        }

        sprintf(buf, "pull-coord%d-dim", i);
        STYPE(buf,              dim_buf,     "Y Y Y");
        sprintf(buf, "pull-coord%d-origin", i);
        STYPE(buf,              origin_buf,  "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-vec", i);
        STYPE(buf,              vec_buf,     "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-start", i);
        EETYPE(buf,             pcrd->bStart, yesno_names);
        sprintf(buf, "pull-coord%d-init", i);
        RTYPE(buf,              pcrd->init,  0.0);
        sprintf(buf, "pull-coord%d-rate", i);
        RTYPE(buf,              pcrd->rate,  0.0);
        sprintf(buf, "pull-coord%d-k", i);
        RTYPE(buf,              pcrd->k,     0.0);
        sprintf(buf, "pull-coord%d-kB", i);
        RTYPE(buf,              pcrd->kB,    pcrd->k);

        /* Initialize the pull coordinate */
        init_pull_coord(pcrd, dim_buf, origin_buf, vec_buf, wi);
    }

    *ninp_p   = ninp;
    *inp_p    = inp;

    return grpbuf;
}

void make_pull_groups(pull_params_t *pull,
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
            gmx_fatal(FARGS, "Pull option pull_group%d required by grompp has not been set.", g);
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

void make_pull_coords(pull_params_t *pull)
{
    int           c, d;
    t_pull_coord *pcrd;

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

        if (pcrd->eGeom == epullgCYL)
        {
            if (pull->group[pcrd->group[0]].nweight > 0)
            {
                gmx_fatal(FARGS, "Weights are not supported for the reference group with cylinder pulling");
            }
        }

        if (pcrd->eGeom != epullgDIST)
        {
            for (d = 0; d < DIM; d++)
            {
                if (pcrd->vec[d] != 0 && pcrd->dim[d] == 0)
                {
                    gmx_fatal(FARGS, "ERROR: pull-group%d-vec has non-zero %c-component while pull_dim for the %c-dimension is N\n", c+1, 'x'+d, 'x'+d);
                }
            }
        }

        if ((pcrd->eGeom == epullgDIR || pcrd->eGeom == epullgCYL) &&
            norm2(pcrd->vec) == 0)
        {
            gmx_fatal(FARGS, "pull-group%d-vec can not be zero with geometry %s",
                      c+1, EPULLGEOM(pcrd->eGeom));
        }
    }
}

void set_pull_init(t_inputrec *ir, gmx_mtop_t *mtop, rvec *x, matrix box, real lambda,
                   const output_env_t oenv)
{
    pull_params_t *pull;
    struct pull_t *pull_work;
    t_mdatoms     *md;
    t_pbc          pbc;
    int            c;
    double         t_start;

    pull      = ir->pull;
    pull_work = init_pull(NULL, pull, ir, 0, NULL, mtop, NULL, oenv, lambda, FALSE, 0);
    md        = init_mdatoms(NULL, mtop, ir->efep);
    atoms2md(mtop, ir, 0, NULL, mtop->natoms, md);
    if (ir->efep)
    {
        update_mdatoms(md, lambda);
    }

    set_pbc(&pbc, ir->ePBC, box);

    t_start = ir->init_t + ir->init_step*ir->delta_t;

    pull_calc_coms(NULL, pull_work, md, &pbc, t_start, x, NULL);

    fprintf(stderr, "Pull group  natoms  pbc atom  distance at start  reference at t=0\n");
    for (c = 0; c < pull->ncoord; c++)
    {
        t_pull_coord *pcrd;
        t_pull_group *pgrp0, *pgrp1;
        double        value;
        real          init = 0;

        pcrd  = &pull->coord[c];

        pgrp0 = &pull->group[pcrd->group[0]];
        pgrp1 = &pull->group[pcrd->group[1]];
        fprintf(stderr, "%8d  %8d  %8d\n",
                pcrd->group[0], pgrp0->nat, pgrp0->pbcatom+1);
        fprintf(stderr, "%8d  %8d  %8d ",
                pcrd->group[1], pgrp1->nat, pgrp1->pbcatom+1);

        if (pcrd->bStart)
        {
            init       = pcrd->init;
            pcrd->init = 0;
        }

        get_pull_coord_value(pull_work, c, &pbc, &value);
        fprintf(stderr, " %10.3f nm", value);

        if (pcrd->bStart)
        {
            pcrd->init = value + init;
        }
        fprintf(stderr, "     %10.3f nm\n", pcrd->init);
    }

    finish_pull(pull_work);
}
