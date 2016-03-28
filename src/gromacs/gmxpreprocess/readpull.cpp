/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/topology/topology.h"
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
    int    n;

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

static void process_pull_dim(char *dim_buf, ivec dim, const t_pull_coord *pcrd)
{
    int           ndim, d, nchar;
    char         *ptr, pulldim1[STRLEN];

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
    if ((pcrd->eGeom == epullgDIHEDRAL) && (ndim < 3))
    {
        gmx_fatal(FARGS, "Pull geometry dihedral is only useful with pull-dim = Y Y Y");
    }
    if ((pcrd->eGeom == epullgANGLE || pcrd->eGeom == epullgANGLEAXIS ) && (ndim < 2))
    {
        gmx_fatal(FARGS, "Pull geometry %s is only useful with pull-dim = Y for at least 2 dimensions",
                  EPULLGEOM(pcrd->eGeom));
    }
}

static void init_pull_coord(t_pull_coord *pcrd, int coord_index_for_output,
                            char *dim_buf,
                            const char *origin_buf, const char *vec_buf,
                            warninp_t wi)
{
    int    m;
    dvec   origin, vec;
    char   buf[STRLEN];

    if (pcrd->eType == epullCONSTRAINT && (pcrd->eGeom == epullgCYL ||
                                           pcrd->eGeom == epullgDIRRELATIVE ||
                                           pcrd->eGeom == epullgANGLE ||
                                           pcrd->eGeom == epullgANGLEAXIS ||
                                           pcrd->eGeom == epullgDIHEDRAL))
    {
        gmx_fatal(FARGS, "Pulling of type %s can not be combined with geometry %s. Consider using pull type %s.",
                  epull_names[pcrd->eType],
                  epullg_names[pcrd->eGeom],
                  epull_names[epullUMBRELLA]);
    }

    if (pcrd->eType == epullEXTERNAL)
    {
        if (pcrd->externalPotentialProvider[0] == '\0')
        {
            sprintf(buf, "The use of pull type '%s' for pull coordinate %d requires that the name of the module providing the potential external is set with the option %s%d%s",
                    epull_names[pcrd->eType], coord_index_for_output,
                    "pull-coord", coord_index_for_output, "-potential-provider");
            warning_error(wi, buf);
        }

        if (pcrd->rate != 0)
        {
            sprintf(buf, "The use of pull type '%s' for pull coordinate %d requires that the pull rate is zero",
                    epull_names[pcrd->eType], coord_index_for_output);
            warning_error(wi, buf);
        }

        if (pcrd->eGeom == epullgCYL)
        {
            /* Warn the user of a PBC restriction, caused by the fact that
             * there is no reference value with an external pull potential.
             */
            sprintf(buf, "With pull type '%s' and geometry '%s', the distance component along the cylinder axis between atoms in the cylinder group and the COM of the pull group should be smaller than half the box length",
                    epull_names[pcrd->eType], epullg_names[pcrd->eGeom]);
            warning_note(wi, buf);
        }
    }

    process_pull_dim(dim_buf, pcrd->dim, pcrd);

    string2dvec(origin_buf, origin);
    if (pcrd->group[0] != 0 && dnorm(origin) > 0)
    {
        gmx_fatal(FARGS, "The pull origin can only be set with an absolute reference");
    }

    /* Check the given initial reference value and warn for dangerous values */
    if (pcrd->eGeom == epullgDIST)
    {
        if (pcrd->bStart && pcrd->init < 0)
        {
            sprintf(buf, "The initial reference distance set by pull-coord-init is set to a negative value (%g) with geometry %s while distances need to be non-negative. "
                    "This may work, since you have set pull-coord-start to 'yes' which modifies this value, but only for certain starting distances. "
                    "If this is a mistake you may want to use geometry %s instead.",
                    pcrd->init, EPULLGEOM(pcrd->eGeom), EPULLGEOM(epullgDIR));
            warning(wi, buf);
        }
    }
    else if (pcrd->eGeom == epullgANGLE || pcrd->eGeom == epullgANGLEAXIS)
    {
        if (pcrd->bStart && (pcrd->init < 0 || pcrd->init > 180))
        {
            /* This value of pcrd->init may be ok depending on pcrd->bStart which modifies pcrd->init later on */
            sprintf(buf, "The initial reference angle set by pull-coord-init (%g) is outside of the allowed range [0, 180] degrees for geometry (%s). "
                    "This may work, since you have set pull-coord-start to 'yes' which modifies this value, but only for certain starting angles.",
                    pcrd->init, EPULLGEOM(pcrd->eGeom));
            warning(wi, buf);
        }
    }
    else if (pcrd->eGeom == epullgDIHEDRAL)
    {
        if (pcrd->bStart && (pcrd->init < -180 || pcrd->init > 180))
        {
            sprintf(buf, "The initial reference angle set by pull-coord-init (%g) is outside of the allowed range [-180, 180] degrees for geometry (%s). "
                    "This may work, since you have set pull-coord-start to 'yes' which modifies this value, but only for certain starting angles.",
                    pcrd->init, EPULLGEOM(pcrd->eGeom));
            warning(wi, buf);
        }
    }

    /* Check and set the pull vector */
    clear_dvec(vec);
    string2dvec(vec_buf, vec);

    if (pcrd->eGeom == epullgDIR || pcrd->eGeom == epullgCYL || pcrd->eGeom == epullgDIRPBC || pcrd->eGeom == epullgANGLEAXIS)
    {
        if (dnorm2(vec) == 0)
        {
            gmx_fatal(FARGS, "With pull geometry %s the pull vector can not be 0,0,0",
                      epullg_names[pcrd->eGeom]);
        }
        for (int d = 0; d < DIM; d++)
        {
            if (vec[d] != 0 && pcrd->dim[d] == 0)
            {
                gmx_fatal(FARGS, "pull-coord-vec has non-zero %c-component while pull_dim for the %c-dimension is set to N", 'x'+d, 'x'+d);
            }
        }

        /* Normalize the direction vector */
        dsvmul(1/dnorm(vec), vec, vec);
    }
    else /* This case is for are all the geometries where the pull vector is not used */
    {
        if (dnorm2(vec) > 0)
        {
            sprintf(buf, "A pull vector is given (%g  %g  %g) but will not be used with geometry %s. If you really want to use this "
                    "vector, consider using geometry %s instead.",
                    vec[0], vec[1], vec[2], EPULLGEOM(pcrd->eGeom),
                    pcrd->eGeom == epullgANGLE ? EPULLGEOM(epullgANGLEAXIS) : EPULLGEOM(epullgDIR));
            warning(wi, buf);
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
    int           ninp, nscan, idum;
    t_inpfile    *inp;
    const char   *tmp;
    char        **grpbuf;
    char          buf[STRLEN];
    char          provider[STRLEN], groups[STRLEN], dim_buf[STRLEN];
    char          wbuf[STRLEN], origin_buf[STRLEN], vec_buf[STRLEN];

    t_pull_group *pgrp;
    t_pull_coord *pcrd;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read pull parameters */
    CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
    RTYPE("pull-cylinder-r",  pull->cylinder_r, 1.5);
    RTYPE("pull-constr-tol",  pull->constr_tol, 1E-6);
    EETYPE("pull-print-com", pull->bPrintCOM, yesno_names);
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
    CTYPE("Group and coordinate parameters");

    /* Read the pull groups */
    snew(grpbuf, pull->ngroup);
    /* Group 0 is the absolute reference, we don't read anything for 0 */
    for (int groupNum = 1; groupNum < pull->ngroup; groupNum++)
    {
        pgrp = &pull->group[groupNum];
        snew(grpbuf[groupNum], STRLEN);
        sprintf(buf, "pull-group%d-name", groupNum);
        STYPE(buf,              grpbuf[groupNum], "");
        sprintf(buf, "pull-group%d-weights", groupNum);
        STYPE(buf,              wbuf, "");
        sprintf(buf, "pull-group%d-pbcatom", groupNum);
        ITYPE(buf,              pgrp->pbcatom, 0);

        /* Initialize the pull group */
        init_pull_group(pgrp, wbuf);
    }

    /* Read the pull coordinates */
    for (int coordNum = 1; coordNum < pull->ncoord + 1; coordNum++)
    {
        pcrd = &pull->coord[coordNum - 1];
        sprintf(buf, "pull-coord%d-type", coordNum);
        EETYPE(buf,             pcrd->eType, epull_names);
        sprintf(buf, "pull-coord%d-potential-provider", coordNum);
        STYPE(buf,              provider, "");
        pcrd->externalPotentialProvider = gmx_strdup(provider);
        sprintf(buf, "pull-coord%d-geometry", coordNum);
        EETYPE(buf,             pcrd->eGeom, epullg_names);
        sprintf(buf, "pull-coord%d-groups", coordNum);
        STYPE(buf,              groups, "");

        switch (pcrd->eGeom)
        {
            case epullgDIHEDRAL:
                pcrd->ngroup = 6; break;
            case epullgDIRRELATIVE:
            case epullgANGLE:
                pcrd->ngroup = 4; break;
            default:
                pcrd->ngroup = 2; break;
        }

        nscan = sscanf(groups, "%d %d %d %d %d %d %d",
                       &pcrd->group[0], &pcrd->group[1], &pcrd->group[2], &pcrd->group[3],
                       &pcrd->group[4], &pcrd->group[5], &idum);
        if (nscan != pcrd->ngroup)
        {
            sprintf(wbuf, "%s should contain %d pull group indices with geometry %s",
                    buf, pcrd->ngroup, epullg_names[pcrd->eGeom]);
            set_warning_line(wi, nullptr, -1);
            warning_error(wi, wbuf);
        }
        for (int g = 0; g < pcrd->ngroup; g++)
        {
            if (pcrd->group[g] < 0 || pcrd->group[g] >= pull->ngroup)
            {
                /* Quit with a fatal error to avoid invalid memory access */
                gmx_fatal(FARGS, "%s contains an invalid pull group %d, you should have %d <= group <= %d",
                          buf, pcrd->group[g], 0, pull->ngroup - 1);
            }
        }

        sprintf(buf, "pull-coord%d-dim", coordNum);
        STYPE(buf,              dim_buf,     "Y Y Y");
        sprintf(buf, "pull-coord%d-origin", coordNum);
        STYPE(buf,              origin_buf,  "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-vec", coordNum);
        STYPE(buf,              vec_buf,     "0.0 0.0 0.0");
        sprintf(buf, "pull-coord%d-start", coordNum);
        EETYPE(buf,             pcrd->bStart, yesno_names);
        sprintf(buf, "pull-coord%d-init", coordNum);
        RTYPE(buf,              pcrd->init,  0.0);
        sprintf(buf, "pull-coord%d-rate", coordNum);
        RTYPE(buf,              pcrd->rate,  0.0);
        sprintf(buf, "pull-coord%d-k", coordNum);
        RTYPE(buf,              pcrd->k,     0.0);
        sprintf(buf, "pull-coord%d-kB", coordNum);
        RTYPE(buf,              pcrd->kB,    pcrd->k);

        /* Initialize the pull coordinate */
        init_pull_coord(pcrd, coordNum, dim_buf, origin_buf, vec_buf, wi);
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
    int           c;
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
    }
}

pull_t *set_pull_init(t_inputrec *ir, const gmx_mtop_t *mtop,
                      rvec *x, matrix box, real lambda,
                      const gmx_output_env_t *oenv)
{
    pull_params_t *pull;
    pull_t        *pull_work;
    t_pbc          pbc;
    int            c;
    double         t_start;

    pull      = ir->pull;
    pull_work = init_pull(nullptr, pull, ir, 0, nullptr, mtop, nullptr, oenv, lambda, FALSE, ContinuationOptions());
    auto mdAtoms = gmx::makeMDAtoms(nullptr, *mtop, *ir, false);
    auto md      = mdAtoms->mdatoms();
    atoms2md(mtop, ir, -1, nullptr, mtop->natoms, mdAtoms.get());
    if (ir->efep)
    {
        update_mdatoms(md, lambda);
    }

    set_pbc(&pbc, ir->ePBC, box);

    t_start = ir->init_t + ir->init_step*ir->delta_t;

    pull_calc_coms(nullptr, pull_work, md, &pbc, t_start, x, nullptr);

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

        value  = get_pull_coord_value(pull_work, c, &pbc);

        value *= pull_conversion_factor_internal2userinput(pcrd);
        fprintf(stderr, " %10.3f %s", value, pull_coordinate_units(pcrd));

        if (pcrd->bStart)
        {
            pcrd->init = value + init;
        }

        if (pcrd->eGeom == epullgDIST)
        {
            if (pcrd->init < 0)
            {
                gmx_fatal(FARGS, "The initial pull distance (%g) needs to be non-negative with geometry %s. If you want a signed distance, use geometry %s instead.",
                          pcrd->init, EPULLGEOM(pcrd->eGeom), EPULLGEOM(epullgDIR));
            }

            /* TODO: With a positive init but a negative rate things could still
             * go wrong, but it might be fine if you don't pull too far.
             * We should give a warning or note when there is only one pull dim
             * active, since that is usually the problematic case when you should
             * be using direction. We will do this later, since an already planned
             * generalization of the pull code makes pull dim available here.
             */
        }
        else if (pcrd->eGeom == epullgANGLE || pcrd->eGeom == epullgANGLEAXIS)
        {
            if (pcrd->init < 0 || pcrd->init > 180)
            {
                gmx_fatal(FARGS,  "The initial pull reference angle (%g) is outside of the allowed range [0, 180] degrees.", pcrd->init);
            }
        }
        else if (pcrd->eGeom == epullgDIHEDRAL)
        {
            if (pcrd->init < -180 || pcrd->init > 180)
            {
                gmx_fatal(FARGS,  "The initial pull reference angle (%g) is outside of the allowed range [-180, 180] degrees.",
                          pcrd->init);
            }
        }


        fprintf(stderr, "     %10.3f %s\n", pcrd->init, pull_coordinate_units(pcrd));
    }

    return pull_work;
}
