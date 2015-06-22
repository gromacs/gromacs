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

#include "pull.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

static std::string append_before_extension(std::string pathname,
                                           std::string to_append)
{
    /* Appends to_append before last '.' in pathname */
    size_t extPos = pathname.find_last_of('.');
    if (extPos == std::string::npos)
    {
        return pathname + to_append;
    }
    else
    {
        return pathname.substr(0, extPos) + to_append +
               pathname.substr(extPos, std::string::npos);
    }
}

static void pull_print_group_x(FILE *out, ivec dim,
                               const pull_group_work_t *pgrp)
{
    int m;

    for (m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", pgrp->x[m]);
        }
    }
}

static void pull_print_coord_dr(FILE *out, const pull_coord_work_t *pcrd,
                                gmx_bool bPrintRefValue,
                                gmx_bool bPrintComponents)
{
    fprintf(out, "\t%g", pcrd->value);

    if (bPrintRefValue)
    {
        fprintf(out, "\t%g", pcrd->value_ref);
    }

    if (bPrintComponents)
    {
        int m;

        for (m = 0; m < DIM; m++)
        {
            if (pcrd->params.dim[m])
            {
                fprintf(out, "\t%g", pcrd->dr[m]);
            }
        }
    }
}

static void pull_print_x(FILE *out, struct pull_t *pull, double t)
{
    int c;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pull->params.bPrintCOM1)
        {
            if (pcrd->params.eGeom == epullgCYL)
            {
                pull_print_group_x(out, pcrd->params.dim, &pull->dyna[c]);
            }
            else
            {
                pull_print_group_x(out, pcrd->params.dim,
                                   &pull->group[pcrd->params.group[0]]);
            }
        }
        if (pull->params.bPrintCOM2)
        {
            pull_print_group_x(out, pcrd->params.dim,
                               &pull->group[pcrd->params.group[1]]);
        }
        pull_print_coord_dr(out, pcrd,
                            pull->params.bPrintRefValue,
                            pull->params.bPrintComp);
    }
    fprintf(out, "\n");
}

static void pull_print_f(FILE *out, struct pull_t *pull, double t)
{
    int c;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        fprintf(out, "\t%g", pull->coord[c].f_scal);
    }
    fprintf(out, "\n");
}

void pull_print_output(struct pull_t *pull, gmx_int64_t step, double time)
{
    if ((pull->params.nstxout != 0) && (step % pull->params.nstxout == 0))
    {
        pull_print_x(pull->out_x, pull, time);
    }

    if ((pull->params.nstfout != 0) && (step % pull->params.nstfout == 0))
    {
        pull_print_f(pull->out_f, pull, time);
    }
}

static FILE *open_pull_out(const char *fn, struct pull_t *pull, const output_env_t oenv,
                           gmx_bool bCoord, unsigned long Flags)
{
    FILE  *fp;
    int    nsets, c, m;
    char **setname, buf[10];

    if (Flags & MD_APPENDFILES)
    {
        fp = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        fp = gmx_fio_fopen(fn, "w+");
        if (bCoord)
        {
            xvgr_header(fp, "Pull COM",  "Time (ps)", "Position (nm)",
                        exvggtXNY, oenv);
        }
        else
        {
            xvgr_header(fp, "Pull force", "Time (ps)", "Force (kJ/mol/nm)",
                        exvggtXNY, oenv);
        }

        /* With default mdp options only the actual distance is printed,
         * but optionally 2 COMs, the reference distance and distance
         * components can also be printed.
         */
        snew(setname, pull->ncoord*(DIM + DIM + 1 + 1 + DIM));
        nsets = 0;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (bCoord)
            {
                /* The order of this legend should match the order of printing
                 * the data in print_pull_x above.
                 */

                if (pull->params.bPrintCOM1)
                {
                    /* Legend for reference group position */
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].params.dim[m])
                        {
                            sprintf(buf, "%d %s%d%c", c+1, "c", 1, 'X'+m);
                            setname[nsets] = gmx_strdup(buf);
                            nsets++;
                        }
                    }
                }
                if (pull->params.bPrintCOM2)
                {
                    /* Legend for reference group position */
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].params.dim[m])
                        {
                            sprintf(buf, "%d %s%d%c", c+1, "c", 2, 'X'+m);
                            setname[nsets] = gmx_strdup(buf);
                            nsets++;
                        }
                    }
                }
                /* The pull coord distance */
                sprintf(buf, "%d", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
                if (pull->params.bPrintRefValue)
                {
                    sprintf(buf, "%c%d", 'r', c+1);
                    setname[nsets] = gmx_strdup(buf);
                    nsets++;
                }
                if (pull->params.bPrintComp)
                {
                    /* The pull coord distance components */
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].params.dim[m])
                        {
                            sprintf(buf, "%d %s%c", c+1, "d", 'X'+m);
                            setname[nsets] = gmx_strdup(buf);
                            nsets++;
                        }
                    }
                }
            }
            else
            {
                /* For the pull force we always only use one scalar */
                sprintf(buf, "%d", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
            }
        }
        if (nsets > 1)
        {
            xvgr_legend(fp, nsets, (const char**)setname, oenv);
        }
        for (c = 0; c < nsets; c++)
        {
            sfree(setname[c]);
        }
        sfree(setname);
    }

    return fp;
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_grp(const pull_group_work_t *pgrp,
                             const t_mdatoms *md,
                             const dvec f_pull, int sign, rvec *f)
{
    int    i, ii, m;
    double wmass, inv_wm;

    inv_wm = pgrp->mwscale;

    if (pgrp->params.nat == 1 && pgrp->nat_loc == 1)
    {
        /* Only one atom and our rank has this atom: we can skip
         * the mass weighting, which means that this code also works
         * for mass=0, e.g. with a virtual site.
         */
        for (m = 0; m < DIM; m++)
        {
            f[pgrp->ind_loc[0]][m] += sign*f_pull[m];
        }
    }
    else
    {
        for (i = 0; i < pgrp->nat_loc; i++)
        {
            ii    = pgrp->ind_loc[i];
            wmass = md->massT[ii];
            if (pgrp->weight_loc)
            {
                wmass *= pgrp->weight_loc[i];
            }

            for (m = 0; m < DIM; m++)
            {
                f[ii][m] += sign * wmass * f_pull[m] * inv_wm;
            }
        }
    }
}

/* Apply forces in a mass weighted fashion to a cylinder group */
static void apply_forces_cyl_grp(const pull_group_work_t *pgrp,
                                 double dv_corr,
                                 const t_mdatoms *md,
                                 const dvec f_pull, double f_scal,
                                 int sign, rvec *f)
{
    int    i, ii, m;
    double mass, weight, inv_wm, dv_com;

    inv_wm = pgrp->mwscale;

    for (i = 0; i < pgrp->nat_loc; i++)
    {
        ii     = pgrp->ind_loc[i];
        mass   = md->massT[ii];
        weight = pgrp->weight_loc[i];
        /* The stored axial distance from the cylinder center (dv) needs
         * to be corrected for an offset (dv_corr), which was unknown when
         * we calculated dv.
         */
        dv_com = pgrp->dv[i] + dv_corr;

        /* Here we not only add the pull force working along vec (f_pull),
         * but also a radial component, due to the dependence of the weights
         * on the radial distance.
         */
        for (m = 0; m < DIM; m++)
        {
            f[ii][m] += sign*inv_wm*(mass*weight*f_pull[m] +
                                     pgrp->mdw[i][m]*dv_com*f_scal);
        }
    }
}

/* Apply torque forces in a mass weighted fashion to the groups that define
 * the pull vector direction for pull coordinate pcrd.
 */
static void apply_forces_vec_torque(const struct pull_t     *pull,
                                    const pull_coord_work_t *pcrd,
                                    const t_mdatoms         *md,
                                    rvec                    *f)
{
    double inpr;
    int    m;
    dvec   f_perp;

    /* The component inpr along the pull vector is accounted for in the usual
     * way. Here we account for the component perpendicular to vec.
     */
    inpr = 0;
    for (m = 0; m < DIM; m++)
    {
        inpr += pcrd->dr[m]*pcrd->vec[m];
    }
    /* The torque force works along the component of the distance vector
     * of between the two "usual" pull groups that is perpendicular to
     * the pull vector. The magnitude of this force is the "usual" scale force
     * multiplied with the ratio of the distance between the two "usual" pull
     * groups and the distance between the two groups that define the vector.
     */
    for (m = 0; m < DIM; m++)
    {
        f_perp[m] = (pcrd->dr[m] - inpr*pcrd->vec[m])/pcrd->vec_len*pcrd->f_scal;
    }

    /* Apply the force to the groups defining the vector using opposite signs */
    apply_forces_grp(&pull->group[pcrd->params.group[2]], md, f_perp, -1, f);
    apply_forces_grp(&pull->group[pcrd->params.group[3]], md, f_perp,  1, f);
}

/* Apply forces in a mass weighted fashion */
static void apply_forces(struct pull_t * pull, t_mdatoms * md, rvec *f)
{
    int c;

    for (c = 0; c < pull->ncoord; c++)
    {
        const pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eGeom == epullgCYL)
        {
            dvec f_tot;
            int  m;

            apply_forces_cyl_grp(&pull->dyna[c], pcrd->cyl_dev, md,
                                 pcrd->f, pcrd->f_scal, -1, f);

            /* Sum the force along the vector and the radial force */
            for (m = 0; m < DIM; m++)
            {
                f_tot[m] = pcrd->f[m] + pcrd->f_scal*pcrd->ffrad[m];
            }
            apply_forces_grp(&pull->group[pcrd->params.group[1]], md, f_tot, 1, f);
        }
        else
        {
            if (pcrd->params.eGeom == epullgDIRRELATIVE)
            {
                /* We need to apply the torque forces to the pull groups
                 * that define the pull vector.
                 */
                apply_forces_vec_torque(pull, pcrd, md, f);
            }

            if (pull->group[pcrd->params.group[0]].params.nat > 0)
            {
                apply_forces_grp(&pull->group[pcrd->params.group[0]], md, pcrd->f, -1, f);
            }
            apply_forces_grp(&pull->group[pcrd->params.group[1]], md, pcrd->f, 1, f);
        }
    }
}

static double max_pull_distance2(const pull_coord_work_t *pcrd,
                                 const t_pbc             *pbc)
{
    double max_d2;
    int    m;

    max_d2 = GMX_DOUBLE_MAX;

    for (m = 0; m < pbc->ndim_ePBC; m++)
    {
        if (pcrd->params.dim[m] != 0)
        {
            max_d2 = std::min(max_d2, static_cast<double>(norm2(pbc->box[m])));
        }
    }

    return 0.25*max_d2;
}

/* This function returns the distance based on coordinates xg and xref.
 * Note that the pull coordinate struct pcrd is not modified.
 */
static void low_get_pull_coord_dr(const struct pull_t *pull,
                                  const pull_coord_work_t *pcrd,
                                  const t_pbc *pbc,
                                  dvec xg, dvec xref, double max_dist2,
                                  dvec dr)
{
    const pull_group_work_t *pgrp0;
    int                      m;
    dvec                     xrefr, dref = {0, 0, 0};
    double                   dr2;

    pgrp0 = &pull->group[pcrd->params.group[0]];

    /* Only the first group can be an absolute reference, in that case nat=0 */
    if (pgrp0->params.nat == 0)
    {
        for (m = 0; m < DIM; m++)
        {
            xref[m] = pcrd->params.origin[m];
        }
    }

    copy_dvec(xref, xrefr);

    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        for (m = 0; m < DIM; m++)
        {
            dref[m] = pcrd->value_ref*pcrd->vec[m];
        }
        /* Add the reference position, so we use the correct periodic image */
        dvec_inc(xrefr, dref);
    }

    pbc_dx_d(pbc, xg, xrefr, dr);
    dr2 = 0;
    for (m = 0; m < DIM; m++)
    {
        dr[m] *= pcrd->params.dim[m];
        dr2   += dr[m]*dr[m];
    }
    if (max_dist2 >= 0 && dr2 > 0.98*0.98*max_dist2)
    {
        gmx_fatal(FARGS, "Distance between pull groups %d and %d (%f nm) is larger than 0.49 times the box size (%f).\nYou might want to consider using \"pull-geometry = direction-periodic\" instead.\n",
                  pcrd->params.group[0], pcrd->params.group[1],
                  sqrt(dr2), sqrt(max_dist2));
    }

    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        dvec_inc(dr, dref);
    }
}

/* This function returns the distance based on the contents of the pull struct.
 * pull->coord[coord_ind].dr, and potentially vector, are updated.
 */
static void get_pull_coord_dr(struct pull_t *pull,
                              int            coord_ind,
                              const t_pbc   *pbc)
{
    double             md2;
    pull_coord_work_t *pcrd;

    pcrd = &pull->coord[coord_ind];

    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        md2 = -1;
    }
    else
    {
        md2 = max_pull_distance2(pcrd, pbc);
    }

    if (pcrd->params.eGeom == epullgDIRRELATIVE)
    {
        /* We need to determine the pull vector */
        const pull_group_work_t *pgrp2, *pgrp3;
        dvec                     vec;
        int                      m;

        pgrp2 = &pull->group[pcrd->params.group[2]];
        pgrp3 = &pull->group[pcrd->params.group[3]];

        pbc_dx_d(pbc, pgrp3->x, pgrp2->x, vec);

        for (m = 0; m < DIM; m++)
        {
            vec[m] *= pcrd->params.dim[m];
        }
        pcrd->vec_len = dnorm(vec);
        for (m = 0; m < DIM; m++)
        {
            pcrd->vec[m] = vec[m]/pcrd->vec_len;
        }
        if (debug)
        {
            fprintf(debug, "pull coord %d vector: %6.3f %6.3f %6.3f normalized: %6.3f %6.3f %6.3f\n",
                    coord_ind,
                    vec[XX], vec[YY], vec[ZZ],
                    pcrd->vec[XX], pcrd->vec[YY], pcrd->vec[ZZ]);
        }
    }

    low_get_pull_coord_dr(pull, pcrd, pbc,
                          pull->group[pcrd->params.group[1]].x,
                          pcrd->params.eGeom == epullgCYL ? pull->dyna[coord_ind].x : pull->group[pcrd->params.group[0]].x,
                          md2,
                          pcrd->dr);
}

static void update_pull_coord_reference_value(pull_coord_work_t *pcrd, double t)
{
    /* With zero rate the reference value is set initially and doesn't change */
    if (pcrd->params.rate != 0)
    {
        pcrd->value_ref = pcrd->params.init + pcrd->params.rate*t;
    }
}

/* Calculates pull->coord[coord_ind].value.
 * This function also updates pull->coord[coord_ind].dr.
 */
static void get_pull_coord_distance(struct pull_t *pull,
                                    int            coord_ind,
                                    const t_pbc   *pbc)
{
    pull_coord_work_t *pcrd;
    int                m;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_dr(pull, coord_ind, pbc);

    switch (pcrd->params.eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            pcrd->value = dnorm(pcrd->dr);
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgDIRRELATIVE:
        case epullgCYL:
            /* Pull along vec */
            pcrd->value = 0;
            for (m = 0; m < DIM; m++)
            {
                pcrd->value += pcrd->vec[m]*pcrd->dr[m];
            }
            break;
    }
}

/* Returns the deviation from the reference value.
 * Updates pull->coord[coord_ind].dr, .value and .value_ref.
 */
static double get_pull_coord_deviation(struct pull_t *pull,
                                       int            coord_ind,
                                       const t_pbc   *pbc,
                                       double         t)
{
    static gmx_bool    bWarned = FALSE; /* TODO: this should be fixed for thread-safety,
                                           but is fairly benign */
    pull_coord_work_t *pcrd;
    double             dev = 0;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_distance(pull, coord_ind, pbc);

    update_pull_coord_reference_value(pcrd, t);

    switch (pcrd->params.eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            if (pcrd->value_ref < 0 && !bWarned)
            {
                fprintf(stderr, "\nPull reference distance for coordinate %d is negative (%f)\n", coord_ind+1, pcrd->value_ref);
                bWarned = TRUE;
            }
            if (pcrd->value == 0)
            {
                /* With no vector we can not determine the direction for the force,
                 * so we set the force to zero.
                 */
                dev = 0;
            }
            else
            {
                /* Determine the deviation */
                dev = pcrd->value - pcrd->value_ref;
            }
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgDIRRELATIVE:
        case epullgCYL:
            /* Pull along vec */
            dev = pcrd->value - pcrd->value_ref;
            break;
    }

    return dev;
}

void get_pull_coord_value(struct pull_t *pull,
                          int            coord_ind,
                          const t_pbc   *pbc,
                          double        *value)
{
    get_pull_coord_distance(pull, coord_ind, pbc);

    *value = pull->coord[coord_ind].value;
}

void clear_pull_forces(struct pull_t *pull)
{
    int i;

    /* Zeroing the forces is only required for constraint pulling.
     * It can happen that multiple constraint steps need to be applied
     * and therefore the constraint forces need to be accumulated.
     */
    for (i = 0; i < pull->ncoord; i++)
    {
        clear_dvec(pull->coord[i].f);
        pull->coord[i].f_scal = 0;
    }
}

/* Apply constraint using SHAKE */
static void do_constraint(struct pull_t *pull, t_pbc *pbc,
                          rvec *x, rvec *v,
                          gmx_bool bMaster, tensor vir,
                          double dt, double t)
{

    dvec         *r_ij;   /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
    dvec          unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */
    dvec         *rnew;   /* current 'new' positions of the groups */
    double       *dr_tot; /* the total update of the coords */
    dvec          vec;
    double        inpr;
    double        lambda, rm, invdt = 0;
    gmx_bool      bConverged_all, bConverged = FALSE;
    int           niter = 0, g, c, ii, j, m, max_iter = 100;
    double        a;
    dvec          tmp, tmp3;

    snew(r_ij,   pull->ncoord);
    snew(dr_tot, pull->ncoord);

    snew(rnew, pull->ngroup);

    /* copy the current unconstrained positions for use in iterations. We
       iterate until rinew[i] and rjnew[j] obey the constraints. Then
       rinew - pull.x_unc[i] is the correction dr to group i */
    for (g = 0; g < pull->ngroup; g++)
    {
        copy_dvec(pull->group[g].xp, rnew[g]);
    }

    /* Determine the constraint directions from the old positions */
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType != epullCONSTRAINT)
        {
            continue;
        }

        /* Note that get_pull_coord_distance also sets pcrd->dr and pcrd->value.
         * We don't modify dr and value anymore, so these values are also used
         * for printing.
         */
        get_pull_coord_distance(pull, c, pbc);

        if (debug)
        {
            fprintf(debug, "Pull coord %d dr %f %f %f\n",
                    c, pcrd->dr[XX], pcrd->dr[YY], pcrd->dr[ZZ]);
        }

        if (pcrd->params.eGeom == epullgDIR ||
            pcrd->params.eGeom == epullgDIRPBC)
        {
            /* Select the component along vec */
            a = 0;
            for (m = 0; m < DIM; m++)
            {
                a += pcrd->vec[m]*pcrd->dr[m];
            }
            for (m = 0; m < DIM; m++)
            {
                r_ij[c][m] = a*pcrd->vec[m];
            }
        }
        else
        {
            copy_dvec(pcrd->dr, r_ij[c]);
        }

        if (dnorm2(r_ij[c]) == 0)
        {
            gmx_fatal(FARGS, "Distance for pull coordinate %d is zero with constraint pulling, which is not allowed.", c + 1);
        }
    }

    bConverged_all = FALSE;
    while (!bConverged_all && niter < max_iter)
    {
        bConverged_all = TRUE;

        /* loop over all constraints */
        for (c = 0; c < pull->ncoord; c++)
        {
            pull_coord_work_t *pcrd;
            pull_group_work_t *pgrp0, *pgrp1;
            dvec               dr0, dr1;

            pcrd  = &pull->coord[c];

            if (pcrd->params.eType != epullCONSTRAINT)
            {
                continue;
            }

            update_pull_coord_reference_value(pcrd, t);

            pgrp0 = &pull->group[pcrd->params.group[0]];
            pgrp1 = &pull->group[pcrd->params.group[1]];

            /* Get the current difference vector */
            low_get_pull_coord_dr(pull, pcrd, pbc,
                                  rnew[pcrd->params.group[1]],
                                  rnew[pcrd->params.group[0]],
                                  -1, unc_ij);

            if (debug)
            {
                fprintf(debug, "Pull coord %d, iteration %d\n", c, niter);
            }

            rm = 1.0/(pgrp0->invtm + pgrp1->invtm);

            switch (pcrd->params.eGeom)
            {
                case epullgDIST:
                    if (pcrd->value_ref <= 0)
                    {
                        gmx_fatal(FARGS, "The pull constraint reference distance for group %d is <= 0 (%f)", c, pcrd->value_ref);
                    }

                    {
                        double q, c_a, c_b, c_c;

                        c_a = diprod(r_ij[c], r_ij[c]);
                        c_b = diprod(unc_ij, r_ij[c])*2;
                        c_c = diprod(unc_ij, unc_ij) - dsqr(pcrd->value_ref);

                        if (c_b < 0)
                        {
                            q      = -0.5*(c_b - sqrt(c_b*c_b - 4*c_a*c_c));
                            lambda = -q/c_a;
                        }
                        else
                        {
                            q      = -0.5*(c_b + sqrt(c_b*c_b - 4*c_a*c_c));
                            lambda = -c_c/q;
                        }

                        if (debug)
                        {
                            fprintf(debug,
                                    "Pull ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n",
                                    c_a, c_b, c_c, lambda);
                        }
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp1->invtm, r_ij[c], dr1);
                    dsvmul( lambda*rm*pgrp0->invtm, r_ij[c], dr0);
                    dr_tot[c] += -lambda*dnorm(r_ij[c]);
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    /* A 1-dimensional constraint along a vector */
                    a = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pcrd->vec[m];
                        a     += unc_ij[m]*vec[m];
                    }
                    /* Select only the component along the vector */
                    dsvmul(a, vec, unc_ij);
                    lambda = a - pcrd->value_ref;
                    if (debug)
                    {
                        fprintf(debug, "Pull inpr %e lambda: %e\n", a, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp1->invtm, vec, dr1);
                    dsvmul( lambda*rm*pgrp0->invtm, vec, dr0);
                    dr_tot[c] += -lambda;
                    break;
                default:
                    gmx_incons("Invalid enumeration value for eGeom");
                    /* Keep static analyzer happy */
                    clear_dvec(dr0);
                    clear_dvec(dr1);
            }

            /* DEBUG */
            if (debug)
            {
                int g0, g1;

                g0 = pcrd->params.group[0];
                g1 = pcrd->params.group[1];
                low_get_pull_coord_dr(pull, pcrd, pbc, rnew[g1], rnew[g0], -1, tmp);
                low_get_pull_coord_dr(pull, pcrd, pbc, dr1, dr0, -1, tmp3);
                fprintf(debug,
                        "Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        rnew[g0][0], rnew[g0][1], rnew[g0][2],
                        rnew[g1][0], rnew[g1][1], rnew[g1][2], dnorm(tmp));
                fprintf(debug,
                        "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f\n",
                        "", "", "", "", "", "", pcrd->value_ref);
                fprintf(debug,
                        "Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        dr0[0], dr0[1], dr0[2],
                        dr1[0], dr1[1], dr1[2],
                        dnorm(tmp3));
            } /* END DEBUG */

            /* Update the COMs with dr */
            dvec_inc(rnew[pcrd->params.group[1]], dr1);
            dvec_inc(rnew[pcrd->params.group[0]], dr0);
        }

        /* Check if all constraints are fullfilled now */
        for (c = 0; c < pull->ncoord; c++)
        {
            pull_coord_work_t *pcrd;

            pcrd = &pull->coord[c];

            if (pcrd->params.eType != epullCONSTRAINT)
            {
                continue;
            }

            low_get_pull_coord_dr(pull, pcrd, pbc,
                                  rnew[pcrd->params.group[1]],
                                  rnew[pcrd->params.group[0]],
                                  -1, unc_ij);

            switch (pcrd->params.eGeom)
            {
                case epullgDIST:
                    bConverged =
                        fabs(dnorm(unc_ij) - pcrd->value_ref) < pull->params.constr_tol;
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pcrd->vec[m];
                    }
                    inpr = diprod(unc_ij, vec);
                    dsvmul(inpr, vec, unc_ij);
                    bConverged =
                        fabs(diprod(unc_ij, vec) - pcrd->value_ref) < pull->params.constr_tol;
                    break;
            }

            if (!bConverged)
            {
                if (debug)
                {
                    fprintf(debug, "NOT CONVERGED YET: Group %d:"
                            "d_ref = %f, current d = %f\n",
                            g, pcrd->value_ref, dnorm(unc_ij));
                }

                bConverged_all = FALSE;
            }
        }

        niter++;
        /* if after all constraints are dealt with and bConverged is still TRUE
           we're finished, if not we do another iteration */
    }
    if (niter > max_iter)
    {
        gmx_fatal(FARGS, "Too many iterations for constraint run: %d", niter);
    }

    /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

    if (v)
    {
        invdt = 1/dt;
    }

    /* update atoms in the groups */
    for (g = 0; g < pull->ngroup; g++)
    {
        const pull_group_work_t *pgrp;
        dvec                     dr;

        pgrp = &pull->group[g];

        /* get the final constraint displacement dr for group g */
        dvec_sub(rnew[g], pgrp->xp, dr);

        if (dnorm2(dr) == 0)
        {
            /* No displacement, no update necessary */
            continue;
        }

        /* update the atom positions */
        copy_dvec(dr, tmp);
        for (j = 0; j < pgrp->nat_loc; j++)
        {
            ii = pgrp->ind_loc[j];
            if (pgrp->weight_loc)
            {
                dsvmul(pgrp->wscale*pgrp->weight_loc[j], dr, tmp);
            }
            for (m = 0; m < DIM; m++)
            {
                x[ii][m] += tmp[m];
            }
            if (v)
            {
                for (m = 0; m < DIM; m++)
                {
                    v[ii][m] += invdt*tmp[m];
                }
            }
        }
    }

    /* calculate the constraint forces, used for output and virial only */
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType != epullCONSTRAINT)
        {
            continue;
        }

        pcrd->f_scal = dr_tot[c]/((pull->group[pcrd->params.group[0]].invtm + pull->group[pcrd->params.group[1]].invtm)*dt*dt);

        if (vir != NULL && pcrd->params.eGeom != epullgDIRPBC && bMaster)
        {
            double f_invr;

            /* Add the pull contribution to the virial */
            /* We have already checked above that r_ij[c] != 0 */
            f_invr = pcrd->f_scal/dnorm(r_ij[c]);

            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5*f_invr*r_ij[c][j]*r_ij[c][m];
                }
            }
        }
    }

    /* finished! I hope. Give back some memory */
    sfree(r_ij);
    sfree(dr_tot);
    sfree(rnew);
}

/* Pulling with a harmonic umbrella potential or constant force */
static void do_pull_pot(struct pull_t *pull, t_pbc *pbc, double t, real lambda,
                        real *V, tensor vir, real *dVdl)
{
    int    c, j, m;
    double dev, ndr, invdr = 0;
    real   k, dkdl;

    /* loop over the pull coordinates */
    *V    = 0;
    *dVdl = 0;
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType == epullCONSTRAINT)
        {
            continue;
        }

        dev = get_pull_coord_deviation(pull, c, pbc, t);

        k    = (1.0 - lambda)*pcrd->params.k + lambda*pcrd->params.kB;
        dkdl = pcrd->params.kB - pcrd->params.k;

        if (pcrd->params.eGeom == epullgDIST)
        {
            ndr   = dnorm(pcrd->dr);
            if (ndr > 0)
            {
                invdr = 1/ndr;
            }
            else
            {
                /* With an harmonic umbrella, the force is 0 at r=0,
                 * so we can set invdr to any value.
                 * With a constant force, the force at r=0 is not defined,
                 * so we zero it (this is anyhow a very rare event).
                 */
                invdr = 0;
            }
        }
        else
        {
            ndr = 0;
            for (m = 0; m < DIM; m++)
            {
                ndr += pcrd->vec[m]*pcrd->dr[m];
            }
        }

        switch (pcrd->params.eType)
        {
            case epullUMBRELLA:
            case epullFLATBOTTOM:
                /* The only difference between an umbrella and a flat-bottom
                 * potential is that a flat-bottom is zero below zero.
                 */
                if (pcrd->params.eType == epullFLATBOTTOM && dev < 0)
                {
                    dev = 0;
                }

                pcrd->f_scal  =       -k*dev;
                *V           += 0.5*   k*dsqr(dev);
                *dVdl        += 0.5*dkdl*dsqr(dev);
                break;
            case epullCONST_F:
                pcrd->f_scal  =   -k;
                *V           +=    k*ndr;
                *dVdl        += dkdl*ndr;
                break;
            default:
                gmx_incons("Unsupported pull type in do_pull_pot");
        }

        if (pcrd->params.eGeom == epullgDIST)
        {
            for (m = 0; m < DIM; m++)
            {
                pcrd->f[m] = pcrd->f_scal*pcrd->dr[m]*invdr;
            }
        }
        else
        {
            for (m = 0; m < DIM; m++)
            {
                pcrd->f[m] = pcrd->f_scal*pcrd->vec[m];
            }
        }

        if (vir != NULL && pcrd->params.eGeom != epullgDIRPBC)
        {
            /* Add the pull contribution to the virial */
            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5*pcrd->f[j]*pcrd->dr[m];
                }
            }
        }
    }
}

real pull_potential(struct pull_t *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, rvec *f, tensor vir, real *dvdlambda)
{
    real V = 0, dVdl;

    assert(pull != NULL);

    if (pull->comm.bParticipate)
    {
        pull_calc_coms(cr, pull, md, pbc, t, x, NULL);

        do_pull_pot(pull, pbc, t, lambda,
                    &V, MASTER(cr) ? vir : NULL, &dVdl);

        /* Distribute forces over the pull groups */
        apply_forces(pull, md, f);

        if (MASTER(cr))
        {
            *dvdlambda += dVdl;
        }
    }

    return (MASTER(cr) ? V : 0.0);
}

void pull_constraint(struct pull_t *pull, t_mdatoms *md, t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir)
{
    assert(pull != NULL);

    if (pull->comm.bParticipate)
    {
        pull_calc_coms(cr, pull, md, pbc, t, x, xp);

        do_constraint(pull, pbc, xp, v, MASTER(cr), vir, dt, t);
    }
}

static void make_local_pull_group(gmx_ga2la_t ga2la,
                                  pull_group_work_t *pg, int start, int end)
{
    int i, ii;

    pg->nat_loc = 0;
    for (i = 0; i < pg->params.nat; i++)
    {
        ii = pg->params.ind[i];
        if (ga2la)
        {
            if (!ga2la_get_home(ga2la, ii, &ii))
            {
                ii = -1;
            }
        }
        if (ii >= start && ii < end)
        {
            /* This is a home atom, add it to the local pull group */
            if (pg->nat_loc >= pg->nalloc_loc)
            {
                pg->nalloc_loc = over_alloc_dd(pg->nat_loc+1);
                srenew(pg->ind_loc, pg->nalloc_loc);
                if (pg->epgrppbc == epgrppbcCOS || pg->params.weight != NULL)
                {
                    srenew(pg->weight_loc, pg->nalloc_loc);
                }
            }
            pg->ind_loc[pg->nat_loc] = ii;
            if (pg->params.weight != NULL)
            {
                pg->weight_loc[pg->nat_loc] = pg->params.weight[i];
            }
            pg->nat_loc++;
        }
    }
}

void dd_make_local_pull_groups(t_commrec *cr, struct pull_t *pull, t_mdatoms *md)
{
    gmx_domdec_t *dd;
    pull_comm_t  *comm;
    gmx_ga2la_t   ga2la;
    gmx_bool      bMustParticipate;
    int           g;

    dd = cr->dd;

    comm = &pull->comm;

    if (dd)
    {
        ga2la = dd->ga2la;
    }
    else
    {
        ga2la = NULL;
    }

    /* We always make the master node participate, such that it can do i/o
     * and to simplify MC type extensions people might have.
     */
    bMustParticipate = (comm->bParticipateAll || dd == NULL || DDMASTER(dd));

    for (g = 0; g < pull->ngroup; g++)
    {
        int a;

        make_local_pull_group(ga2la, &pull->group[g],
                              0, md->homenr);

        /* We should participate if we have pull or pbc atoms */
        if (!bMustParticipate &&
            (pull->group[g].nat_loc > 0 ||
             (pull->group[g].params.pbcatom >= 0 &&
              ga2la_get_home(dd->ga2la, pull->group[g].params.pbcatom, &a))))
        {
            bMustParticipate = TRUE;
        }
    }

    if (!comm->bParticipateAll)
    {
        /* Keep currently not required ranks in the communicator
         * if they needed to participate up to 20 decompositions ago.
         * This avoids frequent rebuilds due to atoms jumping back and forth.
         */
        const gmx_int64_t history_count = 20;
        gmx_bool          bWillParticipate;
        int               count[2];

        /* Increase the decomposition counter for the current call */
        comm->setup_count++;

        if (bMustParticipate)
        {
            comm->must_count = comm->setup_count;
        }

        bWillParticipate =
            bMustParticipate ||
            (comm->bParticipate &&
             comm->must_count >= comm->setup_count - history_count);

        if (debug && dd != NULL)
        {
            fprintf(debug, "Our DD rank (%3d) pull #atoms>0 or master: %d, will be part %d\n",
                    dd->rank, bMustParticipate, bWillParticipate);
        }

        if (bWillParticipate)
        {
            /* Count the number of ranks that we want to have participating */
            count[0] = 1;
            /* Count the number of ranks that need to be added */
            count[1] = comm->bParticipate ? 0 : 1;
        }
        else
        {
            count[0] = 0;
            count[1] = 0;
        }

        /* The cost of this global operation will be less that the cost
         * of the extra MPI_Comm_split calls that we can avoid.
         */
        gmx_sumi(2, count, cr);

        /* If we are missing ranks or if we have 20% more ranks than needed
         * we make a new sub-communicator.
         */
        if (count[1] > 0 || 6*count[0] < 5*comm->nparticipate)
        {
            if (debug)
            {
                fprintf(debug, "Creating new pull subcommunicator of size %d\n",
                        count[0]);
            }

#ifdef GMX_MPI
            if (comm->mpi_comm_com != MPI_COMM_NULL)
            {
                MPI_Comm_free(&comm->mpi_comm_com);
            }

            /* This might be an extremely expensive operation, so we try
             * to avoid this splitting as much as possible.
             */
            assert(dd != NULL);
            MPI_Comm_split(dd->mpi_comm_all, bWillParticipate ? 0 : 1, dd->rank,
                           &comm->mpi_comm_com);
#endif

            comm->bParticipate = bWillParticipate;
            comm->nparticipate = count[0];
        }
    }

    /* Since the PBC of atoms might have changed, we need to update the PBC */
    pull->bSetPBCatoms = TRUE;
}

static void init_pull_group_index(FILE *fplog, t_commrec *cr,
                                  int g, pull_group_work_t *pg,
                                  gmx_bool bConstraint, ivec pulldim_con,
                                  const gmx_mtop_t *mtop,
                                  const t_inputrec *ir, real lambda)
{
    int                   i, ii, d, nfrozen, ndim;
    real                  m, w, mbd;
    double                tmass, wmass, wwmass;
    const gmx_groups_t   *groups;
    gmx_mtop_atomlookup_t alook;
    t_atom               *atom;

    if (EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
    {
        /* There are no masses in the integrator.
         * But we still want to have the correct mass-weighted COMs.
         * So we store the real masses in the weights.
         * We do not set nweight, so these weights do not end up in the tpx file.
         */
        if (pg->params.nweight == 0)
        {
            snew(pg->params.weight, pg->params.nat);
        }
    }

    if (cr && PAR(cr))
    {
        pg->nat_loc    = 0;
        pg->nalloc_loc = 0;
        pg->ind_loc    = NULL;
        pg->weight_loc = NULL;
    }
    else
    {
        pg->nat_loc = pg->params.nat;
        pg->ind_loc = pg->params.ind;
        if (pg->epgrppbc == epgrppbcCOS)
        {
            snew(pg->weight_loc, pg->params.nat);
        }
        else
        {
            pg->weight_loc = pg->params.weight;
        }
    }

    groups = &mtop->groups;

    alook = gmx_mtop_atomlookup_init(mtop);

    nfrozen = 0;
    tmass   = 0;
    wmass   = 0;
    wwmass  = 0;
    for (i = 0; i < pg->params.nat; i++)
    {
        ii = pg->params.ind[i];
        gmx_mtop_atomnr_to_atom(alook, ii, &atom);
        if (bConstraint && ir->opts.nFreeze)
        {
            for (d = 0; d < DIM; d++)
            {
                if (pulldim_con[d] == 1 &&
                    ir->opts.nFreeze[ggrpnr(groups, egcFREEZE, ii)][d])
                {
                    nfrozen++;
                }
            }
        }
        if (ir->efep == efepNO)
        {
            m = atom->m;
        }
        else
        {
            m = (1 - lambda)*atom->m + lambda*atom->mB;
        }
        if (pg->params.nweight > 0)
        {
            w = pg->params.weight[i];
        }
        else
        {
            w = 1;
        }
        if (EI_ENERGY_MINIMIZATION(ir->eI))
        {
            /* Move the mass to the weight */
            w                   *= m;
            m                    = 1;
            pg->params.weight[i] = w;
        }
        else if (ir->eI == eiBD)
        {
            if (ir->bd_fric)
            {
                mbd = ir->bd_fric*ir->delta_t;
            }
            else
            {
                if (groups->grpnr[egcTC] == NULL)
                {
                    mbd = ir->delta_t/ir->opts.tau_t[0];
                }
                else
                {
                    mbd = ir->delta_t/ir->opts.tau_t[groups->grpnr[egcTC][ii]];
                }
            }
            w                   *= m/mbd;
            m                    = mbd;
            pg->params.weight[i] = w;
        }
        tmass  += m;
        wmass  += m*w;
        wwmass += m*w*w;
    }

    gmx_mtop_atomlookup_destroy(alook);

    if (wmass == 0)
    {
        /* We can have single atom groups with zero mass with potential pulling
         * without cosine weighting.
         */
        if (pg->params.nat == 1 && !bConstraint && pg->epgrppbc != epgrppbcCOS)
        {
            /* With one atom the mass doesn't matter */
            wwmass = 1;
        }
        else
        {
            gmx_fatal(FARGS, "The total%s mass of pull group %d is zero",
                      pg->params.weight ? " weighted" : "", g);
        }
    }
    if (fplog)
    {
        fprintf(fplog,
                "Pull group %d: %5d atoms, mass %9.3f",
                g, pg->params.nat, tmass);
        if (pg->params.weight ||
            EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
        {
            fprintf(fplog, ", weighted mass %9.3f", wmass*wmass/wwmass);
        }
        if (pg->epgrppbc == epgrppbcCOS)
        {
            fprintf(fplog, ", cosine weighting will be used");
        }
        fprintf(fplog, "\n");
    }

    if (nfrozen == 0)
    {
        /* A value != 0 signals not frozen, it is updated later */
        pg->invtm  = -1.0;
    }
    else
    {
        ndim = 0;
        for (d = 0; d < DIM; d++)
        {
            ndim += pulldim_con[d]*pg->params.nat;
        }
        if (fplog && nfrozen > 0 && nfrozen < ndim)
        {
            fprintf(fplog,
                    "\nWARNING: In pull group %d some, but not all of the degrees of freedom\n"
                    "         that are subject to pulling are frozen.\n"
                    "         For constraint pulling the whole group will be frozen.\n\n",
                    g);
        }
        pg->invtm  = 0.0;
        pg->wscale = 1.0;
    }
}

struct pull_t *
init_pull(FILE *fplog, const pull_params_t *pull_params, const t_inputrec *ir,
          int nfile, const t_filenm fnm[],
          gmx_mtop_t *mtop, t_commrec *cr, const output_env_t oenv, real lambda,
          gmx_bool bOutFile, unsigned long Flags)
{
    struct pull_t *pull;
    pull_comm_t   *comm;
    int            g, c, m;

    snew(pull, 1);

    /* Copy the pull parameters */
    pull->params       = *pull_params;
    /* Avoid pointer copies */
    pull->params.group = NULL;
    pull->params.coord = NULL;

    pull->ncoord       = pull_params->ncoord;
    pull->ngroup       = pull_params->ngroup;
    snew(pull->coord, pull->ncoord);
    snew(pull->group, pull->ngroup);

    pull->bPotential  = FALSE;
    pull->bConstraint = FALSE;
    pull->bCylinder   = FALSE;

    for (g = 0; g < pull->ngroup; g++)
    {
        pull_group_work_t *pgrp;
        int                i;

        pgrp = &pull->group[g];

        /* Copy the pull group parameters */
        pgrp->params = pull_params->group[g];

        /* Avoid pointer copies by allocating and copying arrays */
        snew(pgrp->params.ind, pgrp->params.nat);
        for (i = 0; i < pgrp->params.nat; i++)
        {
            pgrp->params.ind[i] = pull_params->group[g].ind[i];
        }
        if (pgrp->params.nweight > 0)
        {
            snew(pgrp->params.ind, pgrp->params.nweight);
            for (i = 0; i < pgrp->params.nweight; i++)
            {
                pgrp->params.weight[i] = pull_params->group[g].weight[i];
            }
        }
    }

    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;
        int                calc_com_start, calc_com_end, g;

        pcrd = &pull->coord[c];

        /* Copy all pull coordinate parameters */
        pcrd->params = pull_params->coord[c];

        switch (pcrd->params.eGeom)
        {
            case epullgDIST:
            case epullgDIRRELATIVE:
                /* Direction vector is determined at each step */
                break;
            case epullgDIR:
            case epullgDIRPBC:
            case epullgCYL:
                copy_rvec(pull_params->coord[c].vec, pcrd->vec);
                break;
            default:
                gmx_incons("Pull geometry not handled");
        }

        if (pcrd->params.eType == epullCONSTRAINT)
        {
            /* Check restrictions of the constraint pull code */
            if (pcrd->params.eGeom == epullgCYL ||
                pcrd->params.eGeom == epullgDIRRELATIVE)
            {
                gmx_fatal(FARGS, "Pulling of type %s can not be combined with geometry %s. Consider using pull type %s.",
                          epull_names[pcrd->params.eType],
                          epullg_names[pcrd->params.eGeom],
                          epull_names[epullUMBRELLA]);
            }

            pull->bConstraint = TRUE;
        }
        else
        {
            pull->bPotential = TRUE;
        }

        if (pcrd->params.eGeom == epullgCYL)
        {
            pull->bCylinder = TRUE;
        }
        /* We only need to calculate the plain COM of a group
         * when it is not only used as a cylinder group.
         */
        calc_com_start = (pcrd->params.eGeom == epullgCYL         ? 1 : 0);
        calc_com_end   = (pcrd->params.eGeom == epullgDIRRELATIVE ? 4 : 2);

        for (g = calc_com_start; g <= calc_com_end; g++)
        {
            pull->group[pcrd->params.group[g]].bCalcCOM = TRUE;
        }

        /* With non-zero rate the reference value is set at every step */
        if (pcrd->params.rate == 0)
        {
            /* Initialize the constant reference value */
            pcrd->value_ref = pcrd->params.init;
        }
    }

    pull->ePBC = ir->ePBC;
    switch (pull->ePBC)
    {
        case epbcNONE: pull->npbcdim = 0; break;
        case epbcXY:   pull->npbcdim = 2; break;
        default:       pull->npbcdim = 3; break;
    }

    if (fplog)
    {
        gmx_bool bAbs, bCos;

        bAbs = FALSE;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (pull->group[pull->coord[c].params.group[0]].params.nat == 0 ||
                pull->group[pull->coord[c].params.group[1]].params.nat == 0)
            {
                bAbs = TRUE;
            }
        }

        fprintf(fplog, "\n");
        if (pull->bPotential)
        {
            fprintf(fplog, "Will apply potential COM pulling\n");
        }
        if (pull->bConstraint)
        {
            fprintf(fplog, "Will apply constraint COM pulling\n");
        }
        fprintf(fplog, "with %d pull coordinate%s and %d group%s\n",
                pull->ncoord, pull->ncoord == 1 ? "" : "s",
                pull->ngroup, pull->ngroup == 1 ? "" : "s");
        if (bAbs)
        {
            fprintf(fplog, "with an absolute reference\n");
        }
        bCos = FALSE;
        for (g = 0; g < pull->ngroup; g++)
        {
            if (pull->group[g].params.nat > 1 &&
                pull->group[g].params.pbcatom < 0)
            {
                /* We are using cosine weighting */
                fprintf(fplog, "Cosine weighting is used for group %d\n", g);
                bCos = TRUE;
            }
        }
        if (bCos)
        {
            please_cite(fplog, "Engin2010");
        }
    }

    pull->bRefAt   = FALSE;
    pull->cosdim   = -1;
    for (g = 0; g < pull->ngroup; g++)
    {
        pull_group_work_t *pgrp;

        pgrp           = &pull->group[g];
        pgrp->epgrppbc = epgrppbcNONE;
        if (pgrp->params.nat > 0)
        {
            /* There is an issue when a group is used in multiple coordinates
             * and constraints are applied in different dimensions with atoms
             * frozen in some, but not all dimensions.
             * Since there is only one mass array per group, we can't have
             * frozen/non-frozen atoms for different coords at the same time.
             * But since this is a very exotic case, we don't check for this.
             * A warning is printed in init_pull_group_index.
             */

            gmx_bool bConstraint;
            ivec     pulldim, pulldim_con;

            /* Loop over all pull coordinates to see along which dimensions
             * this group is pulled and if it is involved in constraints.
             */
            bConstraint = FALSE;
            clear_ivec(pulldim);
            clear_ivec(pulldim_con);
            for (c = 0; c < pull->ncoord; c++)
            {
                if (pull->coord[c].params.group[0] == g ||
                    pull->coord[c].params.group[1] == g ||
                    (pull->coord[c].params.eGeom == epullgDIRRELATIVE &&
                     (pull->coord[c].params.group[2] == g ||
                      pull->coord[c].params.group[3] == g)))
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].params.dim[m] == 1)
                        {
                            pulldim[m] = 1;

                            if (pull->coord[c].params.eType == epullCONSTRAINT)
                            {
                                bConstraint    = TRUE;
                                pulldim_con[m] = 1;
                            }
                        }
                    }
                }
            }

            /* Determine if we need to take PBC into account for calculating
             * the COM's of the pull groups.
             */
            GMX_RELEASE_ASSERT(pull->npbcdim <= DIM, "npbcdim must be <= the number of dimensions");
            for (m = 0; m < pull->npbcdim; m++)
            {
                GMX_RELEASE_ASSERT(m <= DIM, "npbcdim must be <= the number of dimensions");
                if (pulldim[m] == 1 && pgrp->params.nat > 1)
                {
                    if (pgrp->params.pbcatom >= 0)
                    {
                        pgrp->epgrppbc = epgrppbcREFAT;
                        pull->bRefAt   = TRUE;
                    }
                    else
                    {
                        if (pgrp->params.weight != NULL)
                        {
                            gmx_fatal(FARGS, "Pull groups can not have relative weights and cosine weighting at same time");
                        }
                        pgrp->epgrppbc = epgrppbcCOS;
                        if (pull->cosdim >= 0 && pull->cosdim != m)
                        {
                            gmx_fatal(FARGS, "Can only use cosine weighting with pulling in one dimension (use mdp option pull-coord?-dim)");
                        }
                        pull->cosdim = m;
                    }
                }
            }
            /* Set the indices */
            init_pull_group_index(fplog, cr, g, pgrp,
                                  bConstraint, pulldim_con,
                                  mtop, ir, lambda);
        }
        else
        {
            /* Absolute reference, set the inverse mass to zero.
             * This is only relevant (and used) with constraint pulling.
             */
            pgrp->invtm  = 0;
            pgrp->wscale = 1;
        }
    }

    /* If we use cylinder coordinates, do some initialising for them */
    if (pull->bCylinder)
    {
        snew(pull->dyna, pull->ncoord);

        for (c = 0; c < pull->ncoord; c++)
        {
            const pull_coord_work_t *pcrd;

            pcrd = &pull->coord[c];

            if (pcrd->params.eGeom == epullgCYL)
            {
                if (pull->group[pcrd->params.group[0]].params.nat == 0)
                {
                    gmx_fatal(FARGS, "A cylinder pull group is not supported when using absolute reference!\n");
                }
            }
        }
    }

    comm = &pull->comm;

#ifdef GMX_MPI
    /* Use a sub-communicator when we have more than 32 ranks */
    comm->bParticipateAll = (cr == NULL || !DOMAINDECOMP(cr) ||
                             cr->dd->nnodes <= 32 ||
                             getenv("GMX_PULL_PARTICIPATE_ALL") != NULL);
    /* This sub-commicator is not used with comm->bParticipateAll,
     * so we can always initialize it to NULL.
     */
    comm->mpi_comm_com    = MPI_COMM_NULL;
    comm->nparticipate    = 0;
#else
    /* No MPI: 1 rank: all ranks pull */
    comm->bParticipateAll = TRUE;
#endif
    comm->bParticipate    = comm->bParticipateAll;
    comm->setup_count     = 0;
    comm->must_count      = 0;

    if (!comm->bParticipateAll && fplog != NULL)
    {
        fprintf(fplog, "Will use a sub-communicator for pull communication\n");
    }

    comm->rbuf     = NULL;
    comm->dbuf     = NULL;
    comm->dbuf_cyl = NULL;

    /* We still need to initialize the PBC reference coordinates */
    pull->bSetPBCatoms = TRUE;

    /* Only do I/O when we are doing dynamics and if we are the MASTER */
    pull->out_x = NULL;
    pull->out_f = NULL;
    if (bOutFile)
    {
        /* Check for px and pf filename collision, if we are writing
           both files */
        std::string px_filename, pf_filename;
        std::string px_appended, pf_appended;
        try
        {
            px_filename  = std::string(opt2fn("-px", nfile, fnm));
            pf_filename  = std::string(opt2fn("-pf", nfile, fnm));
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;


        if ((pull->params.nstxout != 0) &&
            (pull->params.nstfout != 0) &&
            (px_filename == pf_filename))
        {
            if (!opt2bSet("-px", nfile, fnm) && !opt2bSet("-pf", nfile, fnm))
            {
                /* We are writing both pull files but neither set directly. */
                try
                {
                    px_appended   = append_before_extension(px_filename, "_pullx");
                    pf_appended   = append_before_extension(pf_filename, "_pullf");
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                pull->out_x   = open_pull_out(px_appended.c_str(), pull, oenv,
                                              TRUE, Flags);
                pull->out_f = open_pull_out(pf_appended.c_str(), pull, oenv,
                                            FALSE, Flags);
                return pull;
            }
            else
            {
                /* If one of -px and -pf is set but the filenames are identical: */
                gmx_fatal(FARGS, "Identical pull_x and pull_f output filenames %s",
                          px_filename.c_str());
            }
        }
        if (pull->params.nstxout != 0)
        {
            pull->out_x = open_pull_out(opt2fn("-px", nfile, fnm), pull, oenv,
                                        TRUE, Flags);
        }
        if (pull->params.nstfout != 0)
        {
            pull->out_f = open_pull_out(opt2fn("-pf", nfile, fnm), pull, oenv,
                                        FALSE, Flags);
        }
    }

    return pull;
}

static void destroy_pull_group(pull_group_work_t *pgrp)
{
    if (pgrp->ind_loc != pgrp->params.ind)
    {
        sfree(pgrp->ind_loc);
    }
    if (pgrp->weight_loc != pgrp->params.weight)
    {
        sfree(pgrp->weight_loc);
    }
    sfree(pgrp->mdw);
    sfree(pgrp->dv);

    sfree(pgrp->params.ind);
    sfree(pgrp->params.weight);
}

static void destroy_pull(struct pull_t *pull)
{
    int i;

    for (i = 0; i < pull->ngroup; i++)
    {
        destroy_pull_group(&pull->group[i]);
    }
    for (i = 0; i < pull->ncoord; i++)
    {
        if (pull->coord[i].params.eGeom == epullgCYL)
        {
            destroy_pull_group(&(pull->dyna[i]));
        }
    }
    sfree(pull->group);
    sfree(pull->dyna);
    sfree(pull->coord);

#ifdef GMX_MPI
    if (pull->comm.mpi_comm_com != MPI_COMM_NULL)
    {
        MPI_Comm_free(&pull->comm.mpi_comm_com);
    }
#endif
    sfree(pull->comm.rbuf);
    sfree(pull->comm.dbuf);
    sfree(pull->comm.dbuf_cyl);

    sfree(pull);
}

void finish_pull(struct pull_t *pull)
{
    if (pull->out_x)
    {
        gmx_fio_fclose(pull->out_x);
    }
    if (pull->out_f)
    {
        gmx_fio_fclose(pull->out_f);
    }

    destroy_pull(pull);
}

gmx_bool pull_have_potential(const struct pull_t *pull)
{
    return pull->bPotential;
}

gmx_bool pull_have_constraint(const struct pull_t *pull)
{
    return pull->bConstraint;
}
