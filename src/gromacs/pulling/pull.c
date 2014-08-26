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

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/utility/futil.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/network.h"
#include "pull.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

static void pull_print_group_x(FILE *out, ivec dim, const t_pull_group *pgrp)
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

static void pull_print_coord_dr(FILE *out, const t_pull_coord *pcrd)
{
    int m;

    for (m = 0; m < DIM; m++)
    {
        if (pcrd->dim[m])
        {
            fprintf(out, "\t%g", pcrd->dr[m]);
        }
    }
}

static void pull_print_x(FILE *out, t_pull *pull, double t)
{
    int           c;
    t_pull_coord *pcrd;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd = &pull->coord[c];

        if (pull->bPrintRef)
        {
            if (pcrd->eGeom == epullgCYL)
            {
                pull_print_group_x(out, pcrd->dim, &pull->dyna[c]);
            }
            else
            {
                pull_print_group_x(out, pcrd->dim, &pull->group[pcrd->group[0]]);
            }
        }
        pull_print_coord_dr(out, pcrd);
    }
    fprintf(out, "\n");
}

static void pull_print_f(FILE *out, t_pull *pull, double t)
{
    int c, d;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        fprintf(out, "\t%g", pull->coord[c].f_scal);
    }
    fprintf(out, "\n");
}

void pull_print_output(t_pull *pull, gmx_int64_t step, double time)
{
    if ((pull->nstxout != 0) && (step % pull->nstxout == 0))
    {
        pull_print_x(pull->out_x, pull, time);
    }

    if ((pull->nstfout != 0) && (step % pull->nstfout == 0))
    {
        pull_print_f(pull->out_f, pull, time);
    }
}

static FILE *open_pull_out(const char *fn, t_pull *pull, const output_env_t oenv,
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

        snew(setname, 2*pull->ncoord*DIM);
        nsets = 0;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (bCoord)
            {
                if (pull->bPrintRef)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].dim[m])
                        {
                            sprintf(buf, "%d %s%c", c+1, "c", 'X'+m);
                            setname[nsets] = gmx_strdup(buf);
                            nsets++;
                        }
                    }
                }
                for (m = 0; m < DIM; m++)
                {
                    if (pull->coord[c].dim[m])
                    {
                        sprintf(buf, "%d %s%c", c+1, "d", 'X'+m);
                        setname[nsets] = gmx_strdup(buf);
                        nsets++;
                    }
                }
            }
            else
            {
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
static void apply_forces_grp(const t_pull_group *pgrp, const t_mdatoms *md,
                             const dvec f_pull, int sign, rvec *f)
{
    int    i, ii, m;
    double wmass, inv_wm;

    inv_wm = pgrp->mwscale;

    dvec sd;
    clear_dvec(sd);
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

            sd[m]    += sign * wmass * f_pull[m] * inv_wm;
        }
    }
    //printf("sd grp %6.3f %6.3f\n", sd[YY], sd[ZZ]);
}

/* Apply forces in a mass weighted fashion to a cylinder group */
static void apply_forces_cyl_grp(const t_pull_group *pgrp, double dv_corr,
                                 const t_mdatoms *md,
                                 const dvec f_pull, double f_scal,
                                 int sign, rvec *f)
{
    int    i, ii, m;
    double mass, weight, inv_wm, dv_com;

    inv_wm = pgrp->mwscale;

    dvec sd;
    clear_dvec(sd);
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

            sd[m] += sign*inv_wm*(mass*weight*f_pull[m] +
                                  pgrp->mdw[i][m]*dv_com*f_scal);
        }
    }
    //printf("sd cyl %6.3f %6.3f\n", sd[YY], sd[ZZ]);
}

/* Apply forces in a mass weighted fashion */
static void apply_forces(t_pull * pull, t_mdatoms * md, rvec *f)
{
    int                 c;
    const t_pull_coord *pcrd;

    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd = &pull->coord[c];

        if (pcrd->eGeom == epullgCYL)
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
            apply_forces_grp(&pull->group[pcrd->group[1]], md, f_tot, 1, f);
        }
        else
        {
            if (pull->group[pcrd->group[0]].nat > 0)
            {
                apply_forces_grp(&pull->group[pcrd->group[0]], md, pcrd->f, -1, f);
            }
            apply_forces_grp(&pull->group[pcrd->group[1]], md, pcrd->f, 1, f);
        }
    }
}

static double max_pull_distance2(const t_pull_coord *pcrd, const t_pbc *pbc)
{
    double max_d2;
    int    m;

    max_d2 = GMX_DOUBLE_MAX;

    for (m = 0; m < pbc->ndim_ePBC; m++)
    {
        if (pcrd->dim[m] != 0)
        {
            max_d2 = min(max_d2, norm2(pbc->box[m]));
        }
    }

    return 0.25*max_d2;
}

static void low_get_pull_coord_dr(const t_pull *pull,
                                  const t_pull_coord *pcrd,
                                  const t_pbc *pbc, double t,
                                  dvec xg, dvec xref, double max_dist2,
                                  dvec dr)
{
    const t_pull_group *pgrp0, *pgrp1;
    int                 m;
    dvec                xrefr, dref = {0, 0, 0};
    double              dr2;

    pgrp0 = &pull->group[pcrd->group[0]];
    pgrp1 = &pull->group[pcrd->group[1]];

    /* Only the first group can be an absolute reference, in that case nat=0 */
    if (pgrp0->nat == 0)
    {
        for (m = 0; m < DIM; m++)
        {
            xref[m] = pcrd->origin[m];
        }
    }

    copy_dvec(xref, xrefr);

    if (pcrd->eGeom == epullgDIRPBC)
    {
        for (m = 0; m < DIM; m++)
        {
            dref[m] = (pcrd->init + pcrd->rate*t)*pcrd->vec[m];
        }
        /* Add the reference position, so we use the correct periodic image */
        dvec_inc(xrefr, dref);
    }

    pbc_dx_d(pbc, xg, xrefr, dr);
    dr2 = 0;
    for (m = 0; m < DIM; m++)
    {
        dr[m] *= pcrd->dim[m];
        dr2   += dr[m]*dr[m];
    }
    if (max_dist2 >= 0 && dr2 > 0.98*0.98*max_dist2)
    {
        gmx_fatal(FARGS, "Distance between pull groups %d and %d (%f nm) is larger than 0.49 times the box size (%f).\nYou might want to consider using \"pull-geometry = direction-periodic\" instead.\n",
                  pcrd->group[0], pcrd->group[1], sqrt(dr2), sqrt(max_dist2));
    }

    if (pcrd->eGeom == epullgDIRPBC)
    {
        dvec_inc(dr, dref);
    }
}

static void get_pull_coord_dr(const t_pull *pull,
                              int coord_ind,
                              const t_pbc *pbc, double t,
                              dvec dr)
{
    double              md2;
    const t_pull_coord *pcrd;

    pcrd = &pull->coord[coord_ind];

    if (pcrd->eGeom == epullgDIRPBC)
    {
        md2 = -1;
    }
    else
    {
        md2 = max_pull_distance2(pcrd, pbc);
    }

    low_get_pull_coord_dr(pull, pcrd, pbc, t,
                          pull->group[pcrd->group[1]].x,
                          pcrd->eGeom == epullgCYL ? pull->dyna[coord_ind].x : pull->group[pcrd->group[0]].x,
                          md2,
                          dr);
}

void get_pull_coord_distance(const t_pull *pull,
                             int coord_ind,
                             const t_pbc *pbc, double t,
                             dvec dr, double *dev)
{
    static gmx_bool     bWarned = FALSE; /* TODO: this should be fixed for thread-safety,
                                            but is fairly benign */
    const t_pull_coord *pcrd;
    int                 m;
    double              ref, drs, inpr;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_dr(pull, coord_ind, pbc, t, dr);

    ref = pcrd->init + pcrd->rate*t;

    switch (pcrd->eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            if (ref < 0 && !bWarned)
            {
                fprintf(stderr, "\nPull reference distance for coordinate %d is negative (%f)\n", coord_ind+1, ref);
                bWarned = TRUE;
            }
            drs = dnorm(dr);
            if (drs == 0)
            {
                /* With no vector we can not determine the direction for the force,
                 * so we set the force to zero.
                 */
                *dev = 0;
            }
            else
            {
                /* Determine the deviation */
                *dev = drs - ref;
            }
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgCYL:
            /* Pull along vec */
            inpr = 0;
            for (m = 0; m < DIM; m++)
            {
                inpr += pcrd->vec[m]*dr[m];
            }
            *dev = inpr - ref;
            break;
    }
}

void clear_pull_forces(t_pull *pull)
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
static void do_constraint(t_pull *pull, t_pbc *pbc,
                          rvec *x, rvec *v,
                          gmx_bool bMaster, tensor vir,
                          double dt, double t)
{

    dvec         *r_ij;   /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
    dvec          unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */
    dvec         *rnew;   /* current 'new' positions of the groups */
    double       *dr_tot; /* the total update of the coords */
    double        ref;
    dvec          vec;
    double        d0, inpr;
    double        lambda, rm, mass, invdt = 0;
    gmx_bool      bConverged_all, bConverged = FALSE;
    int           niter = 0, g, c, ii, j, m, max_iter = 100;
    double        a;
    dvec          f;       /* the pull force */
    dvec          tmp, tmp3;
    t_pull_group *pdyna, *pgrp0, *pgrp1;
    t_pull_coord *pcrd;
    int           cylinder_coord = -1;

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
        pcrd = &pull->coord[c];

        if (pcrd->eType != epullCONSTRAINT)
        {
            continue;
        }

        if (pcrd->eGeom == epullgCYL)
        {
            /* group[0] is a cylinder group (and is used only once) */
            copy_dvec(pull->dyna[c].xp, rnew[pcrd->group[0]]);
            cylinder_coord = c;
        }

        get_pull_coord_dr(pull, c, pbc, t, r_ij[c]);
        /* Store the difference vector at time t for printing */
        copy_dvec(r_ij[c], pcrd->dr);
        if (debug)
        {
            fprintf(debug, "Pull coord %d dr %f %f %f\n",
                    c, r_ij[c][XX], r_ij[c][YY], r_ij[c][ZZ]);
        }

        if (pcrd->eGeom == epullgDIR || pcrd->eGeom == epullgDIRPBC)
        {
            /* Select the component along vec */
            a = 0;
            for (m = 0; m < DIM; m++)
            {
                a += pcrd->vec[m]*r_ij[c][m];
            }
            for (m = 0; m < DIM; m++)
            {
                r_ij[c][m] = a*pcrd->vec[m];
            }
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
            dvec dr0, dr1;

            pcrd  = &pull->coord[c];

            if (pcrd->eType != epullCONSTRAINT)
            {
                continue;
            }

            pgrp0 = &pull->group[pcrd->group[0]];
            pgrp1 = &pull->group[pcrd->group[1]];

            /* Get the current difference vector */
            low_get_pull_coord_dr(pull, pcrd, pbc, t,
                                  rnew[pcrd->group[1]],
                                  rnew[pcrd->group[0]],
                                  -1, unc_ij);

            ref = pcrd->init + pcrd->rate*t;

            if (debug)
            {
                fprintf(debug, "Pull coord %d, iteration %d\n", c, niter);
            }

            rm = 1.0/(pgrp0->invtm + pgrp1->invtm);

            switch (pcrd->eGeom)
            {
                case epullgDIST:
                    if (ref <= 0)
                    {
                        gmx_fatal(FARGS, "The pull constraint reference distance for group %d is <= 0 (%f)", c, ref);
                    }

                    {
                        double q, c_a, c_b, c_c;

                        c_a = diprod(r_ij[c], r_ij[c]);
                        c_b = diprod(unc_ij, r_ij[c])*2;
                        c_c = diprod(unc_ij, unc_ij) - dsqr(ref);

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
                    lambda = a - ref;
                    if (debug)
                    {
                        fprintf(debug, "Pull inpr %e lambda: %e\n", a, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp1->invtm, vec, dr1);
                    dsvmul( lambda*rm*pgrp0->invtm, vec, dr0);
                    dr_tot[c] += -lambda;
                    break;
            }

            /* DEBUG */
            if (debug)
            {
                int g0, g1;

                g0 = pcrd->group[0];
                g1 = pcrd->group[1];
                low_get_pull_coord_dr(pull, pcrd, pbc, t, rnew[g1], rnew[g0], -1, tmp);
                low_get_pull_coord_dr(pull, pcrd, pbc, t, dr1, dr0, -1, tmp3);
                fprintf(debug,
                        "Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        rnew[g0][0], rnew[g0][1], rnew[g0][2],
                        rnew[g1][0], rnew[g1][1], rnew[g1][2], dnorm(tmp));
                fprintf(debug,
                        "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f\n",
                        "", "", "", "", "", "", ref);
                fprintf(debug,
                        "Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        dr0[0], dr0[1], dr0[2],
                        dr1[0], dr1[1], dr1[2],
                        dnorm(tmp3));
            } /* END DEBUG */

            /* Update the COMs with dr */
            dvec_inc(rnew[pcrd->group[1]], dr1);
            dvec_inc(rnew[pcrd->group[0]], dr0);
        }

        /* Check if all constraints are fullfilled now */
        for (c = 0; c < pull->ncoord; c++)
        {
            pcrd = &pull->coord[c];

            if (pcrd->eType != epullCONSTRAINT)
            {
                continue;
            }

            ref  = pcrd->init + pcrd->rate*t;

            low_get_pull_coord_dr(pull, pcrd, pbc, t,
                                  rnew[pcrd->group[1]],
                                  rnew[pcrd->group[0]],
                                  -1, unc_ij);

            switch (pcrd->eGeom)
            {
                case epullgDIST:
                    bConverged = fabs(dnorm(unc_ij) - ref) < pull->constr_tol;
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
                        fabs(diprod(unc_ij, vec) - ref) < pull->constr_tol;
                    break;
            }

            if (!bConverged)
            {
                if (debug)
                {
                    fprintf(debug, "NOT CONVERGED YET: Group %d:"
                            "d_ref = %f, current d = %f\n",
                            g, ref, dnorm(unc_ij));
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
        const t_pull_group *pgrp;
        dvec                dr;

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
        pcrd = &pull->coord[c];

        if (pcrd->eType != epullCONSTRAINT)
        {
            continue;
        }

        pcrd->f_scal = dr_tot[c]/((pull->group[pcrd->group[0]].invtm + pull->group[pcrd->group[1]].invtm)*dt*dt);

        if (vir != NULL && pcrd->eGeom != epullgDIRPBC && bMaster)
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
static void do_pull_pot(t_pull *pull, t_pbc *pbc, double t, real lambda,
                        real *V, tensor vir, real *dVdl)
{
    int           c, j, m;
    double        dev, ndr, invdr = 0;
    real          k, dkdl;
    t_pull_coord *pcrd;

    /* loop over the pull coordinates */
    *V    = 0;
    *dVdl = 0;
    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd = &pull->coord[c];

        if (pcrd->eType == epullCONSTRAINT)
        {
            continue;
        }

        get_pull_coord_distance(pull, c, pbc, t, pcrd->dr, &dev);

        k    = (1.0 - lambda)*pcrd->k + lambda*pcrd->kB;
        dkdl = pcrd->kB - pcrd->k;

        if (pcrd->eGeom == epullgDIST)
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

        switch (pcrd->eType)
        {
            case epullUMBRELLA:
            case epullFLATBOTTOM:
                /* The only difference between an umbrella and a flat-bottom
                 * potential is that a flat-bottom is zero below zero.
                 */
                if (pcrd->eType == epullFLATBOTTOM && dev < 0)
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

        if (pcrd->eGeom == epullgDIST)
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

        if (vir != NULL && pcrd->eGeom != epullgDIRPBC)
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

real pull_potential(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, rvec *f, tensor vir, real *dvdlambda)
{
    real V, dVdl;

    pull_calc_coms(cr, pull, md, pbc, t, x, NULL);

    do_pull_pot(pull, pbc, t, lambda,
                &V, MASTER(cr) ? vir : NULL, &dVdl);

    /* Distribute forces over pulled groups */
    apply_forces(pull, md, f);

    if (MASTER(cr))
    {
        *dvdlambda += dVdl;
    }

    return (MASTER(cr) ? V : 0.0);
}

void pull_constraint(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir)
{
    pull_calc_coms(cr, pull, md, pbc, t, x, xp);

    do_constraint(pull, pbc, xp, v, MASTER(cr), vir, dt, t);
}

static void make_local_pull_group(gmx_ga2la_t ga2la,
                                  t_pull_group *pg, int start, int end)
{
    int i, ii;

    pg->nat_loc = 0;
    for (i = 0; i < pg->nat; i++)
    {
        ii = pg->ind[i];
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
                if (pg->epgrppbc == epgrppbcCOS || pg->weight)
                {
                    srenew(pg->weight_loc, pg->nalloc_loc);
                }
            }
            pg->ind_loc[pg->nat_loc] = ii;
            if (pg->weight)
            {
                pg->weight_loc[pg->nat_loc] = pg->weight[i];
            }
            pg->nat_loc++;
        }
    }
}

void dd_make_local_pull_groups(gmx_domdec_t *dd, t_pull *pull, t_mdatoms *md)
{
    gmx_ga2la_t ga2la;
    int         g;

    if (dd)
    {
        ga2la = dd->ga2la;
    }
    else
    {
        ga2la = NULL;
    }

    for (g = 0; g < pull->ngroup; g++)
    {
        make_local_pull_group(ga2la, &pull->group[g],
                              0, md->homenr);
    }

    /* Since the PBC of atoms might have changed, we need to update the PBC */
    pull->bSetPBCatoms = TRUE;
}

static void init_pull_group_index(FILE *fplog, t_commrec *cr,
                                  int g, t_pull_group *pg,
                                  gmx_bool bConstraint, ivec pulldim_con,
                                  gmx_mtop_t *mtop, t_inputrec *ir, real lambda)
{
    int                   i, ii, d, nfrozen, ndim;
    real                  m, w, mbd;
    double                tmass, wmass, wwmass;
    gmx_groups_t         *groups;
    gmx_mtop_atomlookup_t alook;
    t_atom               *atom;

    if (EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
    {
        /* There are no masses in the integrator.
         * But we still want to have the correct mass-weighted COMs.
         * So we store the real masses in the weights.
         * We do not set nweight, so these weights do not end up in the tpx file.
         */
        if (pg->nweight == 0)
        {
            snew(pg->weight, pg->nat);
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
        pg->nat_loc = pg->nat;
        pg->ind_loc = pg->ind;
        if (pg->epgrppbc == epgrppbcCOS)
        {
            snew(pg->weight_loc, pg->nat);
        }
        else
        {
            pg->weight_loc = pg->weight;
        }
    }

    groups = &mtop->groups;

    alook = gmx_mtop_atomlookup_init(mtop);

    nfrozen = 0;
    tmass   = 0;
    wmass   = 0;
    wwmass  = 0;
    for (i = 0; i < pg->nat; i++)
    {
        ii = pg->ind[i];
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
        if (pg->nweight > 0)
        {
            w = pg->weight[i];
        }
        else
        {
            w = 1;
        }
        if (EI_ENERGY_MINIMIZATION(ir->eI))
        {
            /* Move the mass to the weight */
            w            *= m;
            m             = 1;
            pg->weight[i] = w;
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
            w            *= m/mbd;
            m             = mbd;
            pg->weight[i] = w;
        }
        tmass  += m;
        wmass  += m*w;
        wwmass += m*w*w;
    }

    gmx_mtop_atomlookup_destroy(alook);

    if (wmass == 0)
    {
        gmx_fatal(FARGS, "The total%s mass of pull group %d is zero",
                  pg->weight ? " weighted" : "", g);
    }
    if (fplog)
    {
        fprintf(fplog,
                "Pull group %d: %5d atoms, mass %9.3f", g, pg->nat, tmass);
        if (pg->weight || EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
        {
            fprintf(fplog, ", weighted mass %9.3f", wmass*wmass/wwmass);
        }
        if (pg->epgrppbc == epgrppbcCOS)
        {
            fprintf(fplog, ", cosine weighting will be used");
        }
        fprintf(fplog, "\n");
    }

    if (nfrozen > 0)
    {
        ndim = 0;
        for (d = 0; d < DIM; d++)
        {
            ndim += pulldim_con[d]*pg->nat;
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

void init_pull(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
               gmx_mtop_t *mtop, t_commrec *cr, const output_env_t oenv, real lambda,
               gmx_bool bOutFile, unsigned long Flags)
{
    t_pull       *pull;
    t_pull_group *pgrp;
    int           c, g, start = 0, end = 0, m;

    pull = ir->pull;

    pull->bPotential  = FALSE;
    pull->bConstraint = FALSE;
    pull->bCylinder   = FALSE;
    for (c = 0; c < pull->ncoord; c++)
    {
        t_pull_coord *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->eType == epullCONSTRAINT)
        {
            pull->bConstraint = TRUE;
        }
        else
        {
            pull->bPotential = TRUE;
        }

        if (pcrd->eGeom == epullgCYL)
        {
            pull->bCylinder = TRUE;

            if (pcrd->eType == epullCONSTRAINT)
            {
                gmx_fatal(FARGS, "Pulling of type %s can not be combined with geometry %s. Consider using pull type %s.",
                          epull_names[pcrd->eType],
                          epullg_names[pcrd->eGeom],
                          epull_names[epullUMBRELLA]);
            }
        }
        else
        {
            /* We only need to calculate the plain COM of a group
             * when it is not only used as a cylinder group.
             */
            if (pull->group[pcrd->group[0]].nat > 0)
            {
                pull->group[pcrd->group[0]].bCalcCOM = TRUE;
            }
        }
        if (pull->group[pcrd->group[1]].nat > 0)
        {
            pull->group[pcrd->group[1]].bCalcCOM = TRUE;
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
            if (pull->group[pull->coord[c].group[0]].nat == 0 ||
                pull->group[pull->coord[c].group[1]].nat == 0)
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
            if (pull->group[g].nat > 1 &&
                pull->group[g].pbcatom < 0)
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

    pull->rbuf     = NULL;
    pull->dbuf     = NULL;
    pull->dbuf_cyl = NULL;
    pull->bRefAt   = FALSE;
    pull->cosdim   = -1;
    for (g = 0; g < pull->ngroup; g++)
    {
        pgrp           = &pull->group[g];
        pgrp->epgrppbc = epgrppbcNONE;
        if (pgrp->nat > 0)
        {
            gmx_bool bConstraint;
            ivec     pulldim, pulldim_con;

            /* There is an issue when a group is used in multiple coordinates
             * and constraints are applied in different dimensions with atoms
             * frozen in some, but not all dimensions.
             * Since there is only one mass array per group, we can't have
             * frozen/non-frozen atoms for different coords at the same time.
             * But since this is a very exotic case, we don't check for this.
             * A warning is printed in init_pull_group_index.
             */
            bConstraint = FALSE;
            clear_ivec(pulldim);
            clear_ivec(pulldim_con);
            for (c = 0; c < pull->ncoord; c++)
            {
                if (pull->coord[c].group[0] == g ||
                    pull->coord[c].group[1] == g)
                {
                    bConstraint = TRUE;
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->coord[c].dim[m] == 1)
                        {
                            pulldim[m] = 1;

                            if (pull->coord[c].eType == epullCONSTRAINT)
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
            for (m = 0; m < pull->npbcdim; m++)
            {
                if (pulldim[m] == 1 && pgrp->nat > 1)
                {
                    if (pgrp->pbcatom >= 0)
                    {
                        pgrp->epgrppbc = epgrppbcREFAT;
                        pull->bRefAt   = TRUE;
                    }
                    else
                    {
                        if (pgrp->weight)
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
            const t_pull_coord *pcrd;

            pcrd = &pull->coord[c];

            if (pcrd->eGeom == epullgCYL)
            {
                if (pull->group[pcrd->group[0]].nat == 0)
                {
                    gmx_fatal(FARGS, "A cylinder pull group is not supported when using absolute reference!\n");
                }
            }
        }
    }

    /* We still need to initialize the PBC reference coordinates */
    pull->bSetPBCatoms = TRUE;

    /* Only do I/O when we are doing dynamics and if we are the MASTER */
    pull->out_x = NULL;
    pull->out_f = NULL;
    if (bOutFile)
    {
        if (pull->nstxout > 0)
        {
            pull->out_x = open_pull_out(opt2fn("-px", nfile, fnm), pull, oenv, TRUE, Flags);
        }
        if (pull->nstfout > 0)
        {
            pull->out_f = open_pull_out(opt2fn("-pf", nfile, fnm), pull, oenv,
                                        FALSE, Flags);
        }
    }
}

void finish_pull(t_pull *pull)
{
    if (pull->out_x)
    {
        gmx_fio_fclose(pull->out_x);
    }
    if (pull->out_f)
    {
        gmx_fio_fclose(pull->out_f);
    }
}
