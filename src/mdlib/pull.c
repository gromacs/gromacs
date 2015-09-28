/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "futil.h"
#include "index.h"
#include "statutil.h"
#include "gmxfio.h"
#include "vec.h"
#include "typedefs.h"
#include "network.h"
#include "filenm.h"
#include <string.h>
#include "smalloc.h"
#include "pull.h"
#include "xvgr.h"
#include "names.h"
#include "partdec.h"
#include "pbc.h"
#include "mtop_util.h"
#include "mdrun.h"
#include "gmx_ga2la.h"
#include "copyrite.h"

static void pull_print_x_grp(FILE *out, gmx_bool bRef, ivec dim, t_pullgrp *pgrp)
{
    int m;

    for (m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", bRef ? pgrp->x[m] : pgrp->dr[m]);
        }
    }
}

static void pull_print_x(FILE *out, t_pull *pull, double t)
{
    int g;

    fprintf(out, "%.4f", t);

    if (PULL_CYL(pull))
    {
        for (g = 1; g < 1+pull->ngrp; g++)
        {
            pull_print_x_grp(out, TRUE, pull->dim, &pull->dyna[g]);
            pull_print_x_grp(out, FALSE, pull->dim, &pull->grp[g]);
        }
    }
    else
    {
        for (g = 0; g < 1+pull->ngrp; g++)
        {
            if (pull->grp[g].nat > 0)
            {
                pull_print_x_grp(out, g == 0, pull->dim, &pull->grp[g]);
            }
        }
    }
    fprintf(out, "\n");
}

static void pull_print_f(FILE *out, t_pull *pull, double t)
{
    int g, d;

    fprintf(out, "%.4f", t);

    for (g = 1; g < 1+pull->ngrp; g++)
    {
        if (pull->eGeom == epullgPOS)
        {
            for (d = 0; d < DIM; d++)
            {
                if (pull->dim[d])
                {
                    fprintf(out, "\t%g", pull->grp[g].f[d]);
                }
            }
        }
        else
        {
            fprintf(out, "\t%g", pull->grp[g].f_scal);
        }
    }
    fprintf(out, "\n");
}

void pull_print_output(t_pull *pull, gmx_large_int_t step, double time)
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
    int    nsets, g, m;
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

        snew(setname, (1+pull->ngrp)*DIM);
        nsets = 0;
        for (g = 0; g < 1+pull->ngrp; g++)
        {
            if (pull->grp[g].nat > 0 &&
                (g > 0 || (bCoord && !PULL_CYL(pull))))
            {
                if (bCoord || pull->eGeom == epullgPOS)
                {
                    if (PULL_CYL(pull))
                    {
                        for (m = 0; m < DIM; m++)
                        {
                            if (pull->dim[m])
                            {
                                sprintf(buf, "%d %s%c", g, "c", 'X'+m);
                                setname[nsets] = strdup(buf);
                                nsets++;
                            }
                        }
                    }
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->dim[m])
                        {
                            sprintf(buf, "%d %s%c",
                                    g, (bCoord && g > 0) ? "d" : "", 'X'+m);
                            setname[nsets] = strdup(buf);
                            nsets++;
                        }
                    }
                }
                else
                {
                    sprintf(buf, "%d", g);
                    setname[nsets] = strdup(buf);
                    nsets++;
                }
            }
        }
        if (bCoord || nsets > 1)
        {
            xvgr_legend(fp, nsets, (const char**)setname, oenv);
        }
        for (g = 0; g < nsets; g++)
        {
            sfree(setname[g]);
        }
        sfree(setname);
    }

    return fp;
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_grp(t_pullgrp *pgrp, t_mdatoms * md,
                             gmx_ga2la_t ga2la,
                             dvec f_pull, int sign, rvec *f)
{
    int    i, ii, m, start, end;
    double wmass, inv_wm;

    start = md->start;
    end   = md->homenr + start;

    inv_wm = pgrp->wscale*pgrp->invtm;

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

/* Apply forces in a mass weighted fashion */
static void apply_forces(t_pull * pull, t_mdatoms * md, gmx_ga2la_t ga2la,
                         rvec *f)
{
    int        i;
    t_pullgrp *pgrp;

    for (i = 1; i < pull->ngrp+1; i++)
    {
        pgrp = &(pull->grp[i]);
        apply_forces_grp(pgrp, md, ga2la, pgrp->f, 1, f);
        if (pull->grp[0].nat)
        {
            if (PULL_CYL(pull))
            {
                apply_forces_grp(&(pull->dyna[i]), md, ga2la, pgrp->f, -1, f);
            }
            else
            {
                apply_forces_grp(&(pull->grp[0]), md, ga2la, pgrp->f, -1, f);
            }
        }
    }
}

static double max_pull_distance2(const t_pull *pull, const t_pbc *pbc)
{
    double max_d2;
    int    m;

    max_d2 = GMX_DOUBLE_MAX;

    if (pull->eGeom != epullgDIRPBC)
    {
        for (m = 0; m < pbc->ndim_ePBC; m++)
        {
            if (pull->dim[m] != 0)
            {
                max_d2 = min(max_d2, norm2(pbc->box[m]));
            }
        }
    }

    return 0.25*max_d2;
}

static void get_pullgrps_dr(const t_pull *pull, const t_pbc *pbc, int g, double t,
                            dvec xg, dvec xref, double max_dist2,
                            dvec dr)
{
    t_pullgrp *pref, *pgrp;
    int        m;
    dvec       xrefr, dref = {0, 0, 0};
    double     dr2;

    pgrp = &pull->grp[g];

    copy_dvec(xref, xrefr);

    if (pull->eGeom == epullgDIRPBC)
    {
        for (m = 0; m < DIM; m++)
        {
            dref[m] = (pgrp->init[0] + pgrp->rate*t)*pull->grp[g].vec[m];
        }
        /* Add the reference position, so we use the correct periodic image */
        dvec_inc(xrefr, dref);
    }

    pbc_dx_d(pbc, xg, xrefr, dr);
    dr2 = 0;
    for (m = 0; m < DIM; m++)
    {
        dr[m] *= pull->dim[m];
        dr2   += dr[m]*dr[m];
    }
    if (max_dist2 >= 0 && dr2 > 0.98*0.98*max_dist2)
    {
        gmx_fatal(FARGS, "Distance of pull group %d (%f nm) is larger than 0.49 times the box size (%f).\nYou might want to consider using \"pull-geometry = direction-periodic\" instead.\n", g, sqrt(dr2), sqrt(max_dist2));

    }

    if (pull->eGeom == epullgDIRPBC)
    {
        dvec_inc(dr, dref);
    }
}

static void get_pullgrp_dr(const t_pull *pull, const t_pbc *pbc, int g, double t,
                           dvec dr)
{
    double md2;

    if (pull->eGeom == epullgDIRPBC)
    {
        md2 = -1;
    }
    else
    {
        md2 = max_pull_distance2(pull, pbc);
    }

    get_pullgrps_dr(pull, pbc, g, t,
                    pull->grp[g].x,
                    PULL_CYL(pull) ? pull->dyna[g].x : pull->grp[0].x,
                    md2,
                    dr);
}

void get_pullgrp_distance(t_pull *pull, t_pbc *pbc, int g, double t,
                          dvec dr, dvec dev)
{
    static gmx_bool bWarned = FALSE; /* TODO: this should be fixed for thread-safety,
                                        but is fairly benign */
    t_pullgrp      *pgrp;
    int             m;
    dvec            ref;
    double          drs, inpr;

    pgrp = &pull->grp[g];

    get_pullgrp_dr(pull, pbc, g, t, dr);

    if (pull->eGeom == epullgPOS)
    {
        for (m = 0; m < DIM; m++)
        {
            ref[m] = pgrp->init[m] + pgrp->rate*t*pgrp->vec[m];
        }
    }
    else
    {
        ref[0] = pgrp->init[0] + pgrp->rate*t;
    }

    switch (pull->eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            if (ref[0] < 0 && !bWarned)
            {
                fprintf(stderr, "\nPull reference distance for group %d is negative (%f)\n", g, ref[0]);
                bWarned = TRUE;
            }
            drs = dnorm(dr);
            if (drs == 0)
            {
                /* With no vector we can not determine the direction for the force,
                 * so we set the force to zero.
                 */
                dev[0] = 0;
            }
            else
            {
                /* Determine the deviation */
                dev[0] = drs - ref[0];
            }
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgCYL:
            /* Pull along vec */
            inpr = 0;
            for (m = 0; m < DIM; m++)
            {
                inpr += pgrp->vec[m]*dr[m];
            }
            dev[0] = inpr - ref[0];
            break;
        case epullgPOS:
            /* Determine the difference of dr and ref along each dimension */
            for (m = 0; m < DIM; m++)
            {
                dev[m] = (dr[m] - ref[m])*pull->dim[m];
            }
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
    for (i = 0; i < 1+pull->ngrp; i++)
    {
        clear_dvec(pull->grp[i].f);
        pull->grp[i].f_scal = 0;
    }
}

/* Apply constraint using SHAKE */
static void do_constraint(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                          rvec *x, rvec *v,
                          gmx_bool bMaster, tensor vir,
                          double dt, double t)
{

    dvec      *r_ij;   /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
    dvec       unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */

    dvec      *rinew;  /* current 'new' position of group i */
    dvec      *rjnew;  /* current 'new' position of group j */
    dvec       ref, vec;
    double     d0, inpr;
    double     lambda, rm, mass, invdt = 0;
    gmx_bool   bConverged_all, bConverged = FALSE;
    int        niter = 0, g, ii, j, m, max_iter = 100;
    double     q, a, b, c; /* for solving the quadratic equation,
                              see Num. Recipes in C ed 2 p. 184 */
    dvec      *dr;         /* correction for group i */
    dvec       ref_dr;     /* correction for group j */
    dvec       f;          /* the pull force */
    dvec       tmp, tmp3;
    t_pullgrp *pdyna, *pgrp, *pref;

    snew(r_ij, pull->ngrp+1);
    if (PULL_CYL(pull))
    {
        snew(rjnew, pull->ngrp+1);
    }
    else
    {
        snew(rjnew, 1);
    }
    snew(dr, pull->ngrp+1);
    snew(rinew, pull->ngrp+1);

    /* copy the current unconstrained positions for use in iterations. We
       iterate until rinew[i] and rjnew[j] obey the constraints. Then
       rinew - pull.x_unc[i] is the correction dr to group i */
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        copy_dvec(pull->grp[g].xp, rinew[g]);
    }
    if (PULL_CYL(pull))
    {
        for (g = 1; g < 1+pull->ngrp; g++)
        {
            copy_dvec(pull->dyna[g].xp, rjnew[g]);
        }
    }
    else
    {
        copy_dvec(pull->grp[0].xp, rjnew[0]);
    }

    /* Determine the constraint directions from the old positions */
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        get_pullgrp_dr(pull, pbc, g, t, r_ij[g]);
        /* Store the difference vector at time t for printing */
        copy_dvec(r_ij[g], pull->grp[g].dr);
        if (debug)
        {
            fprintf(debug, "Pull group %d dr %f %f %f\n",
                    g, r_ij[g][XX], r_ij[g][YY], r_ij[g][ZZ]);
        }

        if (pull->eGeom == epullgDIR || pull->eGeom == epullgDIRPBC)
        {
            /* Select the component along vec */
            a = 0;
            for (m = 0; m < DIM; m++)
            {
                a += pull->grp[g].vec[m]*r_ij[g][m];
            }
            for (m = 0; m < DIM; m++)
            {
                r_ij[g][m] = a*pull->grp[g].vec[m];
            }
        }
    }

    bConverged_all = FALSE;
    while (!bConverged_all && niter < max_iter)
    {
        bConverged_all = TRUE;

        /* loop over all constraints */
        for (g = 1; g < 1+pull->ngrp; g++)
        {
            pgrp = &pull->grp[g];
            if (PULL_CYL(pull))
            {
                pref = &pull->dyna[g];
            }
            else
            {
                pref = &pull->grp[0];
            }

            /* Get the current difference vector */
            get_pullgrps_dr(pull, pbc, g, t, rinew[g], rjnew[PULL_CYL(pull) ? g : 0],
                            -1, unc_ij);

            if (pull->eGeom == epullgPOS)
            {
                for (m = 0; m < DIM; m++)
                {
                    ref[m] = pgrp->init[m] + pgrp->rate*t*pgrp->vec[m];
                }
            }
            else
            {
                ref[0] = pgrp->init[0] + pgrp->rate*t;
                /* Keep the compiler happy */
                ref[1] = 0;
                ref[2] = 0;
            }

            if (debug)
            {
                fprintf(debug, "Pull group %d, iteration %d\n", g, niter);
            }

            rm = 1.0/(pull->grp[g].invtm + pref->invtm);

            switch (pull->eGeom)
            {
                case epullgDIST:
                    if (ref[0] <= 0)
                    {
                        gmx_fatal(FARGS, "The pull constraint reference distance for group %d is <= 0 (%f)", g, ref[0]);
                    }

                    a = diprod(r_ij[g], r_ij[g]);
                    b = diprod(unc_ij, r_ij[g])*2;
                    c = diprod(unc_ij, unc_ij) - dsqr(ref[0]);

                    if (b < 0)
                    {
                        q      = -0.5*(b - sqrt(b*b - 4*a*c));
                        lambda = -q/a;
                    }
                    else
                    {
                        q      = -0.5*(b + sqrt(b*b - 4*a*c));
                        lambda = -c/q;
                    }

                    if (debug)
                    {
                        fprintf(debug,
                                "Pull ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n",
                                a, b, c, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp->invtm, r_ij[g],  dr[g]);
                    dsvmul( lambda*rm*pref->invtm, r_ij[g], ref_dr);
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    /* A 1-dimensional constraint along a vector */
                    a = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pgrp->vec[m];
                        a     += unc_ij[m]*vec[m];
                    }
                    /* Select only the component along the vector */
                    dsvmul(a, vec, unc_ij);
                    lambda = a - ref[0];
                    if (debug)
                    {
                        fprintf(debug, "Pull inpr %e lambda: %e\n", a, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pull->grp[g].invtm, vec, dr[g]);
                    dsvmul( lambda*rm*       pref->invtm, vec, ref_dr);
                    break;
                case epullgPOS:
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->dim[m])
                        {
                            lambda = r_ij[g][m] - ref[m];
                            /* The position corrections dr due to the constraints */
                            dr[g][m]  = -lambda*rm*pull->grp[g].invtm;
                            ref_dr[m] =  lambda*rm*pref->invtm;
                        }
                        else
                        {
                            dr[g][m]  = 0;
                            ref_dr[m] = 0;
                        }
                    }
                    break;
            }

            /* DEBUG */
            if (debug)
            {
                j = (PULL_CYL(pull) ? g : 0);
                get_pullgrps_dr(pull, pbc, g, t, rinew[g], rjnew[j], -1, tmp);
                get_pullgrps_dr(pull, pbc, g, t, dr[g], ref_dr, -1, tmp3);
                fprintf(debug,
                        "Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        rinew[g][0], rinew[g][1], rinew[g][2],
                        rjnew[j][0], rjnew[j][1], rjnew[j][2], dnorm(tmp));
                if (pull->eGeom == epullgPOS)
                {
                    fprintf(debug,
                            "Pull ref %8.5f %8.5f %8.5f\n",
                            pgrp->vec[0], pgrp->vec[1], pgrp->vec[2]);
                }
                else
                {
                    fprintf(debug,
                            "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f %8.5f %8.5f\n",
                            "", "", "", "", "", "", ref[0], ref[1], ref[2]);
                }
                fprintf(debug,
                        "Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        dr[g][0], dr[g][1], dr[g][2],
                        ref_dr[0], ref_dr[1], ref_dr[2],
                        dnorm(tmp3));
                fprintf(debug,
                        "Pull cor %10.7f %10.7f %10.7f\n",
                        dr[g][0], dr[g][1], dr[g][2]);
            } /* END DEBUG */

            /* Update the COMs with dr */
            dvec_inc(rinew[g],                     dr[g]);
            dvec_inc(rjnew[PULL_CYL(pull) ? g : 0], ref_dr);
        }

        /* Check if all constraints are fullfilled now */
        for (g = 1; g < 1+pull->ngrp; g++)
        {
            pgrp = &pull->grp[g];

            if (pull->eGeom == epullgPOS)
            {
                for (m = 0; m < DIM; m++)
                {
                    ref[m] = pgrp->init[m] + pgrp->rate*t*pgrp->vec[m];
                }
            }
            else
            {
                ref[0] = pgrp->init[0] + pgrp->rate*t;
                /* Keep the compiler happy */
                ref[1] = 0;
                ref[2] = 0;
            }

            get_pullgrps_dr(pull, pbc, g, t, rinew[g], rjnew[PULL_CYL(pull) ? g : 0],
                            -1, unc_ij);

            switch (pull->eGeom)
            {
                case epullgDIST:
                    bConverged = fabs(dnorm(unc_ij) - ref[0]) < pull->constr_tol;
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pgrp->vec[m];
                    }
                    inpr = diprod(unc_ij, vec);
                    dsvmul(inpr, vec, unc_ij);
                    bConverged =
                        fabs(diprod(unc_ij, vec) - ref[0]) < pull->constr_tol;
                    break;
                case epullgPOS:
                    bConverged = TRUE;
                    for (m = 0; m < DIM; m++)
                    {
                        if (pull->dim[m] &&
                            fabs(unc_ij[m] - ref[m]) >= pull->constr_tol)
                        {
                            bConverged = FALSE;
                        }
                    }
                    break;
            }

            if (!bConverged)
            {
                if (debug)
                {
                    fprintf(debug, "NOT CONVERGED YET: Group %d:"
                            "d_ref = %f %f %f, current d = %f\n",
                            g, ref[0], ref[1], ref[2], dnorm(unc_ij));
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

    /* update the normal groups */
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        /* get the final dr and constraint force for group i */
        dvec_sub(rinew[g], pgrp->xp, dr[g]);
        /* select components of dr */
        for (m = 0; m < DIM; m++)
        {
            dr[g][m] *= pull->dim[m];
        }
        dsvmul(1.0/(pgrp->invtm*dt*dt), dr[g], f);
        dvec_inc(pgrp->f, f);
        switch (pull->eGeom)
        {
            case epullgDIST:
                for (m = 0; m < DIM; m++)
                {
                    pgrp->f_scal += r_ij[g][m]*f[m]/dnorm(r_ij[g]);
                }
                break;
            case epullgDIR:
            case epullgDIRPBC:
            case epullgCYL:
                for (m = 0; m < DIM; m++)
                {
                    pgrp->f_scal += pgrp->vec[m]*f[m];
                }
                break;
            case epullgPOS:
                break;
        }

        if (vir && bMaster)
        {
            /* Add the pull contribution to the virial */
            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5*f[j]*r_ij[g][m];
                }
            }
        }

        /* update the atom positions */
        copy_dvec(dr[g], tmp);
        for (j = 0; j < pgrp->nat_loc; j++)
        {
            ii = pgrp->ind_loc[j];
            if (pgrp->weight_loc)
            {
                dsvmul(pgrp->wscale*pgrp->weight_loc[j], dr[g], tmp);
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

    /* update the reference groups */
    if (PULL_CYL(pull))
    {
        /* update the dynamic reference groups */
        for (g = 1; g < 1+pull->ngrp; g++)
        {
            pdyna = &pull->dyna[g];
            dvec_sub(rjnew[g], pdyna->xp, ref_dr);
            /* select components of ref_dr */
            for (m = 0; m < DIM; m++)
            {
                ref_dr[m] *= pull->dim[m];
            }

            for (j = 0; j < pdyna->nat_loc; j++)
            {
                /* reset the atoms with dr, weighted by w_i */
                dsvmul(pdyna->wscale*pdyna->weight_loc[j], ref_dr, tmp);
                ii = pdyna->ind_loc[j];
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
    }
    else
    {
        pgrp = &pull->grp[0];
        /* update the reference group */
        dvec_sub(rjnew[0], pgrp->xp, ref_dr);
        /* select components of ref_dr */
        for (m = 0; m < DIM; m++)
        {
            ref_dr[m] *= pull->dim[m];
        }

        copy_dvec(ref_dr, tmp);
        for (j = 0; j < pgrp->nat_loc; j++)
        {
            ii = pgrp->ind_loc[j];
            if (pgrp->weight_loc)
            {
                dsvmul(pgrp->wscale*pgrp->weight_loc[j], ref_dr, tmp);
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

    /* finished! I hope. Give back some memory */
    sfree(r_ij);
    sfree(rinew);
    sfree(rjnew);
    sfree(dr);
}

/* Pulling with a harmonic umbrella potential or constant force */
static void do_pull_pot(int ePull,
                        t_pull *pull, t_pbc *pbc, double t, real lambda,
                        real *V, tensor vir, real *dVdl)
{
    int        g, j, m;
    dvec       dev;
    double     ndr, invdr;
    real       k, dkdl;
    t_pullgrp *pgrp;

    /* loop over the groups that are being pulled */
    *V    = 0;
    *dVdl = 0;
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        get_pullgrp_distance(pull, pbc, g, t, pgrp->dr, dev);

        k    = (1.0 - lambda)*pgrp->k + lambda*pgrp->kB;
        dkdl = pgrp->kB - pgrp->k;

        switch (pull->eGeom)
        {
            case epullgDIST:
                ndr   = dnorm(pgrp->dr);
                invdr = 1/ndr;
                if (ePull == epullUMBRELLA)
                {
                    pgrp->f_scal  =       -k*dev[0];
                    *V           += 0.5*   k*dsqr(dev[0]);
                    *dVdl        += 0.5*dkdl*dsqr(dev[0]);
                }
                else
                {
                    pgrp->f_scal  =   -k;
                    *V           +=    k*ndr;
                    *dVdl        += dkdl*ndr;
                }
                for (m = 0; m < DIM; m++)
                {
                    pgrp->f[m]    = pgrp->f_scal*pgrp->dr[m]*invdr;
                }
                break;
            case epullgDIR:
            case epullgDIRPBC:
            case epullgCYL:
                if (ePull == epullUMBRELLA)
                {
                    pgrp->f_scal  =       -k*dev[0];
                    *V           += 0.5*   k*dsqr(dev[0]);
                    *dVdl        += 0.5*dkdl*dsqr(dev[0]);
                }
                else
                {
                    ndr = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        ndr += pgrp->vec[m]*pgrp->dr[m];
                    }
                    pgrp->f_scal  =   -k;
                    *V           +=    k*ndr;
                    *dVdl        += dkdl*ndr;
                }
                for (m = 0; m < DIM; m++)
                {
                    pgrp->f[m]    = pgrp->f_scal*pgrp->vec[m];
                }
                break;
            case epullgPOS:
                for (m = 0; m < DIM; m++)
                {
                    if (ePull == epullUMBRELLA)
                    {
                        pgrp->f[m]  =       -k*dev[m];
                        *V         += 0.5*   k*dsqr(dev[m]);
                        *dVdl      += 0.5*dkdl*dsqr(dev[m]);
                    }
                    else
                    {
                        pgrp->f[m]  =   -k*pull->dim[m];
                        *V         +=    k*pgrp->dr[m]*pull->dim[m];
                        *dVdl      += dkdl*pgrp->dr[m]*pull->dim[m];
                    }
                }
                break;
        }

        if (vir)
        {
            /* Add the pull contribution to the virial */
            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5*pgrp->f[j]*pgrp->dr[m];
                }
            }
        }
    }
}

real pull_potential(int ePull, t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, rvec *f, tensor vir, real *dvdlambda)
{
    real V, dVdl;

    pull_calc_coms(cr, pull, md, pbc, t, x, NULL);

    do_pull_pot(ePull, pull, pbc, t, lambda,
                &V, pull->bVirial && MASTER(cr) ? vir : NULL, &dVdl);

    /* Distribute forces over pulled groups */
    apply_forces(pull, md, DOMAINDECOMP(cr) ? cr->dd->ga2la : NULL, f);

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

    do_constraint(pull, md, pbc, xp, v, pull->bVirial && MASTER(cr), vir, dt, t);
}

static void make_local_pull_group(gmx_ga2la_t ga2la,
                                  t_pullgrp *pg, int start, int end)
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

    if (pull->grp[0].nat > 0)
    {
        make_local_pull_group(ga2la, &pull->grp[0], md->start, md->start+md->homenr);
    }
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        make_local_pull_group(ga2la, &pull->grp[g], md->start, md->start+md->homenr);
    }
}

static void init_pull_group_index(FILE *fplog, t_commrec *cr,
                                  int start, int end,
                                  int g, t_pullgrp *pg, ivec pulldims,
                                  gmx_mtop_t *mtop, t_inputrec *ir, real lambda)
{
    int                   i, ii, d, nfrozen, ndim;
    real                  m, w, mbd;
    double                tmass, wmass, wwmass;
    gmx_bool              bDomDec;
    gmx_ga2la_t           ga2la = NULL;
    gmx_groups_t         *groups;
    gmx_mtop_atomlookup_t alook;
    t_atom               *atom;

    bDomDec = (cr && DOMAINDECOMP(cr));
    if (bDomDec)
    {
        ga2la = cr->dd->ga2la;
    }

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
        if (cr && PAR(cr) && !bDomDec && ii >= start && ii < end)
        {
            pg->ind_loc[pg->nat_loc++] = ii;
        }
        if (ir->opts.nFreeze)
        {
            for (d = 0; d < DIM; d++)
            {
                if (pulldims[d] && ir->opts.nFreeze[ggrpnr(groups, egcFREEZE, ii)][d])
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

    if (nfrozen == 0)
    {
        /* A value > 0 signals not frozen, it is updated later */
        pg->invtm  = 1.0;
    }
    else
    {
        ndim = 0;
        for (d = 0; d < DIM; d++)
        {
            ndim += pulldims[d]*pg->nat;
        }
        if (fplog && nfrozen > 0 && nfrozen < ndim)
        {
            fprintf(fplog,
                    "\nWARNING: In pull group %d some, but not all of the degrees of freedom\n"
                    "         that are subject to pulling are frozen.\n"
                    "         For pulling the whole group will be frozen.\n\n",
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
    t_pullgrp    *pgrp;
    int           g, start = 0, end = 0, m;
    gmx_bool      bCite;

    pull = ir->pull;

    pull->ePBC = ir->ePBC;
    switch (pull->ePBC)
    {
        case epbcNONE: pull->npbcdim = 0; break;
        case epbcXY:   pull->npbcdim = 2; break;
        default:       pull->npbcdim = 3; break;
    }

    if (fplog)
    {
        fprintf(fplog, "\nWill apply %s COM pulling in geometry '%s'\n",
                EPULLTYPE(ir->ePull), EPULLGEOM(pull->eGeom));
        if (pull->grp[0].nat > 0)
        {
            fprintf(fplog, "between a reference group and %d group%s\n",
                    pull->ngrp, pull->ngrp == 1 ? "" : "s");
        }
        else
        {
            fprintf(fplog, "with an absolute reference on %d group%s\n",
                    pull->ngrp, pull->ngrp == 1 ? "" : "s");
        }
        bCite = FALSE;
        for (g = 0; g < pull->ngrp+1; g++)
        {
            if (pull->grp[g].nat > 1 &&
                pull->grp[g].pbcatom < 0)
            {
                /* We are using cosine weighting */
                fprintf(fplog, "Cosine weighting is used for group %d\n", g);
                bCite = TRUE;
            }
        }
        if (bCite)
        {
            please_cite(fplog, "Engin2010");
        }
    }

    /* We always add the virial contribution,
     * except for geometry = direction_periodic where this is impossible.
     */
    pull->bVirial = (pull->eGeom != epullgDIRPBC);
    if (getenv("GMX_NO_PULLVIR") != NULL)
    {
        if (fplog)
        {
            fprintf(fplog, "Found env. var., will not add the virial contribution of the COM pull forces\n");
        }
        pull->bVirial = FALSE;
    }

    if (cr && PARTDECOMP(cr))
    {
        pd_at_range(cr, &start, &end);
    }
    pull->rbuf     = NULL;
    pull->dbuf     = NULL;
    pull->dbuf_cyl = NULL;
    pull->bRefAt   = FALSE;
    pull->cosdim   = -1;
    for (g = 0; g < pull->ngrp+1; g++)
    {
        pgrp           = &pull->grp[g];
        pgrp->epgrppbc = epgrppbcNONE;
        if (pgrp->nat > 0)
        {
            /* Determine if we need to take PBC into account for calculating
             * the COM's of the pull groups.
             */
            for (m = 0; m < pull->npbcdim; m++)
            {
                if (pull->dim[m] && pgrp->nat > 1)
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
                            gmx_fatal(FARGS, "Can only use cosine weighting with pulling in one dimension (use mdp option pull_dim)");
                        }
                        pull->cosdim = m;
                    }
                }
            }
            /* Set the indices */
            init_pull_group_index(fplog, cr, start, end, g, pgrp, pull->dim, mtop, ir, lambda);
            if (PULL_CYL(pull) && pgrp->invtm == 0)
            {
                gmx_fatal(FARGS, "Can not have frozen atoms in a cylinder pull group");
            }
        }
        else
        {
            /* Absolute reference, set the inverse mass to zero */
            pgrp->invtm  = 0;
            pgrp->wscale = 1;
        }
    }

    /* if we use dynamic reference groups, do some initialising for them */
    if (PULL_CYL(pull))
    {
        if (pull->grp[0].nat == 0)
        {
            gmx_fatal(FARGS, "Dynamic reference groups are not supported when using absolute reference!\n");
        }
        snew(pull->dyna, pull->ngrp+1);
    }

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

void finish_pull(FILE *fplog, t_pull *pull)
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
