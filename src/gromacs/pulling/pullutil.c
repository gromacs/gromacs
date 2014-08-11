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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

#include "sysstuff.h"
#include "princ.h"
#include "gromacs/fileio/futil.h"
#include "vec.h"
#include "gromacs/utility/smalloc.h"
#include "typedefs.h"
#include "types/commrec.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "symtab.h"
#include "index.h"
#include "gromacs/fileio/confio.h"
#include "network.h"
#include "pbc.h"
#include "pull.h"
#include "gmx_ga2la.h"

static void pull_set_pbcatom(t_commrec *cr, t_pull_group *pgrp,
                             rvec *x,
                             rvec x_pbc)
{
    int a, m;

    if (cr != NULL && DOMAINDECOMP(cr))
    {
        if (ga2la_get_home(cr->dd->ga2la, pgrp->pbcatom, &a))
        {
            copy_rvec(x[a], x_pbc);
        }
        else
        {
            clear_rvec(x_pbc);
        }
    }
    else
    {
        copy_rvec(x[pgrp->pbcatom], x_pbc);
    }
}

static void pull_set_pbcatoms(t_commrec *cr, t_pull *pull,
                              rvec *x,
                              rvec *x_pbc)
{
    int g, n, m;

    n = 0;
    for (g = 0; g < pull->ngroup; g++)
    {
        if ((g == 0 && PULL_CYL(pull)) || pull->group[g].pbcatom == -1)
        {
            clear_rvec(x_pbc[g]);
        }
        else
        {
            pull_set_pbcatom(cr, &pull->group[g], x, x_pbc[g]);
            for (m = 0; m < DIM; m++)
            {
                if (pull->dim[m] == 0)
                {
                    x_pbc[g][m] = 0.0;
                }
            }
            n++;
        }
    }

    if (cr && PAR(cr) && n > 0)
    {
        /* Sum over the nodes to get x_pbc from the home node of pbcatom */
        gmx_sum(pull->ngroup*DIM, x_pbc[0], cr);
    }
}

/* switch function, x between r and w */
static real get_weight(real x, real r1, real r0)
{
    real weight;

    if (x >= r0)
    {
        weight = 0;
    }
    else if (x <= r1)
    {
        weight = 1;
    }
    else
    {
        weight = (r0 - x)/(r0 - r1);
    }

    return weight;
}

static void make_cyl_refgrps(t_commrec *cr, t_pull *pull, t_mdatoms *md,
                             t_pbc *pbc, double t, rvec *x, rvec *xp)
{
    int           c, i, ii, m, start, end;
    rvec          g_x, dx, dir;
    double        r0_2, sum_a, sum_ap, dr2, mass, weight, wmass, wwmass, inp;
    t_pull_coord *pcrd;
    t_pull_group *pref, *pgrp, *pdyna;
    gmx_ga2la_t   ga2la = NULL;

    if (pull->dbuf_cyl == NULL)
    {
        snew(pull->dbuf_cyl, pull->ncoord*4);
    }

    if (cr && DOMAINDECOMP(cr))
    {
        ga2la = cr->dd->ga2la;
    }

    start = 0;
    end   = md->homenr;

    r0_2 = dsqr(pull->cyl_r0);

    /* loop over all groups to make a reference group for each*/
    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd  = &pull->coord[c];

        /* pref will be the same group for all pull coordinates */
        pref  = &pull->group[pcrd->group[0]];
        pgrp  = &pull->group[pcrd->group[1]];
        pdyna = &pull->dyna[c];
        copy_rvec(pcrd->vec, dir);
        sum_a          = 0;
        sum_ap         = 0;
        wmass          = 0;
        wwmass         = 0;
        pdyna->nat_loc = 0;

        for (m = 0; m < DIM; m++)
        {
            g_x[m] = pgrp->x[m] - pcrd->vec[m]*(pcrd->init + pcrd->rate*t);
        }

        /* loop over all atoms in the main ref group */
        for (i = 0; i < pref->nat; i++)
        {
            ii = pref->ind[i];
            if (ga2la)
            {
                if (!ga2la_get_home(ga2la, pref->ind[i], &ii))
                {
                    ii = -1;
                }
            }
            if (ii >= start && ii < end)
            {
                pbc_dx_aiuc(pbc, x[ii], g_x, dx);
                inp = iprod(dir, dx);
                dr2 = 0;
                for (m = 0; m < DIM; m++)
                {
                    dr2 += dsqr(dx[m] - inp*dir[m]);
                }

                if (dr2 < r0_2)
                {
                    /* add to index, to sum of COM, to weight array */
                    if (pdyna->nat_loc >= pdyna->nalloc_loc)
                    {
                        pdyna->nalloc_loc = over_alloc_large(pdyna->nat_loc+1);
                        srenew(pdyna->ind_loc, pdyna->nalloc_loc);
                        srenew(pdyna->weight_loc, pdyna->nalloc_loc);
                    }
                    pdyna->ind_loc[pdyna->nat_loc] = ii;
                    mass   = md->massT[ii];
                    weight = get_weight(sqrt(dr2), pull->cyl_r1, pull->cyl_r0);
                    pdyna->weight_loc[pdyna->nat_loc] = weight;
                    sum_a += mass*weight*inp;
                    if (xp)
                    {
                        pbc_dx_aiuc(pbc, xp[ii], g_x, dx);
                        inp     = iprod(dir, dx);
                        sum_ap += mass*weight*inp;
                    }
                    wmass  += mass*weight;
                    wwmass += mass*sqr(weight);
                    pdyna->nat_loc++;
                }
            }
        }
        pull->dbuf_cyl[c*4+0] = wmass;
        pull->dbuf_cyl[c*4+1] = wwmass;
        pull->dbuf_cyl[c*4+2] = sum_a;
        pull->dbuf_cyl[c*4+3] = sum_ap;
    }

    if (cr && PAR(cr))
    {
        /* Sum the contributions over the nodes */
        gmx_sumd(pull->ncoord*4, pull->dbuf_cyl, cr);
    }

    for (c = 0; c < pull->ncoord; c++)
    {
        pcrd  = &pull->coord[c];

        pdyna = &pull->dyna[c];
        pgrp  = &pull->group[pcrd->group[1]];

        wmass         = pull->dbuf_cyl[c*4+0];
        wwmass        = pull->dbuf_cyl[c*4+1];
        pdyna->wscale = wmass/wwmass;
        pdyna->invtm  = 1.0/(pdyna->wscale*wmass);

        for (m = 0; m < DIM; m++)
        {
            g_x[m]      = pgrp->x[m] - pcrd->vec[m]*(pcrd->init + pcrd->rate*t);
            pdyna->x[m] = g_x[m] + pcrd->vec[m]*pull->dbuf_cyl[c*4+2]/wmass;
            if (xp)
            {
                pdyna->xp[m] = g_x[m] + pcrd->vec[m]*pull->dbuf_cyl[c*4+3]/wmass;
            }
        }

        if (debug)
        {
            fprintf(debug, "Pull cylinder group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
                    c, pdyna->x[0], pdyna->x[1],
                    pdyna->x[2], 1.0/pdyna->invtm);
        }
    }
}

static double atan2_0_2pi(double y, double x)
{
    double a;

    a = atan2(y, x);
    if (a < 0)
    {
        a += 2.0*M_PI;
    }
    return a;
}

/* calculates center of mass of selection index from all coordinates x */
void pull_calc_coms(t_commrec *cr,
                    t_pull *pull, t_mdatoms *md, t_pbc *pbc, double t,
                    rvec x[], rvec *xp)
{
    int           g, i, ii, m;
    real          mass, w, wm, twopi_box = 0;
    double        wmass, wwmass, invwmass;
    dvec          com, comp;
    double        cm, sm, cmp, smp, ccm, csm, ssm, csw, snw;
    rvec         *xx[2], x_pbc = {0, 0, 0}, dx;
    t_pull_group *pgrp;

    if (pull->rbuf == NULL)
    {
        snew(pull->rbuf, pull->ngroup);
    }
    if (pull->dbuf == NULL)
    {
        snew(pull->dbuf, 3*pull->ngroup);
    }

    if (pull->bRefAt && pull->bSetPBCatoms)
    {
        pull_set_pbcatoms(cr, pull, x, pull->rbuf);

        if (cr != NULL && DOMAINDECOMP(cr))
        {
            /* We can keep these PBC reference coordinates fixed for nstlist
             * steps, since atoms won't jump over PBC.
             * This avoids a global reduction at the next nstlist-1 steps.
             * Note that the exact values of the pbc reference coordinates
             * are irrelevant, as long all atoms in the group are within
             * half a box distance of the reference coordinate.
             */
            pull->bSetPBCatoms = FALSE;
        }
    }

    if (pull->cosdim >= 0)
    {
        for (m = pull->cosdim+1; m < pull->npbcdim; m++)
        {
            if (pbc->box[m][pull->cosdim] != 0)
            {
                gmx_fatal(FARGS, "Can not do cosine weighting for trilinic dimensions");
            }
        }
        twopi_box = 2.0*M_PI/pbc->box[pull->cosdim][pull->cosdim];
    }

    for (g = 0; g < pull->ngroup; g++)
    {
        pgrp = &pull->group[g];
        clear_dvec(com);
        clear_dvec(comp);
        wmass  = 0;
        wwmass = 0;
        cm     = 0;
        sm     = 0;
        cmp    = 0;
        smp    = 0;
        ccm    = 0;
        csm    = 0;
        ssm    = 0;
        if (!(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc == epgrppbcREFAT)
            {
                /* Set the pbc atom */
                copy_rvec(pull->rbuf[g], x_pbc);
            }
            w = 1;
            for (i = 0; i < pgrp->nat_loc; i++)
            {
                ii   = pgrp->ind_loc[i];
                mass = md->massT[ii];
                if (pgrp->epgrppbc != epgrppbcCOS)
                {
                    if (pgrp->weight_loc)
                    {
                        w = pgrp->weight_loc[i];
                    }
                    wm      = w*mass;
                    wmass  += wm;
                    wwmass += wm*w;
                    if (pgrp->epgrppbc == epgrppbcNONE)
                    {
                        /* Plain COM: sum the coordinates */
                        for (m = 0; m < DIM; m++)
                        {
                            com[m]    += wm*x[ii][m];
                        }
                        if (xp)
                        {
                            for (m = 0; m < DIM; m++)
                            {
                                comp[m] += wm*xp[ii][m];
                            }
                        }
                    }
                    else
                    {
                        /* Sum the difference with the reference atom */
                        pbc_dx(pbc, x[ii], x_pbc, dx);
                        for (m = 0; m < DIM; m++)
                        {
                            com[m]    += wm*dx[m];
                        }
                        if (xp)
                        {
                            /* For xp add the difference between xp and x to dx,
                             * such that we use the same periodic image,
                             * also when xp has a large displacement.
                             */
                            for (m = 0; m < DIM; m++)
                            {
                                comp[m] += wm*(dx[m] + xp[ii][m] - x[ii][m]);
                            }
                        }
                    }
                }
                else
                {
                    /* Determine cos and sin sums */
                    csw  = cos(x[ii][pull->cosdim]*twopi_box);
                    snw  = sin(x[ii][pull->cosdim]*twopi_box);
                    cm  += csw*mass;
                    sm  += snw*mass;
                    ccm += csw*csw*mass;
                    csm += csw*snw*mass;
                    ssm += snw*snw*mass;

                    if (xp)
                    {
                        csw  = cos(xp[ii][pull->cosdim]*twopi_box);
                        snw  = sin(xp[ii][pull->cosdim]*twopi_box);
                        cmp += csw*mass;
                        smp += snw*mass;
                    }
                }
            }
        }

        /* Copy local sums to a buffer for global summing */
        switch (pgrp->epgrppbc)
        {
            case epgrppbcNONE:
            case epgrppbcREFAT:
                copy_dvec(com, pull->dbuf[g*3]);
                copy_dvec(comp, pull->dbuf[g*3+1]);
                pull->dbuf[g*3+2][0] = wmass;
                pull->dbuf[g*3+2][1] = wwmass;
                pull->dbuf[g*3+2][2] = 0;
                break;
            case epgrppbcCOS:
                pull->dbuf[g*3  ][0] = cm;
                pull->dbuf[g*3  ][1] = sm;
                pull->dbuf[g*3  ][2] = 0;
                pull->dbuf[g*3+1][0] = ccm;
                pull->dbuf[g*3+1][1] = csm;
                pull->dbuf[g*3+1][2] = ssm;
                pull->dbuf[g*3+2][0] = cmp;
                pull->dbuf[g*3+2][1] = smp;
                pull->dbuf[g*3+2][2] = 0;
                break;
        }
    }

    if (cr && PAR(cr))
    {
        /* Sum the contributions over the nodes */
        gmx_sumd(pull->ngroup*3*DIM, pull->dbuf[0], cr);
    }

    for (g = 0; g < pull->ngroup; g++)
    {
        pgrp = &pull->group[g];
        if (pgrp->nat > 0 && !(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                /* Determine the inverse mass */
                wmass    = pull->dbuf[g*3+2][0];
                wwmass   = pull->dbuf[g*3+2][1];
                invwmass = 1/wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm > 0)
                {
                    pgrp->wscale = wmass/wwmass;
                    pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
                }
                /* Divide by the total mass */
                for (m = 0; m < DIM; m++)
                {
                    pgrp->x[m]    = pull->dbuf[g*3  ][m]*invwmass;
                    if (xp)
                    {
                        pgrp->xp[m] = pull->dbuf[g*3+1][m]*invwmass;
                    }
                    if (pgrp->epgrppbc == epgrppbcREFAT)
                    {
                        pgrp->x[m]    += pull->rbuf[g][m];
                        if (xp)
                        {
                            pgrp->xp[m] += pull->rbuf[g][m];
                        }
                    }
                }
            }
            else
            {
                /* Determine the optimal location of the cosine weight */
                csw                   = pull->dbuf[g*3][0];
                snw                   = pull->dbuf[g*3][1];
                pgrp->x[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                /* Set the weights for the local atoms */
                wmass  = sqrt(csw*csw + snw*snw);
                wwmass = (pull->dbuf[g*3+1][0]*csw*csw +
                          pull->dbuf[g*3+1][1]*csw*snw +
                          pull->dbuf[g*3+1][2]*snw*snw)/(wmass*wmass);
                pgrp->wscale = wmass/wwmass;
                pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
                /* Set the weights for the local atoms */
                csw *= pgrp->invtm;
                snw *= pgrp->invtm;
                for (i = 0; i < pgrp->nat_loc; i++)
                {
                    ii                  = pgrp->ind_loc[i];
                    pgrp->weight_loc[i] = csw*cos(twopi_box*x[ii][pull->cosdim]) +
                        snw*sin(twopi_box*x[ii][pull->cosdim]);
                }
                if (xp)
                {
                    csw                    = pull->dbuf[g*3+2][0];
                    snw                    = pull->dbuf[g*3+2][1];
                    pgrp->xp[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                }
            }
            if (debug)
            {
                fprintf(debug, "Pull group %d wmass %f wwmass %f invtm %f\n",
                        g, wmass, wwmass, pgrp->invtm);
            }
        }
    }

    if (PULL_CYL(pull))
    {
        /* Calculate the COMs for the cyclinder reference groups */
        make_cyl_refgrps(cr, pull, md, pbc, t, x, xp);
    }
}
