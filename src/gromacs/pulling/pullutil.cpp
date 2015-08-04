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

#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static void pull_reduce_real(t_commrec   *cr,
                             pull_comm_t *comm,
                             int          n,
                             real        *data)
{
    if (cr != NULL && PAR(cr))
    {
        if (comm->bParticipateAll)
        {
            /* Sum the contributions over all DD ranks */
            gmx_sum(n, data, cr);
        }
        else
        {
#ifdef GMX_MPI
#ifdef MPI_IN_PLACE_EXISTS
            MPI_Allreduce(MPI_IN_PLACE, data, n, GMX_MPI_REAL, MPI_SUM,
                          comm->mpi_comm_com);
#else
            real *buf;

            snew(buf, n);

            MPI_Allreduce(data, buf, n, GMX_MPI_REAL, MPI_SUM,
                          comm->mpi_comm_com);

            /* Copy the result from the buffer to the input/output data */
            for (int i = 0; i < n; i++)
            {
                data[i] = buf[i];
            }
            sfree(buf);
#endif
#else
            gmx_incons("comm->bParticipateAll=FALSE without GMX_MPI");
#endif
        }
    }
}

static void pull_reduce_double(t_commrec   *cr,
                               pull_comm_t *comm,
                               int          n,
                               double      *data)
{
    if (cr != NULL && PAR(cr))
    {
        if (comm->bParticipateAll)
        {
            /* Sum the contributions over all DD ranks */
            gmx_sumd(n, data, cr);
        }
        else
        {
#ifdef GMX_MPI
#ifdef MPI_IN_PLACE_EXISTS
            MPI_Allreduce(MPI_IN_PLACE, data, n, MPI_DOUBLE, MPI_SUM,
                          comm->mpi_comm_com);
#else
            double *buf;

            snew(buf, n);

            MPI_Allreduce(data, buf, n, MPI_DOUBLE, MPI_SUM,
                          comm->mpi_comm_com);

            /* Copy the result from the buffer to the input/output data */
            for (int i = 0; i < n; i++)
            {
                data[i] = buf[i];
            }
            sfree(buf);
#endif
#else
            gmx_incons("comm->bParticipateAll=FALSE without GMX_MPI");
#endif
        }
    }
}

static void pull_set_pbcatom(t_commrec *cr, pull_group_work_t *pgrp,
                             rvec *x,
                             rvec x_pbc)
{
    int a;

    if (cr != NULL && DOMAINDECOMP(cr))
    {
        if (ga2la_get_home(cr->dd->ga2la, pgrp->params.pbcatom, &a))
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
        copy_rvec(x[pgrp->params.pbcatom], x_pbc);
    }
}

static void pull_set_pbcatoms(t_commrec *cr, struct pull_t *pull,
                              rvec *x,
                              rvec *x_pbc)
{
    int g, n;

    n = 0;
    for (g = 0; g < pull->ngroup; g++)
    {
        if (!pull->group[g].bCalcCOM || pull->group[g].params.pbcatom == -1)
        {
            clear_rvec(x_pbc[g]);
        }
        else
        {
            pull_set_pbcatom(cr, &pull->group[g], x, x_pbc[g]);
            n++;
        }
    }

    if (cr && PAR(cr) && n > 0)
    {
        /* Sum over participating ranks to get x_pbc from the home ranks.
         * This can be very expensive at high parallelization, so we only
         * do this after each DD repartitioning.
         */
        pull_reduce_real(cr, &pull->comm, pull->ngroup*DIM, x_pbc[0]);
    }
}

static void make_cyl_refgrps(t_commrec *cr, struct pull_t *pull, t_mdatoms *md,
                             t_pbc *pbc, double t, rvec *x)
{
    /* The size and stride per coord for the reduction buffer */
    const int     stride = 9;
    int           c, i, ii, m, start, end;
    rvec          g_x, dx, dir;
    double        inv_cyl_r2;
    pull_comm_t  *comm;
    gmx_ga2la_t   ga2la = NULL;

    comm = &pull->comm;

    if (comm->dbuf_cyl == NULL)
    {
        snew(comm->dbuf_cyl, pull->ncoord*stride);
    }

    if (cr && DOMAINDECOMP(cr))
    {
        ga2la = cr->dd->ga2la;
    }

    start = 0;
    end   = md->homenr;

    inv_cyl_r2 = 1/dsqr(pull->params.cylinder_r);

    /* loop over all groups to make a reference group for each*/
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;
        double             sum_a, wmass, wwmass;
        dvec               radf_fac0, radf_fac1;

        pcrd   = &pull->coord[c];

        sum_a  = 0;
        wmass  = 0;
        wwmass = 0;
        clear_dvec(radf_fac0);
        clear_dvec(radf_fac1);

        if (pcrd->params.eGeom == epullgCYL)
        {
            pull_group_work_t *pref, *pgrp, *pdyna;

            /* pref will be the same group for all pull coordinates */
            pref  = &pull->group[pcrd->params.group[0]];
            pgrp  = &pull->group[pcrd->params.group[1]];
            pdyna = &pull->dyna[c];
            copy_rvec(pcrd->vec, dir);
            pdyna->nat_loc = 0;

            /* We calculate distances with respect to the reference location
             * of this cylinder group (g_x), which we already have now since
             * we reduced the other group COM over the ranks. This resolves
             * any PBC issues and we don't need to use a PBC-atom here.
             */
            if (pcrd->params.rate != 0)
            {
                /* With rate=0, value_ref is set initially */
                pcrd->value_ref = pcrd->params.init + pcrd->params.rate*t;
            }
            for (m = 0; m < DIM; m++)
            {
                g_x[m] = pgrp->x[m] - pcrd->vec[m]*pcrd->value_ref;
            }

            /* loop over all atoms in the main ref group */
            for (i = 0; i < pref->params.nat; i++)
            {
                ii = pref->params.ind[i];
                if (ga2la)
                {
                    if (!ga2la_get_home(ga2la, pref->params.ind[i], &ii))
                    {
                        ii = -1;
                    }
                }
                if (ii >= start && ii < end)
                {
                    double dr2, dr2_rel, inp;
                    dvec   dr;

                    pbc_dx_aiuc(pbc, x[ii], g_x, dx);
                    inp = iprod(dir, dx);
                    dr2 = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        /* Determine the radial components */
                        dr[m] = dx[m] - inp*dir[m];
                        dr2  += dr[m]*dr[m];
                    }
                    dr2_rel = dr2*inv_cyl_r2;

                    if (dr2_rel < 1)
                    {
                        double mass, weight, dweight_r;
                        dvec   mdw;

                        /* add to index, to sum of COM, to weight array */
                        if (pdyna->nat_loc >= pdyna->nalloc_loc)
                        {
                            pdyna->nalloc_loc = over_alloc_large(pdyna->nat_loc+1);
                            srenew(pdyna->ind_loc,    pdyna->nalloc_loc);
                            srenew(pdyna->weight_loc, pdyna->nalloc_loc);
                            srenew(pdyna->mdw,        pdyna->nalloc_loc);
                            srenew(pdyna->dv,         pdyna->nalloc_loc);
                        }
                        pdyna->ind_loc[pdyna->nat_loc] = ii;

                        mass      = md->massT[ii];
                        /* The radial weight function is 1-2x^2+x^4,
                         * where x=r/cylinder_r. Since this function depends
                         * on the radial component, we also get radial forces
                         * on both groups.
                         */
                        weight    = 1 + (-2 + dr2_rel)*dr2_rel;
                        dweight_r = (-4 + 4*dr2_rel)*inv_cyl_r2;
                        pdyna->weight_loc[pdyna->nat_loc] = weight;
                        sum_a    += mass*weight*inp;
                        wmass    += mass*weight;
                        wwmass   += mass*weight*weight;
                        dsvmul(mass*dweight_r, dr, mdw);
                        copy_dvec(mdw, pdyna->mdw[pdyna->nat_loc]);
                        /* Currently we only have the axial component of the
                         * distance (inp) up to an unkown offset. We add this
                         * offset after the reduction needs to determine the
                         * COM of the cylinder group.
                         */
                        pdyna->dv[pdyna->nat_loc] = inp;
                        for (m = 0; m < DIM; m++)
                        {
                            radf_fac0[m] += mdw[m];
                            radf_fac1[m] += mdw[m]*inp;
                        }
                        pdyna->nat_loc++;
                    }
                }
            }
        }
        comm->dbuf_cyl[c*stride+0] = wmass;
        comm->dbuf_cyl[c*stride+1] = wwmass;
        comm->dbuf_cyl[c*stride+2] = sum_a;
        comm->dbuf_cyl[c*stride+3] = radf_fac0[XX];
        comm->dbuf_cyl[c*stride+4] = radf_fac0[YY];
        comm->dbuf_cyl[c*stride+5] = radf_fac0[ZZ];
        comm->dbuf_cyl[c*stride+6] = radf_fac1[XX];
        comm->dbuf_cyl[c*stride+7] = radf_fac1[YY];
        comm->dbuf_cyl[c*stride+8] = radf_fac1[ZZ];
    }

    if (cr != NULL && PAR(cr))
    {
        /* Sum the contributions over the ranks */
        pull_reduce_double(cr, comm, pull->ncoord*stride, comm->dbuf_cyl);
    }

    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd  = &pull->coord[c];

        if (pcrd->params.eGeom == epullgCYL)
        {
            pull_group_work_t *pdyna, *pgrp;
            double             wmass, wwmass, dist;

            pdyna = &pull->dyna[c];
            pgrp  = &pull->group[pcrd->params.group[1]];

            wmass          = comm->dbuf_cyl[c*stride+0];
            wwmass         = comm->dbuf_cyl[c*stride+1];
            pdyna->mwscale = 1.0/wmass;
            /* Cylinder pulling can't be used with constraints, but we set
             * wscale and invtm anyhow, in case someone would like to use them.
             */
            pdyna->wscale  = wmass/wwmass;
            pdyna->invtm   = wwmass/(wmass*wmass);

            /* We store the deviation of the COM from the reference location
             * used above, since we need it when we apply the radial forces
             * to the atoms in the cylinder group.
             */
            pcrd->cyl_dev  = 0;
            for (m = 0; m < DIM; m++)
            {
                g_x[m]         = pgrp->x[m] - pcrd->vec[m]*pcrd->value_ref;
                dist           = -pcrd->vec[m]*comm->dbuf_cyl[c*stride+2]*pdyna->mwscale;
                pdyna->x[m]    = g_x[m] - dist;
                pcrd->cyl_dev += dist;
            }
            /* Now we know the exact COM of the cylinder reference group,
             * we can determine the radial force factor (ffrad) that when
             * multiplied with the axial pull force will give the radial
             * force on the pulled (non-cylinder) group.
             */
            for (m = 0; m < DIM; m++)
            {
                pcrd->ffrad[m] = (comm->dbuf_cyl[c*stride+6+m] +
                                  comm->dbuf_cyl[c*stride+3+m]*pcrd->cyl_dev)/wmass;
            }

            if (debug)
            {
                fprintf(debug, "Pull cylinder group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
                        c, pdyna->x[0], pdyna->x[1],
                        pdyna->x[2], 1.0/pdyna->invtm);
                fprintf(debug, "ffrad %8.3f %8.3f %8.3f\n",
                        pcrd->ffrad[XX], pcrd->ffrad[YY], pcrd->ffrad[ZZ]);
            }
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
                    struct pull_t *pull, t_mdatoms *md, t_pbc *pbc, double t,
                    rvec x[], rvec *xp)
{
    int          g;
    real         twopi_box = 0;
    pull_comm_t *comm;

    comm = &pull->comm;

    if (comm->rbuf == NULL)
    {
        snew(comm->rbuf, pull->ngroup);
    }
    if (comm->dbuf == NULL)
    {
        snew(comm->dbuf, 3*pull->ngroup);
    }

    if (pull->bRefAt && pull->bSetPBCatoms)
    {
        pull_set_pbcatoms(cr, pull, x, comm->rbuf);

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
        int m;

        assert(pull->npbcdim <= DIM);

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
        pull_group_work_t *pgrp;

        pgrp = &pull->group[g];

        if (pgrp->bCalcCOM)
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                dvec   com, comp;
                double wmass, wwmass;
                rvec   x_pbc = { 0, 0, 0 };
                int    i;

                clear_dvec(com);
                clear_dvec(comp);
                wmass  = 0;
                wwmass = 0;

                if (pgrp->epgrppbc == epgrppbcREFAT)
                {
                    /* Set the pbc atom */
                    copy_rvec(comm->rbuf[g], x_pbc);
                }

                for (i = 0; i < pgrp->nat_loc; i++)
                {
                    int  ii, m;
                    real mass, wm;

                    ii   = pgrp->ind_loc[i];
                    mass = md->massT[ii];
                    if (pgrp->weight_loc == NULL)
                    {
                        wm     = mass;
                        wmass += wm;
                    }
                    else
                    {
                        real w;

                        w       = pgrp->weight_loc[i];
                        wm      = w*mass;
                        wmass  += wm;
                        wwmass += wm*w;
                    }
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
                        rvec dx;

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

                /* We do this check after the loop above to avoid more nesting.
                 * If we have a single-atom group the mass is irrelevant, so
                 * we can remove the mass factor to avoid division by zero.
                 * Note that with constraint pulling the mass does matter, but
                 * in that case a check group mass != 0 has been done before.
                 */
                if (pgrp->params.nat == 1 && pgrp->nat_loc == 1 && wmass == 0)
                {
                    int m;

                    /* Copy the single atom coordinate */
                    for (m = 0; m < DIM; m++)
                    {
                        com[m] = x[pgrp->ind_loc[0]][m];
                    }
                    /* Set all mass factors to 1 to get the correct COM */
                    wmass  = 1;
                    wwmass = 1;
                }

                if (pgrp->weight_loc == NULL)
                {
                    wwmass = wmass;
                }

                /* Copy local sums to a buffer for global summing */
                copy_dvec(com,  comm->dbuf[g*3]);
                copy_dvec(comp, comm->dbuf[g*3 + 1]);
                comm->dbuf[g*3 + 2][0] = wmass;
                comm->dbuf[g*3 + 2][1] = wwmass;
                comm->dbuf[g*3 + 2][2] = 0;
            }
            else
            {
                /* Cosine weighting geometry */
                double cm, sm, cmp, smp, ccm, csm, ssm, csw, snw;
                int    i;

                cm  = 0;
                sm  = 0;
                cmp = 0;
                smp = 0;
                ccm = 0;
                csm = 0;
                ssm = 0;

                for (i = 0; i < pgrp->nat_loc; i++)
                {
                    int  ii;
                    real mass;

                    ii   = pgrp->ind_loc[i];
                    mass = md->massT[ii];
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

                /* Copy local sums to a buffer for global summing */
                comm->dbuf[g*3  ][0] = cm;
                comm->dbuf[g*3  ][1] = sm;
                comm->dbuf[g*3  ][2] = 0;
                comm->dbuf[g*3+1][0] = ccm;
                comm->dbuf[g*3+1][1] = csm;
                comm->dbuf[g*3+1][2] = ssm;
                comm->dbuf[g*3+2][0] = cmp;
                comm->dbuf[g*3+2][1] = smp;
                comm->dbuf[g*3+2][2] = 0;
            }
        }
    }

    pull_reduce_double(cr, comm, pull->ngroup*3*DIM, comm->dbuf[0]);

    for (g = 0; g < pull->ngroup; g++)
    {
        pull_group_work_t *pgrp;

        pgrp = &pull->group[g];
        if (pgrp->params.nat > 0 && pgrp->bCalcCOM)
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                double wmass, wwmass;
                int    m;

                /* Determine the inverse mass */
                wmass             = comm->dbuf[g*3+2][0];
                wwmass            = comm->dbuf[g*3+2][1];
                pgrp->mwscale     = 1.0/wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm != 0)
                {
                    pgrp->wscale  = wmass/wwmass;
                    pgrp->invtm   = wwmass/(wmass*wmass);
                }
                /* Divide by the total mass */
                for (m = 0; m < DIM; m++)
                {
                    pgrp->x[m]      = comm->dbuf[g*3  ][m]*pgrp->mwscale;
                    if (xp)
                    {
                        pgrp->xp[m] = comm->dbuf[g*3+1][m]*pgrp->mwscale;
                    }
                    if (pgrp->epgrppbc == epgrppbcREFAT)
                    {
                        pgrp->x[m]      += comm->rbuf[g][m];
                        if (xp)
                        {
                            pgrp->xp[m] += comm->rbuf[g][m];
                        }
                    }
                }
            }
            else
            {
                /* Cosine weighting geometry */
                double csw, snw, wmass, wwmass;
                int    i, ii;

                /* Determine the optimal location of the cosine weight */
                csw                   = comm->dbuf[g*3][0];
                snw                   = comm->dbuf[g*3][1];
                pgrp->x[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                /* Set the weights for the local atoms */
                wmass  = sqrt(csw*csw + snw*snw);
                wwmass = (comm->dbuf[g*3+1][0]*csw*csw +
                          comm->dbuf[g*3+1][1]*csw*snw +
                          comm->dbuf[g*3+1][2]*snw*snw)/(wmass*wmass);

                pgrp->mwscale = 1.0/wmass;
                pgrp->wscale  = wmass/wwmass;
                pgrp->invtm   = wwmass/(wmass*wmass);
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
                    csw                    = comm->dbuf[g*3+2][0];
                    snw                    = comm->dbuf[g*3+2][1];
                    pgrp->xp[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                }
            }
            if (debug)
            {
                fprintf(debug, "Pull group %d wmass %f invtm %f\n",
                        g, 1.0/pgrp->mwscale, pgrp->invtm);
            }
        }
    }

    if (pull->bCylinder)
    {
        /* Calculate the COMs for the cyclinder reference groups */
        make_cyl_refgrps(cr, pull, md, pbc, t, x);
    }
}
