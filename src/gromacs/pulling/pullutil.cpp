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

#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static void pull_reduce_real(t_commrec   *cr,
                             pull_comm_t *comm,
                             int          n,
                             real        *data)
{
    if (cr != nullptr && PAR(cr))
    {
        if (comm->bParticipateAll)
        {
            /* Sum the contributions over all DD ranks */
            gmx_sum(n, data, cr);
        }
        else
        {
#if GMX_MPI
#if MPI_IN_PLACE_EXISTS
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
    if (cr != nullptr && PAR(cr))
    {
        if (comm->bParticipateAll)
        {
            /* Sum the contributions over all DD ranks */
            gmx_sumd(n, data, cr);
        }
        else
        {
#if GMX_MPI
#if MPI_IN_PLACE_EXISTS
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

    if (cr != nullptr && DOMAINDECOMP(cr))
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
    const int       stride = 9;
    int             c, i, ii, m, start, end;
    rvec            g_x, dx, dir;
    double          inv_cyl_r2;
    pull_comm_t    *comm;
    gmx_ga2la_t    *ga2la = nullptr;

    comm = &pull->comm;

    if (comm->dbuf_cyl == nullptr)
    {
        snew(comm->dbuf_cyl, pull->ncoord*stride);
    }

    if (cr && DOMAINDECOMP(cr))
    {
        ga2la = cr->dd->ga2la;
    }

    start = 0;
    end   = md->homenr;

    inv_cyl_r2 = 1.0/gmx::square(pull->params.cylinder_r);

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
            copy_dvec_to_rvec(pcrd->vec, dir);
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

    if (cr != nullptr && PAR(cr))
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

static void sum_com_part(const pull_group_work_t *pgrp,
                         int ind_start, int ind_end,
                         const rvec *x, const rvec *xp,
                         const real *mass,
                         const t_pbc *pbc,
                         const rvec x_pbc,
                         pull_sum_com_t *sum_com)
{
    double sum_wm   = 0;
    double sum_wwm  = 0;
    dvec   sum_wmx  = { 0, 0, 0 };
    dvec   sum_wmxp = { 0, 0, 0 };

    for (int i = ind_start; i < ind_end; i++)
    {
        int  ii = pgrp->ind_loc[i];
        real wm;
        if (pgrp->weight_loc == nullptr)
        {
            wm      = mass[ii];
            sum_wm += wm;
        }
        else
        {
            real w;

            w        = pgrp->weight_loc[i];
            wm       = w*mass[ii];
            sum_wm  += wm;
            sum_wwm += wm*w;
        }
        if (pgrp->epgrppbc == epgrppbcNONE)
        {
            /* Plain COM: sum the coordinates */
            for (int d = 0; d < DIM; d++)
            {
                sum_wmx[d]      += wm*x[ii][d];
            }
            if (xp)
            {
                for (int d = 0; d < DIM; d++)
                {
                    sum_wmxp[d] += wm*xp[ii][d];
                }
            }
        }
        else
        {
            rvec dx;

            /* Sum the difference with the reference atom */
            pbc_dx(pbc, x[ii], x_pbc, dx);
            for (int d = 0; d < DIM; d++)
            {
                sum_wmx[d]     += wm*dx[d];
            }
            if (xp)
            {
                /* For xp add the difference between xp and x to dx,
                 * such that we use the same periodic image,
                 * also when xp has a large displacement.
                 */
                for (int d = 0; d < DIM; d++)
                {
                    sum_wmxp[d] += wm*(dx[d] + xp[ii][d] - x[ii][d]);
                }
            }
        }
    }

    sum_com->sum_wm  = sum_wm;
    sum_com->sum_wwm = sum_wwm;
    copy_dvec(sum_wmx, sum_com->sum_wmx);
    if (xp)
    {
        copy_dvec(sum_wmxp, sum_com->sum_wmxp);
    }
}

static void sum_com_part_cosweight(const pull_group_work_t *pgrp,
                                   int ind_start, int ind_end,
                                   int cosdim, real twopi_box,
                                   const rvec *x, const rvec *xp,
                                   const real *mass,
                                   pull_sum_com_t *sum_com)
{
    /* Cosine weighting geometry */
    double sum_cm  = 0;
    double sum_sm  = 0;
    double sum_ccm = 0;
    double sum_csm = 0;
    double sum_ssm = 0;
    double sum_cmp = 0;
    double sum_smp = 0;

    for (int i = ind_start; i < ind_end; i++)
    {
        int  ii  = pgrp->ind_loc[i];
        real m   = mass[ii];
        /* Determine cos and sin sums */
        real cw  = std::cos(x[ii][cosdim]*twopi_box);
        real sw  = std::sin(x[ii][cosdim]*twopi_box);
        sum_cm  += static_cast<double>(cw*m);
        sum_sm  += static_cast<double>(sw*m);
        sum_ccm += static_cast<double>(cw*cw*m);
        sum_csm += static_cast<double>(cw*sw*m);
        sum_ssm += static_cast<double>(sw*sw*m);

        if (xp != nullptr)
        {
            real cw  = std::cos(xp[ii][cosdim]*twopi_box);
            real sw  = std::sin(xp[ii][cosdim]*twopi_box);
            sum_cmp += static_cast<double>(cw*m);
            sum_smp += static_cast<double>(sw*m);
        }
    }

    sum_com->sum_cm  = sum_cm;
    sum_com->sum_sm  = sum_sm;
    sum_com->sum_ccm = sum_ccm;
    sum_com->sum_csm = sum_csm;
    sum_com->sum_ssm = sum_ssm;
    sum_com->sum_cmp = sum_cmp;
    sum_com->sum_smp = sum_smp;
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

    if (comm->rbuf == nullptr)
    {
        snew(comm->rbuf, pull->ngroup);
    }
    if (comm->dbuf == nullptr)
    {
        snew(comm->dbuf, 3*pull->ngroup);
    }

    if (pull->bRefAt && pull->bSetPBCatoms)
    {
        pull_set_pbcatoms(cr, pull, x, comm->rbuf);

        if (cr != nullptr && DOMAINDECOMP(cr))
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
                rvec   x_pbc = { 0, 0, 0 };

                if (pgrp->epgrppbc == epgrppbcREFAT)
                {
                    /* Set the pbc atom */
                    copy_rvec(comm->rbuf[g], x_pbc);
                }

                /* The final sums should end up in sum_com[0] */
                pull_sum_com_t *sum_com = &pull->sum_com[0];

                /* If we have a single-atom group the mass is irrelevant, so
                 * we can remove the mass factor to avoid division by zero.
                 * Note that with constraint pulling the mass does matter, but
                 * in that case a check group mass != 0 has been done before.
                 */
                if (pgrp->params.nat == 1 &&
                    pgrp->nat_loc == 1 &&
                    md->massT[pgrp->ind_loc[0]] == 0)
                {
                    GMX_ASSERT(xp == NULL, "We should not have groups with zero mass with constraints, i.e. xp!=NULL");

                    /* Copy the single atom coordinate */
                    for (int d = 0; d < DIM; d++)
                    {
                        sum_com->sum_wmx[d] = x[pgrp->ind_loc[0]][d];
                    }
                    /* Set all mass factors to 1 to get the correct COM */
                    sum_com->sum_wm  = 1;
                    sum_com->sum_wwm = 1;
                }
                else if (pgrp->nat_loc <= c_pullMaxNumLocalAtomsSingleThreaded)
                {
                    sum_com_part(pgrp, 0, pgrp->nat_loc,
                                 x, xp, md->massT,
                                 pbc, x_pbc,
                                 sum_com);
                }
                else
                {
#pragma omp parallel for num_threads(pull->nthreads) schedule(static)
                    for (int t = 0; t < pull->nthreads; t++)
                    {
                        int ind_start = (pgrp->nat_loc*(t + 0))/pull->nthreads;
                        int ind_end   = (pgrp->nat_loc*(t + 1))/pull->nthreads;
                        sum_com_part(pgrp, ind_start, ind_end,
                                     x, xp, md->massT,
                                     pbc, x_pbc,
                                     &pull->sum_com[t]);
                    }

                    /* Reduce the thread contributions to sum_com[0] */
                    for (int t = 1; t < pull->nthreads; t++)
                    {
                        sum_com->sum_wm  += pull->sum_com[t].sum_wm;
                        sum_com->sum_wwm += pull->sum_com[t].sum_wwm;
                        dvec_inc(sum_com->sum_wmx, pull->sum_com[t].sum_wmx);
                        dvec_inc(sum_com->sum_wmxp, pull->sum_com[t].sum_wmxp);
                    }
                }

                if (pgrp->weight_loc == nullptr)
                {
                    sum_com->sum_wwm = sum_com->sum_wm;
                }

                /* Copy local sums to a buffer for global summing */
                copy_dvec(sum_com->sum_wmx,  comm->dbuf[g*3]);
                copy_dvec(sum_com->sum_wmxp, comm->dbuf[g*3 + 1]);
                comm->dbuf[g*3 + 2][0] = sum_com->sum_wm;
                comm->dbuf[g*3 + 2][1] = sum_com->sum_wwm;
                comm->dbuf[g*3 + 2][2] = 0;
            }
            else
            {
                /* Cosine weighting geometry.
                 * This uses a slab of the system, thus we always have many
                 * atoms in the pull groups. Therefore, always use threads.
                 */
#pragma omp parallel for num_threads(pull->nthreads) schedule(static)
                for (int t = 0; t < pull->nthreads; t++)
                {
                    int ind_start = (pgrp->nat_loc*(t + 0))/pull->nthreads;
                    int ind_end   = (pgrp->nat_loc*(t + 1))/pull->nthreads;
                    sum_com_part_cosweight(pgrp, ind_start, ind_end,
                                           pull->cosdim, twopi_box,
                                           x, xp, md->massT,
                                           &pull->sum_com[t]);
                }

                /* Reduce the thread contributions to sum_com[0] */
                pull_sum_com_t *sum_com = &pull->sum_com[0];
                for (int t = 1; t < pull->nthreads; t++)
                {
                    sum_com->sum_cm  += pull->sum_com[t].sum_cm;
                    sum_com->sum_sm  += pull->sum_com[t].sum_sm;
                    sum_com->sum_ccm += pull->sum_com[t].sum_ccm;
                    sum_com->sum_csm += pull->sum_com[t].sum_csm;
                    sum_com->sum_ssm += pull->sum_com[t].sum_ssm;
                    sum_com->sum_cmp += pull->sum_com[t].sum_cmp;
                    sum_com->sum_smp += pull->sum_com[t].sum_smp;
                }

                /* Copy local sums to a buffer for global summing */
                comm->dbuf[g*3    ][0] = sum_com->sum_cm;
                comm->dbuf[g*3    ][1] = sum_com->sum_sm;
                comm->dbuf[g*3    ][2] = 0;
                comm->dbuf[g*3 + 1][0] = sum_com->sum_ccm;
                comm->dbuf[g*3 + 1][1] = sum_com->sum_csm;
                comm->dbuf[g*3 + 1][2] = sum_com->sum_ssm;
                comm->dbuf[g*3 + 2][0] = sum_com->sum_cmp;
                comm->dbuf[g*3 + 2][1] = sum_com->sum_smp;
                comm->dbuf[g*3 + 2][2] = 0;
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
