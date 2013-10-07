/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <math.h>
#include "macros.h"
#include "main.h"
#include "smalloc.h"
#include "futil.h"
#include "tgroup.h"
#include "vec.h"
#include "network.h"
#include "smalloc.h"
#include "update.h"
#include "rbin.h"
#include "mtop_util.h"
#include "gmx_omp_nthreads.h"

static void init_temperature_coupling_group_outputs(int ngtc, gmx_temperature_coupling_group_outputs_t tc[])
{
    int i, j;

    for (i = 0; (i < ngtc); i++)
    {
        tc[i].T = 0;
        clear_mat(tc[i].ekinh);
        clear_mat(tc[i].ekinh_old);
        clear_mat(tc[i].ekinf);
    }
}

void init_temperature_coupling_outputs(const t_grpopts *opts,
                                       gmx_temperature_coupling_outputs_t *temperature_coupling_outputs)
{
    int i;

    if (debug)
    {
        fprintf(debug, "ngtc: %d\n", opts->ngtc);
    }

    temperature_coupling_outputs->ngroups = opts->ngtc;
    snew(temperature_coupling_outputs->group_data, opts->ngtc);
    init_temperature_coupling_group_outputs(opts->ngtc, temperature_coupling_outputs->group_data);
    /* Set Berendsen tcoupl lambda's to 1,
     * so runs without Berendsen coupling are not affected.
     */
    for (i = 0; i < opts->ngtc; i++)
    {
        temperature_coupling_outputs->group_data[i].lambda         = 1.0;
        temperature_coupling_outputs->group_data[i].vscale_nhc     = 1.0;
        temperature_coupling_outputs->group_data[i].ekinscaleh_nhc = 1.0;
        temperature_coupling_outputs->group_data[i].ekinscalef_nhc = 1.0;
    }

}

static void init_constant_acceleration_groups(const gmx_mtop_t *mtop, int ngacc, t_grp_acc *group_data)
{
    const gmx_groups_t     *groups;
    gmx_mtop_atomloop_all_t aloop;
    int                     i, grp;
    t_atom                 *atom;

    if (ngacc > 0)
    {
        groups = &mtop->groups;
        aloop  = gmx_mtop_atomloop_all_init(mtop);
        while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
        {
            grp = ggrpnr(groups, egcACC, i);
            if ((grp < 0) && (grp >= ngacc))
            {
                gmx_incons("Input for acceleration groups wrong");
            }
            group_data[grp].nat++;
            /* This will not work for integrator BD */
            group_data[grp].mA += atom->m;
            group_data[grp].mB += atom->mB;
        }
    }
}

void
init_constant_acceleration(const gmx_mtop_t *mtop,
                           const t_grpopts *opts,
                           gmx_constant_acceleration_t *const_acc)
{
    if (debug)
    {
        fprintf(debug, "Number of acceleration groups: %d\n", opts->ngacc);
    }

    /* bDoAcceleration tells if we should remove remove the COM velocity
     * from the velocities during velocity scaling in T-coupling.
     * Turn this on when we have multiple acceleration groups
     * or one accelerated group.
     *
     * TODO Why should this be on when there are multiple acceleration
     * groups and perhaps all of their norms are zero?
     */
    const_acc->bDoAcceleration = (opts->ngacc > 1 || norm(opts->acc[0]) > 0);

    const_acc->ngroups = opts->ngacc;
    snew(const_acc->group_data, opts->ngacc);
    init_constant_acceleration_groups(mtop, opts->ngacc, const_acc->group_data);
}

void init_ekindata(const t_grpopts *opts,
                   gmx_ekindata_t  *ekind)
{
    int i;
    int nthread, thread;

    nthread = gmx_omp_nthreads_get(emntUpdate);

    snew(ekind->ekin_work_alloc, nthread);
    snew(ekind->ekin_work, nthread);
    snew(ekind->dekindl_work, nthread);
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
#define EKIN_WORK_BUFFER_SIZE 2
        /* Allocate 2 extra elements on both sides, so in single
         * precision we have
         * EKIN_WORK_BUFFER_SIZE*DIM*DIM*sizeof(real) = 72/144 bytes
         * buffer on both sides to avoid cache pollution.
         */
        snew(ekind->ekin_work_alloc[thread], opts->ngtc+2*EKIN_WORK_BUFFER_SIZE);
        ekind->ekin_work[thread] = ekind->ekin_work_alloc[thread] + EKIN_WORK_BUFFER_SIZE;
        /* Nasty hack so we can have the per-thread accumulation
         * variable for dekindl in the same thread-local cache lines
         * as the per-thread accumulation tensors for ekin[fh],
         * because they are accumulated in the same loop. */
        ekind->dekindl_work[thread] = &(ekind->ekin_work[thread][opts->ngtc][0][0]);
#undef EKIN_WORK_BUFFER_SIZE
    }
}

void accumulate_u(t_commrec *cr,
                  gmx_constant_acceleration_t *constant_acceleration)
{
    t_bin *rb;
    int    g;

    if (!PAR(cr) || !constant_acceleration->bDoAcceleration)
    {
        return;
    }

    rb = mk_bin();

    for (g = 0; (g < constant_acceleration->ngroups); g++)
    {
        add_binr(rb, DIM, constant_acceleration->group_data[g].u);
    }
    sum_bin(rb, cr);

    for (g = 0; (g < constant_acceleration->ngroups); g++)
    {
        extract_binr(rb, DIM*g, DIM, constant_acceleration->group_data[g].u);
    }
    destroy_bin(rb);
}

// TODO find out if this is even called anywhere!!!
void update_ekindata(int start, int homenr,
                     gmx_temperature_coupling_outputs_t *temperature_coupling_outputs,
                     gmx_constant_acceleration_t *constant_acceleration,
                     t_grpopts *opts, rvec v[], t_mdatoms *md, real lambda)
{
    int  d, g, n;
    real mv;

    /* calculate mean velocities at whole timestep */
    for (g = 0; (g < temperature_coupling_outputs->ngroups); g++)
    {
        temperature_coupling_outputs->group_data[g].T = 0;
    }

    if (constant_acceleration->bDoAcceleration)
    {
        for (g = 0; (g < constant_acceleration->ngroups); g++)
        {
            clear_rvec(constant_acceleration->group_data[g].u);
        }

        g = 0;
        for (n = start; (n < start+homenr); n++)
        {
            if (md->cACC)
            {
                g = md->cACC[n];
            }
            for (d = 0; (d < DIM); d++)
            {
                mv                                         = md->massT[n]*v[n][d];
                constant_acceleration->group_data[g].u[d] += mv;
            }
        }

        for (g = 0; (g < constant_acceleration->ngroups); g++)
        {
            for (d = 0; (d < DIM); d++)
            {
                constant_acceleration->group_data[g].u[d] /=
                    (1-lambda)*constant_acceleration->group_data[g].mA +
                    lambda*constant_acceleration->group_data[g].mB;
            }
        }
    }
}

real sum_ekin(t_grpopts *opts,
              gmx_ekindata_t *ekind,
              gmx_temperature_coupling_outputs_t *temperature_coupling_outputs,
              real *dekindlambda,
              gmx_bool bEkinAveVel, gmx_bool bScaleEkin)
{
    int           i, j, m, ngtc;
    real          T, ek;
    gmx_temperature_coupling_group_outputs_t *tcstat;
    real          nrdf, nd, *ndf;

    ngtc = opts->ngtc;
    ndf  = opts->nrdf;

    T    = 0;
    nrdf = 0;

    clear_mat(ekind->ekin);

    for (i = 0; (i < ngtc); i++)
    {

        nd     = ndf[i];
        tcstat = &temperature_coupling_outputs->group_data[i];
        /* Sometimes a group does not have degrees of freedom, e.g.
         * when it consists of shells and virtual sites, then we just
         * set the temperatue to 0 and also neglect the kinetic
         * energy, which should be  zero anyway.
         */

        if (nd > 0)
        {
            if (bEkinAveVel)
            {
                if (!bScaleEkin)
                {
                    /* in this case, kinetic energy is from the current velocities already */
                    msmul(tcstat->ekinf, tcstat->ekinscalef_nhc, tcstat->ekinf);
                }
            }
            else
            {
                /* Calculate the full step Ekin as the average of the half steps */
                for (j = 0; (j < DIM); j++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        tcstat->ekinf[j][m] =
                            0.5*(tcstat->ekinh[j][m]*tcstat->ekinscaleh_nhc + tcstat->ekinh_old[j][m]);
                    }
                }
            }
            m_add(tcstat->ekinf, ekind->ekin, ekind->ekin);

            tcstat->Th = calc_temp(trace(tcstat->ekinh), nd);
            tcstat->T  = calc_temp(trace(tcstat->ekinf), nd);

            /* after the scaling factors have been multiplied in, we can remove them */
            if (bEkinAveVel)
            {
                tcstat->ekinscalef_nhc = 1.0;
            }
            else
            {
                tcstat->ekinscaleh_nhc = 1.0;
            }
        }
        else
        {
            tcstat->T  = 0;
            tcstat->Th = 0;
        }
        T    += nd*tcstat->T;
        nrdf += nd;
    }
    if (nrdf > 0)
    {
        T /= nrdf;
    }
    if (dekindlambda)
    {
        if (bEkinAveVel)
        {
            *dekindlambda = ekind->dekindl;
        }
        else
        {
            *dekindlambda = 0.5*(ekind->dekindl + ekind->dekindl_old);
        }
    }
    return T;
}
