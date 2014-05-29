/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

#include "state.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/math/vec.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

static void exchange_reals(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, real *v, int n)
{
    real *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}


static void exchange_ints(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, int *v, int n)
{
    int *buf;
    int  i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
             buf,n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
             ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_doubles(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, double *v, int n)
{
    double *buf;
    int     i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_rvecs(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, rvec *v, int n)
{
    rvec *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            copy_rvec(buf[i], v[i]);
        }
        sfree(buf);
    }
}

void exchange_state(const gmx_multisim_t *ms, int b, t_state *state)
{
    /* When t_state changes, this code should be updated. */
    int ngtc, nnhpres;
    ngtc    = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    exchange_rvecs(ms, b, state->box, DIM);
    exchange_rvecs(ms, b, state->box_rel, DIM);
    exchange_rvecs(ms, b, state->boxv, DIM);
    exchange_reals(ms, b, &(state->veta), 1);
    exchange_reals(ms, b, &(state->vol0), 1);
    exchange_rvecs(ms, b, state->svir_prev, DIM);
    exchange_rvecs(ms, b, state->fvir_prev, DIM);
    exchange_rvecs(ms, b, state->pres_prev, DIM);
    exchange_doubles(ms, b, state->nosehoover_xi, ngtc);
    exchange_doubles(ms, b, state->nosehoover_vxi, ngtc);
    exchange_doubles(ms, b, state->nhpres_xi, nnhpres);
    exchange_doubles(ms, b, state->nhpres_vxi, nnhpres);
    exchange_doubles(ms, b, state->therm_integral, state->ngtc);
    exchange_rvecs(ms, b, state->x, state->natoms);
    exchange_rvecs(ms, b, state->v, state->natoms);
    exchange_rvecs(ms, b, state->sd_X, state->natoms);
}

static void copy_rvecs(const rvec *source, rvec *sink, int n)
{
    int i;

    if (sink != NULL)
    {
        for (i = 0; i < n; i++)
        {
            copy_rvec(source[i], sink[i]);
        }
    }
}

static void copy_doubles(const double *source, double *sink, int n)
{
    int i;

    if (sink != NULL)
    {
        for (i = 0; i < n; i++)
        {
            sink[i] = source[i];
        }
    }
}

static void copy_reals(const real *source, real *sink, int n)
{
    int i;

    if (sink != NULL)
    {
        for (i = 0; i < n; i++)
        {
            sink[i] = source[i];
        }
    }
}

static void copy_ints(const int *source, int *sink, int n)
{
    int i;

    if (sink != NULL)
    {
        for (i = 0; i < n; i++)
        {
            sink[i] = source[i];
        }
    }
}

void copy_state_nonatomdata(const t_state *source, t_state *sink)
{
    /* When t_state changes, this code should be updated. */
    int ngtc, nnhpres;
    ngtc    = source->ngtc * source->nhchainlength;
    nnhpres = source->nnhpres* source->nhchainlength;
    copy_rvecs(source->box, sink->box, DIM);
    copy_rvecs(source->box_rel, sink->box_rel, DIM);
    copy_rvecs(source->boxv, sink->boxv, DIM);
    sink->veta = source->veta;
    sink->vol0 = source->vol0;
    copy_rvecs(source->svir_prev, sink->svir_prev, DIM);
    copy_rvecs(source->fvir_prev, sink->fvir_prev, DIM);
    copy_rvecs(source->pres_prev, sink->pres_prev, DIM);
    copy_doubles(source->nosehoover_xi, sink->nosehoover_xi, ngtc);
    copy_doubles(source->nosehoover_vxi, sink->nosehoover_vxi, ngtc);
    copy_doubles(source->nhpres_xi, sink->nhpres_xi, nnhpres);
    copy_doubles(source->nhpres_vxi, sink->nhpres_vxi, nnhpres);
    copy_doubles(source->therm_integral, sink->therm_integral, source->ngtc);
    copy_rvecs(source->x, sink->x, source->natoms);
    copy_rvecs(source->v, sink->v, source->natoms);
    copy_rvecs(source->sd_X, sink->sd_X, source->natoms);
    copy_ints(&(source->fep_state), &(sink->fep_state), 1);
    copy_reals(source->lambda, sink->lambda, efptNR);
}
