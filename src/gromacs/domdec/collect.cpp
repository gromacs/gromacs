/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/* \internal \file
 *
 * \brief Implements functions to collect state data to the master rank.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "collect.h"

#include "config.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "distribute.h"
#include "domdec_internal.h"

static void dd_collect_cg(gmx_domdec_t  *dd,
                          const t_state *state_local)
{
    gmx_domdec_master_t *ma = nullptr;
    int                  buf2[2], *ibuf, i, ncg_home = 0, nat_home = 0;

    if (state_local->ddp_count == dd->comm->master_cg_ddp_count)
    {
        /* The master has the correct distribution */
        return;
    }

    const int *cg;

    if (state_local->ddp_count == dd->ddp_count)
    {
        /* The local state and DD are in sync, use the DD indices */
        ncg_home = dd->ncg_home;
        cg       = dd->index_gl;
        nat_home = dd->nat_home;
    }
    else if (state_local->ddp_count_cg_gl == state_local->ddp_count)
    {
        /* The DD is out of sync with the local state, but we have stored
         * the cg indices with the local state, so we can use those.
         */
        t_block *cgs_gl;

        cgs_gl = &dd->comm->cgs_gl;

        ncg_home = state_local->cg_gl.size();
        cg       = state_local->cg_gl.data();
        nat_home = 0;
        for (i = 0; i < ncg_home; i++)
        {
            nat_home += cgs_gl->index[cg[i]+1] - cgs_gl->index[cg[i]];
        }
    }
    else
    {
        gmx_incons("Attempted to collect a vector for a state for which the charge group distribution is unknown");
    }

    buf2[0] = ncg_home;
    buf2[1] = nat_home;
    if (DDMASTER(dd))
    {
        ma   = dd->ma;
        ibuf = ma->ibuf;
    }
    else
    {
        ibuf = nullptr;
    }
    /* Collect the charge group and atom counts on the master */
    dd_gather(dd, 2*sizeof(int), buf2, ibuf);

    if (DDMASTER(dd))
    {
        ma->index[0] = 0;
        for (i = 0; i < dd->nnodes; i++)
        {
            ma->ncg[i]     = ma->ibuf[2*i];
            ma->nat[i]     = ma->ibuf[2*i+1];
            ma->index[i+1] = ma->index[i] + ma->ncg[i];

        }
        /* Make byte counts and indices */
        for (i = 0; i < dd->nnodes; i++)
        {
            ma->ibuf[i]            = ma->ncg[i]*sizeof(int);
            ma->ibuf[dd->nnodes+i] = ma->index[i]*sizeof(int);
        }
        if (debug)
        {
            fprintf(debug, "Initial charge group distribution: ");
            for (i = 0; i < dd->nnodes; i++)
            {
                fprintf(debug, " %d", ma->ncg[i]);
            }
            fprintf(debug, "\n");
        }
    }

    /* Collect the charge group indices on the master */
    dd_gatherv(dd,
               ncg_home*sizeof(int), cg,
               DDMASTER(dd) ? ma->ibuf : nullptr,
               DDMASTER(dd) ? ma->ibuf+dd->nnodes : nullptr,
               DDMASTER(dd) ? ma->cg : nullptr);

    dd->comm->master_cg_ddp_count = state_local->ddp_count;
}

static void dd_collect_vec_sendrecv(gmx_domdec_t                  *dd,
                                    gmx::ArrayRef<const gmx::RVec> lv,
                                    gmx::ArrayRef<gmx::RVec>       v)
{
    gmx_domdec_master_t *ma;
    int                  n, i, c, a, nalloc = 0;
    rvec                *buf = nullptr;
    t_block             *cgs_gl;

    ma = dd->ma;

    if (!DDMASTER(dd))
    {
#if GMX_MPI
        MPI_Send(const_cast<void *>(static_cast<const void *>(lv.data())), dd->nat_home*sizeof(rvec), MPI_BYTE,
                 dd->masterrank, dd->rank, dd->mpi_comm_all);
#endif
    }
    else
    {
        /* Copy the master coordinates to the global array */
        cgs_gl = &dd->comm->cgs_gl;

        n = dd->masterrank;
        a = 0;
        for (i = ma->index[n]; i < ma->index[n+1]; i++)
        {
            for (c = cgs_gl->index[ma->cg[i]]; c < cgs_gl->index[ma->cg[i]+1]; c++)
            {
                copy_rvec(lv[a++], v[c]);
            }
        }

        for (n = 0; n < dd->nnodes; n++)
        {
            if (n != dd->rank)
            {
                if (ma->nat[n] > nalloc)
                {
                    nalloc = over_alloc_dd(ma->nat[n]);
                    srenew(buf, nalloc);
                }
#if GMX_MPI
                MPI_Recv(buf, ma->nat[n]*sizeof(rvec), MPI_BYTE, n,
                         n, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
                a = 0;
                for (i = ma->index[n]; i < ma->index[n+1]; i++)
                {
                    for (c = cgs_gl->index[ma->cg[i]]; c < cgs_gl->index[ma->cg[i]+1]; c++)
                    {
                        copy_rvec(buf[a++], v[c]);
                    }
                }
            }
        }
        sfree(buf);
    }
}

static void dd_collect_vec_gatherv(gmx_domdec_t                  *dd,
                                   gmx::ArrayRef<const gmx::RVec> lv,
                                   gmx::ArrayRef<gmx::RVec>       v)
{
    gmx_domdec_master_t *ma;
    int                 *rcounts = nullptr, *disps = nullptr;
    int                  n, i, c, a;
    rvec                *buf = nullptr;
    t_block             *cgs_gl;

    ma = dd->ma;

    if (DDMASTER(dd))
    {
        get_commbuffer_counts(dd, &rcounts, &disps);

        buf = ma->vbuf;
    }

    dd_gatherv(dd, dd->nat_home*sizeof(rvec), lv.data(), rcounts, disps, buf);

    if (DDMASTER(dd))
    {
        cgs_gl = &dd->comm->cgs_gl;

        a = 0;
        for (n = 0; n < dd->nnodes; n++)
        {
            for (i = ma->index[n]; i < ma->index[n+1]; i++)
            {
                for (c = cgs_gl->index[ma->cg[i]]; c < cgs_gl->index[ma->cg[i]+1]; c++)
                {
                    copy_rvec(buf[a++], v[c]);
                }
            }
        }
    }
}

void dd_collect_vec(gmx_domdec_t                  *dd,
                    const t_state                 *state_local,
                    gmx::ArrayRef<const gmx::RVec> lv,
                    gmx::ArrayRef<gmx::RVec>       v)
{
    dd_collect_cg(dd, state_local);

    if (dd->nnodes <= c_maxNumRanksUseSendRecvForScatterAndGather)
    {
        dd_collect_vec_sendrecv(dd, lv, v);
    }
    else
    {
        dd_collect_vec_gatherv(dd, lv, v);
    }
}


void dd_collect_state(gmx_domdec_t *dd,
                      const t_state *state_local, t_state *state)
{
    int nh = state_local->nhchainlength;

    if (DDMASTER(dd))
    {
        GMX_RELEASE_ASSERT(state->nhchainlength == nh, "The global and local Nose-Hoover chain lengths should match");

        for (int i = 0; i < efptNR; i++)
        {
            state->lambda[i] = state_local->lambda[i];
        }
        state->fep_state = state_local->fep_state;
        state->veta      = state_local->veta;
        state->vol0      = state_local->vol0;
        copy_mat(state_local->box, state->box);
        copy_mat(state_local->boxv, state->boxv);
        copy_mat(state_local->svir_prev, state->svir_prev);
        copy_mat(state_local->fvir_prev, state->fvir_prev);
        copy_mat(state_local->pres_prev, state->pres_prev);

        for (int i = 0; i < state_local->ngtc; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state->nosehoover_xi[i*nh+j]        = state_local->nosehoover_xi[i*nh+j];
                state->nosehoover_vxi[i*nh+j]       = state_local->nosehoover_vxi[i*nh+j];
            }
            state->therm_integral[i] = state_local->therm_integral[i];
        }
        for (int i = 0; i < state_local->nnhpres; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state->nhpres_xi[i*nh+j]        = state_local->nhpres_xi[i*nh+j];
                state->nhpres_vxi[i*nh+j]       = state_local->nhpres_vxi[i*nh+j];
            }
        }
        state->baros_integral = state_local->baros_integral;
    }
    if (state_local->flags & (1 << estX))
    {
        gmx::ArrayRef<gmx::RVec> globalXRef = state ? gmx::makeArrayRef(state->x) : gmx::EmptyArrayRef();
        dd_collect_vec(dd, state_local, state_local->x, globalXRef);
    }
    if (state_local->flags & (1 << estV))
    {
        gmx::ArrayRef<gmx::RVec> globalVRef = state ? gmx::makeArrayRef(state->v) : gmx::EmptyArrayRef();
        dd_collect_vec(dd, state_local, state_local->v, globalVRef);
    }
    if (state_local->flags & (1 << estCGP))
    {
        gmx::ArrayRef<gmx::RVec> globalCgpRef = state ? gmx::makeArrayRef(state->cg_p) : gmx::EmptyArrayRef();
        dd_collect_vec(dd, state_local, state_local->cg_p, globalCgpRef);
    }
}
