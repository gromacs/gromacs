/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#include "atomdistribution.h"
#include "distribute.h"
#include "domdec_internal.h"

static void dd_collect_cg(gmx_domdec_t* dd, const t_state* state_local)
{
    if (state_local->ddp_count == dd->comm->master_cg_ddp_count)
    {
        /* The master has the correct distribution */
        return;
    }

    gmx::ArrayRef<const int> atomGroups;
    int                      nat_home = 0;

    if (state_local->ddp_count == dd->ddp_count)
    {
        /* The local state and DD are in sync, use the DD indices */
        atomGroups = gmx::constArrayRefFromArray(dd->globalAtomGroupIndices.data(), dd->ncg_home);
        nat_home   = dd->comm->atomRanges.numHomeAtoms();
    }
    else if (state_local->ddp_count_cg_gl == state_local->ddp_count)
    {
        /* The DD is out of sync with the local state, but we have stored
         * the cg indices with the local state, so we can use those.
         */
        atomGroups = state_local->cg_gl;
        nat_home   = atomGroups.size();
    }
    else
    {
        gmx_incons(
                "Attempted to collect a vector for a state for which the charge group distribution "
                "is unknown");
    }

    AtomDistribution* ma = dd->ma.get();

    /* Collect the charge group and atom counts on the master */
    int localBuffer[2] = { static_cast<int>(atomGroups.size()), nat_home };
    dd_gather(dd, 2 * sizeof(int), localBuffer, DDMASTER(dd) ? ma->intBuffer.data() : nullptr);

    if (DDMASTER(dd))
    {
        int groupOffset = 0;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            auto& domainGroups = ma->domainGroups[rank];
            int   numGroups    = ma->intBuffer[2 * rank];

            domainGroups.atomGroups =
                    gmx::constArrayRefFromArray(ma->atomGroups.data() + groupOffset, numGroups);

            domainGroups.numAtoms = ma->intBuffer[2 * rank + 1];

            groupOffset += numGroups;
        }

        if (debug)
        {
            fprintf(debug, "Initial charge group distribution: ");
            for (int rank = 0; rank < dd->nnodes; rank++)
            {
                fprintf(debug, " %td", ma->domainGroups[rank].atomGroups.ssize());
            }
            fprintf(debug, "\n");
        }

        /* Make byte counts and indices */
        int offset = 0;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            int numGroups                    = ma->domainGroups[rank].atomGroups.size();
            ma->intBuffer[rank]              = numGroups * sizeof(int);
            ma->intBuffer[dd->nnodes + rank] = offset * sizeof(int);
            offset += numGroups;
        }
    }

    /* Collect the charge group indices on the master */
    dd_gatherv(dd, atomGroups.size() * sizeof(int), atomGroups.data(),
               DDMASTER(dd) ? ma->intBuffer.data() : nullptr,
               DDMASTER(dd) ? ma->intBuffer.data() + dd->nnodes : nullptr,
               DDMASTER(dd) ? ma->atomGroups.data() : nullptr);

    dd->comm->master_cg_ddp_count = state_local->ddp_count;
}

static void dd_collect_vec_sendrecv(gmx_domdec_t*                  dd,
                                    gmx::ArrayRef<const gmx::RVec> lv,
                                    gmx::ArrayRef<gmx::RVec>       v)
{
    if (!DDMASTER(dd))
    {
#if GMX_MPI
        const int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
        MPI_Send(const_cast<void*>(static_cast<const void*>(lv.data())),
                 numHomeAtoms * sizeof(rvec), MPI_BYTE, dd->masterrank, dd->rank, dd->mpi_comm_all);
#endif
    }
    else
    {
        AtomDistribution& ma = *dd->ma;

        int rank      = dd->masterrank;
        int localAtom = 0;
        for (const int& globalAtom : ma.domainGroups[rank].atomGroups)
        {
            copy_rvec(lv[localAtom++], v[globalAtom]);
        }

        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            if (rank != dd->rank)
            {
                const auto& domainGroups = ma.domainGroups[rank];

                GMX_RELEASE_ASSERT(v.data() != ma.rvecBuffer.data(),
                                   "We need different communication and return buffers");

                /* When we send/recv instead of scatter/gather, we might need
                 * to increase the communication buffer size here.
                 */
                if (static_cast<size_t>(domainGroups.numAtoms) > ma.rvecBuffer.size())
                {
                    ma.rvecBuffer.resize(domainGroups.numAtoms);
                }

#if GMX_MPI
                MPI_Recv(ma.rvecBuffer.data(), domainGroups.numAtoms * sizeof(rvec), MPI_BYTE, rank,
                         rank, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
                int localAtom = 0;
                for (const int& globalAtom : domainGroups.atomGroups)
                {
                    copy_rvec(ma.rvecBuffer[localAtom++], v[globalAtom]);
                }
            }
        }
    }
}

static void dd_collect_vec_gatherv(gmx_domdec_t*                  dd,
                                   gmx::ArrayRef<const gmx::RVec> lv,
                                   gmx::ArrayRef<gmx::RVec>       v)
{
    int* recvCounts    = nullptr;
    int* displacements = nullptr;

    if (DDMASTER(dd))
    {
        get_commbuffer_counts(dd->ma.get(), &recvCounts, &displacements);
    }

    const int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
    dd_gatherv(dd, numHomeAtoms * sizeof(rvec), lv.data(), recvCounts, displacements,
               DDMASTER(dd) ? dd->ma->rvecBuffer.data() : nullptr);

    if (DDMASTER(dd))
    {
        const AtomDistribution& ma = *dd->ma;

        int bufferAtom = 0;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            const auto& domainGroups = ma.domainGroups[rank];
            for (const int& globalAtom : domainGroups.atomGroups)
            {
                copy_rvec(ma.rvecBuffer[bufferAtom++], v[globalAtom]);
            }
        }
    }
}

void dd_collect_vec(gmx_domdec_t*                  dd,
                    const t_state*                 state_local,
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


void dd_collect_state(gmx_domdec_t* dd, const t_state* state_local, t_state* state)
{
    int nh = state_local->nhchainlength;

    if (DDMASTER(dd))
    {
        GMX_RELEASE_ASSERT(state->nhchainlength == nh,
                           "The global and local Nose-Hoover chain lengths should match");

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
                state->nosehoover_xi[i * nh + j]  = state_local->nosehoover_xi[i * nh + j];
                state->nosehoover_vxi[i * nh + j] = state_local->nosehoover_vxi[i * nh + j];
            }
            state->therm_integral[i] = state_local->therm_integral[i];
        }
        for (int i = 0; i < state_local->nnhpres; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state->nhpres_xi[i * nh + j]  = state_local->nhpres_xi[i * nh + j];
                state->nhpres_vxi[i * nh + j] = state_local->nhpres_vxi[i * nh + j];
            }
        }
        state->baros_integral     = state_local->baros_integral;
        state->pull_com_prev_step = state_local->pull_com_prev_step;
    }
    if (state_local->flags & (1 << estX))
    {
        auto globalXRef = state ? state->x : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd, state_local, state_local->x, globalXRef);
    }
    if (state_local->flags & (1 << estV))
    {
        auto globalVRef = state ? state->v : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd, state_local, state_local->v, globalVRef);
    }
    if (state_local->flags & (1 << estCGP))
    {
        auto globalCgpRef = state ? state->cg_p : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd, state_local, state_local->cg_p, globalCgpRef);
    }
}
