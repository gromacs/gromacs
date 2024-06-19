/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/* \internal \file
 *
 * \brief Implements functions to collect state data to the main rank.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "collect.h"

#include "config.h"

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "atomdistribution.h"
#include "distribute.h"
#include "domdec_internal.h"

enum class FreeEnergyPerturbationCouplingType : int;

static void dd_collect_cg(gmx_domdec_t*            dd,
                          const int                ddpCount,
                          const int                ddpCountCgGl,
                          gmx::ArrayRef<const int> localCGNumbers)
{
    if (ddpCount == dd->comm->main_cg_ddp_count)
    {
        /* The main has the correct distribution */
        return;
    }

    gmx::ArrayRef<const int> atomGroups;
    int                      nat_home = 0;

    if (ddpCount == dd->ddp_count)
    {
        /* The local state and DD are in sync, use the DD indices */
        atomGroups = gmx::constArrayRefFromArray(dd->globalAtomIndices.data(), dd->numHomeAtoms);
        nat_home   = dd->comm->atomRanges.numHomeAtoms();
    }
    else if (ddpCountCgGl == ddpCount)
    {
        /* The DD is out of sync with the local state, but we have stored
         * the cg indices with the local state, so we can use those.
         */
        atomGroups = localCGNumbers;
        nat_home   = atomGroups.size();
    }
    else
    {
        gmx_incons(
                "Attempted to collect a vector for a state for which the charge group distribution "
                "is unknown");
    }

    AtomDistribution* ma = dd->ma.get();

    /* Collect the charge group and atom counts on the main */
    int localBuffer[2] = { static_cast<int>(atomGroups.size()), nat_home };
    dd_gather(dd, 2 * sizeof(int), localBuffer, DDMAIN(dd) ? ma->intBuffer.data() : nullptr);

    if (DDMAIN(dd))
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
            ma->intBuffer[rank]              = numGroups;
            ma->intBuffer[dd->nnodes + rank] = offset;
            offset += numGroups;
        }
    }

    /* Collect the charge group indices on the main */
    dd_gatherv(*dd,
               atomGroups,
               DDMAIN(dd) ? gmx::makeArrayRef(ma->intBuffer).subArray(0, dd->nnodes)
                          : gmx::ArrayRef<int>(),
               DDMAIN(dd) ? gmx::makeArrayRef(ma->intBuffer).subArray(dd->nnodes, dd->nnodes)
                          : gmx::ArrayRef<int>(),
               DDMAIN(dd) ? ma->atomGroups : gmx::ArrayRef<int>());

    dd->comm->main_cg_ddp_count = ddpCount;
}

static void dd_collect_vec_sendrecv(gmx_domdec_t*                  dd,
                                    gmx::ArrayRef<const gmx::RVec> lv,
                                    gmx::ArrayRef<gmx::RVec>       v)
{
    if (!DDMAIN(dd))
    {
#if GMX_MPI
        const int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
        MPI_Send(const_cast<void*>(static_cast<const void*>(lv.data())),
                 numHomeAtoms * sizeof(rvec),
                 MPI_BYTE,
                 dd->mainrank,
                 dd->rank,
                 dd->mpi_comm_all);
#endif
    }
    else
    {
        AtomDistribution& ma = *dd->ma;

        int rank      = dd->mainrank;
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
                MPI_Recv(ma.rvecBuffer.data(),
                         domainGroups.numAtoms * sizeof(rvec),
                         MPI_BYTE,
                         rank,
                         rank,
                         dd->mpi_comm_all,
                         MPI_STATUS_IGNORE);
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
    gmx::ArrayRef<const int> recvCounts;
    gmx::ArrayRef<const int> displacements;

    if (DDMAIN(dd))
    {
        get_commbuffer_counts(dd->ma.get(), &recvCounts, &displacements);
    }

    const int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
    dd_gatherv(*dd,
               lv.subArray(0, numHomeAtoms),
               recvCounts,
               displacements,
               DDMAIN(dd) ? dd->ma->rvecBuffer : gmx::ArrayRef<gmx::RVec>());

    if (DDMAIN(dd))
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
                    const int                      ddpCount,
                    const int                      ddpCountCgGl,
                    gmx::ArrayRef<const int>       localCGNumbers,
                    gmx::ArrayRef<const gmx::RVec> localVector,
                    gmx::ArrayRef<gmx::RVec>       globalVector)
{
    dd_collect_cg(dd, ddpCount, ddpCountCgGl, localCGNumbers);

    if (dd->nnodes <= c_maxNumRanksUseSendRecvForScatterAndGather)
    {
        dd_collect_vec_sendrecv(dd, localVector, globalVector);
    }
    else
    {
        dd_collect_vec_gatherv(dd, localVector, globalVector);
    }
}


void dd_collect_state(gmx_domdec_t* dd, const t_state* state_local, t_state* state)
{
    int nh = state_local->nhchainlength;

    if (DDMAIN(dd))
    {
        GMX_RELEASE_ASSERT(state->nhchainlength == nh,
                           "The global and local Nose-Hoover chain lengths should match");

        for (auto i : gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real>::keys())
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
    if (state_local->hasEntry(StateEntry::X))
    {
        auto globalXRef = state ? state->x : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd,
                       state_local->ddp_count,
                       state_local->ddp_count_cg_gl,
                       state_local->cg_gl,
                       state_local->x,
                       globalXRef);
    }
    if (state_local->hasEntry(StateEntry::V))
    {
        auto globalVRef = state ? state->v : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd,
                       state_local->ddp_count,
                       state_local->ddp_count_cg_gl,
                       state_local->cg_gl,
                       state_local->v,
                       globalVRef);
    }
    if (state_local->hasEntry(StateEntry::Cgp))
    {
        auto globalCgpRef = state ? state->cg_p : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(dd,
                       state_local->ddp_count,
                       state_local->ddp_count_cg_gl,
                       state_local->cg_gl,
                       state_local->cg_p,
                       globalCgpRef);
    }
}
