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
 * \brief Implements atom distribution functions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "distribute.h"

#include "config.h"

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "atomdistribution.h"
#include "cellsizes.h"
#include "domdec_internal.h"
#include "utility.h"

static void distributeVecSendrecv(gmx_domdec_t*                  dd,
                                  gmx::ArrayRef<const gmx::RVec> globalVec,
                                  gmx::ArrayRef<gmx::RVec>       localVec)
{
    if (DDMAIN(dd))
    {
        std::vector<gmx::RVec> buffer;

        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            if (rank != dd->rank)
            {
                const auto& domainGroups = dd->ma->domainGroups[rank];

                buffer.resize(domainGroups.numAtoms);

                int localAtom = 0;
                for (const int& globalAtom : domainGroups.atomGroups)
                {
                    buffer[localAtom++] = globalVec[globalAtom];
                }
                GMX_RELEASE_ASSERT(localAtom == domainGroups.numAtoms,
                                   "The index count and number of indices should match");

#if GMX_MPI
                MPI_Send(buffer.data(), domainGroups.numAtoms * sizeof(gmx::RVec), MPI_BYTE, rank, rank, dd->mpi_comm_all);
#endif
            }
        }

        const auto& domainGroups = dd->ma->domainGroups[dd->mainrank];
        int         localAtom    = 0;
        for (const int& globalAtom : domainGroups.atomGroups)
        {
            localVec[localAtom++] = globalVec[globalAtom];
        }
    }
    else
    {
#if GMX_MPI
        int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
        MPI_Recv(localVec.data(),
                 numHomeAtoms * sizeof(gmx::RVec),
                 MPI_BYTE,
                 dd->mainrank,
                 MPI_ANY_TAG,
                 dd->mpi_comm_all,
                 MPI_STATUS_IGNORE);
#endif
    }
}

static void distributeVecScatterv(gmx_domdec_t*                  dd,
                                  gmx::ArrayRef<const gmx::RVec> globalVec,
                                  gmx::ArrayRef<gmx::RVec>       localVec)
{
    gmx::ArrayRef<const int> sendCounts;
    gmx::ArrayRef<const int> displacements;

    if (DDMAIN(dd))
    {
        AtomDistribution& ma = *dd->ma;

        get_commbuffer_counts(&ma, &sendCounts, &displacements);

        gmx::ArrayRef<gmx::RVec> buffer    = ma.rvecBuffer;
        int                      localAtom = 0;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            const auto& domainGroups = ma.domainGroups[rank];
            for (const int& globalAtom : domainGroups.atomGroups)
            {
                buffer[localAtom++] = globalVec[globalAtom];
            }
        }
    }

    int numHomeAtoms = dd->comm->atomRanges.numHomeAtoms();
    dd_scatterv(dd,
                sendCounts,
                displacements,
                DDMAIN(dd) ? dd->ma->rvecBuffer.data() : nullptr,
                numHomeAtoms,
                localVec.data());
}

static void distributeVec(gmx_domdec_t*                  dd,
                          gmx::ArrayRef<const gmx::RVec> globalVec,
                          gmx::ArrayRef<gmx::RVec>       localVec)
{
    if (dd->nnodes <= c_maxNumRanksUseSendRecvForScatterAndGather)
    {
        distributeVecSendrecv(dd, globalVec, localVec);
    }
    else
    {
        distributeVecScatterv(dd, globalVec, localVec);
    }
}

void dd_distribute_dfhist(gmx_domdec_t* dd, df_history_t* dfhist)
{
    if (dfhist == nullptr)
    {
        return;
    }

    dd_bcast(dd, sizeof(int), &dfhist->bEquil);
    dd_bcast(dd, sizeof(int), &dfhist->nlambda);
    dd_bcast(dd, sizeof(real), &dfhist->wl_delta);

    if (dfhist->nlambda > 0)
    {
        int nlam = dfhist->nlambda;
        dd_bcast(dd, sizeof(int) * nlam, dfhist->n_at_lam);
        dd_bcast(dd, sizeof(real) * nlam, dfhist->wl_histo);
        dd_bcast(dd, sizeof(real) * nlam, dfhist->sum_weights);
        dd_bcast(dd, sizeof(real) * nlam, dfhist->sum_dg);
        dd_bcast(dd, sizeof(real) * nlam, dfhist->sum_minvar);
        dd_bcast(dd, sizeof(real) * nlam, dfhist->sum_variance);

        for (int i = 0; i < nlam; i++)
        {
            dd_bcast(dd, sizeof(real) * nlam, dfhist->accum_p[i]);
            dd_bcast(dd, sizeof(real) * nlam, dfhist->accum_m[i]);
            dd_bcast(dd, sizeof(real) * nlam, dfhist->accum_p2[i]);
            dd_bcast(dd, sizeof(real) * nlam, dfhist->accum_m2[i]);
            dd_bcast(dd, sizeof(real) * nlam, dfhist->Tij[i]);
            dd_bcast(dd, sizeof(real) * nlam, dfhist->Tij_empirical[i]);
        }
    }
}

static void dd_distribute_state(gmx_domdec_t* dd, const t_state* state, t_state* state_local)
{
    int nh = state_local->nhchainlength;

    if (DDMAIN(dd))
    {
        GMX_RELEASE_ASSERT(state->nhchainlength == nh,
                           "The global and local Nose-Hoover chain lengths should match");

        for (auto i : gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real>::keys())
        {
            state_local->lambda[i] = state->lambda[i];
        }
        state_local->fep_state = state->fep_state;
        state_local->veta      = state->veta;
        state_local->vol0      = state->vol0;
        copy_mat(state->box, state_local->box);
        copy_mat(state->box_rel, state_local->box_rel);
        copy_mat(state->boxv, state_local->boxv);
        copy_mat(state->svir_prev, state_local->svir_prev);
        copy_mat(state->fvir_prev, state_local->fvir_prev);
        if (state->dfhist != nullptr)
        {
            copy_df_history(state_local->dfhist, state->dfhist);
        }
        for (int i = 0; i < state_local->ngtc; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state_local->nosehoover_xi[i * nh + j]  = state->nosehoover_xi[i * nh + j];
                state_local->nosehoover_vxi[i * nh + j] = state->nosehoover_vxi[i * nh + j];
            }
            state_local->therm_integral[i] = state->therm_integral[i];
        }
        for (int i = 0; i < state_local->nnhpres; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state_local->nhpres_xi[i * nh + j]  = state->nhpres_xi[i * nh + j];
                state_local->nhpres_vxi[i * nh + j] = state->nhpres_vxi[i * nh + j];
            }
        }
        state_local->baros_integral = state->baros_integral;
    }
    dd_bcast(dd,
             (static_cast<int>(FreeEnergyPerturbationCouplingType::Count) * sizeof(real)),
             state_local->lambda.data());
    dd_bcast(dd, sizeof(int), &state_local->fep_state);
    dd_bcast(dd, sizeof(real), &state_local->veta);
    dd_bcast(dd, sizeof(real), &state_local->vol0);
    dd_bcast(dd, sizeof(state_local->box), state_local->box);
    dd_bcast(dd, sizeof(state_local->box_rel), state_local->box_rel);
    dd_bcast(dd, sizeof(state_local->boxv), state_local->boxv);
    dd_bcast(dd, sizeof(state_local->svir_prev), state_local->svir_prev);
    dd_bcast(dd, sizeof(state_local->fvir_prev), state_local->fvir_prev);
    dd_bcast(dd, ((state_local->ngtc * nh) * sizeof(double)), state_local->nosehoover_xi.data());
    dd_bcast(dd, ((state_local->ngtc * nh) * sizeof(double)), state_local->nosehoover_vxi.data());
    dd_bcast(dd, state_local->ngtc * sizeof(double), state_local->therm_integral.data());
    dd_bcast(dd, ((state_local->nnhpres * nh) * sizeof(double)), state_local->nhpres_xi.data());
    dd_bcast(dd, ((state_local->nnhpres * nh) * sizeof(double)), state_local->nhpres_vxi.data());

    /* communicate df_history -- required for restarting from checkpoint */
    dd_distribute_dfhist(dd, state_local->dfhist);

    state_local->changeNumAtoms(dd->comm->atomRanges.numHomeAtoms());

    if (state_local->hasEntry(StateEntry::X))
    {
        distributeVec(dd, DDMAIN(dd) ? state->x : gmx::ArrayRef<const gmx::RVec>(), state_local->x);
    }
    if (state_local->hasEntry(StateEntry::V))
    {
        distributeVec(dd, DDMAIN(dd) ? state->v : gmx::ArrayRef<const gmx::RVec>(), state_local->v);
    }
    if (state_local->hasEntry(StateEntry::Cgp))
    {
        distributeVec(dd, DDMAIN(dd) ? state->cg_p : gmx::ArrayRef<const gmx::RVec>(), state_local->cg_p);
    }
}

/* Computes and returns the domain index for the given atom group.
 *
 * Also updates the coordinates in pos for PBC, when necessary.
 */
static inline int computeAtomGroupDomainIndex(const gmx_domdec_t& dd,
                                              const gmx_ddbox_t&  ddbox,
                                              const matrix&       triclinicCorrectionMatrix,
                                              gmx::ArrayRef<const std::vector<real>> cellBoundaries,
                                              int                                    atomBegin,
                                              int                                    atomEnd,
                                              const matrix                           box,
                                              rvec*                                  pos,
                                              rvec*                                  velocities)
{
    const bool haveBoxDeformation = dd.comm->systemInfo.haveBoxDeformation;

    /* Set the reference location for assigning the group */
    rvec cog;
    int  numAtoms = atomEnd - atomBegin;
    if (numAtoms == 1)
    {
        copy_rvec(pos[atomBegin], cog);
    }
    else
    {
        real invNumAtoms = 1 / static_cast<real>(numAtoms);

        clear_rvec(cog);
        for (int a = atomBegin; a < atomEnd; a++)
        {
            rvec_inc(cog, pos[a]);
        }
        for (int d = 0; d < DIM; d++)
        {
            cog[d] *= invNumAtoms;
        }
    }
    /* Put the charge group in the box and determine the cell index ind */
    ivec ind;
    for (int d = DIM - 1; d >= 0; d--)
    {
        real pos_d = cog[d];
        if (d < dd.unitCellInfo.npbcdim)
        {
            bool bScrew = (dd.unitCellInfo.haveScrewPBC && d == XX);
            if (ddbox.tric_dir[d] && dd.numCells[d] > 1)
            {
                /* Use triclinic coordinates for this dimension */
                for (int j = d + 1; j < DIM; j++)
                {
                    pos_d += cog[j] * triclinicCorrectionMatrix[j][d];
                }
            }
            while (pos_d >= box[d][d])
            {
                pos_d -= box[d][d];
                rvec_dec(cog, box[d]);
                if (bScrew)
                {
                    cog[YY] = box[YY][YY] - cog[YY];
                    cog[ZZ] = box[ZZ][ZZ] - cog[ZZ];
                }
                for (int a = atomBegin; a < atomEnd; a++)
                {
                    rvec_dec(pos[a], box[d]);
                    if (haveBoxDeformation)
                    {
                        correctVelocityForDisplacement<true>(
                                dd.comm->systemInfo.boxDeformationRate, velocities[a], box[d]);
                    }
                    if (bScrew)
                    {
                        pos[a][YY] = box[YY][YY] - pos[a][YY];
                        pos[a][ZZ] = box[ZZ][ZZ] - pos[a][ZZ];
                    }
                }
            }
            while (pos_d < 0)
            {
                pos_d += box[d][d];
                rvec_inc(cog, box[d]);
                if (bScrew)
                {
                    cog[YY] = box[YY][YY] - cog[YY];
                    cog[ZZ] = box[ZZ][ZZ] - cog[ZZ];
                }
                for (int a = atomBegin; a < atomEnd; a++)
                {
                    rvec_inc(pos[a], box[d]);
                    if (haveBoxDeformation)
                    {
                        correctVelocityForDisplacement<false>(
                                dd.comm->systemInfo.boxDeformationRate, velocities[a], box[d]);
                    }
                    if (bScrew)
                    {
                        pos[a][YY] = box[YY][YY] - pos[a][YY];
                        pos[a][ZZ] = box[ZZ][ZZ] - pos[a][ZZ];
                    }
                }
            }
        }
        /* This could be done more efficiently */
        ind[d] = 0;
        while (ind[d] + 1 < dd.numCells[d] && pos_d >= cellBoundaries[d][ind[d] + 1])
        {
            ind[d]++;
        }
    }

    return dd_index(dd.numCells, ind);
}


static std::vector<std::vector<int>> getAtomGroupDistribution(const gmx::MDLogger& mdlog,
                                                              const gmx_mtop_t&    mtop,
                                                              const matrix         box,
                                                              const gmx_ddbox_t&   ddbox,
                                                              rvec                 pos[],
                                                              rvec                 velocities[],
                                                              gmx_domdec_t*        dd)
{
    AtomDistribution& ma = *dd->ma;

    /* Clear the count */
    for (int rank = 0; rank < dd->nnodes; rank++)
    {
        ma.domainGroups[rank].numAtoms = 0;
    }

    matrix triclinicCorrectionMatrix;
    make_tric_corr_matrix(dd->unitCellInfo.npbcdim, box, triclinicCorrectionMatrix);

    ivec       npulse;
    const auto cellBoundaries = set_dd_cell_sizes_slb(dd, &ddbox, setcellsizeslbMAIN, npulse);

    std::vector<std::vector<int>> indices(dd->nnodes);

    if (dd->comm->systemInfo.useUpdateGroups)
    {
        int atomOffset = 0;
        for (const gmx_molblock_t& molblock : mtop.molblock)
        {
            const auto& updateGrouping =
                    dd->comm->systemInfo.updateGroupingsPerMoleculeType[molblock.type];

            for (int mol = 0; mol < molblock.nmol; mol++)
            {
                for (int g = 0; g < updateGrouping.numBlocks(); g++)
                {
                    const auto& block       = updateGrouping.block(g);
                    const int   atomBegin   = atomOffset + block.begin();
                    const int   atomEnd     = atomOffset + block.end();
                    const int   domainIndex = computeAtomGroupDomainIndex(
                            *dd, ddbox, triclinicCorrectionMatrix, cellBoundaries, atomBegin, atomEnd, box, pos, velocities);

                    for (int atomIndex : block)
                    {
                        indices[domainIndex].push_back(atomOffset + atomIndex);
                    }
                    ma.domainGroups[domainIndex].numAtoms += block.size();
                }

                atomOffset += updateGrouping.fullRange().end();
            }
        }

        GMX_RELEASE_ASSERT(atomOffset == mtop.natoms, "Should distribute all atoms");
    }
    else
    {
        /* Compute the center of geometry for all atoms */
        for (int atom = 0; atom < mtop.natoms; atom++)
        {
            int domainIndex = computeAtomGroupDomainIndex(
                    *dd, ddbox, triclinicCorrectionMatrix, cellBoundaries, atom, atom + 1, box, pos, velocities);

            indices[domainIndex].push_back(atom);
            ma.domainGroups[domainIndex].numAtoms += 1;
        }
    }

    {
        // Use double for the sums to avoid natoms^2 overflowing
        // (65537^2 > 2^32)
        int    nat_sum  = 0;
        double nat2_sum = 0;
        int    nat_min  = ma.domainGroups[0].numAtoms;
        int    nat_max  = ma.domainGroups[0].numAtoms;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            int numAtoms = ma.domainGroups[rank].numAtoms;
            nat_sum += numAtoms;
            // convert to double to avoid integer overflows when squaring
            nat2_sum += gmx::square(double(numAtoms));
            nat_min = std::min(nat_min, numAtoms);
            nat_max = std::max(nat_max, numAtoms);
        }
        nat_sum /= dd->nnodes;
        nat2_sum /= dd->nnodes;

        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Atom distribution over %d domains: av %d stddev %d min %d max %d",
                        dd->nnodes,
                        nat_sum,
                        gmx::roundToInt(std::sqrt(nat2_sum - gmx::square(static_cast<double>(nat_sum)))),
                        nat_min,
                        nat_max);
    }

    return indices;
}

static void distributeAtomGroups(const gmx::MDLogger& mdlog,
                                 gmx_domdec_t*        dd,
                                 const gmx_mtop_t&    mtop,
                                 const matrix         box,
                                 const gmx_ddbox_t*   ddbox,
                                 rvec                 pos[],
                                 rvec                 velocities[])
{
    AtomDistribution* ma   = dd->ma.get();
    int *             ibuf = nullptr, buf2[2] = { 0, 0 };
    gmx_bool          bMain = DDMAIN(dd);

    std::vector<std::vector<int>> groupIndices;

    if (bMain)
    {
        GMX_ASSERT(box && pos, "box or pos not set on main");

        if (dd->unitCellInfo.haveScrewPBC)
        {
            check_screw_box(box);
        }

        groupIndices = getAtomGroupDistribution(mdlog, mtop, box, *ddbox, pos, velocities, dd);

        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            ma->intBuffer[rank * 2]     = groupIndices[rank].size();
            ma->intBuffer[rank * 2 + 1] = ma->domainGroups[rank].numAtoms;
        }
        ibuf = ma->intBuffer.data();
    }
    else
    {
        ibuf = nullptr;
    }
    dd_scatter(dd, 2 * sizeof(int), ibuf, buf2);

    dd->numHomeAtoms = buf2[0];
    dd->comm->atomRanges.setEnd(DDAtomRanges::Type::Home, buf2[1]);
    dd->globalAtomIndices.resize(dd->numHomeAtoms);

    if (bMain)
    {
        ma->atomGroups.clear();

        int groupOffset = 0;
        for (int rank = 0; rank < dd->nnodes; rank++)
        {
            ma->intBuffer[rank]              = groupIndices[rank].size();
            ma->intBuffer[dd->nnodes + rank] = groupOffset;

            ma->atomGroups.insert(
                    ma->atomGroups.end(), groupIndices[rank].begin(), groupIndices[rank].end());

            ma->domainGroups[rank].atomGroups = gmx::constArrayRefFromArray(
                    ma->atomGroups.data() + groupOffset, groupIndices[rank].size());

            groupOffset += groupIndices[rank].size();
        }
    }

    dd_scatterv(dd,
                bMain ? gmx::makeArrayRef(ma->intBuffer).subArray(0, dd->nnodes) : gmx::ArrayRef<int>(),
                bMain ? gmx::makeArrayRef(ma->intBuffer).subArray(dd->nnodes, dd->nnodes)
                      : gmx::ArrayRef<int>(),
                bMain ? ma->atomGroups.data() : nullptr,
                dd->numHomeAtoms,
                dd->globalAtomIndices.data());

    if (gmx_debug_at)
    {
        fprintf(debug, "Home charge groups:\n");
        for (int i = 0; i < dd->numHomeAtoms; i++)
        {
            fprintf(debug, " %d", dd->globalAtomIndices[i]);
            if (i % 10 == 9)
            {
                fprintf(debug, "\n");
            }
        }
        fprintf(debug, "\n");
    }
}

void distributeState(const gmx::MDLogger& mdlog,
                     gmx_domdec_t*        dd,
                     const gmx_mtop_t&    mtop,
                     t_state*             state_global,
                     const gmx_ddbox_t&   ddbox,
                     t_state*             state_local)
{
    rvec* xGlobal = (DDMAIN(dd) ? state_global->x.rvec_array() : nullptr);
    rvec* vGlobal = (DDMAIN(dd) ? state_global->v.rvec_array() : nullptr);

    distributeAtomGroups(mdlog, dd, mtop, DDMAIN(dd) ? state_global->box : nullptr, &ddbox, xGlobal, vGlobal);

    dd_distribute_state(dd, state_global, state_local);
}
