/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2005- The GROMACS Authors
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
 * \brief Implements atom redistribution functions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "redistribute.h"

#include <cinttypes>
#include <cstring>

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/template_mp.h"

#include "domdec_internal.h"
#include "utility.h"

/* The size per charge group of the cggl_flag buffer in gmx_domdec_comm_t */
static constexpr int DD_CGIBS = 2;

/* The flags for the cggl_flag buffer in gmx_domdec_comm_t */

/* The lower 16 bits are reserved for the charge group size */
static constexpr int DD_FLAG_NRCG = 65535;

/* Returns which bit tells whether to move a group forward along dimension d */
static inline int DD_FLAG_FW(int d)
{
    return 1 << (16 + d * 2);
}

/* Returns which bit tells whether to move a group backward along dimension d */
static inline int DD_FLAG_BW(int d)
{
    return 1 << (16 + d * 2 + 1);
}

static void copyMovedAtomsToBufferPerAtom(gmx::ArrayRef<const int> move,
                                          int                      nvec,
                                          int                      vec,
                                          rvec*                    src,
                                          gmx_domdec_comm_t*       comm)
{
    int pos_vec[DIM * 2] = { 0 };

    for (gmx::Index i = 0; i < move.ssize(); i++)
    {
        /* Skip moved atoms */
        const int m = move[i];
        if (m >= 0)
        {
            /* Copy to the communication buffer */
            pos_vec[m] += 1 + vec;
            copy_rvec(src[i], comm->cgcm_state[m][pos_vec[m]++]);
            pos_vec[m] += nvec - vec - 1;
        }
    }
}

static void copyMovedUpdateGroupCogs(gmx::ArrayRef<const int>       move,
                                     int                            nvec,
                                     gmx::ArrayRef<const gmx::RVec> coordinates,
                                     gmx_domdec_comm_t*             comm)
{
    int pos_vec[DIM * 2] = { 0 };

    for (gmx::Index g = 0; g < move.ssize(); g++)
    {
        /* Skip moved atoms */
        const int m = move[g];
        if (m >= 0)
        {
            /* Copy to the communication buffer */
            const gmx::RVec& cog =
                    (comm->systemInfo.useUpdateGroups ? comm->updateGroupsCog->cogForAtom(g)
                                                      : coordinates[g]);
            copy_rvec(cog, comm->cgcm_state[m][pos_vec[m]]);
            pos_vec[m] += 1 + nvec;
        }
    }
}

static void clear_and_mark_ind(gmx::ArrayRef<const int> move,
                               gmx::ArrayRef<const int> globalAtomIndices,
                               gmx_ga2la_t*             ga2la,
                               int*                     cell_index)
{
    for (gmx::Index a = 0; a < move.ssize(); a++)
    {
        if (move[a] >= 0)
        {
            /* Clear the global indices */
            ga2la->erase(globalAtomIndices[a]);
            /* Signal that this atom has moved using the ns cell index.
             * Here we set it to -1. fill_grid will change it
             * from -1 to NSGRID_SIGNAL_MOVED_FAC*grid->ncells.
             */
            cell_index[a] = -1;
        }
    }
}

static void print_cg_move(FILE*               fplog,
                          const gmx_domdec_t* dd,
                          int64_t             step,
                          int                 cg,
                          int                 dim,
                          int                 dir,
                          gmx_bool            bHaveCgcmOld,
                          real                limitd,
                          rvec                cm_old,
                          rvec                cm_new,
                          real                pos_d)
{
    const gmx_domdec_comm_t* comm = dd->comm.get();
    std::string              mesg;

    fprintf(fplog, "\nStep %" PRId64 ":\n", step);

    if (comm->systemInfo.useUpdateGroups)
    {
        mesg += "The update group starting at atom";
    }
    else
    {
        mesg += "Atom";
    }
    mesg += gmx::formatString(
            " %d moved more than the distance allowed by the domain decomposition", ddglatnr(dd, cg));
    if (limitd > 0)
    {
        mesg += gmx::formatString(" (%f)", limitd);
    }
    mesg += gmx::formatString(" in direction %c\n", dim2char(dim));
    fprintf(fplog, "%s", mesg.c_str());

    fprintf(fplog,
            "distance out of cell %f\n",
            dir == 1 ? pos_d - comm->cell_x1[dim] : pos_d - comm->cell_x0[dim]);
    if (bHaveCgcmOld)
    {
        fprintf(fplog, "Old coordinates: %8.3f %8.3f %8.3f\n", cm_old[XX], cm_old[YY], cm_old[ZZ]);
    }
    fprintf(fplog, "New coordinates: %8.3f %8.3f %8.3f\n", cm_new[XX], cm_new[YY], cm_new[ZZ]);
    fprintf(fplog,
            "Old cell boundaries in direction %c: %8.3f %8.3f\n",
            dim2char(dim),
            comm->old_cell_x0[dim],
            comm->old_cell_x1[dim]);
    fprintf(fplog,
            "New cell boundaries in direction %c: %8.3f %8.3f\n",
            dim2char(dim),
            comm->cell_x0[dim],
            comm->cell_x1[dim]);
}

[[noreturn]] static void cg_move_error(FILE*               fplog,
                                       const gmx_domdec_t* dd,
                                       int64_t             step,
                                       int                 cg,
                                       int                 dim,
                                       int                 dir,
                                       gmx_bool            bHaveCgcmOld,
                                       real                limitd,
                                       rvec                cm_old,
                                       rvec                cm_new,
                                       real                pos_d)
{
    if (fplog)
    {
        print_cg_move(fplog, dd, step, cg, dim, dir, bHaveCgcmOld, limitd, cm_old, cm_new, pos_d);
    }
    print_cg_move(stderr, dd, step, cg, dim, dir, bHaveCgcmOld, limitd, cm_old, cm_new, pos_d);
    gmx_fatal(FARGS,
              "One or more atoms moved too far between two domain decomposition steps.\n"
              "This usually means that your system is not well equilibrated");
}

static void rotate_state_atom(t_state* state, int a)
{
    if (state->hasEntry(StateEntry::X))
    {
        auto x = makeArrayRef(state->x);
        /* Rotate the complete state; for a rectangular box only */
        x[a][YY] = state->box[YY][YY] - x[a][YY];
        x[a][ZZ] = state->box[ZZ][ZZ] - x[a][ZZ];
    }
    if (state->hasEntry(StateEntry::V))
    {
        auto v   = makeArrayRef(state->v);
        v[a][YY] = -v[a][YY];
        v[a][ZZ] = -v[a][ZZ];
    }
    if (state->hasEntry(StateEntry::Cgp))
    {
        auto cg_p   = makeArrayRef(state->cg_p);
        cg_p[a][YY] = -cg_p[a][YY];
        cg_p[a][ZZ] = -cg_p[a][ZZ];
    }
}

/* Returns a pointer to a buffer of size \p numAtomsNew with contents 0 from numAtomsOld upwards
 *
 * Note: numAtomsOld should either be 0 or match the current buffer size.
 */
static int* getMovedBuffer(gmx_domdec_comm_t* comm, size_t numAtomsOld, size_t numAtomsNew)
{
    std::vector<int>& movedBuffer = comm->movedBuffer;

    GMX_RELEASE_ASSERT(numAtomsOld == 0 || movedBuffer.size() == numAtomsOld,
                       "numAtomsOld should either be 0 or match the current size");

    /* Contents up to numAtomsOld should be preserved, clear from numAtomsOld */
    if (numAtomsOld == 0)
    {
        movedBuffer.clear();
    }
    movedBuffer.resize(numAtomsNew);

    return movedBuffer.data();
}

/* Bounds to determine whether an atom group moved to far between DD steps */
struct MoveLimits
{
    rvec distance; /* Maximum distance over which a group can move */
    rvec lower;    /* Lower bound for group location */
    rvec upper;    /* Upper bound for group location */
};

/* Compute flag that tells where to move an atom group
 *
 * The input dev tells whether the group moves fw,remains,bw (-1,0,1) along
 * dimensions.
 * The return value is -1 when the group does not move.
 * The return has move flags set when the group does move and the lower 4 bits
 * are (mis)used to tell which is the first dimension (bit 1,2,3) the group
 * needs to be moved along and in which direction (bit 0 not set for fw
 * and set for bw), however only values in [0,6) are used.
 */
static int computeMoveFlag(const gmx_domdec_t& dd, const ivec& dev)
{
    int flag              = 0;
    int firstMoveDimValue = -1;
    for (int d = 0; d < dd.ndim; d++)
    {
        const int dim = dd.dim[d];
        if (dev[dim] == 1)
        {
            flag |= DD_FLAG_FW(d);
            if (firstMoveDimValue == -1)
            {
                firstMoveDimValue = d * 2;
            }
        }
        else if (dev[dim] == -1)
        {
            flag |= DD_FLAG_BW(d);
            if (firstMoveDimValue == -1)
            {
                if (dd.numCells[dim] > 2)
                {
                    firstMoveDimValue = d * 2 + 1;
                }
                else
                {
                    firstMoveDimValue = d * 2;
                }
            }
        }
    }

    return firstMoveDimValue + flag;
}

/* Determine to which atoms in the range \p cg_start, \p cg_end should go.
 *
 * Returns in the move array where the atoms should go.
 */
static void calc_cg_move(FILE*                  fplog,
                         int64_t                step,
                         gmx_domdec_t*          dd,
                         t_state*               state,
                         const ivec             tric_dir,
                         matrix                 tcm,
                         const rvec             cell_x0,
                         const rvec             cell_x1,
                         const MoveLimits&      moveLimits,
                         const gmx::Range<int>& atomRange,
                         gmx::ArrayRef<int>     move)
{
    const int npbcdim = dd->unitCellInfo.npbcdim;
    auto      x       = makeArrayRef(state->x);

    for (int a : atomRange)
    {
        // TODO: Rename this center of geometry variable to cogNew
        rvec cm_new;
        copy_rvec(x[a], cm_new);

        ivec dev = { 0 };
        /* Do pbc and check DD cell boundary crossings */
        for (int d = DIM - 1; d >= 0; d--)
        {
            if (dd->numCells[d] > 1)
            {
                bool bScrew = (dd->unitCellInfo.haveScrewPBC && d == XX);
                /* Determine the location of this cg in lattice coordinates */
                real pos_d = cm_new[d];
                if (tric_dir[d])
                {
                    for (int d2 = d + 1; d2 < DIM; d2++)
                    {
                        pos_d += cm_new[d2] * tcm[d2][d];
                    }
                }
                /* Put the charge group in the triclinic unit-cell */
                if (pos_d >= cell_x1[d])
                {
                    if (pos_d >= moveLimits.upper[d])
                    {
                        cg_move_error(
                                fplog, dd, step, a, d, 1, false, moveLimits.distance[d], cm_new, cm_new, pos_d);
                    }
                    dev[d] = 1;
                    if (dd->ci[d] == dd->numCells[d] - 1)
                    {
                        rvec_dec(cm_new, state->box[d]);
                        if (bScrew)
                        {
                            cm_new[YY] = state->box[YY][YY] - cm_new[YY];
                            cm_new[ZZ] = state->box[ZZ][ZZ] - cm_new[ZZ];
                        }
                        rvec_dec(x[a], state->box[d]);
                        if (bScrew)
                        {
                            rotate_state_atom(state, a);
                        }
                    }
                }
                else if (pos_d < cell_x0[d])
                {
                    if (pos_d < moveLimits.lower[d])
                    {
                        cg_move_error(
                                fplog, dd, step, a, d, -1, false, moveLimits.distance[d], cm_new, cm_new, pos_d);
                    }
                    dev[d] = -1;
                    if (dd->ci[d] == 0)
                    {
                        rvec_inc(cm_new, state->box[d]);
                        if (bScrew)
                        {
                            cm_new[YY] = state->box[YY][YY] - cm_new[YY];
                            cm_new[ZZ] = state->box[ZZ][ZZ] - cm_new[ZZ];
                        }
                        rvec_inc(x[a], state->box[d]);
                        if (bScrew)
                        {
                            rotate_state_atom(state, a);
                        }
                    }
                }
            }
            else if (d < npbcdim)
            {
                /* Put the charge group in the rectangular unit-cell */
                while (cm_new[d] >= state->box[d][d])
                {
                    rvec_dec(cm_new, state->box[d]);
                    rvec_dec(x[a], state->box[d]);
                }
                while (cm_new[d] < 0)
                {
                    rvec_inc(cm_new, state->box[d]);
                    rvec_inc(x[a], state->box[d]);
                }
            }
        }

        /* Temporarily store the flag in move */
        move[a] = computeMoveFlag(*dd, dev);
    }
}

struct PbcAndFlag
{
    /* Constructor that purposely does not initialize anything */
    PbcAndFlag() {}

    gmx::RVec pbcShift;
    int       moveFlag;
};

/* Determine to which domains update groups in the range \p groupBegin, \p groupEnd should go.
 *
 * Returns in the move array where the groups should go.
 * Also updates the COGs and coordinates for jumps over periodic boundaries.
 */
static void calcGroupMove(FILE*                     fplog,
                          int64_t                   step,
                          const gmx_domdec_t*       dd,
                          const t_state*            state,
                          const ivec                tric_dir,
                          matrix                    tcm,
                          const rvec                cell_x0,
                          const rvec                cell_x1,
                          const MoveLimits&         moveLimits,
                          const gmx::Range<int>&    groupRange,
                          gmx::ArrayRef<PbcAndFlag> pbcAndFlags)
{
    GMX_RELEASE_ASSERT(!dd->unitCellInfo.haveScrewPBC, "Screw PBC is not supported here");

    const int npbcdim = dd->unitCellInfo.npbcdim;

    gmx::UpdateGroupsCog* updateGroupsCog = dd->comm->updateGroupsCog.get();

    for (int g : groupRange)
    {

        gmx::RVec& cog    = updateGroupsCog->cog(g);
        gmx::RVec  cogOld = cog;

        ivec dev = { 0 };
        /* Do pbc and check DD cell boundary crossings */
        for (int d = DIM - 1; d >= 0; d--)
        {
            if (dd->numCells[d] > 1)
            {
                /* Determine the location of this COG in lattice coordinates */
                real pos_d = cog[d];
                if (tric_dir[d])
                {
                    for (int d2 = d + 1; d2 < DIM; d2++)
                    {
                        pos_d += cog[d2] * tcm[d2][d];
                    }
                }
                /* Put the COG in the triclinic unit-cell */
                if (pos_d >= cell_x1[d])
                {
                    if (pos_d >= moveLimits.upper[d])
                    {
                        cg_move_error(
                                fplog, dd, step, g, d, 1, true, moveLimits.distance[d], cogOld, cog, pos_d);
                    }
                    dev[d] = 1;
                    if (dd->ci[d] == dd->numCells[d] - 1)
                    {
                        rvec_dec(cog, state->box[d]);
                    }
                }
                else if (pos_d < cell_x0[d])
                {
                    if (pos_d < moveLimits.lower[d])
                    {
                        cg_move_error(
                                fplog, dd, step, g, d, -1, true, moveLimits.distance[d], cogOld, cog, pos_d);
                    }
                    dev[d] = -1;
                    if (dd->ci[d] == 0)
                    {
                        rvec_inc(cog, state->box[d]);
                    }
                }
            }
            else if (d < npbcdim)
            {
                /* Put the COG in the rectangular unit-cell */
                while (cog[d] >= state->box[d][d])
                {
                    rvec_dec(cog, state->box[d]);
                }
                while (cog[d] < 0)
                {
                    rvec_inc(cog, state->box[d]);
                }
            }
        }

        /* Store the PBC and move flag, so we can later apply them to the atoms */
        PbcAndFlag& pbcAndFlag = pbcAndFlags[g];

        rvec_sub(cog, cogOld, pbcAndFlag.pbcShift);
        pbcAndFlag.moveFlag = computeMoveFlag(*dd, dev);
    }
}

template<bool haveBoxDeformation>
static void applyPbcAndSetMoveFlags(const gmx::UpdateGroupsCog&         updateGroupsCog,
                                    gmx::ArrayRef<const PbcAndFlag>     pbcAndFlags,
                                    const matrix gmx_unused             boxDeformationRate,
                                    const gmx::Range<int>&              atomRange,
                                    gmx::ArrayRef<gmx::RVec>            atomCoords,
                                    gmx::ArrayRef<gmx::RVec> gmx_unused atomVelocities,
                                    gmx::ArrayRef<int>                  move)
{
    for (int a : atomRange)
    {
        const PbcAndFlag& pbcAndFlag = pbcAndFlags[updateGroupsCog.cogIndex(a)];
        rvec_inc(atomCoords[a], pbcAndFlag.pbcShift);
        if constexpr (haveBoxDeformation)
        {
            // Correct the velocity for the position change along the flow profile
            correctVelocityForDisplacement<false>(boxDeformationRate, atomVelocities[a], pbcAndFlag.pbcShift);
        }
        /* Temporarily store the flag in move */
        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        move[a] = pbcAndFlag.moveFlag;
    }
}

// Return a item range for thread \p thread when dividing work equally over \p numThreads threads
static gmx::Range<int> getThreadLocalRange(const int numItems, const int numThreads, const int thread)
{
    GMX_ASSERT(numThreads > 0, "There should be at least one thread");
    GMX_ASSERT(thread < numThreads, "The thread index should be smaller than the number of threads");

    return { (numItems * thread) / numThreads, (numItems * (thread + 1)) / numThreads };
}

void dd_redistribute_cg(FILE*         fplog,
                        int64_t       step,
                        gmx_domdec_t* dd,
                        ivec          tric_dir,
                        t_state*      state,
                        t_forcerec*   fr,
                        t_nrnb*       nrnb,
                        int*          ncg_moved)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    if (dd->unitCellInfo.haveScrewPBC)
    {
        check_screw_box(state->box);
    }

    // Positions are always present, so there's nothing to flag
    const bool bV   = state->hasEntry(StateEntry::V);
    const bool bCGP = state->hasEntry(StateEntry::Cgp);

    DDBufferAccess<int> moveBuffer(comm->intBuffer, dd->numHomeAtoms);
    gmx::ArrayRef<int>  move = moveBuffer.buffer;

    const int npbcdim = dd->unitCellInfo.npbcdim;

    rvec       cell_x0, cell_x1;
    MoveLimits moveLimits;
    for (int d = 0; (d < DIM); d++)
    {
        moveLimits.distance[d] = dd->comm->cellsize_min[d];
        if (d >= npbcdim && dd->ci[d] == 0)
        {
            cell_x0[d] = -GMX_FLOAT_MAX;
        }
        else
        {
            cell_x0[d] = comm->cell_x0[d];
        }
        if (d >= npbcdim && dd->ci[d] == dd->numCells[d] - 1)
        {
            cell_x1[d] = GMX_FLOAT_MAX;
        }
        else
        {
            cell_x1[d] = comm->cell_x1[d];
        }
        if (d < npbcdim)
        {
            moveLimits.lower[d] = comm->old_cell_x0[d] - moveLimits.distance[d];
            moveLimits.upper[d] = comm->old_cell_x1[d] + moveLimits.distance[d];
        }
        else
        {
            /* We check after communication if a charge group moved
             * more than one cell. Set the pre-comm check limit to float_max.
             */
            moveLimits.lower[d] = -GMX_FLOAT_MAX;
            moveLimits.upper[d] = GMX_FLOAT_MAX;
        }
    }

    matrix tcm;
    make_tric_corr_matrix(npbcdim, state->box, tcm);

    const int nthread = gmx_omp_nthreads_get(ModuleMultiThread::Domdec);

    /* Compute the center of geometry for all home charge groups
     * and put them in the box and determine where they should go.
     */
    std::vector<PbcAndFlag> pbcAndFlags(
            comm->systemInfo.useUpdateGroups ? comm->updateGroupsCog->numCogs() : 0);

#pragma omp parallel num_threads(nthread)
    {
        try
        {
            const int thread = gmx_omp_get_thread_num();

            const int             numHomeAtoms    = comm->atomRanges.numHomeAtoms();
            const gmx::Range<int> threadAtomRange = getThreadLocalRange(numHomeAtoms, nthread, thread);

            if (comm->systemInfo.useUpdateGroups)
            {
                const auto&           updateGroupsCog = *comm->updateGroupsCog;
                const int             numGroups       = updateGroupsCog.numCogs();
                const gmx::Range<int> threadGroupRange = getThreadLocalRange(numGroups, nthread, thread);
                calcGroupMove(
                        fplog, step, dd, state, tric_dir, tcm, cell_x0, cell_x1, moveLimits, threadGroupRange, pbcAndFlags);
                /* We need a barrier as atoms below can be in a COG of a different thread */
#pragma omp barrier
                gmx::dispatchTemplatedFunction(
                        [&](auto haveBoxDeformation) {
                            applyPbcAndSetMoveFlags<haveBoxDeformation>(updateGroupsCog,
                                                                        pbcAndFlags,
                                                                        comm->systemInfo.boxDeformationRate,
                                                                        threadAtomRange,
                                                                        state->x,
                                                                        state->v,
                                                                        move);
                        },
                        comm->systemInfo.haveBoxDeformation);
            }
            else
            {
                /* Here we handle single atoms or charge groups */
                calc_cg_move(fplog, step, dd, state, tric_dir, tcm, cell_x0, cell_x1, moveLimits, threadAtomRange, move);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    // The counts of atoms to move, forward or backward, over the
    // possible DIM dimensions.
    int nat[DIM * 2] = { 0 };
    for (int cg = 0; cg < dd->numHomeAtoms; cg++)
    {
        if (move[cg] >= 0)
        {
            // The value in move[cg] was computed by computeMoveFlags
            // and describes how this atom should move between domains.
            const int flag = move[cg] & ~DD_FLAG_NRCG;
            // mc contains 4 bits that tell which is the first
            // dimension (bit 1,2,3) that the group needs to be moved
            // along, and in which direction (bit 0; not set for fw
            // and set for bw). However the value is always in
            // the range [0,6)
            const int mc = move[cg] & DD_FLAG_NRCG;
            move[cg]     = mc;

            std::vector<int>& cggl_flag = comm->cggl_flag[mc];

            /* TODO: See if we can use push_back instead */
            if ((nat[mc] + 1) * DD_CGIBS > gmx::Index(cggl_flag.size()))
            {
                cggl_flag.resize((nat[mc] + 1) * DD_CGIBS);
            }
            cggl_flag[nat[mc] * DD_CGIBS]     = dd->globalAtomIndices[cg];
            cggl_flag[nat[mc] * DD_CGIBS + 1] = flag;
            nat[mc]++;
        }
    }

    inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());
    inc_nrnb(nrnb, eNR_RESETX, dd->numHomeAtoms);

    *ncg_moved = 0;
    for (int i = 0; i < dd->ndim * 2; i++)
    {
        *ncg_moved += nat[i];
    }

    int nvec = 1;
    if (bV)
    {
        nvec++;
    }
    if (bCGP)
    {
        nvec++;
    }

    /* Make sure the communication buffers are large enough */
    for (int mc = 0; mc < dd->ndim * 2; mc++)
    {
        size_t nvr = nat[mc] * (1 + nvec);
        if (nvr > comm->cgcm_state[mc].size())
        {
            comm->cgcm_state[mc].resize(nvr);
        }
    }

    /* With update groups we send over their COGs.
     * Without update groups we send the moved atom coordinates
     * over twice. This is so the code further down can be used
     * without many conditionals both with and without update groups.
     */
    copyMovedUpdateGroupCogs(move, nvec, state->x, comm);

    int vectorIndex = 0;
    copyMovedAtomsToBufferPerAtom(move, nvec, vectorIndex++, state->x.rvec_array(), comm);
    if (bV)
    {
        copyMovedAtomsToBufferPerAtom(move, nvec, vectorIndex++, state->v.rvec_array(), comm);
    }
    if (bCGP)
    {
        copyMovedAtomsToBufferPerAtom(move, nvec, vectorIndex++, state->cg_p.rvec_array(), comm);
    }

    int* moved = getMovedBuffer(comm, 0, dd->numHomeAtoms);

    clear_and_mark_ind(move, dd->globalAtomIndices, dd->ga2la.get(), moved);

    /* Now we can remove the excess global atom indices from the list */
    dd->globalAtomIndices.resize(dd->numHomeAtoms);

    /* We reuse the intBuffer without reacquiring since we are in the same scope */
    DDBufferAccess<int>& flagBuffer = moveBuffer;

    gmx::ArrayRef<const gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock =
            fr->atomInfoForEachMoleculeBlock;

    /* Temporarily store atoms passed to our rank at the end of the range */
    int home_pos_at = dd->numHomeAtoms;
    for (int d = 0; d < dd->ndim; d++)
    {
        DDBufferAccess<gmx::RVec> rvecBuffer(comm->rvecBuffer, 0);

        const int dim                = dd->dim[d];
        int       totalAtomsReceived = 0;
        int       nvr                = 0;
        for (int dir = 0; dir < (dd->numCells[dim] == 2 ? 1 : 2); dir++)
        {
            const int cdd = d * 2 + dir;
            /* Communicate the atom counts */
            if (debug)
            {
                fprintf(debug, "Sending ddim %d dir %d: nat %d\n", d, dir, nat[cdd]);
            }
            int atomCountToReceive;
            ddSendrecv(dd,
                       d,
                       dir,
                       gmx::arrayRefFromArray(&nat[cdd], 1),
                       gmx::arrayRefFromArray(&atomCountToReceive, 1));

            flagBuffer.resize((totalAtomsReceived + atomCountToReceive) * DD_CGIBS);

            /* Communicate the charge group indices, sizes and flags */
            ddSendrecv(dd,
                       d,
                       dir,
                       gmx::arrayRefFromArray(comm->cggl_flag[cdd].data(), nat[cdd] * DD_CGIBS),
                       gmx::arrayRefFromArray(flagBuffer.buffer.data() + totalAtomsReceived * DD_CGIBS,
                                              atomCountToReceive * DD_CGIBS));

            const int nvs = nat[cdd] * (1 + nvec);
            const int i   = atomCountToReceive * (1 + nvec);
            rvecBuffer.resize(nvr + i);

            /* Communicate cgcm and state */
            ddSendrecv(dd,
                       d,
                       dir,
                       gmx::arrayRefFromArray(comm->cgcm_state[cdd].data(), nvs),
                       gmx::arrayRefFromArray(rvecBuffer.buffer.data() + nvr, i));
            totalAtomsReceived += atomCountToReceive;
            nvr += i;
        }

        dd_resize_atominfo_and_state(fr, state, home_pos_at + totalAtomsReceived);

        /* Process the received charge or update groups */
        int buf_pos = 0;
        for (int cg = 0; cg < totalAtomsReceived; cg++)
        {
            /* Extract the move flags and COG for the charge or update group */
            int              flag = flagBuffer.buffer[cg * DD_CGIBS + 1];
            const gmx::RVec& cog  = rvecBuffer.buffer[buf_pos];

            if (dim >= npbcdim && dd->numCells[dim] > 2)
            {
                /* No pbc in this dim and more than one domain boundary.
                 * We do a separate check if a charge group didn't move too far.
                 */
                if (((flag & DD_FLAG_FW(d)) && cog[dim] > cell_x1[dim])
                    || ((flag & DD_FLAG_BW(d)) && cog[dim] < cell_x0[dim]))
                {
                    rvec pos = { cog[0], cog[1], cog[2] };
                    cg_move_error(
                            fplog, dd, step, cg, dim, (flag & DD_FLAG_FW(d)) ? 1 : 0, false, 0, pos, pos, pos[dim]);
                }
            }

            int mc = -1;
            if (d < dd->ndim - 1)
            {
                /* Check which direction this cg should go */
                for (int d2 = d + 1; (d2 < dd->ndim && mc == -1); d2++)
                {
                    if (isDlbOn(dd->comm->dlbState))
                    {
                        /* The cell boundaries for dimension d2 are not equal
                         * for each cell row of the lower dimension(s),
                         * therefore we might need to redetermine where
                         * this cg should go.
                         */
                        const int dim2 = dd->dim[d2];
                        /* The DD grid is not staggered at the box boundaries,
                         * so we do not need to handle boundary crossings.
                         * This also means we do not have to handle PBC here.
                         */
                        if (!((dd->ci[dim2] == dd->numCells[dim2] - 1 && (flag & DD_FLAG_FW(d2)))
                              || (dd->ci[dim2] == 0 && (flag & DD_FLAG_BW(d2)))))
                        {
                            /* Clear the two flags for this dimension */
                            flag &= ~(DD_FLAG_FW(d2) | DD_FLAG_BW(d2));
                            /* Determine the location of this cg
                             * in lattice coordinates
                             */
                            real pos_d = cog[dim2];
                            if (tric_dir[dim2])
                            {
                                for (int d3 = dim2 + 1; d3 < DIM; d3++)
                                {
                                    pos_d += cog[d3] * tcm[d3][dim2];
                                }
                            }

                            GMX_ASSERT(dim2 >= 0 && dim2 < DIM, "Keep the static analyzer happy");

                            /* Check if we need to move this group
                             * to an adjacent cell because of the
                             * staggering.
                             */
                            if (pos_d >= cell_x1[dim2] && dd->ci[dim2] != dd->numCells[dim2] - 1)
                            {
                                flag |= DD_FLAG_FW(d2);
                            }
                            else if (pos_d < cell_x0[dim2] && dd->ci[dim2] != 0)
                            {
                                flag |= DD_FLAG_BW(d2);
                            }

                            flagBuffer.buffer[cg * DD_CGIBS + 1] = flag;
                        }
                    }
                    /* Set to which neighboring cell this cg should go */
                    if (flag & DD_FLAG_FW(d2))
                    {
                        mc = d2 * 2;
                    }
                    else if (flag & DD_FLAG_BW(d2))
                    {
                        if (dd->numCells[dd->dim[d2]] > 2)
                        {
                            mc = d2 * 2 + 1;
                        }
                        else
                        {
                            mc = d2 * 2;
                        }
                    }
                }
            }

            if (mc == -1)
            {
                /* Set the global charge group index and size */
                const int globalAtomIndex = flagBuffer.buffer[cg * DD_CGIBS];
                dd->globalAtomIndices.push_back(globalAtomIndex);
                /* Skip the COG entry in the buffer */
                buf_pos++;

                /* Set the cginfo */
                fr->atomInfo[home_pos_at] = ddGetAtomInfo(atomInfoForEachMoleculeBlock, globalAtomIndex);

                auto  x       = makeArrayRef(state->x);
                auto  v       = makeArrayRef(state->v);
                auto  cg_p    = makeArrayRef(state->cg_p);
                rvec* rvecPtr = as_rvec_array(rvecBuffer.buffer.data());
                copy_rvec(rvecPtr[buf_pos++], x[home_pos_at]);
                if (bV)
                {
                    copy_rvec(rvecPtr[buf_pos++], v[home_pos_at]);
                }
                if (bCGP)
                {
                    copy_rvec(rvecPtr[buf_pos++], cg_p[home_pos_at]);
                }
                home_pos_at++;
            }
            else
            {
                /* Reallocate the buffers if necessary  */
                if ((nat[mc] + 1) * DD_CGIBS > gmx::Index(comm->cggl_flag[mc].size()))
                {
                    comm->cggl_flag[mc].resize((nat[mc] + 1) * DD_CGIBS);
                }
                size_t nvr = nat[mc] * (1 + nvec);
                if (nvr + 1 + nvec > comm->cgcm_state[mc].size())
                {
                    comm->cgcm_state[mc].resize(nvr + 1 + nvec);
                }
                /* Copy from the receive to the send buffers */
                memcpy(comm->cggl_flag[mc].data() + nat[mc] * DD_CGIBS,
                       flagBuffer.buffer.data() + cg * DD_CGIBS,
                       DD_CGIBS * sizeof(int));
                memcpy(comm->cgcm_state[mc][nvr],
                       rvecBuffer.buffer.data() + buf_pos,
                       (1 + nvec) * sizeof(rvec));
                buf_pos += 1 + nvec;
                nat[mc]++;
            }
        }
    }

    /* Note that the indices are now only partially up to date
     * and numHomeAtoms and nat_home are not the real count, since there are
     * "holes" in the arrays for the charge groups that moved to neighbors.
     */

    /* We need to clear the moved flags for the received atoms,
     * because the moved buffer will be passed to the nbnxm gridding call.
     */
    moved = getMovedBuffer(comm, dd->numHomeAtoms, home_pos_at);

    for (int i = dd->numHomeAtoms; i < home_pos_at; i++)
    {
        moved[i] = 0;
    }

    dd->numHomeAtoms = home_pos_at;
    comm->atomRanges.setEnd(DDAtomRanges::Type::Home, home_pos_at);

    if (debug)
    {
        fprintf(debug,
                "Finished repartitioning: cgs moved out %d, new home %d\n",
                *ncg_moved,
                dd->numHomeAtoms - *ncg_moved);
    }
}
