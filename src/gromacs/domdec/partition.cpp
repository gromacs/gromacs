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
/*! \internal \file
 *
 * \brief This file defines functions for mdrun to call to make a new
 * domain decomposition, and check it.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "partition.h"

#include "config.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/domdec/localtopology.h"
#include "gromacs/domdec/localtopologychecker.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/nsgrid.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/imd/imd.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "box.h"
#include "cellsizes.h"
#include "distribute.h"
#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_vsite.h"
#include "dump.h"
#include "redistribute.h"
#include "utility.h"

/*! \brief Turn on DLB when the load imbalance causes this amount of total loss.
 *
 * There is a bit of overhead with DLB and it's difficult to achieve
 * a load imbalance of less than 2% with DLB.
 */
#define DD_PERF_LOSS_DLB_ON 0.02

//! Warn about imbalance due to PP or PP/PME load imbalance at this loss.
#define DD_PERF_LOSS_WARN 0.05


//! Debug helper printing a DD zone
static void print_ddzone(FILE* fp, int d, int i, int j, gmx_ddzone_t* zone)
{
    fprintf(fp,
            "zone d0 %d d1 %d d2 %d  min0 %6.3f max1 %6.3f mch0 %6.3f mch1 %6.3f p1_0 %6.3f p1_1 "
            "%6.3f\n",
            d,
            i,
            j,
            zone->min0,
            zone->max1,
            zone->mch0,
            zone->mch0,
            zone->p1_0,
            zone->p1_1);
}

/*! \brief Using the home grid size as input in cell_ns_x0 and cell_ns_x1
 * takes the extremes over all home and remote zones in the halo
 * and returns the results in cell_ns_x0 and cell_ns_x1.
 * Note: only used with the group cut-off scheme.
 */
static void dd_move_cellx(gmx_domdec_t* dd, const gmx_ddbox_t* ddbox, rvec cell_ns_x0, rvec cell_ns_x1)
{
    constexpr int      c_ddZoneCommMaxNumZones = 5;
    gmx_ddzone_t       buf_s[c_ddZoneCommMaxNumZones];
    gmx_ddzone_t       buf_r[c_ddZoneCommMaxNumZones];
    gmx_ddzone_t       buf_e[c_ddZoneCommMaxNumZones];
    gmx_domdec_comm_t* comm = dd->comm.get();

    std::array<gmx::RVec, 2> extr_s;
    std::array<gmx::RVec, 2> extr_r;
    for (int d = 1; d < dd->ndim; d++)
    {
        int           dim = dd->dim[d];
        gmx_ddzone_t& zp  = (d == 1) ? comm->zone_d1[0] : comm->zone_d2[0][0];

        /* Copy the base sizes of the home zone */
        zp.min0    = cell_ns_x0[dim];
        zp.max1    = cell_ns_x1[dim];
        zp.min1    = cell_ns_x1[dim];
        zp.mch0    = cell_ns_x0[dim];
        zp.mch1    = cell_ns_x1[dim];
        zp.p1_0    = cell_ns_x0[dim];
        zp.p1_1    = cell_ns_x1[dim];
        zp.dataSet = 1;
    }

    gmx::ArrayRef<DDCellsizesWithDlb> cellsizes = comm->cellsizesWithDlb;

    /* Loop backward over the dimensions and aggregate the extremes
     * of the cell sizes.
     */
    for (int d = dd->ndim - 2; d >= 0; d--)
    {
        const int  dim      = dd->dim[d];
        const bool applyPbc = (dim < ddbox->npbcdim);

        /* Use an rvec to store two reals */
        extr_s[d][0] = cellsizes[d + 1].fracLower;
        extr_s[d][1] = cellsizes[d + 1].fracUpper;
        extr_s[d][2] = cellsizes[d + 1].fracUpper;

        int pos = 0;
        GMX_ASSERT(pos < c_ddZoneCommMaxNumZones, "The buffers should be sufficiently large");
        /* Store the extremes in the backward sending buffer,
         * so they get updated separately from the forward communication.
         */
        for (int d1 = d; d1 < dd->ndim - 1; d1++)
        {
            gmx_ddzone_t& buf = buf_s[pos];

            /* We invert the order to be able to use the same loop for buf_e */
            buf.min0 = extr_s[d1][1];
            buf.max1 = extr_s[d1][0];
            buf.min1 = extr_s[d1][2];
            buf.mch0 = 0;
            buf.mch1 = 0;
            /* Store the cell corner of the dimension we communicate along */
            buf.p1_0    = comm->cell_x0[dim];
            buf.p1_1    = 0;
            buf.dataSet = 1;
            pos++;
        }

        buf_s[pos] = (dd->ndim == 2) ? comm->zone_d1[0] : comm->zone_d2[0][0];
        pos++;

        if (dd->ndim == 3 && d == 0)
        {
            buf_s[pos] = comm->zone_d2[0][1];
            pos++;
            buf_s[pos] = comm->zone_d1[0];
            pos++;
        }

        /* We only need to communicate the extremes
         * in the forward direction
         */
        int numPulses = comm->cd[d].numPulses();
        int numPulsesMin;
        if (applyPbc)
        {
            /* Take the minimum to avoid double communication */
            numPulsesMin = std::min(numPulses, dd->numCells[dim] - 1 - numPulses);
        }
        else
        {
            /* Without PBC we should really not communicate over
             * the boundaries, but implementing that complicates
             * the communication setup and therefore we simply
             * do all communication, but ignore some data.
             */
            numPulsesMin = numPulses;
        }
        for (int pulse = 0; pulse < numPulsesMin; pulse++)
        {
            /* Communicate the extremes forward */
            bool receiveValidData = (applyPbc || dd->ci[dim] > 0);

            int numElements = dd->ndim - d - 1;
            ddSendrecv(dd,
                       d,
                       dddirForward,
                       gmx::arrayRefFromArray(extr_s.data() + d, numElements),
                       gmx::arrayRefFromArray(extr_r.data() + d, numElements));

            if (receiveValidData)
            {
                for (int d1 = d; d1 < dd->ndim - 1; d1++)
                {
                    extr_s[d1][0] = std::max(extr_s[d1][0], extr_r[d1][0]);
                    extr_s[d1][1] = std::min(extr_s[d1][1], extr_r[d1][1]);
                    extr_s[d1][2] = std::min(extr_s[d1][2], extr_r[d1][2]);
                }
            }
        }

        const int numElementsInBuffer = pos;
        for (int pulse = 0; pulse < numPulses; pulse++)
        {
            /* Communicate all the zone information backward */
            bool receiveValidData = (applyPbc || dd->ci[dim] < dd->numCells[dim] - 1);

            static_assert(
                    sizeof(gmx_ddzone_t) == c_ddzoneNumReals * sizeof(real),
                    "Here we expect gmx_ddzone_t to consist of c_ddzoneNumReals reals (only)");

            int numReals = numElementsInBuffer * c_ddzoneNumReals;
            ddSendrecv(dd,
                       d,
                       dddirBackward,
                       gmx::arrayRefFromArray(&buf_s[0].min0, numReals),
                       gmx::arrayRefFromArray(&buf_r[0].min0, numReals));

            rvec dh = { 0 };
            if (pulse > 0)
            {
                for (int d1 = d + 1; d1 < dd->ndim; d1++)
                {
                    /* Determine the decrease of maximum required
                     * communication height along d1 due to the distance along d,
                     * this avoids a lot of useless atom communication.
                     */
                    real dist_d = comm->cell_x1[dim] - buf_r[0].p1_0;

                    real c;
                    if (ddbox->tric_dir[dim])
                    {
                        /* c is the off-diagonal coupling between the cell planes
                         * along directions d and d1.
                         */
                        c = ddbox->v[dim][dd->dim[d1]][dim];
                    }
                    else
                    {
                        c = 0;
                    }
                    real det = (1 + c * c) * gmx::square(comm->systemInfo.cutoff) - dist_d * dist_d;
                    if (det > 0)
                    {
                        dh[d1] = comm->systemInfo.cutoff - (c * dist_d + std::sqrt(det)) / (1 + c * c);
                    }
                    else
                    {
                        /* A negative value signals out of range */
                        dh[d1] = -1;
                    }
                }
            }

            /* Accumulate the extremes over all pulses */
            for (int i = 0; i < numElementsInBuffer; i++)
            {
                if (pulse == 0)
                {
                    buf_e[i] = buf_r[i];
                }
                else
                {
                    if (receiveValidData)
                    {
                        buf_e[i].min0 = std::min(buf_e[i].min0, buf_r[i].min0);
                        buf_e[i].max1 = std::max(buf_e[i].max1, buf_r[i].max1);
                        buf_e[i].min1 = std::min(buf_e[i].min1, buf_r[i].min1);
                    }

                    int d1;
                    if (dd->ndim == 3 && d == 0 && i == numElementsInBuffer - 1)
                    {
                        d1 = 1;
                    }
                    else
                    {
                        d1 = d + 1;
                    }
                    if (receiveValidData && dh[d1] >= 0)
                    {
                        buf_e[i].mch0 = std::max(buf_e[i].mch0, buf_r[i].mch0 - dh[d1]);
                        buf_e[i].mch1 = std::max(buf_e[i].mch1, buf_r[i].mch1 - dh[d1]);
                    }
                }
                /* Copy the received buffer to the send buffer,
                 * to pass the data through with the next pulse.
                 */
                buf_s[i] = buf_r[i];
            }
            if (((applyPbc || dd->ci[dim] + numPulses < dd->numCells[dim]) && pulse == numPulses - 1)
                || (!applyPbc && dd->ci[dim] + 1 + pulse == dd->numCells[dim] - 1))
            {
                /* Store the extremes */
                int pos = 0;

                for (int d1 = d; d1 < dd->ndim - 1; d1++)
                {
                    extr_s[d1][1] = std::min(extr_s[d1][1], buf_e[pos].min0);
                    extr_s[d1][0] = std::max(extr_s[d1][0], buf_e[pos].max1);
                    extr_s[d1][2] = std::min(extr_s[d1][2], buf_e[pos].min1);
                    pos++;
                }

                if (d == 1 || (d == 0 && dd->ndim == 3))
                {
                    for (int i = d; i < 2; i++)
                    {
                        comm->zone_d2[1 - d][i] = buf_e[pos];
                        pos++;
                    }
                }
                if (d == 0)
                {
                    comm->zone_d1[1] = buf_e[pos];
                    pos++;
                }
            }
            else
            {
                if (d == 1 || (d == 0 && dd->ndim == 3))
                {
                    for (int i = d; i < 2; i++)
                    {
                        comm->zone_d2[1 - d][i].dataSet = 0;
                    }
                }
                if (d == 0)
                {
                    comm->zone_d1[1].dataSet = 0;
                }
            }
        }
    }

    if (dd->ndim >= 2)
    {
        int dim = dd->dim[1];
        for (int i = 0; i < 2; i++)
        {
            if (comm->zone_d1[i].dataSet != 0)
            {
                if (debug)
                {
                    print_ddzone(debug, 1, i, 0, &comm->zone_d1[i]);
                }
                cell_ns_x0[dim] = std::min(cell_ns_x0[dim], comm->zone_d1[i].min0);
                cell_ns_x1[dim] = std::max(cell_ns_x1[dim], comm->zone_d1[i].max1);
            }
        }
    }
    if (dd->ndim >= 3)
    {
        int dim = dd->dim[2];
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (comm->zone_d2[i][j].dataSet != 0)
                {
                    if (debug)
                    {
                        print_ddzone(debug, 2, i, j, &comm->zone_d2[i][j]);
                    }
                    cell_ns_x0[dim] = std::min(cell_ns_x0[dim], comm->zone_d2[i][j].min0);
                    cell_ns_x1[dim] = std::max(cell_ns_x1[dim], comm->zone_d2[i][j].max1);
                }
            }
        }
    }
    for (int d = 1; d < dd->ndim; d++)
    {
        cellsizes[d].fracLowerMax = extr_s[d - 1][0];
        cellsizes[d].fracUpperMin = extr_s[d - 1][1];
        if (debug)
        {
            fprintf(debug,
                    "Cell fraction d %d, max0 %f, min1 %f\n",
                    d,
                    cellsizes[d].fracLowerMax,
                    cellsizes[d].fracUpperMin);
        }
    }
}

//! Sets the atom-range for the home zone to \p dd->numHomeAtoms and all other zones empty
static void set_zones_numHomeAtoms(gmx_domdec_t* dd)
{
    for (int zone = 0; zone < dd->zones.numZones(); zone++)
    {
        dd->zones.setAtomRangeEnd(zone, dd->numHomeAtoms, true);
    }
}

//! Restore atom groups for the charge groups.
static void restoreAtomGroups(gmx_domdec_t* dd, const t_state* state)
{
    gmx::ArrayRef<const int> atomsState = state->cg_gl;

    std::vector<int>& globalAtomIndices = dd->globalAtomIndices;

    globalAtomIndices.resize(atomsState.size());

    /* Copy back the global charge group indices from state
     * and rebuild the local charge group to atom index.
     */
    for (gmx::Index i = 0; i < atomsState.ssize(); i++)
    {
        globalAtomIndices[i] = atomsState[i];
    }

    dd->numHomeAtoms = atomsState.size();
    dd->comm->atomRanges.setEnd(DDAtomRanges::Type::Home, atomsState.ssize());

    set_zones_numHomeAtoms(dd);
}

//! Sets the atom info structures.
static void dd_set_atominfo(gmx::ArrayRef<const int> index_gl, int atomStart, int atomEnd, t_forcerec* fr)
{
    if (fr != nullptr)
    {
        gmx::ArrayRef<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock =
                fr->atomInfoForEachMoleculeBlock;
        gmx::ArrayRef<int32_t> atomInfo = fr->atomInfo;

        const int gmx_unused numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Domdec);
#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int a = atomStart; a < atomEnd; a++)
        {
            atomInfo[a] = ddGetAtomInfo(atomInfoForEachMoleculeBlock, index_gl[a]);
        }
    }
}

//! Makes the mappings between global and local atom indices during DD repartioning.
static void make_dd_indices(gmx_domdec_t* dd, const int atomStart)
{
    const gmx::DomdecZones&  zones             = dd->zones;
    const int                numZones          = zones.numZones();
    gmx::ArrayRef<const int> globalAtomIndices = dd->globalAtomIndices;

    gmx_ga2la_t& ga2la = *dd->ga2la;

    GMX_ASSERT(*zones.atomRange(0).end() == dd->numHomeAtoms, "zones should be up to date");

    /* Make the local to global and global to local atom index */
    int a = atomStart;
    for (int zone = 0; zone < numZones; zone++)
    {
        int cg0;
        if (zone == 0)
        {
            cg0 = atomStart;
        }
        else
        {
            cg0 = *zones.atomRange(zone).begin();
        }
        int cg1    = *zones.atomRange(zone).end();
        int cg1_p1 = zones.directNeighborAtomRangeEnd(zone);

        for (int cg = cg0; cg < cg1; cg++)
        {
            int zone1 = zone;
            if (cg >= cg1_p1)
            {
                /* Signal that this cg is from more than one pulse away */
                zone1 += numZones;
            }
            int globalAtomIndex = globalAtomIndices[cg];
            ga2la.insert(globalAtomIndex, { a, zone1 });
            a++;
        }
    }
}

//! Checks whether global and local atom indices are consistent.
static void check_index_consistency(const gmx_domdec_t* dd, int natoms_sys, const char* where)
{
    int nerr = 0;

    const int numAtomsInZones = dd->comm->atomRanges.end(DDAtomRanges::Type::Zones);

    if (dd->comm->ddSettings.DD_debug > 1)
    {
        std::vector<int> have(natoms_sys);
        for (int a = 0; a < numAtomsInZones; a++)
        {
            int globalAtomIndex = dd->globalAtomIndices[a];
            if (have[globalAtomIndex] > 0)
            {
                fprintf(stderr,
                        "DD rank %d: global atom %d occurs twice: index %d and %d\n",
                        dd->rank,
                        globalAtomIndex + 1,
                        have[globalAtomIndex],
                        a + 1);
            }
            else
            {
                have[globalAtomIndex] = a + 1;
            }
        }
    }

    std::vector<int> have(numAtomsInZones);

    int ngl = 0;
    for (int i = 0; i < natoms_sys; i++)
    {
        if (const auto* entry = dd->ga2la->find(i))
        {
            const int a = entry->la;
            if (a >= numAtomsInZones)
            {
                fprintf(stderr,
                        "DD rank %d: global atom %d marked as local atom %d, which is larger than "
                        "nat_tot (%d)\n",
                        dd->rank,
                        i + 1,
                        a + 1,
                        numAtomsInZones);
                nerr++;
            }
            else
            {
                have[a] = 1;
                if (dd->globalAtomIndices[a] != i)
                {
                    fprintf(stderr,
                            "DD rank %d: global atom %d marked as local atom %d, which has global "
                            "atom index %d\n",
                            dd->rank,
                            i + 1,
                            a + 1,
                            dd->globalAtomIndices[a] + 1);
                    nerr++;
                }
            }
            ngl++;
        }
    }
    if (ngl != numAtomsInZones)
    {
        fprintf(stderr, "DD rank %d, %s: %d global atom indices, %d local atoms\n", dd->rank, where, ngl, numAtomsInZones);
    }
    for (int a = 0; a < numAtomsInZones; a++)
    {
        if (have[a] == 0)
        {
            fprintf(stderr,
                    "DD rank %d, %s: local atom %d, global %d has no global index\n",
                    dd->rank,
                    where,
                    a + 1,
                    dd->globalAtomIndices[a] + 1);
        }
    }

    if (nerr > 0)
    {
        gmx_fatal(FARGS, "DD rank %d, %s: %d atom(group) index inconsistencies", dd->rank, where, nerr);
    }
}

//! Clear all DD global state indices
static void clearDDStateIndices(gmx_domdec_t* dd, const bool keepLocalAtomIndices)
{
    gmx_ga2la_t& ga2la = *dd->ga2la;

    if (!keepLocalAtomIndices)
    {
        /* Clear the whole list without the overhead of searching */
        ga2la.clear(true);
    }
    else
    {
        const int numAtomsInZones = dd->comm->atomRanges.end(DDAtomRanges::Type::Zones);
        for (int i = 0; i < numAtomsInZones; i++)
        {
            ga2la.erase(dd->globalAtomIndices[i]);
        }
    }

    dd_clear_local_vsite_indices(dd);

    if (dd->constraints)
    {
        dd_clear_local_constraint_indices(dd);
    }
}

//! Return the duration of force calculations on this rank.
static float dd_force_load(gmx_domdec_comm_t* comm)
{
    float load;

    if (comm->ddSettings.eFlop)
    {
        load = comm->flop;
        if (comm->ddSettings.eFlop > 1)
        {
            load *= 1.0 + (comm->ddSettings.eFlop - 1) * (0.1 * rand() / RAND_MAX - 0.05);
        }
    }
    else
    {
        load = comm->cycl[ddCyclF];
        if (comm->cycl_n[ddCyclF] > 1)
        {
            /* Subtract the maximum of the last n cycle counts
             * to get rid of possible high counts due to other sources,
             * for instance system activity, that would otherwise
             * affect the dynamic load balancing.
             */
            load -= comm->cycl_max[ddCyclF];
        }

#if GMX_MPI
        if (comm->cycl_n[ddCyclWaitGPU] && comm->nrank_gpu_shared > 1)
        {
            float gpu_wait, gpu_wait_sum;

            gpu_wait = comm->cycl[ddCyclWaitGPU];
            if (comm->cycl_n[ddCyclF] > 1)
            {
                /* We should remove the WaitGPU time of the same MD step
                 * as the one with the maximum F time, since the F time
                 * and the wait time are not independent.
                 * Furthermore, the step for the max F time should be chosen
                 * the same on all ranks that share the same GPU.
                 * But to keep the code simple, we remove the average instead.
                 * The main reason for artificially long times at some steps
                 * is spurious CPU activity or MPI time, so we don't expect
                 * that changes in the GPU wait time matter a lot here.
                 */
                gpu_wait *= (comm->cycl_n[ddCyclF] - 1) / static_cast<float>(comm->cycl_n[ddCyclF]);
            }
            /* Sum the wait times over the ranks that share the same GPU */
            MPI_Allreduce(&gpu_wait, &gpu_wait_sum, 1, MPI_FLOAT, MPI_SUM, comm->mpi_comm_gpu_shared);
            /* Replace the wait time by the average over the ranks */
            load += -gpu_wait + gpu_wait_sum / comm->nrank_gpu_shared;
        }
#endif
    }

    return load;
}

//! Runs cell size checks and communicates the boundaries.
static void comm_dd_ns_cell_sizes(gmx_domdec_t* dd, gmx_ddbox_t* ddbox, rvec cell_ns_x0, rvec cell_ns_x1, int64_t step)
{
    gmx_domdec_comm_t* comm;
    int                dim_ind, dim;

    comm = dd->comm.get();

    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        dim = dd->dim[dim_ind];

        /* Without PBC we don't have restrictions on the outer cells */
        if (!(dim >= ddbox->npbcdim && (dd->ci[dim] == 0 || dd->ci[dim] == dd->numCells[dim] - 1))
            && isDlbOn(comm->dlbState)
            && (comm->cell_x1[dim] - comm->cell_x0[dim]) * ddbox->skew_fac[dim] < comm->cellsize_min[dim])
        {
            char buf[22];
            gmx_fatal(FARGS,
                      "step %s: The %c-size (%f) times the triclinic skew factor (%f) is smaller "
                      "than the smallest allowed cell size (%f) for domain decomposition grid cell "
                      "%d %d %d",
                      gmx_step_str(step, buf),
                      dim2char(dim),
                      comm->cell_x1[dim] - comm->cell_x0[dim],
                      ddbox->skew_fac[dim],
                      dd->comm->cellsize_min[dim],
                      dd->ci[XX],
                      dd->ci[YY],
                      dd->ci[ZZ]);
        }
    }

    if ((isDlbOn(dd->comm->dlbState) && dd->ndim > 1) || ddbox->nboundeddim < DIM)
    {
        /* Communicate the boundaries and update cell_ns_x0/1 */
        dd_move_cellx(dd, ddbox, cell_ns_x0, cell_ns_x1);
        if (isDlbOn(dd->comm->dlbState) && dd->ndim > 1)
        {
            gmx::check_grid_jump(step, dd, dd->comm->systemInfo.cutoff, ddbox, TRUE);
        }
    }
}

//! Compute and communicate to determine the load distribution across PP ranks.
static void get_load_distribution(gmx_domdec_t* dd, gmx_wallcycle* wcycle)
{
    gmx_domdec_comm_t* comm;
    domdec_load_t*     load;
    float              cell_frac = 0, sbuf[DD_NLOAD_MAX];
    gmx_bool           bSepPME;

    if (debug)
    {
        fprintf(debug, "get_load_distribution start\n");
    }

    wallcycle_start(wcycle, WallCycleCounter::DDCommLoad);

    comm = dd->comm.get();

    bSepPME = (dd->pme_nodeid >= 0);

    if (dd->ndim == 0 && bSepPME)
    {
        /* Without decomposition, but with PME nodes, we need the load */
        comm->load[0].mdf = comm->cycl[ddCyclPPduringPME];
        comm->load[0].pme = comm->cycl[ddCyclPME];
    }

    // Either we have DLB off, or we have it on and the array is large enough
    GMX_ASSERT(!isDlbOn(dd->comm->dlbState)
                       || static_cast<int>(dd->comm->cellsizesWithDlb.size()) == dd->ndim,
               "DLB cell sizes data not set up properly ");
    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        const int dim = dd->dim[d];
        /* Check if we participate in the communication in this dimension */
        if (d == dd->ndim - 1 || (dd->ci[dd->dim[d + 1]] == 0 && dd->ci[dd->dim[dd->ndim - 1]] == 0))
        {
            load = &comm->load[d];
            if (isDlbOn(dd->comm->dlbState))
            {
                cell_frac = comm->cellsizesWithDlb[d].fracUpper - comm->cellsizesWithDlb[d].fracLower;
            }
            int pos = 0;
            if (d == dd->ndim - 1)
            {
                sbuf[pos++] = dd_force_load(comm);
                sbuf[pos++] = sbuf[0];
                if (isDlbOn(dd->comm->dlbState))
                {
                    sbuf[pos++] = sbuf[0];
                    sbuf[pos++] = cell_frac;
                    if (d > 0)
                    {
                        sbuf[pos++] = comm->cellsizesWithDlb[d].fracLowerMax;
                        sbuf[pos++] = comm->cellsizesWithDlb[d].fracUpperMin;
                    }
                }
                if (bSepPME)
                {
                    sbuf[pos++] = comm->cycl[ddCyclPPduringPME];
                    sbuf[pos++] = comm->cycl[ddCyclPME];
                }
            }
            else
            {
                sbuf[pos++] = comm->load[d + 1].sum;
                sbuf[pos++] = comm->load[d + 1].max;
                if (isDlbOn(dd->comm->dlbState))
                {
                    sbuf[pos++] = comm->load[d + 1].sum_m;
                    sbuf[pos++] = comm->load[d + 1].cvol_min * cell_frac;
                    sbuf[pos++] = comm->load[d + 1].flags;
                    if (d > 0)
                    {
                        sbuf[pos++] = comm->cellsizesWithDlb[d].fracLowerMax;
                        sbuf[pos++] = comm->cellsizesWithDlb[d].fracUpperMin;
                    }
                }
                if (bSepPME)
                {
                    sbuf[pos++] = comm->load[d + 1].mdf;
                    sbuf[pos++] = comm->load[d + 1].pme;
                }
            }
            load->nload = pos;
            /* Communicate a row in DD direction d.
             * The communicators are setup such that the root always has rank 0.
             */
#if GMX_MPI
            MPI_Gather(sbuf,
                       load->nload * sizeof(float),
                       MPI_BYTE,
                       load->load.data(),
                       load->nload * sizeof(float),
                       MPI_BYTE,
                       0,
                       comm->mpi_comm_load[d]);
#endif
            if (dd->ci[dim] == dd->main_ci[dim])
            {
                /* We are the main along this row, process this row */
                RowCoordinator* rowCoordinator = nullptr;

                if (isDlbOn(comm->dlbState))
                {
                    rowCoordinator = comm->cellsizesWithDlb[d].rowCoordinator.get();
                }
                load->sum      = 0;
                load->max      = 0;
                load->sum_m    = 0;
                load->cvol_min = 1;
                load->flags    = 0;
                load->mdf      = 0;
                load->pme      = 0;
                int pos        = 0;
                for (int i = 0; i < dd->numCells[dim]; i++)
                {
                    load->sum += load->load[pos++];
                    load->max = std::max(load->max, load->load[pos]);
                    pos++;
                    if (isDlbOn(dd->comm->dlbState))
                    {
                        if (rowCoordinator->dlbIsLimited)
                        {
                            /* This direction could not be load balanced properly,
                             * therefore we need to use the maximum iso the average load.
                             */
                            load->sum_m = std::max(load->sum_m, load->load[pos]);
                        }
                        else
                        {
                            load->sum_m += load->load[pos];
                        }
                        pos++;
                        load->cvol_min = std::min(load->cvol_min, load->load[pos]);
                        pos++;
                        if (d < dd->ndim - 1)
                        {
                            load->flags = gmx::roundToInt(load->load[pos++]);
                        }
                        if (d > 0)
                        {
                            rowCoordinator->bounds[i].cellFracLowerMax = load->load[pos++];
                            rowCoordinator->bounds[i].cellFracUpperMin = load->load[pos++];
                        }
                    }
                    if (bSepPME)
                    {
                        load->mdf = std::max(load->mdf, load->load[pos]);
                        pos++;
                        load->pme = std::max(load->pme, load->load[pos]);
                        pos++;
                    }
                }
                if (isDlbOn(comm->dlbState) && rowCoordinator->dlbIsLimited)
                {
                    load->sum_m *= dd->numCells[dim];
                    load->flags |= (1 << d);
                }
            }
        }
    }

    if (DDMAIN(dd))
    {
        comm->nload += dd_load_count(comm);
        comm->load_step += comm->cycl[ddCyclStep];
        comm->load_sum += comm->load[0].sum;
        comm->load_max += comm->load[0].max;
        if (isDlbOn(comm->dlbState))
        {
            for (int d = 0; d < dd->ndim; d++)
            {
                if (comm->load[0].flags & (1 << d))
                {
                    comm->load_lim[d]++;
                }
            }
        }
        if (bSepPME)
        {
            comm->load_mdf += comm->load[0].mdf;
            comm->load_pme += comm->load[0].pme;
        }
    }

    wallcycle_stop(wcycle, WallCycleCounter::DDCommLoad);

    if (debug)
    {
        fprintf(debug, "get_load_distribution finished\n");
    }
}

/*! \brief Return the relative performance loss on the total run time
 * due to the force calculation load imbalance. */
static float dd_force_load_fraction(gmx_domdec_t* dd)
{
    if (dd->comm->nload > 0 && dd->comm->load_step > 0)
    {
        return dd->comm->load_sum / (dd->comm->load_step * dd->nnodes);
    }
    else
    {
        return 0;
    }
}

/*! \brief Return the relative performance loss on the total run time
 * due to the force calculation load imbalance. */
static float dd_force_imb_perf_loss(gmx_domdec_t* dd)
{
    if (dd->comm->nload > 0 && dd->comm->load_step > 0)
    {
        return (dd->comm->load_max * dd->nnodes - dd->comm->load_sum) / (dd->comm->load_step * dd->nnodes);
    }
    else
    {
        return 0;
    }
}

//! Print load-balance report e.g. at the end of a run.
static void print_dd_load_av(FILE* fplog, gmx_domdec_t* dd)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    /* Only the main rank prints loads and only if we measured loads */
    if (!DDMAIN(dd) || comm->nload == 0)
    {
        return;
    }

    char buf[STRLEN];
    int  numPpRanks  = dd->nnodes;
    int  numPmeRanks = (comm->ddRankSetup.usePmeOnlyRanks ? comm->ddRankSetup.numRanksDoingPme : 0);
    int  numRanks    = numPpRanks + numPmeRanks;
    float lossFraction = 0;

    /* Print the average load imbalance and performance loss */
    if (dd->nnodes > 1 && comm->load_sum > 0)
    {
        float imbalance = comm->load_max * numPpRanks / comm->load_sum - 1;
        lossFraction    = dd_force_imb_perf_loss(dd);

        std::string msg = "\nDynamic load balancing report:\n";
        std::string dlbStateStr;

        switch (dd->comm->dlbState)
        {
            case DlbState::offUser:
                dlbStateStr = "DLB was off during the run per user request.";
                break;
            case DlbState::offForever:
                /* Currectly this can happen due to performance loss observed, cell size
                 * limitations or incompatibility with other settings observed during
                 * determineInitialDlbState(). */
                dlbStateStr = "DLB got disabled because it was unsuitable to use.";
                break;
            case DlbState::offCanTurnOn:
                dlbStateStr = "DLB was off during the run due to low measured imbalance.";
                break;
            case DlbState::offTemporarilyLocked:
                dlbStateStr =
                        "DLB was locked at the end of the run due to unfinished PP-PME "
                        "balancing.";
                break;
            case DlbState::onCanTurnOff:
                dlbStateStr = "DLB was turned on during the run due to measured imbalance.";
                break;
            case DlbState::onUser:
                dlbStateStr = "DLB was permanently on during the run per user request.";
                break;
            default: GMX_ASSERT(false, "Undocumented DLB state");
        }

        msg += " " + dlbStateStr + "\n";
        msg += gmx::formatString(" Average load imbalance: %.1f%%.\n", imbalance * 100);
        msg += gmx::formatString(
                " The balanceable part of the MD step is %d%%, load imbalance is computed from "
                "this.\n",
                gmx::roundToInt(dd_force_load_fraction(dd) * 100));
        msg += gmx::formatString(
                " Part of the total run time spent waiting due to load imbalance: %.1f%%.\n",
                lossFraction * 100);
        fprintf(fplog, "%s", msg.c_str());
        fprintf(stderr, "\n%s", msg.c_str());
    }

    /* Print during what percentage of steps the  load balancing was limited */
    bool dlbWasLimited = false;
    if (isDlbOn(comm->dlbState))
    {
        sprintf(buf, " Steps where the load balancing was limited by -rdd, -rcon and/or -dds:");
        for (int d = 0; d < dd->ndim; d++)
        {
            int limitPercentage = (200 * comm->load_lim[d] + 1) / (2 * comm->nload);
            sprintf(buf + strlen(buf), " %c %d %%", dim2char(dd->dim[d]), limitPercentage);
            if (limitPercentage >= 50)
            {
                dlbWasLimited = true;
            }
        }
        sprintf(buf + strlen(buf), "\n");
        fprintf(fplog, "%s", buf);
        fprintf(stderr, "%s", buf);
    }

    /* Print the performance loss due to separate PME - PP rank imbalance */
    float lossFractionPme = 0;
    if (numPmeRanks > 0 && comm->load_mdf > 0 && comm->load_step > 0)
    {
        float pmeForceRatio = comm->load_pme / comm->load_mdf;
        lossFractionPme     = (comm->load_pme - comm->load_mdf) / comm->load_step;
        if (lossFractionPme <= 0)
        {
            lossFractionPme *= numPmeRanks / static_cast<float>(numRanks);
        }
        else
        {
            lossFractionPme *= numPpRanks / static_cast<float>(numRanks);
        }
        sprintf(buf, " Average PME mesh/force load: %5.3f\n", pmeForceRatio);
        fprintf(fplog, "%s", buf);
        fprintf(stderr, "%s", buf);
        sprintf(buf,
                " Part of the total run time spent waiting due to PP/PME imbalance: %.1f %%\n",
                std::fabs(lossFractionPme) * 100);
        fprintf(fplog, "%s", buf);
        fprintf(stderr, "%s", buf);
    }
    fprintf(fplog, "\n");
    fprintf(stderr, "\n");

    if ((lossFraction >= DD_PERF_LOSS_WARN) && (dd->comm->dlbState != DlbState::offTemporarilyLocked))
    {
        std::string message = gmx::formatString(
                "NOTE: %.1f %% of the available CPU time was lost due to load imbalance\n"
                "      in the domain decomposition.\n",
                lossFraction * 100);

        bool hadSuggestion = false;
        if (dd->comm->dlbState == DlbState::offUser)
        {
            message += "      You might want to allow dynamic load balancing (option -dlb auto.)\n";
            hadSuggestion = true;
        }
        else if (dd->comm->dlbState == DlbState::offCanTurnOn)
        {
            message +=
                    "      Dynamic load balancing was automatically disabled, but it might be "
                    "beneficial to manually turn it on (option -dlb yes.)\n";
            hadSuggestion = true;
        }
        else if (dlbWasLimited)
        {
            message +=
                    "      You might want to decrease the cell size limit (options -rdd, -rcon "
                    "and/or -dds).\n";
            hadSuggestion = true;
        }
        message += gmx::formatString(
                "      You can %sconsider manually changing the decomposition (option -dd);\n"
                "      e.g. by using fewer domains along the box dimension in which there is\n"
                "      considerable inhomogeneity in the simulated system.",
                hadSuggestion ? "also " : "");

        fprintf(fplog, "%s\n", message.c_str());
        fprintf(stderr, "%s\n", message.c_str());
    }
    if (numPmeRanks > 0 && std::fabs(lossFractionPme) >= DD_PERF_LOSS_WARN)
    {
        sprintf(buf,
                "NOTE: %.1f %% performance was lost because the PME ranks\n"
                "      had %s work to do than the PP ranks.\n"
                "      You might want to %s the number of PME ranks\n"
                "      or %s the cut-off and the grid spacing.\n",
                std::fabs(lossFractionPme * 100),
                (lossFractionPme < 0) ? "less" : "more",
                (lossFractionPme < 0) ? "decrease" : "increase",
                (lossFractionPme < 0) ? "decrease" : "increase");
        fprintf(fplog, "%s\n", buf);
        fprintf(stderr, "%s\n", buf);
    }
}

//! Return the minimum communication volume.
static float dd_vol_min(gmx_domdec_t* dd)
{
    return dd->comm->load[0].cvol_min * dd->nnodes;
}

//! Return the DD load flags.
static int dd_load_flags(gmx_domdec_t* dd)
{
    return dd->comm->load[0].flags;
}

//! Return the reported load imbalance in force calculations.
static float dd_f_imbal(gmx_domdec_t* dd)
{
    if (dd->comm->load[0].sum > 0)
    {
        return dd->comm->load[0].max * dd->nnodes / dd->comm->load[0].sum - 1.0F;
    }
    else
    {
        /* Something is wrong in the cycle counting, report no load imbalance */
        return 0.0F;
    }
}

//! Returns DD load balance report.
static std::string dd_print_load(gmx_domdec_t* dd, int64_t step)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         log(&stream);

    int flags = dd_load_flags(dd);
    if (flags)
    {
        log.writeString("DD  load balancing is limited by minimum cell size in dimension");
        for (int d = 0; d < dd->ndim; d++)
        {
            if (flags & (1 << d))
            {
                log.writeStringFormatted(" %c", dim2char(dd->dim[d]));
            }
        }
        log.ensureLineBreak();
    }
    log.writeString("DD  step " + gmx::toString(step));
    if (isDlbOn(dd->comm->dlbState))
    {
        log.writeStringFormatted("  vol min/aver %5.3f%c", dd_vol_min(dd), flags ? '!' : ' ');
    }
    if (dd->nnodes > 1)
    {
        log.writeStringFormatted(" load imb.: force %4.1f%%", dd_f_imbal(dd) * 100);
    }
    if (dd->comm->cycl_n[ddCyclPME])
    {
        log.writeStringFormatted("  pme mesh/force %5.3f", dd_pme_f_ratio(dd));
    }
    log.ensureLineBreak();
    return stream.toString();
}

//! Prints DD load balance report in mdrun verbose mode.
static void dd_print_load_verbose(gmx_domdec_t* dd)
{
    if (isDlbOn(dd->comm->dlbState))
    {
        fprintf(stderr, "vol %4.2f%c ", dd_vol_min(dd), dd_load_flags(dd) ? '!' : ' ');
    }
    if (dd->nnodes > 1)
    {
        fprintf(stderr, "imb F %2d%% ", gmx::roundToInt(dd_f_imbal(dd) * 100));
    }
    if (dd->comm->cycl_n[ddCyclPME])
    {
        fprintf(stderr, "pme/F %4.2f ", dd_pme_f_ratio(dd));
    }
}

//! Turns on dynamic load balancing if possible and needed.
static void turn_on_dlb(const gmx::MDLogger& mdlog, gmx_domdec_t* dd, int64_t step)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    real cellsize_min = comm->cellsize_min[dd->dim[0]];
    for (int d = 1; d < dd->ndim; d++)
    {
        cellsize_min = std::min(cellsize_min, comm->cellsize_min[dd->dim[d]]);
    }

    /* Turn off DLB if we're too close to the cell size limit. */
    if (cellsize_min < comm->cellsize_limit * 1.05)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "step %s Measured %.1f %% performance loss due to load imbalance, "
                        "but the minimum cell size is smaller than 1.05 times the cell size limit. "
                        "Will no longer try dynamic load balancing.",
                        gmx::toString(step).c_str(),
                        dd_force_imb_perf_loss(dd) * 100);

        comm->dlbState = DlbState::offForever;
        return;
    }

    GMX_LOG(mdlog.info)
            .appendTextFormatted(
                    "step %s Turning on dynamic load balancing, because the performance loss due "
                    "to load imbalance is %.1f %%.",
                    gmx::toString(step).c_str(),
                    dd_force_imb_perf_loss(dd) * 100);
    comm->dlbState = DlbState::onCanTurnOff;

    /* Store the non-DLB performance, so we can check if DLB actually
     * improves performance.
     */
    GMX_RELEASE_ASSERT(comm->cycl_n[ddCyclStep] > 0,
                       "When we turned on DLB, we should have measured cycles");
    comm->cyclesPerStepBeforeDLB = comm->cycl[ddCyclStep] / comm->cycl_n[ddCyclStep];

    set_dlb_limits(dd);

    /* We can set the required cell size info here,
     * so we do not need to communicate this.
     * The grid is completely uniform.
     */
    for (int d = 0; d < dd->ndim; d++)
    {
        RowCoordinator* rowCoordinator = comm->cellsizesWithDlb[d].rowCoordinator.get();

        if (rowCoordinator)
        {
            comm->load[d].sum_m = comm->load[d].sum;

            int nc = dd->numCells[dd->dim[d]];
            for (int i = 0; i < nc; i++)
            {
                rowCoordinator->cellFrac[i] = i / static_cast<real>(nc);
                if (d > 0)
                {
                    rowCoordinator->bounds[i].cellFracLowerMax = i / static_cast<real>(nc);
                    rowCoordinator->bounds[i].cellFracUpperMin = (i + 1) / static_cast<real>(nc);
                }
            }
            rowCoordinator->cellFrac[nc] = 1.0;
        }
    }
}

//! Turns off dynamic load balancing (but leave it able to turn back on).
static void turn_off_dlb(const gmx::MDLogger& mdlog, gmx_domdec_t* dd, int64_t step)
{
    GMX_LOG(mdlog.info)
            .appendText(
                    "step " + gmx::toString(step)
                    + " Turning off dynamic load balancing, because it is degrading performance.");
    dd->comm->dlbState                     = DlbState::offCanTurnOn;
    dd->comm->haveTurnedOffDlb             = true;
    dd->comm->ddPartioningCountFirstDlbOff = dd->ddp_count;
}

//! Turns off dynamic load balancing permanently.
static void turn_off_dlb_forever(const gmx::MDLogger& mdlog, gmx_domdec_t* dd, int64_t step)
{
    GMX_RELEASE_ASSERT(dd->comm->dlbState == DlbState::offCanTurnOn,
                       "Can only turn off DLB forever when it was in the can-turn-on state");
    GMX_LOG(mdlog.info)
            .appendText(
                    "step " + gmx::toString(step)
                    + " Will no longer try dynamic load balancing, as it degraded performance.");
    dd->comm->dlbState = DlbState::offForever;
}

void set_dd_dlb_max_cutoff(t_commrec* cr, real cutoff)
{
    gmx_domdec_comm_t* comm;

    comm = cr->dd->comm.get();

    /* Turn on the DLB limiting (might have been on already) */
    comm->bPMELoadBalDLBLimits = TRUE;

    /* Change the cut-off limit */
    comm->PMELoadBal_max_cutoff = cutoff;

    if (debug)
    {
        fprintf(debug,
                "PME load balancing set a limit to the DLB staggering such that a %f cut-off will "
                "continue to fit\n",
                comm->PMELoadBal_max_cutoff);
    }
}

/*! \brief Merge received atoms for one pulse and zone into the atom buffers
 *
 * \param[in]     numZones  The number of zones to apply the merging to
 * \param[in,out] cd      The communication setup for the pulses along the current dimension
 * \param[in]     pulse   The index of the current pulse
 * \param[in,out] zones   The DD zone information in which the atom ranges will be updated
 * \param[in,out] index_gl        The global atom indices
 * \param[in]     recv_i  List of received atom indices for this pulse
 * \param[in,out] x       The home + halo coordinate buffer
 * \param[in]     recv_vr Buffer with received coordinates
 * \param[in]     atomInfoForEachMoleculeBlock  List of atom information for molecule blocks
 * \param[in,out] atomInfo  List of home + halo atom information
 */
static void mergeAtomBuffers(const int                                       numZones,
                             gmx_domdec_comm_dim_t*                          cd,
                             const int                                       pulse,
                             gmx::DomdecZones*                               zones,
                             gmx::ArrayRef<int>                              index_gl,
                             const int*                                      recv_i,
                             gmx::ArrayRef<gmx::RVec>                        x,
                             gmx::ArrayRef<const gmx::RVec>                  recv_vr,
                             gmx::ArrayRef<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock,
                             gmx::ArrayRef<int32_t>                          atomInfo)
{
    GMX_ASSERT(zones->numZones() >= 2 * numZones, "zones should contain at least 2*numZones zones");

    const gmx_domdec_ind_t& ind = cd->ind[pulse];

    /* First correct the already stored data */
    int shift = ind.nrecv[numZones];
    for (int zone = numZones - 1; zone >= 0; zone--)
    {
        shift -= ind.nrecv[zone];
        if (shift > 0)
        {
            /* Move the atoms present from previous grid pulses */
            const int atomStart = *zones->atomRange(numZones + zone).begin();
            const int atomEnd   = *zones->atomRange(numZones + zone).end();
            for (int a = atomEnd - 1; a >= atomStart; a--)
            {
                index_gl[a + shift] = index_gl[a];
                x[a + shift]        = x[a];
                atomInfo[a + shift] = atomInfo[a];
            }
            /* Correct the already stored send indices for the shift */
            for (int p = 1; p <= pulse; p++)
            {
                gmx_domdec_ind_t& ind_p          = cd->ind[p];
                int               pulseAtomStart = 0;
                for (int z = 0; z < zone; z++)
                {
                    pulseAtomStart += ind_p.nsend[z];
                }
                int pulseAtomEnd = pulseAtomStart + ind_p.nsend[zone];
                for (int a = pulseAtomStart; a < pulseAtomEnd; a++)
                {
                    ind_p.index[a] += shift;
                }
            }
        }
    }

    /* Merge in the communicated buffers */
    shift       = 0;
    int atomSrc = 0;
    for (int zone = 0; zone < numZones; zone++)
    {
        int atomDest = *zones->atomRange(numZones + zone).end() + shift;
        for (int a = 0; a < ind.nrecv[zone]; a++)
        {
            /* Copy this atom from the buffer */
            index_gl[atomDest] = recv_i[atomSrc];
            x[atomDest]        = recv_vr[atomSrc];
            /* Copy information */
            atomInfo[atomDest] = ddGetAtomInfo(atomInfoForEachMoleculeBlock, index_gl[atomDest]);
            atomSrc++;
            atomDest++;
        }
        shift += ind.nrecv[zone];
        zones->setAtomRangeEnd(numZones + zone, atomDest, false);
    }
}

//! Makes a range partitioning for the atom groups wthin a cell
static void make_cell2at_index(gmx_domdec_comm_dim_t* cd, int nzone, int atomGroupStart)
{
    /* Store the atom block boundaries for easy copying of communication buffers
     */
    int g = atomGroupStart;
    for (int zone = 0; zone < nzone; zone++)
    {
        for (gmx_domdec_ind_t& ind : cd->ind)
        {
            ind.cell2at0[zone] = g;
            g += ind.nrecv[zone];
            ind.cell2at1[zone] = g;
        }
    }
}

//! Returns whether a link is missing.
static bool missing_link(const gmx::ListOfLists<int>& link, const int globalAtomIndex, const gmx_ga2la_t& ga2la)
{
    return std::any_of(link[globalAtomIndex].begin(), link[globalAtomIndex].end(), [&](const int a) {
        return ga2la.findHome(a) == nullptr;
    });
}

//! Domain corners for communication, a maximum of 4 i-zones see a j domain
typedef struct
{
    //! The corners for the non-bonded communication.
    real c[DIM][4];
    //! Corner for rounding.
    real cr0;
    //! Corners for rounding.
    real cr1[4];
    //! Corners for bounded communication.
    real bc[DIM];
    //! Corner for rounding for bonded communication.
    real bcr1;
} dd_corners_t;

//! Determine the corners of the domain(s) we are communicating with.
static void set_dd_corners(const gmx_domdec_t* dd, int dim0, int dim1, int dim2, gmx_bool bDistMB, dd_corners_t* c)
{
    const gmx_domdec_comm_t* comm = dd->comm.get();

    const gmx::DomdecZones& zones = dd->zones;

    /* Keep the compiler happy */
    c->cr0  = 0;
    c->bcr1 = 0;

    /* The first dimension is equal for all cells */
    c->c[0][0] = comm->cell_x0[dim0];
    if (bDistMB)
    {
        c->bc[0] = c->c[0][0];
    }
    if (dd->ndim >= 2)
    {
        dim1 = dd->dim[1];
        /* This cell row is only seen from the first row */
        c->c[1][0] = comm->cell_x0[dim1];
        /* All rows can see this row */
        c->c[1][1] = comm->cell_x0[dim1];
        if (isDlbOn(dd->comm->dlbState))
        {
            c->c[1][1] = std::max(comm->cell_x0[dim1], comm->zone_d1[1].mch0);
            if (bDistMB)
            {
                /* For the multi-body distance we need the maximum */
                c->bc[1] = std::max(comm->cell_x0[dim1], comm->zone_d1[1].p1_0);
            }
        }
        /* Set the upper-right corner for rounding */
        c->cr0 = comm->cell_x1[dim0];

        if (dd->ndim >= 3)
        {
            dim2 = dd->dim[2];
            for (int j = 0; j < 4; j++)
            {
                c->c[2][j] = comm->cell_x0[dim2];
            }
            if (isDlbOn(dd->comm->dlbState))
            {
                /* Use the maximum of the i-cells that see a j-cell */
                for (int iz = 0; iz < zones.numIZones(); iz++)
                {
                    const auto& jZoneRange = zones.jZoneRange(iz);

                    for (int jZone : jZoneRange)
                    {
                        if (jZone >= 4)
                        {
                            c->c[2][jZone - 4] = std::max(
                                    c->c[2][jZone - 4],
                                    comm->zone_d2[zones.shift(iz)[dim0]][zones.shift(iz)[dim1]].mch0);
                        }
                    }
                }
                if (bDistMB)
                {
                    /* For the multi-body distance we need the maximum */
                    c->bc[2] = comm->cell_x0[dim2];
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            c->bc[2] = std::max(c->bc[2], comm->zone_d2[i][j].p1_0);
                        }
                    }
                }
            }

            /* Set the upper-right corner for rounding */
            /* Cell (0,0,0) and cell (1,0,0) can see cell 4 (0,1,1)
             * Only cell (0,0,0) can see cell 7 (1,1,1)
             */
            c->cr1[0] = comm->cell_x1[dim1];
            c->cr1[3] = comm->cell_x1[dim1];
            if (isDlbOn(dd->comm->dlbState))
            {
                c->cr1[0] = std::max(comm->cell_x1[dim1], comm->zone_d1[1].mch1);
                if (bDistMB)
                {
                    /* For the multi-body distance we need the maximum */
                    c->bcr1 = std::max(comm->cell_x1[dim1], comm->zone_d1[1].p1_1);
                }
            }
        }
    }
}

/*! \brief Add the atom groups and coordinates we need to send in this
 * pulse from this zone to \p localAtomGroups and \p work. */
static void get_zone_pulse_groups(gmx_domdec_t*                  dd,
                                  int                            zonei,
                                  int                            zone,
                                  int                            cg0,
                                  int                            cg1,
                                  gmx::ArrayRef<const int>       globalAtomIndices,
                                  int                            dim,
                                  int                            dim_ind,
                                  int                            dim0,
                                  int                            dim1,
                                  int                            dim2,
                                  real                           r_comm2,
                                  real                           r_bcomm2,
                                  matrix                         box,
                                  bool                           distanceIsTriclinic,
                                  rvec*                          normal,
                                  real                           skew_fac2_d,
                                  real                           skew_fac_01,
                                  rvec*                          v_d,
                                  rvec*                          v_0,
                                  rvec*                          v_1,
                                  const dd_corners_t*            c,
                                  const rvec                     sf2_round,
                                  gmx_bool                       bDistBonded,
                                  gmx_bool                       bBondComm,
                                  gmx_bool                       bDist2B,
                                  gmx_bool                       bDistMB,
                                  gmx::ArrayRef<const gmx::RVec> coordinates,
                                  gmx::ArrayRef<const int32_t>   atomInfo,
                                  std::vector<int>*              localAtomGroups,
                                  dd_comm_setup_work_t*          work)
{
    gmx_domdec_comm_t* comm;
    gmx_bool           bScrew;
    gmx_bool           bDistMB_pulse;
    int                cg, i;
    real               r2, rb2, r, tric_sh;
    rvec               rn, rb;
    int                dimd;
    int                nsend_z, nat;

    comm = dd->comm.get();

    bScrew = (dd->unitCellInfo.haveScrewPBC && dim == XX);

    bDistMB_pulse = (bDistMB && bDistBonded);

    /* Unpack the work data */
    std::vector<int>&       ibuf = work->atomGroupBuffer;
    std::vector<gmx::RVec>& vbuf = work->positionBuffer;
    nsend_z                      = 0;
    nat                          = work->nat;

    for (cg = cg0; cg < cg1; cg++)
    {
        r2  = 0;
        rb2 = 0;
        if (!distanceIsTriclinic)
        {
            /* Rectangular direction, easy */
            r = coordinates[cg][dim] - c->c[dim_ind][zone];
            if (r > 0)
            {
                r2 += r * r;
            }
            if (bDistMB_pulse)
            {
                r = coordinates[cg][dim] - c->bc[dim_ind];
                if (r > 0)
                {
                    rb2 += r * r;
                }
            }
            /* Rounding gives at most a 16% reduction
             * in communicated atoms
             */
            if (dim_ind >= 1 && (zonei == 1 || zonei == 2))
            {
                r = coordinates[cg][dim0] - c->cr0;
                /* This is the first dimension, so always r >= 0 */
                r2 += r * r;
                if (bDistMB_pulse)
                {
                    rb2 += r * r;
                }
            }
            if (dim_ind == 2 && (zonei == 2 || zonei == 3))
            {
                r = coordinates[cg][dim1] - c->cr1[zone];
                if (r > 0)
                {
                    r2 += r * r;
                }
                if (bDistMB_pulse)
                {
                    r = coordinates[cg][dim1] - c->bcr1;
                    if (r > 0)
                    {
                        rb2 += r * r;
                    }
                }
            }
        }
        else
        {
            /* Triclinic direction, more complicated */
            clear_rvec(rn);
            clear_rvec(rb);
            /* Rounding, conservative as the skew_fac multiplication
             * will slightly underestimate the distance.
             */
            if (dim_ind >= 1 && (zonei == 1 || zonei == 2))
            {
                rn[dim0] = coordinates[cg][dim0] - c->cr0;
                for (i = dim0 + 1; i < DIM; i++)
                {
                    rn[dim0] -= coordinates[cg][i] * v_0[i][dim0];
                }
                r2 = rn[dim0] * rn[dim0] * sf2_round[dim0];
                if (bDistMB_pulse)
                {
                    rb[dim0] = rn[dim0];
                    rb2      = r2;
                }
                /* Take care that the cell planes along dim0 might not
                 * be orthogonal to those along dim1 and dim2.
                 */
                for (i = 1; i <= dim_ind; i++)
                {
                    dimd = dd->dim[i];
                    if (normal[dim0][dimd] > 0)
                    {
                        rn[dimd] -= rn[dim0] * normal[dim0][dimd];
                        if (bDistMB_pulse)
                        {
                            rb[dimd] -= rb[dim0] * normal[dim0][dimd];
                        }
                    }
                }
            }
            if (dim_ind == 2 && (zonei == 2 || zonei == 3))
            {
                GMX_ASSERT(dim1 >= 0 && dim1 < DIM, "Must have a valid dimension index");
                rn[dim1] += coordinates[cg][dim1] - c->cr1[zone];
                tric_sh = 0;
                for (i = dim1 + 1; i < DIM; i++)
                {
                    tric_sh -= coordinates[cg][i] * v_1[i][dim1];
                }
                rn[dim1] += tric_sh;
                if (rn[dim1] > 0)
                {
                    r2 += rn[dim1] * rn[dim1] * sf2_round[dim1];
                    /* Take care of coupling of the distances
                     * to the planes along dim0 and dim1 through dim2.
                     */
                    r2 -= rn[dim0] * rn[dim1] * skew_fac_01;
                    /* Take care that the cell planes along dim1
                     * might not be orthogonal to that along dim2.
                     */
                    if (normal[dim1][dim2] > 0)
                    {
                        rn[dim2] -= rn[dim1] * normal[dim1][dim2];
                    }
                }
                if (bDistMB_pulse)
                {
                    rb[dim1] += coordinates[cg][dim1] - c->bcr1 + tric_sh;
                    if (rb[dim1] > 0)
                    {
                        rb2 += rb[dim1] * rb[dim1] * sf2_round[dim1];
                        /* Take care of coupling of the distances
                         * to the planes along dim0 and dim1 through dim2.
                         */
                        rb2 -= rb[dim0] * rb[dim1] * skew_fac_01;
                        /* Take care that the cell planes along dim1
                         * might not be orthogonal to that along dim2.
                         */
                        if (normal[dim1][dim2] > 0)
                        {
                            rb[dim2] -= rb[dim1] * normal[dim1][dim2];
                        }
                    }
                }
            }
            /* The distance along the communication direction */
            rn[dim] += coordinates[cg][dim] - c->c[dim_ind][zone];
            tric_sh = 0;
            for (i = dim + 1; i < DIM; i++)
            {
                tric_sh -= coordinates[cg][i] * v_d[i][dim];
            }
            rn[dim] += tric_sh;
            if (rn[dim] > 0)
            {
                r2 += rn[dim] * rn[dim] * skew_fac2_d;
                /* Take care of coupling of the distances
                 * to the planes along dim0 and dim1 through dim2.
                 */
                if (dim_ind == 1 && zonei == 1)
                {
                    r2 -= rn[dim0] * rn[dim] * skew_fac_01;
                }
            }
            if (bDistMB_pulse)
            {
                clear_rvec(rb);
                GMX_ASSERT(dim >= 0 && dim < DIM, "Must have a valid dimension index");
                rb[dim] += coordinates[cg][dim] - c->bc[dim_ind] + tric_sh;
                if (rb[dim] > 0)
                {
                    rb2 += rb[dim] * rb[dim] * skew_fac2_d;
                    /* Take care of coupling of the distances
                     * to the planes along dim0 and dim1 through dim2.
                     */
                    if (dim_ind == 1 && zonei == 1)
                    {
                        rb2 -= rb[dim0] * rb[dim] * skew_fac_01;
                    }
                }
            }
        }

        if (r2 < r_comm2
            || (bDistBonded && ((bDistMB && rb2 < r_bcomm2) || (bDist2B && r2 < r_bcomm2))
                && (!bBondComm
                    || ((atomInfo[cg] & gmx::sc_atomInfo_BondCommunication)
                        && missing_link(*comm->bondedLinks, globalAtomIndices[cg], *dd->ga2la)))))
        {
            /* Store the local and global atom group indices and position */
            localAtomGroups->push_back(cg);
            ibuf.push_back(globalAtomIndices[cg]);
            nsend_z++;

            rvec posPbc;
            if (dd->ci[dim] == 0)
            {
                /* Correct coordinates for pbc */
                rvec_add(coordinates[cg], box[dim], posPbc);
                if (bScrew)
                {
                    posPbc[YY] = box[YY][YY] - posPbc[YY];
                    posPbc[ZZ] = box[ZZ][ZZ] - posPbc[ZZ];
                }
            }
            else
            {
                copy_rvec(coordinates[cg], posPbc);
            }
            vbuf.emplace_back(posPbc[XX], posPbc[YY], posPbc[ZZ]);

            nat += 1;
        }
    }

    work->nat        = nat;
    work->nsend_zone = nsend_z;
}

//! Clear data.
static void clearCommSetupData(dd_comm_setup_work_t* work)
{
    work->localAtomGroupBuffer.clear();
    work->atomGroupBuffer.clear();
    work->positionBuffer.clear();
    work->nat        = 0;
    work->nsend_zone = 0;
}

//! Prepare DD communication.
static void setup_dd_communication(gmx_domdec_t* dd, matrix box, gmx_ddbox_t* ddbox, t_forcerec* fr, t_state* state)
{
    int                    dim_ind, dim, dim0, dim1, dim2, dimd, nat_tot;
    int                    nzone, nzone_send, zone, zonei;
    int                    c;
    int                    pos_cg;
    gmx_domdec_comm_t*     comm;
    gmx_domdec_comm_dim_t* cd;
    gmx_bool               bBondComm, bDist2B, bDistMB, bDistBonded;
    dd_corners_t           corners;
    rvec *                 normal, *v_d, *v_0 = nullptr, *v_1 = nullptr;
    real                   skew_fac2_d, skew_fac_01;
    rvec                   sf2_round;

    if (debug)
    {
        fprintf(debug, "Setting up DD communication\n");
    }

    comm = dd->comm.get();

    if (comm->dth.empty())
    {
        /* Initialize the thread data.
         * This can not be done in init_domain_decomposition,
         * as the numbers of threads is determined later.
         */
        int numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Domdec);
        comm->dth.resize(numThreads);
    }

    bBondComm = comm->systemInfo.filterBondedCommunication;

    /* Do we need to determine extra distances for multi-body bondeds? */
    bDistMB = (comm->systemInfo.haveInterDomainMultiBodyBondeds && isDlbOn(dd->comm->dlbState)
               && dd->ndim > 1);

    /* Do we need to determine extra distances for only two-body bondeds? */
    bDist2B = (bBondComm && !bDistMB);

    const real r_comm2 =
            gmx::square(domainToDomainIntoAtomToDomainCutoff(comm->systemInfo, comm->systemInfo.cutoff));
    const real r_bcomm2 =
            gmx::square(domainToDomainIntoAtomToDomainCutoff(comm->systemInfo, comm->cutoff_mbody));

    if (debug)
    {
        fprintf(debug, "bBondComm %s, r_bc %f\n", gmx::boolToString(bBondComm), std::sqrt(r_bcomm2));
    }

    dim0 = dd->dim[0];
    dim1 = (dd->ndim >= 2 ? dd->dim[1] : -1);
    dim2 = (dd->ndim >= 3 ? dd->dim[2] : -1);

    set_dd_corners(dd, dim0, dim1, dim2, bDistMB, &corners);

    /* Triclinic stuff */
    normal      = ddbox->normal;
    skew_fac_01 = 0;
    if (dd->ndim >= 2)
    {
        v_0 = ddbox->v[dim0];
        if (ddbox->tric_dir[dim0] && ddbox->tric_dir[dim1])
        {
            /* Determine the coupling coefficient for the distances
             * to the cell planes along dim0 and dim1 through dim2.
             * This is required for correct rounding.
             */
            skew_fac_01 = ddbox->v[dim0][dim1 + 1][dim0] * ddbox->v[dim1][dim1 + 1][dim1];
            if (debug)
            {
                fprintf(debug, "\nskew_fac_01 %f\n", skew_fac_01);
            }
        }
    }
    if (dd->ndim >= 3)
    {
        v_1 = ddbox->v[dim1];
    }

    gmx::ArrayRef<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock =
            fr->atomInfoForEachMoleculeBlock;

    gmx::DomdecZones& zones = dd->zones;

    zones.setAtomRangeEnd(0, dd->numHomeAtoms, true);
    pos_cg = dd->numHomeAtoms;

    nat_tot = comm->atomRanges.numHomeAtoms();
    nzone   = 1;
    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        dim = dd->dim[dim_ind];
        cd  = &comm->cd[dim_ind];

        /* Check if we need to compute triclinic distances along this dim */
        bool distanceIsTriclinic = false;
        for (int i = 0; i <= dim_ind; i++)
        {
            if (ddbox->tric_dir[dd->dim[i]])
            {
                distanceIsTriclinic = true;
            }
        }

        if (dim >= ddbox->npbcdim && dd->ci[dim] == 0)
        {
            /* No pbc in this dimension, the first node should not comm. */
            nzone_send = 0;
        }
        else
        {
            nzone_send = nzone;
        }

        v_d         = ddbox->v[dim];
        skew_fac2_d = gmx::square(ddbox->skew_fac[dim]);

        cd->receiveInPlace = true;
        for (int p = 0; p < cd->numPulses(); p++)
        {
            /* Only atoms communicated in the first pulse are used
             * for multi-body bonded interactions or for bBondComm.
             */
            bDistBonded = ((bDistMB || bDist2B) && p == 0);

            gmx_domdec_ind_t* ind = &cd->ind[p];

            /* Thread 0 writes in the global index array */
            ind->index.clear();
            clearCommSetupData(&comm->dth[0]);

            for (zone = 0; zone < nzone_send; zone++)
            {
                if (dim_ind > 0 && distanceIsTriclinic)
                {
                    /* Determine slightly more optimized skew_fac's
                     * for rounding.
                     * This reduces the number of communicated atoms
                     * by about 10% for 3D DD of rhombic dodecahedra.
                     */
                    for (dimd = 0; dimd < dim; dimd++)
                    {
                        sf2_round[dimd] = 1;
                        if (ddbox->tric_dir[dimd])
                        {
                            for (int i = dd->dim[dimd] + 1; i < DIM; i++)
                            {
                                /* If we are shifted in dimension i
                                 * and the cell plane is tilted forward
                                 * in dimension i, skip this coupling.
                                 */
                                if (!(zones.shift(nzone + zone)[i] && ddbox->v[dimd][i][dimd] >= 0))
                                {
                                    sf2_round[dimd] += gmx::square(ddbox->v[dimd][i][dimd]);
                                }
                            }
                            sf2_round[dimd] = 1 / sf2_round[dimd];
                        }
                    }
                }

                zonei = zone_perm[dim_ind][zone];

                int atomStart;
                int atomEnd;
                if (p == 0)
                {
                    /* Here we permutate the zones to obtain a convenient order
                     * for neighbor searching
                     */
                    atomStart = *zones.atomRange(zonei).begin();
                    atomEnd   = *zones.atomRange(zonei).end();
                }
                else
                {
                    /* Look only at the atoms received in the previous grid pulse
                     */
                    atomEnd   = *zones.atomRange(nzone + zone).end();
                    atomStart = atomEnd - cd->ind[p - 1].nrecv[zone];
                }

                const int numThreads = gmx::ssize(comm->dth);
#pragma omp parallel for num_threads(numThreads) schedule(static)
                for (int th = 0; th < numThreads; th++)
                {
                    try
                    {
                        dd_comm_setup_work_t& work = comm->dth[th];

                        /* Retain data accumulated into buffers of thread 0 */
                        if (th > 0)
                        {
                            clearCommSetupData(&work);
                        }

                        const int taskAtomStart = atomStart + ((atomEnd - atomStart) * th) / numThreads;
                        const int taskAtomEnd = atomStart + ((atomEnd - atomStart) * (th + 1)) / numThreads;

                        /* Get the atom groups and coordinates for this pulse in this zone */
                        get_zone_pulse_groups(dd,
                                              zonei,
                                              zone,
                                              taskAtomStart,
                                              taskAtomEnd,
                                              dd->globalAtomIndices,
                                              dim,
                                              dim_ind,
                                              dim0,
                                              dim1,
                                              dim2,
                                              r_comm2,
                                              r_bcomm2,
                                              box,
                                              distanceIsTriclinic,
                                              normal,
                                              skew_fac2_d,
                                              skew_fac_01,
                                              v_d,
                                              v_0,
                                              v_1,
                                              &corners,
                                              sf2_round,
                                              bDistBonded,
                                              bBondComm,
                                              bDist2B,
                                              bDistMB,
                                              state->x,
                                              fr->atomInfo,
                                              th == 0 ? &ind->index : &work.localAtomGroupBuffer,
                                              &work);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                } // END

                std::vector<int>&       atomGroups = comm->dth[0].atomGroupBuffer;
                std::vector<gmx::RVec>& positions  = comm->dth[0].positionBuffer;
                ind->nsend[zone]                   = comm->dth[0].nsend_zone;
                /* Append data of threads>=1 to the communication buffers */
                for (int th = 1; th < numThreads; th++)
                {
                    const dd_comm_setup_work_t& dth = comm->dth[th];

                    ind->index.insert(ind->index.end(),
                                      dth.localAtomGroupBuffer.begin(),
                                      dth.localAtomGroupBuffer.end());
                    atomGroups.insert(
                            atomGroups.end(), dth.atomGroupBuffer.begin(), dth.atomGroupBuffer.end());
                    positions.insert(
                            positions.end(), dth.positionBuffer.begin(), dth.positionBuffer.end());
                    comm->dth[0].nat += dth.nat;
                    ind->nsend[zone] += dth.nsend_zone;
                }
            }
            /* Clear the counts in case we do not have pbc */
            for (zone = nzone_send; zone < nzone; zone++)
            {
                ind->nsend[zone] = 0;
            }
            ind->nsend[nzone]     = ind->index.size();
            ind->nsend[nzone + 1] = comm->dth[0].nat;
            /* Communicate the number of atoms to receive */
            ddSendrecv(dd,
                       dim_ind,
                       dddirBackward,
                       gmx::arrayRefFromArray(ind->nsend, nzone + 2),
                       gmx::arrayRefFromArray(ind->nrecv, nzone + 2));

            if (p > 0)
            {
                /* We can receive in place if only the last zone is not empty */
                for (zone = 0; zone < nzone - 1; zone++)
                {
                    if (ind->nrecv[zone] > 0)
                    {
                        cd->receiveInPlace = false;
                    }
                }
            }

            int receiveBufferSize = 0;
            if (!cd->receiveInPlace)
            {
                receiveBufferSize = ind->nrecv[nzone];
            }
            /* These buffer are actually only needed with in-place */
            DDBufferAccess<int>       globalAtomBuffer(comm->intBuffer, receiveBufferSize);
            DDBufferAccess<gmx::RVec> rvecBuffer(comm->rvecBuffer, receiveBufferSize);

            dd_comm_setup_work_t& work = comm->dth[0];

            /* Make space for the global atom indices */
            int numAtomGroupsNew = pos_cg + ind->nrecv[nzone];
            dd->globalAtomIndices.resize(numAtomGroupsNew);
            /* Communicate the global atom indices */
            gmx::ArrayRef<int> integerBufferRef;
            if (cd->receiveInPlace)
            {
                integerBufferRef = gmx::arrayRefFromArray(dd->globalAtomIndices.data() + pos_cg,
                                                          ind->nrecv[nzone]);
            }
            else
            {
                integerBufferRef = globalAtomBuffer.buffer;
            }
            ddSendrecv<int>(dd, dim_ind, dddirBackward, work.atomGroupBuffer, integerBufferRef);

            /* Make space for atominfo */
            dd_resize_atominfo_and_state(fr, state, pos_cg + ind->nrecv[nzone]);

            /* Communicate the coordinates */
            gmx::ArrayRef<gmx::RVec> rvecBufferRef;
            if (cd->receiveInPlace)
            {
                rvecBufferRef = gmx::makeArrayRef(state->x).subArray(pos_cg, ind->nrecv[nzone]);
            }
            else
            {
                rvecBufferRef = rvecBuffer.buffer;
            }
            ddSendrecv<gmx::RVec>(dd, dim_ind, dddirBackward, work.positionBuffer, rvecBufferRef);

            /* Make the atom group index */
            if (cd->receiveInPlace)
            {
                zone = (p == 0 ? 0 : nzone - 1);
                while (zone < nzone)
                {
                    const int gmx_unused numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Domdec);
#pragma omp parallel for num_threads(numThreads) schedule(static)
                    for (int i = 0; i < ind->nrecv[zone]; i++)
                    {
                        int globalAtomIndex = dd->globalAtomIndices[pos_cg + i];
                        fr->atomInfo[pos_cg + i] =
                                ddGetAtomInfo(atomInfoForEachMoleculeBlock, globalAtomIndex);
                    }
                    pos_cg += ind->nrecv[zone];
                    zones.setAtomRangeEnd(nzone + zone, pos_cg, p == 0);

                    zone++;
                }
            }
            else
            {
                /* This part of the code is never executed with bBondComm. */
                mergeAtomBuffers(nzone,
                                 cd,
                                 p,
                                 &zones,
                                 dd->globalAtomIndices,
                                 integerBufferRef.data(),
                                 state->x,
                                 rvecBufferRef,
                                 fr->atomInfoForEachMoleculeBlock,
                                 fr->atomInfo);
                pos_cg += ind->nrecv[nzone];
            }
            nat_tot += ind->nrecv[nzone + 1];
        }
        if (!cd->receiveInPlace)
        {
            /* Store the atom block for easy copying of communication buffers */
            make_cell2at_index(cd, nzone, *zones.atomRange(nzone - 1).end());
        }
        nzone += nzone;
    }

    GMX_ASSERT(*zones.atomRange(zones.numZones() - 1).end() == nat_tot,
               "The zone atom counts should cover the whole atom range");

    comm->atomRanges.setEnd(DDAtomRanges::Type::Zones, nat_tot);

    if (!bBondComm)
    {
        /* We don't need to update atominfo, since that was already done above.
         * So we pass NULL for the forcerec.
         */
        dd_set_atominfo(dd->globalAtomIndices, dd->numHomeAtoms, dd->globalAtomIndices.size(), nullptr);
    }

    if (debug)
    {
        fprintf(debug, "Finished setting up DD communication, zones:");
        for (c = 0; c < zones.numZones(); c++)
        {
            fprintf(debug, " %d", zones.atomRange(c).size());
        }
        fprintf(debug, "\n");
    }
}

/*! \brief Order data in \p dataToSort according to \p sort
 *
 * Note: both buffers should have at least \p sort.size() elements.
 */
template<typename T>
static void orderVector(gmx::ArrayRef<const gmx_cgsort_t> sort,
                        gmx::ArrayRef<T>                  dataToSort,
                        gmx::ArrayRef<T>                  sortBuffer)
{
    GMX_ASSERT(dataToSort.size() >= sort.size(), "The vector needs to be sufficiently large");
    GMX_ASSERT(sortBuffer.size() >= sort.size(),
               "The sorting buffer needs to be sufficiently large");

    /* Order the data into the temporary buffer */
    size_t i = 0;
    for (const gmx_cgsort_t& entry : sort)
    {
        sortBuffer[i++] = dataToSort[entry.ind];
    }

    /* Copy back to the original array */
    std::copy(sortBuffer.begin(), sortBuffer.begin() + sort.size(), dataToSort.begin());
}

/*! \brief Order data in \p dataToSort according to \p sort
 *
 * Note: \p vectorToSort should have at least \p sort.size() elements,
 *       \p workVector is resized when it is too small.
 */
template<typename T>
static void orderVector(gmx::ArrayRef<const gmx_cgsort_t> sort,
                        gmx::ArrayRef<T>                  vectorToSort,
                        std::vector<T>*                   workVector)
{
    if (gmx::Index(workVector->size()) < sort.ssize())
    {
        workVector->resize(sort.size());
    }
    orderVector<T>(sort, vectorToSort, *workVector);
}

//! Returns the sorting order for atoms based on the nbnxn grid order in sort
static void dd_sort_order_nbnxn(const t_forcerec* fr, std::vector<gmx_cgsort_t>* sort)
{
    gmx::ArrayRef<const int> atomOrder = fr->nbv->getLocalAtomOrder();

    /* Using push_back() instead of this resize results in much slower code */
    sort->resize(atomOrder.size());
    gmx::ArrayRef<gmx_cgsort_t> buffer    = *sort;
    size_t                      numSorted = 0;
    for (int i : atomOrder)
    {
        if (i >= 0)
        {
            buffer[numSorted++].ind = i;
        }
    }
    sort->resize(numSorted);
}

//! Returns the sorting state for DD.
static void dd_sort_state(gmx_domdec_t* dd, t_forcerec* fr, t_state* state)
{
    gmx_domdec_sort_t* sort = dd->comm->sort.get();

    dd_sort_order_nbnxn(fr, &sort->sorted);

    /* We alloc with the old size, since cgindex is still old */
    DDBufferAccess<gmx::RVec> rvecBuffer(dd->comm->rvecBuffer, dd->numHomeAtoms);

    /* Set the new home atom/charge group count */
    dd->numHomeAtoms = sort->sorted.size();
    if (debug)
    {
        fprintf(debug, "Set the new home atom count to %d\n", dd->numHomeAtoms);
    }

    /* Reorder the state */
    gmx::ArrayRef<const gmx_cgsort_t> cgsort = sort->sorted;
    GMX_RELEASE_ASSERT(cgsort.ssize() == dd->numHomeAtoms,
                       "We should sort all the home atom groups");

    if (state->hasEntry(StateEntry::X))
    {
        orderVector(cgsort, makeArrayRef(state->x), rvecBuffer.buffer);
    }
    if (state->hasEntry(StateEntry::V))
    {
        orderVector(cgsort, makeArrayRef(state->v), rvecBuffer.buffer);
    }
    if (state->hasEntry(StateEntry::Cgp))
    {
        orderVector(cgsort, makeArrayRef(state->cg_p), rvecBuffer.buffer);
    }

    /* Reorder the global cg index */
    orderVector<int>(cgsort, dd->globalAtomIndices, &sort->intBuffer);
    /* Reorder the atom info */
    orderVector<int>(cgsort, fr->atomInfo, &sort->intBuffer);
    /* Set the home atom number */
    dd->comm->atomRanges.setEnd(DDAtomRanges::Type::Home, dd->numHomeAtoms);

    /* The atoms are now exactly in grid order, update the grid order */
    fr->nbv->setLocalAtomOrder();
}

//! Accumulates load statistics.
static void add_dd_statistics(gmx_domdec_t* dd)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    for (int i = 0; i < static_cast<int>(DDAtomRanges::Type::Number); i++)
    {
        auto range = static_cast<DDAtomRanges::Type>(i);
        comm->sum_nat[i] += comm->atomRanges.end(range) - comm->atomRanges.start(range);
    }
    comm->ndecomp++;
}

void reset_dd_statistics_counters(gmx_domdec_t* dd)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    /* Reset all the statistics and counters for total run counting */
    for (int i = 0; i < static_cast<int>(DDAtomRanges::Type::Number); i++)
    {
        comm->sum_nat[i] = 0;
    }
    comm->ndecomp   = 0;
    comm->nload     = 0;
    comm->load_step = 0;
    comm->load_sum  = 0;
    comm->load_max  = 0;
    clear_ivec(comm->load_lim);
    comm->load_mdf = 0;
    comm->load_pme = 0;
}

namespace gmx
{

bool check_grid_jump(int64_t step, const gmx_domdec_t* dd, real cutoff, const gmx_ddbox_t* ddbox, bool bFatal)
{
    gmx_domdec_comm_t* comm    = dd->comm.get();
    bool               invalid = false;

    for (int d = 1; d < dd->ndim; d++)
    {
        const DDCellsizesWithDlb& cellsizes = comm->cellsizesWithDlb[d];
        const int                 dim       = dd->dim[d];
        const real                limit     = grid_jump_limit(comm, cutoff, d);
        real                      bfac      = ddbox->box_size[dim];
        if (ddbox->tric_dir[dim])
        {
            bfac *= ddbox->skew_fac[dim];
        }
        if ((cellsizes.fracUpper - cellsizes.fracLowerMax) * bfac < limit
            || (cellsizes.fracLower - cellsizes.fracUpperMin) * bfac > -limit)
        {
            invalid = true;

            if (bFatal)
            {
                char buf[22];

                /* This error should never be triggered under normal
                 * circumstances, but you never know ...
                 */
                gmx_fatal(FARGS,
                          "step %s: The domain decomposition grid has shifted too much in the "
                          "%c-direction around cell %d %d %d. This should not have happened. "
                          "Running with fewer ranks might avoid this issue.",
                          gmx_step_str(step, buf),
                          dim2char(dim),
                          dd->ci[XX],
                          dd->ci[YY],
                          dd->ci[ZZ]);
            }
        }
    }

    return invalid;
}

void print_dd_statistics(const t_commrec* cr, const t_inputrec& inputrec, FILE* fplog)
{
    gmx_domdec_comm_t* comm = cr->dd->comm.get();

    const int numRanges = static_cast<int>(DDAtomRanges::Type::Number);
    gmx_sumd(numRanges, comm->sum_nat, cr);

    if (fplog == nullptr)
    {
        return;
    }

    fprintf(fplog, "\n    D O M A I N   D E C O M P O S I T I O N   S T A T I S T I C S\n\n");

    for (int i = static_cast<int>(DDAtomRanges::Type::Zones); i < numRanges; i++)
    {
        auto   range = static_cast<DDAtomRanges::Type>(i);
        double av    = comm->sum_nat[i] / comm->ndecomp;
        switch (range)
        {
            case DDAtomRanges::Type::Zones:
                fprintf(fplog, " av. #atoms communicated per step for force:  %d x %.1f\n", 2, av);
                break;
            case DDAtomRanges::Type::Vsites:
                if (cr->dd->vsite_comm)
                {
                    fprintf(fplog,
                            " av. #atoms communicated per step for vsites: %d x %.1f\n",
                            (usingPme(inputrec.coulombtype)
                             || inputrec.coulombtype == CoulombInteractionType::Ewald)
                                    ? 3
                                    : 2,
                            av);
                }
                break;
            case DDAtomRanges::Type::Constraints:
                if (cr->dd->constraint_comm)
                {
                    fprintf(fplog,
                            " av. #atoms communicated per step for LINCS:  %d x %.1f\n",
                            1 + inputrec.nLincsIter,
                            av);
                }
                break;
            default: gmx_incons(" Unknown type for DD statistics");
        }
    }
    fprintf(fplog, "\n");

    if (comm->ddSettings.recordLoad && EI_DYNAMICS(inputrec.eI))
    {
        print_dd_load_av(fplog, cr->dd);
    }
}

//!\brief TODO Remove fplog when group scheme and charge groups are gone
void dd_partition_system(FILE*                     fplog,
                         const gmx::MDLogger&      mdlog,
                         int64_t                   step,
                         const t_commrec*          cr,
                         bool                      bMainState,
                         t_state*                  state_global,
                         const gmx_mtop_t&         top_global,
                         const t_inputrec&         inputrec,
                         const MDModulesNotifiers& mdModulesNotifiers,
                         gmx::ImdSession*          imdSession,
                         pull_t*                   pull_work,
                         t_state*                  state_local,
                         gmx::ForceBuffers*        f,
                         gmx::MDAtoms*             mdAtoms,
                         gmx_localtop_t*           top_local,
                         t_forcerec*               fr,
                         gmx::VirtualSitesHandler* vsite,
                         gmx::Constraints*         constr,
                         t_nrnb*                   nrnb,
                         gmx_wallcycle*            wcycle,
                         bool                      bVerbose)
{
    gmx_ddbox_t ddbox = { 0 };
    int         ncgindex_set;
    char        sbuf[22];

    wallcycle_start(wcycle, WallCycleCounter::Domdec);

    gmx_domdec_t*      dd   = cr->dd;
    gmx_domdec_comm_t* comm = dd->comm.get();

    // TODO if the update code becomes accessible here, use
    // upd->deform for this logic.
    bool bBoxChanged = (bMainState || inputrecDeform(&inputrec));
    if (inputrec.pressureCouplingOptions.epc != PressureCoupling::No)
    {
        /* With nstpcouple > 1 pressure coupling happens.
         * one step after calculating the pressure.
         * Box scaling happens at the end of the MD step,
         * after the DD partitioning.
         * We therefore have to do DLB in the first partitioning
         * after an MD step where P-coupling occurred.
         * We need to determine the last step in which p-coupling occurred.
         * MRS -- need to validate this for vv?
         */
        int     n = inputrec.pressureCouplingOptions.nstpcouple;
        int64_t step_pcoupl;
        if (n == 1)
        {
            step_pcoupl = step - 1;
        }
        else
        {
            step_pcoupl = ((step - 1) / n) * n + 1;
        }
        if (step_pcoupl >= comm->partition_step)
        {
            bBoxChanged = true;
        }
    }

    bool bDoDLB;
    if (!isDlbOn(comm->dlbState))
    {
        bDoDLB = false;
    }
    else
    {
        /* Should we do dynamic load balacing this step?
         * Since it requires (possibly expensive) global communication,
         * we might want to do DLB less frequently.
         */
        if (bBoxChanged || inputrec.pressureCouplingOptions.epc != PressureCoupling::No)
        {
            bDoDLB = bBoxChanged;
        }
        else
        {
            bDoDLB = (step % dd->comm->nstDDGlobalComm == 0);
        }
    }

    /* Check if we have recorded loads on the nodes */
    if (comm->ddSettings.recordLoad && dd_load_count(comm) > 0)
    {
        bool bCheckWhetherToTurnDlbOn = dd_dlb_get_should_check_whether_to_turn_dlb_on(dd);

        /* Print load every nstlog, first and last step to the log file */
        bool bLogLoad = ((inputrec.nstlog > 0 && step % inputrec.nstlog == 0) || comm->n_load_collect == 0
                         || (inputrec.nsteps >= 0
                             && (step + inputrec.nstlist > inputrec.init_step + inputrec.nsteps)));

        if (bDoDLB || bLogLoad || bCheckWhetherToTurnDlbOn || bVerbose)
        {
            get_load_distribution(dd, wcycle);
            if (DDMAIN(dd))
            {
                if (bLogLoad)
                {
                    GMX_LOG(mdlog.info).asParagraph().appendText(dd_print_load(dd, step - 1));
                }
                if (bVerbose)
                {
                    dd_print_load_verbose(dd);
                }
            }
            comm->n_load_collect++;

            if (isDlbOn(comm->dlbState))
            {
                if (DDMAIN(dd))
                {
                    /* Add the measured cycles to the running average */
                    const float averageFactor = 0.1F;
                    comm->cyclesPerStepDlbExpAverage =
                            (1 - averageFactor) * comm->cyclesPerStepDlbExpAverage
                            + averageFactor * comm->cycl[ddCyclStep] / comm->cycl_n[ddCyclStep];
                }
                if (comm->dlbState == DlbState::onCanTurnOff
                    && dd->comm->n_load_have % c_checkTurnDlbOffInterval == c_checkTurnDlbOffInterval - 1)
                {
                    bool turnOffDlb;
                    if (DDMAIN(dd))
                    {
                        /* If the running averaged cycles with DLB are more
                         * than before we turned on DLB, turn off DLB.
                         * We will again run and check the cycles without DLB
                         * and we can then decide if to turn off DLB forever.
                         */
                        turnOffDlb = (comm->cyclesPerStepDlbExpAverage > comm->cyclesPerStepBeforeDLB);
                    }
                    dd_bcast(dd, sizeof(turnOffDlb), &turnOffDlb);
                    if (turnOffDlb)
                    {
                        /* To turn off DLB, we need to redistribute the atoms */
                        dd_collect_state(dd, state_local, state_global);
                        bMainState = true;
                        turn_off_dlb(mdlog, dd, step);
                    }
                }
            }
            else if (bCheckWhetherToTurnDlbOn)
            {
                bool turnOffDlbForever = false;
                bool turnOnDlb         = false;

                /* Since the timings are node dependent, the main decides */
                if (DDMAIN(dd))
                {
                    /* If we recently turned off DLB, we want to check if
                     * performance is better without DLB. We want to do this
                     * ASAP to minimize the chance that external factors
                     * slowed down the DLB step are gone here and we
                     * incorrectly conclude that DLB was causing the slowdown.
                     * So we measure one nstlist block, no running average.
                     */
                    if (comm->haveTurnedOffDlb
                        && comm->cycl[ddCyclStep] / comm->cycl_n[ddCyclStep] < comm->cyclesPerStepDlbExpAverage)
                    {
                        /* After turning off DLB we ran nstlist steps in fewer
                         * cycles than with DLB. This likely means that DLB
                         * in not benefical, but this could be due to a one
                         * time unlucky fluctuation, so we require two such
                         * observations in close succession to turn off DLB
                         * forever.
                         */
                        if (comm->dlbSlowerPartitioningCount > 0
                            && dd->ddp_count < comm->dlbSlowerPartitioningCount + 10 * c_checkTurnDlbOnInterval)
                        {
                            turnOffDlbForever = true;
                        }
                        comm->haveTurnedOffDlb = false;
                        /* Register when we last measured DLB slowdown */
                        comm->dlbSlowerPartitioningCount = dd->ddp_count;
                    }
                    else
                    {
                        /* Here we check if the max PME rank load is more than 0.98
                         * the max PP force load. If so, PP DLB will not help,
                         * since we are (almost) limited by PME. Furthermore,
                         * DLB will cause a significant extra x/f redistribution
                         * cost on the PME ranks, which will then surely result
                         * in lower total performance.
                         */
                        if (comm->ddRankSetup.usePmeOnlyRanks && dd_pme_f_ratio(dd) > 1 - DD_PERF_LOSS_DLB_ON)
                        {
                            turnOnDlb = false;
                        }
                        else
                        {
                            turnOnDlb = (dd_force_imb_perf_loss(dd) >= DD_PERF_LOSS_DLB_ON);
                        }
                    }
                }
                struct
                {
                    bool turnOffDlbForever;
                    bool turnOnDlb;
                } bools{ turnOffDlbForever, turnOnDlb };
                dd_bcast(dd, sizeof(bools), &bools);
                if (bools.turnOffDlbForever)
                {
                    turn_off_dlb_forever(mdlog, dd, step);
                }
                else if (bools.turnOnDlb)
                {
                    turn_on_dlb(mdlog, dd, step);
                    bDoDLB = true;
                }
            }
        }
        comm->n_load_have++;
    }

    bool bRedist = false;
    if (bMainState)
    {
        /* Clear the old state */
        clearDDStateIndices(dd, false);
        ncgindex_set = 0;

        auto xGlobal = positionsFromStatePointer(state_global);

        set_ddbox(*dd, true, DDMAIN(dd) ? state_global->box : nullptr, true, xGlobal, &ddbox);

        distributeState(mdlog, dd, top_global, state_global, ddbox, state_local);

        /* Ensure that we have space for the new distribution */
        dd_resize_atominfo_and_state(fr, state_local, dd->numHomeAtoms);

        inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());

        dd_set_atominfo(dd->globalAtomIndices, 0, dd->numHomeAtoms, fr);
    }
    else if (state_local->ddp_count != dd->ddp_count)
    {
        if (state_local->ddp_count > dd->ddp_count)
        {
            gmx_fatal(FARGS,
                      "Internal inconsistency state_local->ddp_count (%d) > dd->ddp_count (%" PRId64
                      ")",
                      state_local->ddp_count,
                      dd->ddp_count);
        }

        if (state_local->ddp_count_cg_gl != state_local->ddp_count)
        {
            gmx_fatal(FARGS,
                      "Internal inconsistency state_local->ddp_count_cg_gl (%d) != "
                      "state_local->ddp_count (%d)",
                      state_local->ddp_count_cg_gl,
                      state_local->ddp_count);
        }

        /* Clear the old state */
        clearDDStateIndices(dd, false);

        /* Restore the atom group indices from state_local */
        restoreAtomGroups(dd, state_local);
        make_dd_indices(dd, 0);
        ncgindex_set = dd->numHomeAtoms;

        inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());

        dd_set_atominfo(dd->globalAtomIndices, 0, dd->numHomeAtoms, fr);

        set_ddbox(*dd, bMainState, state_local->box, true, state_local->x, &ddbox);

        bRedist = isDlbOn(comm->dlbState);
    }
    else
    {
        /* We have the full state, only redistribute the cgs */

        /* Clear the non-home indices */
        clearDDStateIndices(dd, true);
        ncgindex_set = 0;

        /* To avoid global communication, we do not recompute the extent
         * of the system for dims without pbc. Therefore we need to copy
         * the previously computed values when we do not communicate.
         */
        const bool doGlobalComm = (step % dd->comm->nstDDGlobalComm == 0);
        if (!doGlobalComm)
        {
            copy_rvec(comm->box0, ddbox.box0);
            copy_rvec(comm->box_size, ddbox.box_size);
        }
        set_ddbox(*dd, bMainState, state_local->box, doGlobalComm, state_local->x, &ddbox);

        bBoxChanged = true;
        bRedist     = true;
    }
    /* Copy needed for dim's without pbc when avoiding communication */
    copy_rvec(ddbox.box0, comm->box0);
    copy_rvec(ddbox.box_size, comm->box_size);

    set_dd_cell_sizes(dd, &ddbox, dd->unitCellInfo.ddBoxIsDynamic, bMainState, bDoDLB, step, wcycle);

    if (comm->ddSettings.nstDDDumpGrid > 0 && step % comm->ddSettings.nstDDDumpGrid == 0)
    {
        write_dd_grid_pdb("dd_grid", step, dd, state_local->box, &ddbox);
    }

    if (comm->systemInfo.useUpdateGroups)
    {
        comm->updateGroupsCog->addCogs(
                gmx::arrayRefFromArray(dd->globalAtomIndices.data(), dd->numHomeAtoms), state_local->x);
    }

    /* Check if we should sort the charge groups */
    const bool bSortCG = (bMainState || bRedist);

    /* When repartitioning we mark atom groups that will move to neighboring
     * DD cells, but we do not move them right away for performance reasons.
     * Thus we need to keep track of how many charge groups will move for
     * obtaining correct local charge group / atom counts.
     */
    int ncg_moved = 0;
    if (bRedist)
    {
        wallcycle_sub_start(wcycle, WallCycleSubCounter::DDRedist);

        ncgindex_set = dd->numHomeAtoms;
        dd_redistribute_cg(fplog, step, dd, ddbox.tric_dir, state_local, fr, nrnb, &ncg_moved);

        GMX_RELEASE_ASSERT(bSortCG, "Sorting is required after redistribution");

        if (comm->systemInfo.useUpdateGroups)
        {
            comm->updateGroupsCog->addCogs(
                    gmx::arrayRefFromArray(dd->globalAtomIndices.data(), dd->numHomeAtoms),
                    state_local->x);
        }

        wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDRedist);
    }

    RVec cell_ns_x0, cell_ns_x1;
    get_nsgrid_boundaries(ddbox.nboundeddim,
                          state_local->box,
                          dd,
                          &ddbox,
                          &comm->cell_x0,
                          &comm->cell_x1,
                          dd->numHomeAtoms,
                          as_rvec_array(state_local->x.data()),
                          cell_ns_x0,
                          cell_ns_x1);

    if (bBoxChanged)
    {
        comm_dd_ns_cell_sizes(dd, &ddbox, cell_ns_x0, cell_ns_x1, step);
    }

    if (bSortCG)
    {
        wallcycle_sub_start(wcycle, WallCycleSubCounter::DDGrid);

        /* Sort the state on charge group position.
         * This enables exact restarts from this step.
         * It also improves performance by about 15% with larger numbers
         * of atoms per node.
         */

        /* Fill the ns grid with the home cell,
         * so we can sort with the indices.
         */
        set_zones_numHomeAtoms(dd);

        dd->zones.setSizes(*dd, state_local->box, &ddbox, { 0, 1 });

        real homeZoneVolume = 1;
        for (int dim = 0; dim < DIM; dim++)
        {
            homeZoneVolume *= dd->zones.sizes(0).x1[dim] - dd->zones.sizes(0).x0[dim];
        }
        // The home atom list still contains moved atoms, compute the new atom count
        const int  newNumHomeAtoms = dd->numHomeAtoms - ncg_moved;
        const real atomDensity     = newNumHomeAtoms / homeZoneVolume;

        fr->nbv->putAtomsOnGrid(state_local->box,
                                0,
                                dd->zones.sizes(0).bb_x0,
                                dd->zones.sizes(0).bb_x1,
                                comm->updateGroupsCog.get(),
                                { 0, dd->numHomeAtoms },
                                newNumHomeAtoms,
                                atomDensity,
                                fr->atomInfo,
                                state_local->x,
                                bRedist ? comm->movedBuffer.data() : nullptr);

        if (debug)
        {
            fprintf(debug, "Step %s, sorting the %d home charge groups\n", gmx_step_str(step, sbuf), dd->numHomeAtoms);
        }
        dd_sort_state(dd, fr, state_local);

        /* After sorting and compacting we set the correct size */
        state_local->changeNumAtoms(comm->atomRanges.numHomeAtoms());

        /* Rebuild all the indices */
        dd->ga2la->clear(false);
        ncgindex_set = 0;

        wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDGrid);
    }
    else
    {
        /* With the group scheme the sorting array is part of the DD state,
         * but it just got out of sync, so mark as invalid by emptying it.
         */
        if (inputrec.cutoff_scheme == CutoffScheme::Group)
        {
            comm->sort->sorted.clear();
        }
    }

    if (comm->systemInfo.useUpdateGroups)
    {
        /* The update groups cog's are invalid after sorting
         * and need to be cleared before the next partitioning anyhow.
         */
        comm->updateGroupsCog->clear();
    }

    wallcycle_sub_start(wcycle, WallCycleSubCounter::DDSetupComm);

    /* Set the induces for the home atoms */
    set_zones_numHomeAtoms(dd);
    make_dd_indices(dd, ncgindex_set);

    /* Setup up the communication and communicate the coordinates */
    setup_dd_communication(dd, state_local->box, &ddbox, fr, state_local);

    /* Set the indices for the halo atoms */
    make_dd_indices(dd, dd->numHomeAtoms);

    /* When bSortCG=true, we have already set the size for zone 0 */
    dd->zones.setSizes(*dd, state_local->box, &ddbox, { bSortCG ? 1 : 0, dd->zones.numZones() });

    wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDSetupComm);

    /*
       write_dd_pdb("dd_home",step,"dump",top_global,cr,
                 -1,state_local->x.rvec_array(),state_local->box);
     */

    wallcycle_sub_start(wcycle, WallCycleSubCounter::DDMakeTop);

    /* Extract a local topology from the global topology */
    IVec numPulses;
    for (int i = 0; i < dd->ndim; i++)
    {
        numPulses[dd->dim[i]] = comm->cd[i].numPulses();
    }
    int numBondedInteractionsToReduce = dd_make_local_top(*dd,
                                                          dd->zones,
                                                          dd->unitCellInfo.npbcdim,
                                                          state_local->box,
                                                          comm->cellsize_min,
                                                          numPulses,
                                                          fr,
                                                          state_local->x,
                                                          top_global,
                                                          fr->atomInfo,
                                                          top_local);
    dd->localTopologyChecker->scheduleCheckOfLocalTopology(numBondedInteractionsToReduce);

    wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDMakeTop);

    wallcycle_sub_start(wcycle, WallCycleSubCounter::DDMakeConstr);

    /* Set up the special atom communication */
    int n = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int i = static_cast<int>(DDAtomRanges::Type::Zones) + 1;
         i < static_cast<int>(DDAtomRanges::Type::Number);
         i++)
    {
        auto range = static_cast<DDAtomRanges::Type>(i);
        switch (range)
        {
            case DDAtomRanges::Type::Vsites:
                if (vsite && vsite->numInterUpdategroupVirtualSites())
                {
                    n = dd_make_local_vsites(dd, n, top_local->idef.il);
                }
                break;
            case DDAtomRanges::Type::Constraints:
                if (dd->comm->systemInfo.mayHaveSplitConstraints || dd->comm->systemInfo.mayHaveSplitSettles)
                {
                    /* Only for inter-atom-group constraints we need special code */
                    n = dd_make_local_constraints(dd,
                                                  n,
                                                  top_global,
                                                  fr->atomInfo,
                                                  constr,
                                                  inputrec.nProjOrder,
                                                  top_local->idef.il);
                }
                break;
            default: gmx_incons("Unknown special atom type setup");
        }
        comm->atomRanges.setEnd(range, n);
    }

    wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDMakeConstr);

    wallcycle_sub_start(wcycle, WallCycleSubCounter::DDTopOther);

    /* Make space for the extra coordinates for virtual site
     * or constraint communication.
     */
    state_local->changeNumAtoms(comm->atomRanges.numAtomsTotal());

    int nat_f_novirsum;
    if (vsite && vsite->numInterUpdategroupVirtualSites())
    {
        nat_f_novirsum = comm->atomRanges.end(DDAtomRanges::Type::Vsites);
    }
    else
    {
        if (usingFullElectrostatics(inputrec.coulombtype) && dd->haveExclusions)
        {
            nat_f_novirsum = comm->atomRanges.end(DDAtomRanges::Type::Zones);
        }
        else
        {
            nat_f_novirsum = comm->atomRanges.numHomeAtoms();
        }
    }

    /* Set the number of atoms required for the force calculation.
     * Forces need to be constrained when doing energy
     * minimization. For simple simulations we could avoid some
     * allocation, zeroing and copying, but this is probably not worth
     * the complications and checking.
     */
    forcerec_set_ranges(fr,
                        comm->atomRanges.end(DDAtomRanges::Type::Zones),
                        comm->atomRanges.end(DDAtomRanges::Type::Constraints),
                        nat_f_novirsum);

    /* Update atom data for mdatoms and several algorithms */
    mdAlgorithmsSetupAtomData(cr, inputrec, top_global, top_local, fr, f, mdAtoms, constr, vsite, nullptr);

    auto* mdatoms = mdAtoms->mdatoms();
    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Send the charges and/or c6/sigmas to our PME only node */
        gmx_pme_send_parameters(cr,
                                *fr->ic,
                                mdatoms->nChargePerturbed != 0,
                                mdatoms->nTypePerturbed != 0,
                                mdatoms->chargeA,
                                mdatoms->chargeB,
                                mdatoms->sqrt_c6A,
                                mdatoms->sqrt_c6B,
                                mdatoms->sigmaA,
                                mdatoms->sigmaB,
                                dd_pme_maxshift_x(*dd),
                                dd_pme_maxshift_y(*dd));
    }

    if (dd->atomSets != nullptr)
    {
        /* Update the local atom sets */
        dd->atomSets->setIndicesInDomainDecomposition(*(dd->ga2la));
    }

    // The pull group construction can need the atom sets updated above
    if (inputrec.bPull)
    {
        /* Update the local pull groups */
        dd_make_local_pull_groups(cr, pull_work);
    }

    /* Update the local atoms to be communicated via the IMD protocol if bIMD is true. */
    imdSession->dd_make_local_IMD_atoms(dd);

    add_dd_statistics(dd);

    /* Make sure we only count the cycles for this DD partitioning */
    clear_dd_cycle_counts(dd);

    /* Because the order of the atoms might have changed since
     * the last vsite construction, we need to communicate the constructing
     * atom coordinates again (for spreading the forces this MD step).
     */
    dd_move_x_vsites(*dd, state_local->box, state_local->x);

    wallcycle_sub_stop(wcycle, WallCycleSubCounter::DDTopOther);

    if (comm->ddSettings.nstDDDump > 0 && step % comm->ddSettings.nstDDDump == 0)
    {
        dd_move_x(dd, state_local->box, state_local->x, nullptr);
        write_dd_pdb("dd_dump",
                     step,
                     "dump",
                     top_global,
                     cr,
                     -1,
                     state_local->x.rvec_array(),
                     state_local->box);
    }

    /* Store the partitioning step */
    comm->partition_step = step;

    /* Increase the DD partitioning counter */
    dd->ddp_count++;
    /* The state currently matches this DD partitioning count, store it */
    state_local->ddp_count = dd->ddp_count;
    if (bMainState)
    {
        /* The DD main node knows the complete atom distribution,
         * store the count so we can possibly skip the atom info communication.
         */
        comm->main_cg_ddp_count = (bSortCG ? 0 : dd->ddp_count);
    }

    if (comm->ddSettings.DD_debug > 0)
    {
        /* Set the env var GMX_DD_DEBUG if you suspect corrupted indices */
        check_index_consistency(dd, top_global.natoms, "after partitioning");
    }

    // Now we have made the local atom sets and x is up to date, MDModules can be signaled
    MDModulesAtomsRedistributedSignal mdModulesAtomsRedistributedSignal(
            state_local->box,
            gmx::makeConstArrayRef(state_local->x).subArray(0, comm->atomRanges.numHomeAtoms()));
    mdModulesNotifiers.simulationSetupNotifier_.notify(mdModulesAtomsRedistributedSignal);

    wallcycle_stop(wcycle, WallCycleCounter::Domdec);
}

} // namespace gmx
