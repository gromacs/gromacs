/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "domdec.h"

#include "config.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/chargegroup.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/constraintrange.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/mdsetup.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_grid.h"
#include "gromacs/mdlib/nsgrid.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/qsort_threadsafe.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "atomdistribution.h"
#include "cellsizes.h"
#include "distribute.h"
#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_vsite.h"
#include "redistribute.h"
#include "utility.h"

#define DD_NLOAD_MAX 9

static const char *edlbs_names[edlbsNR] = { "off", "auto", "locked", "on", "on" };

/* The size per charge group of the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_CGIBS 2

/* The flags for the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_FLAG_NRCG  65535
#define DD_FLAG_FW(d) (1<<(16+(d)*2))
#define DD_FLAG_BW(d) (1<<(16+(d)*2+1))

/* The DD zone order */
static const ivec dd_zo[DD_MAXZONE] =
{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}};

/* The non-bonded zone-pair setup for domain decomposition
 * The first number is the i-zone, the second number the first j-zone seen by
 * this i-zone, the third number the last+1 j-zone seen by this i-zone.
 * As is, this is for 3D decomposition, where there are 4 i-zones.
 * With 2D decomposition use only the first 2 i-zones and a last+1 j-zone of 4.
 * With 1D decomposition use only the first i-zone and a last+1 j-zone of 2.
 */
static const int
    ddNonbondedZonePairRanges[DD_MAXIZONE][3] = {{0, 0, 8},
                                                 {1, 3, 6},
                                                 {2, 5, 6},
                                                 {3, 5, 7}};

/* Turn on DLB when the load imbalance causes this amount of total loss.
 * There is a bit of overhead with DLB and it's difficult to achieve
 * a load imbalance of less than 2% with DLB.
 */
#define DD_PERF_LOSS_DLB_ON  0.02

/* Warn about imbalance due to PP or PP/PME load imbalance at this loss */
#define DD_PERF_LOSS_WARN    0.05


/* We check if to turn on DLB at the first and every 100 DD partitionings.
 * With large imbalance DLB will turn on at the first step, so we can
 * make the interval so large that the MPI overhead of the check is negligible.
 */
static const int c_checkTurnDlbOnInterval  = 100;
/* We need to check if DLB results in worse performance and then turn it off.
 * We check this more often then for turning DLB on, because the DLB can scale
 * the domains very rapidly, so if unlucky the load imbalance can go up quickly
 * and furthermore, we are already synchronizing often with DLB, so
 * the overhead of the MPI Bcast is not that high.
 */
static const int c_checkTurnDlbOffInterval =  20;

/* Forward declaration */
static void dd_dlb_set_should_check_whether_to_turn_dlb_on(gmx_domdec_t *dd, gmx_bool bValue);


/*
   #define dd_index(n,i) ((((i)[ZZ]*(n)[YY] + (i)[YY])*(n)[XX]) + (i)[XX])

   static void index2xyz(ivec nc,int ind,ivec xyz)
   {
   xyz[XX] = ind % nc[XX];
   xyz[YY] = (ind / nc[XX]) % nc[YY];
   xyz[ZZ] = ind / (nc[YY]*nc[XX]);
   }
 */

static void ddindex2xyz(const ivec nc, int ind, ivec xyz)
{
    xyz[XX] = ind / (nc[YY]*nc[ZZ]);
    xyz[YY] = (ind / nc[ZZ]) % nc[YY];
    xyz[ZZ] = ind % nc[ZZ];
}

static int ddcoord2ddnodeid(gmx_domdec_t *dd, ivec c)
{
    int ddindex;
    int ddnodeid = -1;

    ddindex = dd_index(dd->nc, c);
    if (dd->comm->bCartesianPP_PME)
    {
        ddnodeid = dd->comm->ddindex2ddnodeid[ddindex];
    }
    else if (dd->comm->bCartesianPP)
    {
#if GMX_MPI
        MPI_Cart_rank(dd->mpi_comm_all, c, &ddnodeid);
#endif
    }
    else
    {
        ddnodeid = ddindex;
    }

    return ddnodeid;
}

static gmx_bool dynamic_dd_box(const gmx_ddbox_t *ddbox, const t_inputrec *ir)
{
    return (ddbox->nboundeddim < DIM || inputrecDynamicBox(ir));
}

int ddglatnr(const gmx_domdec_t *dd, int i)
{
    int atnr;

    if (dd == nullptr)
    {
        atnr = i + 1;
    }
    else
    {
        if (i >= dd->comm->atomRanges.numAtomsTotal())
        {
            gmx_fatal(FARGS, "glatnr called with %d, which is larger than the local number of atoms (%d)", i, dd->comm->atomRanges.numAtomsTotal());
        }
        atnr = dd->globalAtomIndices[i] + 1;
    }

    return atnr;
}

t_block *dd_charge_groups_global(gmx_domdec_t *dd)
{
    return &dd->comm->cgs_gl;
}

void dd_store_state(gmx_domdec_t *dd, t_state *state)
{
    int i;

    if (state->ddp_count != dd->ddp_count)
    {
        gmx_incons("The MD state does not match the domain decomposition state");
    }

    state->cg_gl.resize(dd->ncg_home);
    for (i = 0; i < dd->ncg_home; i++)
    {
        state->cg_gl[i] = dd->globalAtomGroupIndices[i];
    }

    state->ddp_count_cg_gl = dd->ddp_count;
}

gmx_domdec_zones_t *domdec_zones(gmx_domdec_t *dd)
{
    return &dd->comm->zones;
}

void dd_get_ns_ranges(const gmx_domdec_t *dd, int icg,
                      int *jcg0, int *jcg1, ivec shift0, ivec shift1)
{
    gmx_domdec_zones_t *zones;
    int                 izone, d, dim;

    zones = &dd->comm->zones;

    izone = 0;
    while (icg >= zones->izone[izone].cg1)
    {
        izone++;
    }

    if (izone == 0)
    {
        *jcg0 = icg;
    }
    else if (izone < zones->nizone)
    {
        *jcg0 = zones->izone[izone].jcg0;
    }
    else
    {
        gmx_fatal(FARGS, "DD icg %d out of range: izone (%d) >= nizone (%d)",
                  icg, izone, zones->nizone);
    }

    *jcg1 = zones->izone[izone].jcg1;

    for (d = 0; d < dd->ndim; d++)
    {
        dim         = dd->dim[d];
        shift0[dim] = zones->izone[izone].shift0[dim];
        shift1[dim] = zones->izone[izone].shift1[dim];
        if (dd->comm->tric_dir[dim] || (isDlbOn(dd->comm) && d > 0))
        {
            /* A conservative approach, this can be optimized */
            shift0[dim] -= 1;
            shift1[dim] += 1;
        }
    }
}

int dd_numHomeAtoms(const gmx_domdec_t &dd)
{
    return dd.comm->atomRanges.numHomeAtoms();
}

int dd_natoms_mdatoms(const gmx_domdec_t *dd)
{
    /* We currently set mdatoms entries for all atoms:
     * local + non-local + communicated for vsite + constraints
     */

    return dd->comm->atomRanges.numAtomsTotal();
}

int dd_natoms_vsite(const gmx_domdec_t *dd)
{
    return dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites);
}

void dd_get_constraint_range(const gmx_domdec_t *dd, int *at_start, int *at_end)
{
    *at_start = dd->comm->atomRanges.start(DDAtomRanges::Type::Constraints);
    *at_end   = dd->comm->atomRanges.end(DDAtomRanges::Type::Constraints);
}

void dd_move_x(gmx_domdec_t             *dd,
               matrix                    box,
               gmx::ArrayRef<gmx::RVec>  x,
               gmx_wallcycle            *wcycle)
{
    wallcycle_start(wcycle, ewcMOVEX);

    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    rvec                   shift = {0, 0, 0};
    gmx_bool               bPBC, bScrew;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = 1;
    nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        bPBC   = (dd->ci[dd->dim[d]] == 0);
        bScrew = (bPBC && dd->bScrewPBC && dd->dim[d] == XX);
        if (bPBC)
        {
            copy_rvec(box[dd->dim[d]], shift);
        }
        cd = &comm->cd[d];
        for (const gmx_domdec_ind_t &ind : cd->ind)
        {
            DDBufferAccess<gmx::RVec>  sendBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>  &sendBuffer = sendBufferAccess.buffer;
            int                        n          = 0;
            if (!bPBC)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        sendBuffer[n] = x[j];
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* We need to shift the coordinates */
                        for (int d = 0; d < DIM; d++)
                        {
                            sendBuffer[n][d] = x[j][d] + shift[d];
                        }
                        n++;
                    }
                }
            }
            else
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* Shift x */
                        sendBuffer[n][XX] = x[j][XX] + shift[XX];
                        /* Rotate y and z.
                         * This operation requires a special shift force
                         * treatment, which is performed in calc_vir.
                         */
                        sendBuffer[n][YY] = box[YY][YY] - x[j][YY];
                        sendBuffer[n][ZZ] = box[ZZ][ZZ] - x[j][ZZ];
                        n++;
                    }
                }
            }

            DDBufferAccess<gmx::RVec>  receiveBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec>   receiveBuffer;
            if (cd->receiveInPlace)
            {
                receiveBuffer = gmx::arrayRefFromArray(x.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                receiveBuffer = receiveBufferAccess.buffer;
            }
            /* Send and receive the coordinates */
            ddSendrecv(dd, d, dddirBackward,
                       sendBuffer, receiveBuffer);

            if (!cd->receiveInPlace)
            {
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        x[i] = receiveBuffer[j++];
                    }
                }
            }
            nat_tot += ind.nrecv[nzone+1];
        }
        nzone += nzone;
    }

    wallcycle_stop(wcycle, ewcMOVEX);
}

void dd_move_f(gmx_domdec_t             *dd,
               gmx::ArrayRef<gmx::RVec>  f,
               rvec                     *fshift,
               gmx_wallcycle            *wcycle)
{
    wallcycle_start(wcycle, ewcMOVEF);

    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    ivec                   vis;
    int                    is;
    gmx_bool               bShiftForcesNeedPbc, bScrew;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = comm->zones.n/2;
    nat_tot = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int d = dd->ndim-1; d >= 0; d--)
    {
        /* Only forces in domains near the PBC boundaries need to
           consider PBC in the treatment of fshift */
        bShiftForcesNeedPbc   = (dd->ci[dd->dim[d]] == 0);
        bScrew                = (bShiftForcesNeedPbc && dd->bScrewPBC && dd->dim[d] == XX);
        if (fshift == nullptr && !bScrew)
        {
            bShiftForcesNeedPbc = FALSE;
        }
        /* Determine which shift vector we need */
        clear_ivec(vis);
        vis[dd->dim[d]] = 1;
        is              = IVEC2IS(vis);

        cd = &comm->cd[d];
        for (int p = cd->numPulses() - 1; p >= 0; p--)
        {
            const gmx_domdec_ind_t    &ind  = cd->ind[p];
            DDBufferAccess<gmx::RVec>  receiveBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>  &receiveBuffer = receiveBufferAccess.buffer;

            nat_tot                        -= ind.nrecv[nzone+1];

            DDBufferAccess<gmx::RVec>  sendBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec>   sendBuffer;
            if (cd->receiveInPlace)
            {
                sendBuffer = gmx::arrayRefFromArray(f.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                sendBuffer = sendBufferAccess.buffer;
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        sendBuffer[j++] = f[i];
                    }
                }
            }
            /* Communicate the forces */
            ddSendrecv(dd, d, dddirForward,
                       sendBuffer, receiveBuffer);
            /* Add the received forces */
            int n = 0;
            if (!bShiftForcesNeedPbc)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        for (int d = 0; d < DIM; d++)
                        {
                            f[j][d] += receiveBuffer[n][d];
                        }
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                /* fshift should always be defined if this function is
                 * called when bShiftForcesNeedPbc is true */
                assert(nullptr != fshift);
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        for (int d = 0; d < DIM; d++)
                        {
                            f[j][d] += receiveBuffer[n][d];
                        }
                        /* Add this force to the shift force */
                        for (int d = 0; d < DIM; d++)
                        {
                            fshift[is][d] += receiveBuffer[n][d];
                        }
                        n++;
                    }
                }
            }
            else
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* Rotate the force */
                        f[j][XX] += receiveBuffer[n][XX];
                        f[j][YY] -= receiveBuffer[n][YY];
                        f[j][ZZ] -= receiveBuffer[n][ZZ];
                        if (fshift)
                        {
                            /* Add this force to the shift force */
                            for (int d = 0; d < DIM; d++)
                            {
                                fshift[is][d] += receiveBuffer[n][d];
                            }
                        }
                        n++;
                    }
                }
            }
        }
        nzone /= 2;
    }
    wallcycle_stop(wcycle, ewcMOVEF);
}

/* Convenience function for extracting a real buffer from an rvec buffer
 *
 * To reduce the number of temporary communication buffers and avoid
 * cache polution, we reuse gmx::RVec buffers for storing reals.
 * This functions return a real buffer reference with the same number
 * of elements as the gmx::RVec buffer (so 1/3 of the size in bytes).
 */
static gmx::ArrayRef<real>
realArrayRefFromRvecArrayRef(gmx::ArrayRef<gmx::RVec> arrayRef)
{
    return gmx::arrayRefFromArray(as_rvec_array(arrayRef.data())[0],
                                  arrayRef.size());
}

void dd_atom_spread_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = 1;
    nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        cd = &comm->cd[d];
        for (const gmx_domdec_ind_t &ind : cd->ind)
        {
            /* Note: We provision for RVec instead of real, so a factor of 3
             * more than needed. The buffer actually already has this size
             * and we pass a plain pointer below, so this does not matter.
             */
            DDBufferAccess<gmx::RVec> sendBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<real>       sendBuffer = realArrayRefFromRvecArrayRef(sendBufferAccess.buffer);
            int                       n          = 0;
            for (int g : ind.index)
            {
                for (int j : atomGrouping.block(g))
                {
                    sendBuffer[n++] = v[j];
                }
            }

            DDBufferAccess<gmx::RVec> receiveBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<real>       receiveBuffer;
            if (cd->receiveInPlace)
            {
                receiveBuffer = gmx::arrayRefFromArray(v + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                receiveBuffer = realArrayRefFromRvecArrayRef(receiveBufferAccess.buffer);
            }
            /* Send and receive the data */
            ddSendrecv(dd, d, dddirBackward,
                       sendBuffer, receiveBuffer);
            if (!cd->receiveInPlace)
            {
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        v[i] = receiveBuffer[j++];
                    }
                }
            }
            nat_tot += ind.nrecv[nzone+1];
        }
        nzone += nzone;
    }
}

void dd_atom_sum_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = comm->zones.n/2;
    nat_tot = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int d = dd->ndim-1; d >= 0; d--)
    {
        cd = &comm->cd[d];
        for (int p = cd->numPulses() - 1; p >= 0; p--)
        {
            const gmx_domdec_ind_t &ind = cd->ind[p];

            /* Note: We provision for RVec instead of real, so a factor of 3
             * more than needed. The buffer actually already has this size
             * and we typecast, so this works as intended.
             */
            DDBufferAccess<gmx::RVec> receiveBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<real>       receiveBuffer = realArrayRefFromRvecArrayRef(receiveBufferAccess.buffer);
            nat_tot -= ind.nrecv[nzone + 1];

            DDBufferAccess<gmx::RVec> sendBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<real>       sendBuffer;
            if (cd->receiveInPlace)
            {
                sendBuffer = gmx::arrayRefFromArray(v + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                sendBuffer = realArrayRefFromRvecArrayRef(sendBufferAccess.buffer);
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        sendBuffer[j++] = v[i];
                    }
                }
            }
            /* Communicate the forces */
            ddSendrecv(dd, d, dddirForward,
                       sendBuffer, receiveBuffer);
            /* Add the received forces */
            int n = 0;
            for (int g : ind.index)
            {
                for (int j : atomGrouping.block(g))
                {
                    v[j] += receiveBuffer[n];
                    n++;
                }
            }
        }
        nzone /= 2;
    }
}

static void print_ddzone(FILE *fp, int d, int i, int j, gmx_ddzone_t *zone)
{
    fprintf(fp, "zone d0 %d d1 %d d2 %d  min0 %6.3f max1 %6.3f mch0 %6.3f mch1 %6.3f p1_0 %6.3f p1_1 %6.3f\n",
            d, i, j,
            zone->min0, zone->max1,
            zone->mch0, zone->mch0,
            zone->p1_0, zone->p1_1);
}

/* Using the home grid size as input in cell_ns_x0 and cell_ns_x1
 * takes the extremes over all home and remote zones in the halo
 * and returns the results in cell_ns_x0 and cell_ns_x1.
 * Note: only used with the group cut-off scheme.
 */
static void dd_move_cellx(gmx_domdec_t      *dd,
                          const gmx_ddbox_t *ddbox,
                          rvec               cell_ns_x0,
                          rvec               cell_ns_x1)
{
    constexpr int      c_ddZoneCommMaxNumZones = 5;
    gmx_ddzone_t       buf_s[c_ddZoneCommMaxNumZones];
    gmx_ddzone_t       buf_r[c_ddZoneCommMaxNumZones];
    gmx_ddzone_t       buf_e[c_ddZoneCommMaxNumZones];
    gmx_domdec_comm_t *comm = dd->comm;

    rvec               extr_s[2];
    rvec               extr_r[2];
    for (int d = 1; d < dd->ndim; d++)
    {
        int           dim = dd->dim[d];
        gmx_ddzone_t &zp  = (d == 1) ? comm->zone_d1[0] : comm->zone_d2[0][0];

        /* Copy the base sizes of the home zone */
        zp.min0 = cell_ns_x0[dim];
        zp.max1 = cell_ns_x1[dim];
        zp.min1 = cell_ns_x1[dim];
        zp.mch0 = cell_ns_x0[dim];
        zp.mch1 = cell_ns_x1[dim];
        zp.p1_0 = cell_ns_x0[dim];
        zp.p1_1 = cell_ns_x1[dim];
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
        for (int d1 = d; d1 < dd->ndim-1; d1++)
        {
            /* We invert the order to be able to use the same loop for buf_e */
            buf_s[pos].min0 = extr_s[d1][1];
            buf_s[pos].max1 = extr_s[d1][0];
            buf_s[pos].min1 = extr_s[d1][2];
            buf_s[pos].mch0 = 0;
            buf_s[pos].mch1 = 0;
            /* Store the cell corner of the dimension we communicate along */
            buf_s[pos].p1_0 = comm->cell_x0[dim];
            buf_s[pos].p1_1 = 0;
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
            numPulsesMin = std::min(numPulses, dd->nc[dim] - 1 - numPulses);
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

            int  numElements      = dd->ndim - d - 1;
            ddSendrecv(dd, d, dddirForward,
                       extr_s + d, numElements,
                       extr_r + d, numElements);

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
            bool receiveValidData = (applyPbc || dd->ci[dim] < dd->nc[dim] - 1);

            static_assert(sizeof(gmx_ddzone_t) == c_ddzoneNumReals*sizeof(real), "Here we expect gmx_ddzone_t to consist of c_ddzoneNumReals reals (only)");

            int numReals = numElementsInBuffer*c_ddzoneNumReals;
            ddSendrecv(dd, d, dddirBackward,
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

                    int  c;
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
                    real det = (1 + c*c)*comm->cutoff*comm->cutoff - dist_d*dist_d;
                    if (det > 0)
                    {
                        dh[d1] = comm->cutoff - (c*dist_d + std::sqrt(det))/(1 + c*c);
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
                        buf_e[i].mch0 = std::max(buf_e[i].mch0, buf_r[i].mch0-dh[d1]);
                        buf_e[i].mch1 = std::max(buf_e[i].mch1, buf_r[i].mch1-dh[d1]);
                    }
                }
                /* Copy the received buffer to the send buffer,
                 * to pass the data through with the next pulse.
                 */
                buf_s[i] = buf_r[i];
            }
            if (((applyPbc || dd->ci[dim] + numPulses < dd->nc[dim]) && pulse == numPulses - 1) ||
                (!applyPbc && dd->ci[dim] + 1 + pulse == dd->nc[dim] - 1))
            {
                /* Store the extremes */
                int pos = 0;

                for (int d1 = d; d1 < dd->ndim-1; d1++)
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
                        comm->zone_d2[1-d][i] = buf_e[pos];
                        pos++;
                    }
                }
                if (d == 0)
                {
                    comm->zone_d1[1] = buf_e[pos];
                    pos++;
                }
            }
        }
    }

    if (dd->ndim >= 2)
    {
        int dim = dd->dim[1];
        for (int i = 0; i < 2; i++)
        {
            if (debug)
            {
                print_ddzone(debug, 1, i, 0, &comm->zone_d1[i]);
            }
            cell_ns_x0[dim] = std::min(cell_ns_x0[dim], comm->zone_d1[i].min0);
            cell_ns_x1[dim] = std::max(cell_ns_x1[dim], comm->zone_d1[i].max1);
        }
    }
    if (dd->ndim >= 3)
    {
        int dim = dd->dim[2];
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
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
    for (int d = 1; d < dd->ndim; d++)
    {
        cellsizes[d].fracLowerMax = extr_s[d-1][0];
        cellsizes[d].fracUpperMin = extr_s[d-1][1];
        if (debug)
        {
            fprintf(debug, "Cell fraction d %d, max0 %f, min1 %f\n",
                    d, cellsizes[d].fracLowerMax, cellsizes[d].fracUpperMin);
        }
    }
}

static void write_dd_grid_pdb(const char *fn, int64_t step,
                              gmx_domdec_t *dd, matrix box, gmx_ddbox_t *ddbox)
{
    rvec   grid_s[2], *grid_r = nullptr, cx, r;
    char   fname[STRLEN], buf[22];
    FILE  *out;
    int    a, i, d, z, y, x;
    matrix tric;
    real   vol;

    copy_rvec(dd->comm->cell_x0, grid_s[0]);
    copy_rvec(dd->comm->cell_x1, grid_s[1]);

    if (DDMASTER(dd))
    {
        snew(grid_r, 2*dd->nnodes);
    }

    dd_gather(dd, 2*sizeof(rvec), grid_s, DDMASTER(dd) ? grid_r : nullptr);

    if (DDMASTER(dd))
    {
        for (d = 0; d < DIM; d++)
        {
            for (i = 0; i < DIM; i++)
            {
                if (d == i)
                {
                    tric[d][i] = 1;
                }
                else
                {
                    if (d < ddbox->npbcdim && dd->nc[d] > 1)
                    {
                        tric[d][i] = box[i][d]/box[i][i];
                    }
                    else
                    {
                        tric[d][i] = 0;
                    }
                }
            }
        }
        sprintf(fname, "%s_%s.pdb", fn, gmx_step_str(step, buf));
        out = gmx_fio_fopen(fname, "w");
        gmx_write_pdb_box(out, dd->bScrewPBC ? epbcSCREW : epbcXYZ, box);
        a = 1;
        for (i = 0; i < dd->nnodes; i++)
        {
            vol = dd->nnodes/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
            for (d = 0; d < DIM; d++)
            {
                vol *= grid_r[i*2+1][d] - grid_r[i*2][d];
            }
            for (z = 0; z < 2; z++)
            {
                for (y = 0; y < 2; y++)
                {
                    for (x = 0; x < 2; x++)
                    {
                        cx[XX] = grid_r[i*2+x][XX];
                        cx[YY] = grid_r[i*2+y][YY];
                        cx[ZZ] = grid_r[i*2+z][ZZ];
                        mvmul(tric, cx, r);
                        gmx_fprintf_pdb_atomline(out, epdbATOM, a++, "CA", ' ', "GLY", ' ', i+1, ' ',
                                                 10*r[XX], 10*r[YY], 10*r[ZZ], 1.0, vol, "");
                    }
                }
            }
            for (d = 0; d < DIM; d++)
            {
                for (x = 0; x < 4; x++)
                {
                    switch (d)
                    {
                        case 0: y = 1 + i*8 + 2*x; break;
                        case 1: y = 1 + i*8 + 2*x - (x % 2); break;
                        case 2: y = 1 + i*8 + x; break;
                    }
                    fprintf(out, "%6s%5d%5d\n", "CONECT", y, y+(1<<d));
                }
            }
        }
        gmx_fio_fclose(out);
        sfree(grid_r);
    }
}

void write_dd_pdb(const char *fn, int64_t step, const char *title,
                  const gmx_mtop_t *mtop, const t_commrec *cr,
                  int natoms, const rvec x[], const matrix box)
{
    char          fname[STRLEN], buf[22];
    FILE         *out;
    int           resnr;
    const char   *atomname, *resname;
    gmx_domdec_t *dd;

    dd = cr->dd;
    if (natoms == -1)
    {
        natoms = dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites);
    }

    sprintf(fname, "%s_%s_n%d.pdb", fn, gmx_step_str(step, buf), cr->sim_nodeid);

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, dd->bScrewPBC ? epbcSCREW : epbcXYZ, box);
    int molb = 0;
    for (int i = 0; i < natoms; i++)
    {
        int  ii = dd->globalAtomIndices[i];
        mtopGetAtomAndResidueName(mtop, ii, &molb, &atomname, &resnr, &resname, nullptr);
        int  c;
        real b;
        if (i < dd->comm->atomRanges.end(DDAtomRanges::Type::Zones))
        {
            c = 0;
            while (i >= dd->atomGrouping().subRange(0, dd->comm->zones.cg_range[c + 1]).end())
            {
                c++;
            }
            b = c;
        }
        else if (i < dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites))
        {
            b = dd->comm->zones.n;
        }
        else
        {
            b = dd->comm->zones.n + 1;
        }
        gmx_fprintf_pdb_atomline(out, epdbATOM, ii+1, atomname, ' ', resname, ' ', resnr, ' ',
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ], 1.0, b, "");
    }
    fprintf(out, "TER\n");

    gmx_fio_fclose(out);
}

real dd_cutoff_multibody(const gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm;
    int                di;
    real               r;

    comm = dd->comm;

    r = -1;
    if (comm->bInterCGBondeds)
    {
        if (comm->cutoff_mbody > 0)
        {
            r = comm->cutoff_mbody;
        }
        else
        {
            /* cutoff_mbody=0 means we do not have DLB */
            r = comm->cellsize_min[dd->dim[0]];
            for (di = 1; di < dd->ndim; di++)
            {
                r = std::min(r, comm->cellsize_min[dd->dim[di]]);
            }
            if (comm->bBondComm)
            {
                r = std::max(r, comm->cutoff_mbody);
            }
            else
            {
                r = std::min(r, comm->cutoff);
            }
        }
    }

    return r;
}

real dd_cutoff_twobody(const gmx_domdec_t *dd)
{
    real r_mb;

    r_mb = dd_cutoff_multibody(dd);

    return std::max(dd->comm->cutoff, r_mb);
}


static void dd_cart_coord2pmecoord(const gmx_domdec_t *dd, const ivec coord,
                                   ivec coord_pme)
{
    int nc, ntot;

    nc   = dd->nc[dd->comm->cartpmedim];
    ntot = dd->comm->ntot[dd->comm->cartpmedim];
    copy_ivec(coord, coord_pme);
    coord_pme[dd->comm->cartpmedim] =
        nc + (coord[dd->comm->cartpmedim]*(ntot - nc) + (ntot - nc)/2)/nc;
}

static int ddindex2pmeindex(const gmx_domdec_t *dd, int ddindex)
{
    int npp, npme;

    npp  = dd->nnodes;
    npme = dd->comm->npmenodes;

    /* Here we assign a PME node to communicate with this DD node
     * by assuming that the major index of both is x.
     * We add cr->npmenodes/2 to obtain an even distribution.
     */
    return (ddindex*npme + npme/2)/npp;
}

static int *dd_interleaved_pme_ranks(const gmx_domdec_t *dd)
{
    int *pme_rank;
    int  n, i, p0, p1;

    snew(pme_rank, dd->comm->npmenodes);
    n = 0;
    for (i = 0; i < dd->nnodes; i++)
    {
        p0 = ddindex2pmeindex(dd, i);
        p1 = ddindex2pmeindex(dd, i+1);
        if (i+1 == dd->nnodes || p1 > p0)
        {
            if (debug)
            {
                fprintf(debug, "pme_rank[%d] = %d\n", n, i+1+n);
            }
            pme_rank[n] = i + 1 + n;
            n++;
        }
    }

    return pme_rank;
}

static int gmx_ddcoord2pmeindex(const t_commrec *cr, int x, int y, int z)
{
    gmx_domdec_t *dd;
    ivec          coords;
    int           slab;

    dd = cr->dd;
    /*
       if (dd->comm->bCartesian) {
       gmx_ddindex2xyz(dd->nc,ddindex,coords);
       dd_coords2pmecoords(dd,coords,coords_pme);
       copy_ivec(dd->ntot,nc);
       nc[dd->cartpmedim]         -= dd->nc[dd->cartpmedim];
       coords_pme[dd->cartpmedim] -= dd->nc[dd->cartpmedim];

       slab = (coords_pme[XX]*nc[YY] + coords_pme[YY])*nc[ZZ] + coords_pme[ZZ];
       } else {
       slab = (ddindex*cr->npmenodes + cr->npmenodes/2)/dd->nnodes;
       }
     */
    coords[XX] = x;
    coords[YY] = y;
    coords[ZZ] = z;
    slab       = ddindex2pmeindex(dd, dd_index(dd->nc, coords));

    return slab;
}

static int ddcoord2simnodeid(const t_commrec *cr, int x, int y, int z)
{
    gmx_domdec_comm_t *comm;
    ivec               coords;
    int                ddindex, nodeid = -1;

    comm = cr->dd->comm;

    coords[XX] = x;
    coords[YY] = y;
    coords[ZZ] = z;
    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        MPI_Cart_rank(cr->mpi_comm_mysim, coords, &nodeid);
#endif
    }
    else
    {
        ddindex = dd_index(cr->dd->nc, coords);
        if (comm->bCartesianPP)
        {
            nodeid = comm->ddindex2simnodeid[ddindex];
        }
        else
        {
            if (comm->pmenodes)
            {
                nodeid = ddindex + gmx_ddcoord2pmeindex(cr, x, y, z);
            }
            else
            {
                nodeid = ddindex;
            }
        }
    }

    return nodeid;
}

static int dd_simnode2pmenode(const gmx_domdec_t         *dd,
                              const t_commrec gmx_unused *cr,
                              int                         sim_nodeid)
{
    int pmenode = -1;

    const gmx_domdec_comm_t *comm = dd->comm;

    /* This assumes a uniform x domain decomposition grid cell size */
    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        ivec coord, coord_pme;
        MPI_Cart_coords(cr->mpi_comm_mysim, sim_nodeid, DIM, coord);
        if (coord[comm->cartpmedim] < dd->nc[comm->cartpmedim])
        {
            /* This is a PP node */
            dd_cart_coord2pmecoord(dd, coord, coord_pme);
            MPI_Cart_rank(cr->mpi_comm_mysim, coord_pme, &pmenode);
        }
#endif
    }
    else if (comm->bCartesianPP)
    {
        if (sim_nodeid < dd->nnodes)
        {
            pmenode = dd->nnodes + ddindex2pmeindex(dd, sim_nodeid);
        }
    }
    else
    {
        /* This assumes DD cells with identical x coordinates
         * are numbered sequentially.
         */
        if (dd->comm->pmenodes == nullptr)
        {
            if (sim_nodeid < dd->nnodes)
            {
                /* The DD index equals the nodeid */
                pmenode = dd->nnodes + ddindex2pmeindex(dd, sim_nodeid);
            }
        }
        else
        {
            int i = 0;
            while (sim_nodeid > dd->comm->pmenodes[i])
            {
                i++;
            }
            if (sim_nodeid < dd->comm->pmenodes[i])
            {
                pmenode = dd->comm->pmenodes[i];
            }
        }
    }

    return pmenode;
}

NumPmeDomains getNumPmeDomains(const gmx_domdec_t *dd)
{
    if (dd != nullptr)
    {
        return {
                   dd->comm->npmenodes_x, dd->comm->npmenodes_y
        };
    }
    else
    {
        return {
                   1, 1
        };
    }
}

std::vector<int> get_pme_ddranks(const t_commrec *cr, int pmenodeid)
{
    gmx_domdec_t *dd;
    int           x, y, z;
    ivec          coord, coord_pme;

    dd = cr->dd;

    std::vector<int> ddranks;
    ddranks.reserve((dd->nnodes+cr->npmenodes-1)/cr->npmenodes);

    for (x = 0; x < dd->nc[XX]; x++)
    {
        for (y = 0; y < dd->nc[YY]; y++)
        {
            for (z = 0; z < dd->nc[ZZ]; z++)
            {
                if (dd->comm->bCartesianPP_PME)
                {
                    coord[XX] = x;
                    coord[YY] = y;
                    coord[ZZ] = z;
                    dd_cart_coord2pmecoord(dd, coord, coord_pme);
                    if (dd->ci[XX] == coord_pme[XX] &&
                        dd->ci[YY] == coord_pme[YY] &&
                        dd->ci[ZZ] == coord_pme[ZZ])
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
                else
                {
                    /* The slab corresponds to the nodeid in the PME group */
                    if (gmx_ddcoord2pmeindex(cr, x, y, z) == pmenodeid)
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
            }
        }
    }
    return ddranks;
}

static gmx_bool receive_vir_ener(const gmx_domdec_t *dd, const t_commrec *cr)
{
    gmx_bool bReceive = TRUE;

    if (cr->npmenodes < dd->nnodes)
    {
        gmx_domdec_comm_t *comm = dd->comm;
        if (comm->bCartesianPP_PME)
        {
#if GMX_MPI
            int  pmenode = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
            ivec coords;
            MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, coords);
            coords[comm->cartpmedim]++;
            if (coords[comm->cartpmedim] < dd->nc[comm->cartpmedim])
            {
                int rank;
                MPI_Cart_rank(cr->mpi_comm_mysim, coords, &rank);
                if (dd_simnode2pmenode(dd, cr, rank) == pmenode)
                {
                    /* This is not the last PP node for pmenode */
                    bReceive = FALSE;
                }
            }
#else
            GMX_RELEASE_ASSERT(false, "Without MPI we should not have Cartesian PP-PME with #PMEnodes < #DDnodes");
#endif
        }
        else
        {
            int pmenode = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
            if (cr->sim_nodeid+1 < cr->nnodes &&
                dd_simnode2pmenode(dd, cr, cr->sim_nodeid+1) == pmenode)
            {
                /* This is not the last PP node for pmenode */
                bReceive = FALSE;
            }
        }
    }

    return bReceive;
}

static void set_zones_ncg_home(gmx_domdec_t *dd)
{
    gmx_domdec_zones_t *zones;
    int                 i;

    zones = &dd->comm->zones;

    zones->cg_range[0] = 0;
    for (i = 1; i < zones->n+1; i++)
    {
        zones->cg_range[i] = dd->ncg_home;
    }
    /* zone_ncg1[0] should always be equal to ncg_home */
    dd->comm->zone_ncg1[0] = dd->ncg_home;
}

static void restoreAtomGroups(gmx_domdec_t *dd,
                              const int *gcgs_index, const t_state *state)
{
    gmx::ArrayRef<const int>  atomGroupsState        = state->cg_gl;

    std::vector<int>         &globalAtomGroupIndices = dd->globalAtomGroupIndices;
    gmx::RangePartitioning   &atomGrouping           = dd->atomGrouping_;

    globalAtomGroupIndices.resize(atomGroupsState.size());
    atomGrouping.clear();

    /* Copy back the global charge group indices from state
     * and rebuild the local charge group to atom index.
     */
    for (gmx::index i = 0; i < atomGroupsState.size(); i++)
    {
        const int atomGroupGlobal  = atomGroupsState[i];
        const int groupSize        = gcgs_index[atomGroupGlobal + 1] - gcgs_index[atomGroupGlobal];
        globalAtomGroupIndices[i]  = atomGroupGlobal;
        atomGrouping.appendBlock(groupSize);
    }

    dd->ncg_home = atomGroupsState.size();
    dd->comm->atomRanges.setEnd(DDAtomRanges::Type::Home, atomGrouping.fullRange().end());

    set_zones_ncg_home(dd);
}

static void dd_set_cginfo(gmx::ArrayRef<const int> index_gl, int cg0, int cg1,
                          t_forcerec *fr, char *bLocalCG)
{
    cginfo_mb_t *cginfo_mb;
    int         *cginfo;
    int          cg;

    if (fr != nullptr)
    {
        cginfo_mb = fr->cginfo_mb;
        cginfo    = fr->cginfo;

        for (cg = cg0; cg < cg1; cg++)
        {
            cginfo[cg] = ddcginfo(cginfo_mb, index_gl[cg]);
        }
    }

    if (bLocalCG != nullptr)
    {
        for (cg = cg0; cg < cg1; cg++)
        {
            bLocalCG[index_gl[cg]] = TRUE;
        }
    }
}

static void make_dd_indices(gmx_domdec_t *dd,
                            const int *gcgs_index, int cg_start)
{
    const int                 numZones               = dd->comm->zones.n;
    const int                *zone2cg                = dd->comm->zones.cg_range;
    const int                *zone_ncg1              = dd->comm->zone_ncg1;
    gmx::ArrayRef<const int>  globalAtomGroupIndices = dd->globalAtomGroupIndices;
    const gmx_bool            bCGs                   = dd->comm->bCGs;

    std::vector<int>         &globalAtomIndices      = dd->globalAtomIndices;

    if (zone2cg[1] != dd->ncg_home)
    {
        gmx_incons("dd->ncg_zone is not up to date");
    }

    /* Make the local to global and global to local atom index */
    int a = dd->atomGrouping().subRange(cg_start, cg_start).begin();
    globalAtomIndices.resize(a);
    for (int zone = 0; zone < numZones; zone++)
    {
        int cg0;
        if (zone == 0)
        {
            cg0 = cg_start;
        }
        else
        {
            cg0 = zone2cg[zone];
        }
        int cg1    = zone2cg[zone+1];
        int cg1_p1 = cg0 + zone_ncg1[zone];

        for (int cg = cg0; cg < cg1; cg++)
        {
            int zone1 = zone;
            if (cg >= cg1_p1)
            {
                /* Signal that this cg is from more than one pulse away */
                zone1 += numZones;
            }
            int cg_gl = globalAtomGroupIndices[cg];
            if (bCGs)
            {
                for (int a_gl = gcgs_index[cg_gl]; a_gl < gcgs_index[cg_gl+1]; a_gl++)
                {
                    globalAtomIndices.push_back(a_gl);
                    ga2la_set(dd->ga2la, a_gl, a, zone1);
                    a++;
                }
            }
            else
            {
                globalAtomIndices.push_back(cg_gl);
                ga2la_set(dd->ga2la, cg_gl, a, zone1);
                a++;
            }
        }
    }
}

static int check_bLocalCG(gmx_domdec_t *dd, int ncg_sys, const char *bLocalCG,
                          const char *where)
{
    int nerr = 0;
    if (bLocalCG == nullptr)
    {
        return nerr;
    }
    for (size_t i = 0; i < dd->globalAtomGroupIndices.size(); i++)
    {
        if (!bLocalCG[dd->globalAtomGroupIndices[i]])
        {
            fprintf(stderr,
                    "DD rank %d, %s: atom group %zu, global atom group %d is not marked in bLocalCG (ncg_home %d)\n", dd->rank, where, i + 1, dd->globalAtomGroupIndices[i] + 1, dd->ncg_home);
            nerr++;
        }
    }
    size_t ngl = 0;
    for (int i = 0; i < ncg_sys; i++)
    {
        if (bLocalCG[i])
        {
            ngl++;
        }
    }
    if (ngl != dd->globalAtomGroupIndices.size())
    {
        fprintf(stderr, "DD rank %d, %s: In bLocalCG %zu atom groups are marked as local, whereas there are %zu\n", dd->rank, where, ngl, dd->globalAtomGroupIndices.size());
        nerr++;
    }

    return nerr;
}

static void check_index_consistency(gmx_domdec_t *dd,
                                    int natoms_sys, int ncg_sys,
                                    const char *where)
{
    int       nerr = 0;

    const int numAtomsInZones = dd->comm->atomRanges.end(DDAtomRanges::Type::Zones);

    if (dd->comm->DD_debug > 1)
    {
        std::vector<int> have(natoms_sys);
        for (int a = 0; a < numAtomsInZones; a++)
        {
            int globalAtomIndex = dd->globalAtomIndices[a];
            if (have[globalAtomIndex] > 0)
            {
                fprintf(stderr, "DD rank %d: global atom %d occurs twice: index %d and %d\n", dd->rank, globalAtomIndex + 1, have[globalAtomIndex], a+1);
            }
            else
            {
                have[globalAtomIndex] = a + 1;
            }
        }
    }

    std::vector<int> have(numAtomsInZones);

    int              ngl = 0;
    for (int i = 0; i < natoms_sys; i++)
    {
        int a;
        int cell;
        if (ga2la_get(dd->ga2la, i, &a, &cell))
        {
            if (a >= numAtomsInZones)
            {
                fprintf(stderr, "DD rank %d: global atom %d marked as local atom %d, which is larger than nat_tot (%d)\n", dd->rank, i+1, a+1, numAtomsInZones);
                nerr++;
            }
            else
            {
                have[a] = 1;
                if (dd->globalAtomIndices[a] != i)
                {
                    fprintf(stderr, "DD rank %d: global atom %d marked as local atom %d, which has global atom index %d\n", dd->rank, i+1, a+1, dd->globalAtomIndices[a]+1);
                    nerr++;
                }
            }
            ngl++;
        }
    }
    if (ngl != numAtomsInZones)
    {
        fprintf(stderr,
                "DD rank %d, %s: %d global atom indices, %d local atoms\n",
                dd->rank, where, ngl, numAtomsInZones);
    }
    for (int a = 0; a < numAtomsInZones; a++)
    {
        if (have[a] == 0)
        {
            fprintf(stderr,
                    "DD rank %d, %s: local atom %d, global %d has no global index\n",
                    dd->rank, where, a + 1, dd->globalAtomIndices[a] + 1);
        }
    }

    nerr += check_bLocalCG(dd, ncg_sys, dd->comm->bLocalCG, where);

    if (nerr > 0)
    {
        gmx_fatal(FARGS, "DD rank %d, %s: %d atom(group) index inconsistencies",
                  dd->rank, where, nerr);
    }
}

/* Clear all DD global state indices, starting from \p atomGroupStart and \p atomStart */
static void clearDDStateIndices(gmx_domdec_t *dd,
                                int           atomGroupStart,
                                int           atomStart)
{
    if (atomStart == 0)
    {
        /* Clear the whole list without searching */
        ga2la_clear(dd->ga2la);
    }
    else
    {
        const int numAtomsInZones = dd->comm->atomRanges.end(DDAtomRanges::Type::Zones);
        for (int i = 0; i < numAtomsInZones; i++)
        {
            ga2la_del(dd->ga2la, dd->globalAtomIndices[i]);
        }
    }

    char *bLocalCG = dd->comm->bLocalCG;
    if (bLocalCG)
    {
        for (size_t atomGroup = atomGroupStart; atomGroup < dd->globalAtomGroupIndices.size(); atomGroup++)
        {
            bLocalCG[dd->globalAtomGroupIndices[atomGroup]] = FALSE;
        }
    }

    dd_clear_local_vsite_indices(dd);

    if (dd->constraints)
    {
        dd_clear_local_constraint_indices(dd);
    }
}

static bool check_grid_jump(int64_t             step,
                            const gmx_domdec_t *dd,
                            real                cutoff,
                            const gmx_ddbox_t  *ddbox,
                            gmx_bool            bFatal)
{
    gmx_domdec_comm_t *comm    = dd->comm;
    bool               invalid = false;

    for (int d = 1; d < dd->ndim; d++)
    {
        const DDCellsizesWithDlb &cellsizes = comm->cellsizesWithDlb[d];
        const int                 dim       = dd->dim[d];
        const real                limit     = grid_jump_limit(comm, cutoff, d);
        real                      bfac      = ddbox->box_size[dim];
        if (ddbox->tric_dir[dim])
        {
            bfac *= ddbox->skew_fac[dim];
        }
        if ((cellsizes.fracUpper - cellsizes.fracLowerMax)*bfac <  limit ||
            (cellsizes.fracLower - cellsizes.fracUpperMin)*bfac > -limit)
        {
            invalid = true;

            if (bFatal)
            {
                char buf[22];

                /* This error should never be triggered under normal
                 * circumstances, but you never know ...
                 */
                gmx_fatal(FARGS, "step %s: The domain decomposition grid has shifted too much in the %c-direction around cell %d %d %d. This should not have happened. Running with fewer ranks might avoid this issue.",
                          gmx_step_str(step, buf),
                          dim2char(dim), dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
            }
        }
    }

    return invalid;
}

static float dd_force_load(gmx_domdec_comm_t *comm)
{
    float load;

    if (comm->eFlop)
    {
        load = comm->flop;
        if (comm->eFlop > 1)
        {
            load *= 1.0 + (comm->eFlop - 1)*(0.1*rand()/RAND_MAX - 0.05);
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
                gpu_wait *= (comm->cycl_n[ddCyclF] - 1)/(float)comm->cycl_n[ddCyclF];
            }
            /* Sum the wait times over the ranks that share the same GPU */
            MPI_Allreduce(&gpu_wait, &gpu_wait_sum, 1, MPI_FLOAT, MPI_SUM,
                          comm->mpi_comm_gpu_shared);
            /* Replace the wait time by the average over the ranks */
            load += -gpu_wait + gpu_wait_sum/comm->nrank_gpu_shared;
        }
#endif
    }

    return load;
}

static void set_slb_pme_dim_f(gmx_domdec_t *dd, int dim, real **dim_f)
{
    gmx_domdec_comm_t *comm;
    int                i;

    comm = dd->comm;

    snew(*dim_f, dd->nc[dim]+1);
    (*dim_f)[0] = 0;
    for (i = 1; i < dd->nc[dim]; i++)
    {
        if (comm->slb_frac[dim])
        {
            (*dim_f)[i] = (*dim_f)[i-1] + comm->slb_frac[dim][i-1];
        }
        else
        {
            (*dim_f)[i] = (real)i/(real)dd->nc[dim];
        }
    }
    (*dim_f)[dd->nc[dim]] = 1;
}

static void init_ddpme(gmx_domdec_t *dd, gmx_ddpme_t *ddpme, int dimind)
{
    int  pmeindex, slab, nso, i;
    ivec xyz;

    if (dimind == 0 && dd->dim[0] == YY && dd->comm->npmenodes_x == 1)
    {
        ddpme->dim = YY;
    }
    else
    {
        ddpme->dim = dimind;
    }
    ddpme->dim_match = (ddpme->dim == dd->dim[dimind]);

    ddpme->nslab = (ddpme->dim == 0 ?
                    dd->comm->npmenodes_x :
                    dd->comm->npmenodes_y);

    if (ddpme->nslab <= 1)
    {
        return;
    }

    nso = dd->comm->npmenodes/ddpme->nslab;
    /* Determine for each PME slab the PP location range for dimension dim */
    snew(ddpme->pp_min, ddpme->nslab);
    snew(ddpme->pp_max, ddpme->nslab);
    for (slab = 0; slab < ddpme->nslab; slab++)
    {
        ddpme->pp_min[slab] = dd->nc[dd->dim[dimind]] - 1;
        ddpme->pp_max[slab] = 0;
    }
    for (i = 0; i < dd->nnodes; i++)
    {
        ddindex2xyz(dd->nc, i, xyz);
        /* For y only use our y/z slab.
         * This assumes that the PME x grid size matches the DD grid size.
         */
        if (dimind == 0 || xyz[XX] == dd->ci[XX])
        {
            pmeindex = ddindex2pmeindex(dd, i);
            if (dimind == 0)
            {
                slab = pmeindex/nso;
            }
            else
            {
                slab = pmeindex % ddpme->nslab;
            }
            ddpme->pp_min[slab] = std::min(ddpme->pp_min[slab], xyz[dimind]);
            ddpme->pp_max[slab] = std::max(ddpme->pp_max[slab], xyz[dimind]);
        }
    }

    set_slb_pme_dim_f(dd, ddpme->dim, &ddpme->slb_dim_f);
}

int dd_pme_maxshift_x(const gmx_domdec_t *dd)
{
    if (dd->comm->ddpme[0].dim == XX)
    {
        return dd->comm->ddpme[0].maxshift;
    }
    else
    {
        return 0;
    }
}

int dd_pme_maxshift_y(const gmx_domdec_t *dd)
{
    if (dd->comm->ddpme[0].dim == YY)
    {
        return dd->comm->ddpme[0].maxshift;
    }
    else if (dd->comm->npmedecompdim >= 2 && dd->comm->ddpme[1].dim == YY)
    {
        return dd->comm->ddpme[1].maxshift;
    }
    else
    {
        return 0;
    }
}

static void comm_dd_ns_cell_sizes(gmx_domdec_t *dd,
                                  gmx_ddbox_t *ddbox,
                                  rvec cell_ns_x0, rvec cell_ns_x1,
                                  int64_t step)
{
    gmx_domdec_comm_t *comm;
    int                dim_ind, dim;

    comm = dd->comm;

    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        dim = dd->dim[dim_ind];

        /* Without PBC we don't have restrictions on the outer cells */
        if (!(dim >= ddbox->npbcdim &&
              (dd->ci[dim] == 0 || dd->ci[dim] == dd->nc[dim] - 1)) &&
            isDlbOn(comm) &&
            (comm->cell_x1[dim] - comm->cell_x0[dim])*ddbox->skew_fac[dim] <
            comm->cellsize_min[dim])
        {
            char buf[22];
            gmx_fatal(FARGS, "step %s: The %c-size (%f) times the triclinic skew factor (%f) is smaller than the smallest allowed cell size (%f) for domain decomposition grid cell %d %d %d",
                      gmx_step_str(step, buf), dim2char(dim),
                      comm->cell_x1[dim] - comm->cell_x0[dim],
                      ddbox->skew_fac[dim],
                      dd->comm->cellsize_min[dim],
                      dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
        }
    }

    if ((isDlbOn(dd->comm) && dd->ndim > 1) || ddbox->nboundeddim < DIM)
    {
        /* Communicate the boundaries and update cell_ns_x0/1 */
        dd_move_cellx(dd, ddbox, cell_ns_x0, cell_ns_x1);
        if (isDlbOn(dd->comm) && dd->ndim > 1)
        {
            check_grid_jump(step, dd, dd->comm->cutoff, ddbox, TRUE);
        }
    }
}

void dd_cycles_add(const gmx_domdec_t *dd, float cycles, int ddCycl)
{
    /* Note that the cycles value can be incorrect, either 0 or some
     * extremely large value, when our thread migrated to another core
     * with an unsynchronized cycle counter. If this happens less often
     * that once per nstlist steps, this will not cause issues, since
     * we later subtract the maximum value from the sum over nstlist steps.
     * A zero count will slightly lower the total, but that's a small effect.
     * Note that the main purpose of the subtraction of the maximum value
     * is to avoid throwing off the load balancing when stalls occur due
     * e.g. system activity or network congestion.
     */
    dd->comm->cycl[ddCycl] += cycles;
    dd->comm->cycl_n[ddCycl]++;
    if (cycles > dd->comm->cycl_max[ddCycl])
    {
        dd->comm->cycl_max[ddCycl] = cycles;
    }
}

static double force_flop_count(t_nrnb *nrnb)
{
    int         i;
    double      sum;
    const char *name;

    sum = 0;
    for (i = 0; i < eNR_NBKERNEL_FREE_ENERGY; i++)
    {
        /* To get closer to the real timings, we half the count
         * for the normal loops and again half it for water loops.
         */
        name = nrnb_str(i);
        if (strstr(name, "W3") != nullptr || strstr(name, "W4") != nullptr)
        {
            sum += nrnb->n[i]*0.25*cost_nrnb(i);
        }
        else
        {
            sum += nrnb->n[i]*0.50*cost_nrnb(i);
        }
    }
    for (i = eNR_NBKERNEL_FREE_ENERGY; i <= eNR_NB14; i++)
    {
        name = nrnb_str(i);
        if (strstr(name, "W3") != nullptr || strstr(name, "W4") != nullptr)
        {
            sum += nrnb->n[i]*cost_nrnb(i);
        }
    }
    for (i = eNR_BONDS; i <= eNR_WALLS; i++)
    {
        sum += nrnb->n[i]*cost_nrnb(i);
    }

    return sum;
}

void dd_force_flop_start(gmx_domdec_t *dd, t_nrnb *nrnb)
{
    if (dd->comm->eFlop)
    {
        dd->comm->flop -= force_flop_count(nrnb);
    }
}
void dd_force_flop_stop(gmx_domdec_t *dd, t_nrnb *nrnb)
{
    if (dd->comm->eFlop)
    {
        dd->comm->flop += force_flop_count(nrnb);
        dd->comm->flop_n++;
    }
}

static void clear_dd_cycle_counts(gmx_domdec_t *dd)
{
    int i;

    for (i = 0; i < ddCyclNr; i++)
    {
        dd->comm->cycl[i]     = 0;
        dd->comm->cycl_n[i]   = 0;
        dd->comm->cycl_max[i] = 0;
    }
    dd->comm->flop   = 0;
    dd->comm->flop_n = 0;
}

static void get_load_distribution(gmx_domdec_t *dd, gmx_wallcycle_t wcycle)
{
    gmx_domdec_comm_t *comm;
    domdec_load_t     *load;
    float              cell_frac = 0, sbuf[DD_NLOAD_MAX];
    gmx_bool           bSepPME;

    if (debug)
    {
        fprintf(debug, "get_load_distribution start\n");
    }

    wallcycle_start(wcycle, ewcDDCOMMLOAD);

    comm = dd->comm;

    bSepPME = (dd->pme_nodeid >= 0);

    if (dd->ndim == 0 && bSepPME)
    {
        /* Without decomposition, but with PME nodes, we need the load */
        comm->load[0].mdf = comm->cycl[ddCyclPPduringPME];
        comm->load[0].pme = comm->cycl[ddCyclPME];
    }

    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        const DDCellsizesWithDlb &cellsizes = comm->cellsizesWithDlb[d];
        const int                 dim       = dd->dim[d];
        /* Check if we participate in the communication in this dimension */
        if (d == dd->ndim-1 ||
            (dd->ci[dd->dim[d+1]] == 0 && dd->ci[dd->dim[dd->ndim-1]] == 0))
        {
            load = &comm->load[d];
            if (isDlbOn(dd->comm))
            {
                cell_frac = cellsizes.fracUpper - cellsizes.fracLower;
            }
            int pos = 0;
            if (d == dd->ndim-1)
            {
                sbuf[pos++] = dd_force_load(comm);
                sbuf[pos++] = sbuf[0];
                if (isDlbOn(dd->comm))
                {
                    sbuf[pos++] = sbuf[0];
                    sbuf[pos++] = cell_frac;
                    if (d > 0)
                    {
                        sbuf[pos++] = cellsizes.fracLowerMax;
                        sbuf[pos++] = cellsizes.fracUpperMin;
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
                sbuf[pos++] = comm->load[d+1].sum;
                sbuf[pos++] = comm->load[d+1].max;
                if (isDlbOn(dd->comm))
                {
                    sbuf[pos++] = comm->load[d+1].sum_m;
                    sbuf[pos++] = comm->load[d+1].cvol_min*cell_frac;
                    sbuf[pos++] = comm->load[d+1].flags;
                    if (d > 0)
                    {
                        sbuf[pos++] = cellsizes.fracLowerMax;
                        sbuf[pos++] = cellsizes.fracUpperMin;
                    }
                }
                if (bSepPME)
                {
                    sbuf[pos++] = comm->load[d+1].mdf;
                    sbuf[pos++] = comm->load[d+1].pme;
                }
            }
            load->nload = pos;
            /* Communicate a row in DD direction d.
             * The communicators are setup such that the root always has rank 0.
             */
#if GMX_MPI
            MPI_Gather(sbuf, load->nload*sizeof(float), MPI_BYTE,
                       load->load, load->nload*sizeof(float), MPI_BYTE,
                       0, comm->mpi_comm_load[d]);
#endif
            if (dd->ci[dim] == dd->master_ci[dim])
            {
                /* We are the master along this row, process this row */
                RowMaster *rowMaster = nullptr;

                if (isDlbOn(comm))
                {
                    rowMaster = cellsizes.rowMaster.get();
                }
                load->sum      = 0;
                load->max      = 0;
                load->sum_m    = 0;
                load->cvol_min = 1;
                load->flags    = 0;
                load->mdf      = 0;
                load->pme      = 0;
                int pos        = 0;
                for (int i = 0; i < dd->nc[dim]; i++)
                {
                    load->sum += load->load[pos++];
                    load->max  = std::max(load->max, load->load[pos]);
                    pos++;
                    if (isDlbOn(dd->comm))
                    {
                        if (rowMaster->dlbIsLimited)
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
                        if (d < dd->ndim-1)
                        {
                            load->flags = (int)(load->load[pos++] + 0.5);
                        }
                        if (d > 0)
                        {
                            rowMaster->bounds[i].cellFracLowerMax = load->load[pos++];
                            rowMaster->bounds[i].cellFracUpperMin = load->load[pos++];
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
                if (isDlbOn(comm) && rowMaster->dlbIsLimited)
                {
                    load->sum_m *= dd->nc[dim];
                    load->flags |= (1<<d);
                }
            }
        }
    }

    if (DDMASTER(dd))
    {
        comm->nload      += dd_load_count(comm);
        comm->load_step  += comm->cycl[ddCyclStep];
        comm->load_sum   += comm->load[0].sum;
        comm->load_max   += comm->load[0].max;
        if (isDlbOn(comm))
        {
            for (int d = 0; d < dd->ndim; d++)
            {
                if (comm->load[0].flags & (1<<d))
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

    wallcycle_stop(wcycle, ewcDDCOMMLOAD);

    if (debug)
    {
        fprintf(debug, "get_load_distribution finished\n");
    }
}

static float dd_force_load_fraction(gmx_domdec_t *dd)
{
    /* Return the relative performance loss on the total run time
     * due to the force calculation load imbalance.
     */
    if (dd->comm->nload > 0 && dd->comm->load_step > 0)
    {
        return dd->comm->load_sum/(dd->comm->load_step*dd->nnodes);
    }
    else
    {
        return 0;
    }
}

static float dd_force_imb_perf_loss(gmx_domdec_t *dd)
{
    /* Return the relative performance loss on the total run time
     * due to the force calculation load imbalance.
     */
    if (dd->comm->nload > 0 && dd->comm->load_step > 0)
    {
        return
            (dd->comm->load_max*dd->nnodes - dd->comm->load_sum)/
            (dd->comm->load_step*dd->nnodes);
    }
    else
    {
        return 0;
    }
}

static void print_dd_load_av(FILE *fplog, gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm = dd->comm;

    /* Only the master rank prints loads and only if we measured loads */
    if (!DDMASTER(dd) || comm->nload == 0)
    {
        return;
    }

    char  buf[STRLEN];
    int   numPpRanks   = dd->nnodes;
    int   numPmeRanks  = (dd->pme_nodeid >= 0) ? comm->npmenodes : 0;
    int   numRanks     = numPpRanks + numPmeRanks;
    float lossFraction = 0;

    /* Print the average load imbalance and performance loss */
    if (dd->nnodes > 1 && comm->load_sum > 0)
    {
        float imbalance = comm->load_max*numPpRanks/comm->load_sum - 1;
        lossFraction    = dd_force_imb_perf_loss(dd);

        std::string msg         = "\n Dynamic load balancing report:\n";
        std::string dlbStateStr;

        switch (dd->comm->dlbState)
        {
            case edlbsOffUser:
                dlbStateStr = "DLB was off during the run per user request.";
                break;
            case edlbsOffForever:
                /* Currectly this can happen due to performance loss observed, cell size
                 * limitations or incompatibility with other settings observed during
                 * determineInitialDlbState(). */
                dlbStateStr = "DLB got disabled because it was unsuitable to use.";
                break;
            case edlbsOffCanTurnOn:
                dlbStateStr = "DLB was off during the run due to low measured imbalance.";
                break;
            case edlbsOffTemporarilyLocked:
                dlbStateStr = "DLB was locked at the end of the run due to unfinished PP-PME balancing.";
                break;
            case edlbsOnCanTurnOff:
                dlbStateStr = "DLB was turned on during the run due to measured imbalance.";
                break;
            case edlbsOnUser:
                dlbStateStr = "DLB was permanently on during the run per user request.";
                break;
            default:
                GMX_ASSERT(false, "Undocumented DLB state");
        }

        msg += " " + dlbStateStr + "\n";
        msg += gmx::formatString(" Average load imbalance: %.1f%%.\n", imbalance*100);
        msg += gmx::formatString(" The balanceable part of the MD step is %d%%, load imbalance is computed from this.\n",
                                 static_cast<int>(dd_force_load_fraction(dd)*100 + 0.5));
        msg += gmx::formatString(" Part of the total run time spent waiting due to load imbalance: %.1f%%.\n",
                                 lossFraction*100);
        fprintf(fplog, "%s", msg.c_str());
        fprintf(stderr, "%s", msg.c_str());
    }

    /* Print during what percentage of steps the  load balancing was limited */
    bool dlbWasLimited = false;
    if (isDlbOn(comm))
    {
        sprintf(buf, " Steps where the load balancing was limited by -rdd, -rcon and/or -dds:");
        for (int d = 0; d < dd->ndim; d++)
        {
            int limitPercentage = (200*comm->load_lim[d] + 1)/(2*comm->nload);
            sprintf(buf+strlen(buf), " %c %d %%",
                    dim2char(dd->dim[d]), limitPercentage);
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
        float pmeForceRatio = comm->load_pme/comm->load_mdf;
        lossFractionPme     = (comm->load_pme - comm->load_mdf)/comm->load_step;
        if (lossFractionPme <= 0)
        {
            lossFractionPme *= numPmeRanks/static_cast<float>(numRanks);
        }
        else
        {
            lossFractionPme *= numPpRanks/static_cast<float>(numRanks);
        }
        sprintf(buf, " Average PME mesh/force load: %5.3f\n", pmeForceRatio);
        fprintf(fplog, "%s", buf);
        fprintf(stderr, "%s", buf);
        sprintf(buf, " Part of the total run time spent waiting due to PP/PME imbalance: %.1f %%\n", std::fabs(lossFractionPme)*100);
        fprintf(fplog, "%s", buf);
        fprintf(stderr, "%s", buf);
    }
    fprintf(fplog, "\n");
    fprintf(stderr, "\n");

    if (lossFraction >= DD_PERF_LOSS_WARN)
    {
        sprintf(buf,
                "NOTE: %.1f %% of the available CPU time was lost due to load imbalance\n"
                "      in the domain decomposition.\n", lossFraction*100);
        if (!isDlbOn(comm))
        {
            sprintf(buf+strlen(buf), "      You might want to use dynamic load balancing (option -dlb.)\n");
        }
        else if (dlbWasLimited)
        {
            sprintf(buf+strlen(buf), "      You might want to decrease the cell size limit (options -rdd, -rcon and/or -dds).\n");
        }
        fprintf(fplog, "%s\n", buf);
        fprintf(stderr, "%s\n", buf);
    }
    if (numPmeRanks > 0 && std::fabs(lossFractionPme) >= DD_PERF_LOSS_WARN)
    {
        sprintf(buf,
                "NOTE: %.1f %% performance was lost because the PME ranks\n"
                "      had %s work to do than the PP ranks.\n"
                "      You might want to %s the number of PME ranks\n"
                "      or %s the cut-off and the grid spacing.\n",
                std::fabs(lossFractionPme*100),
                (lossFractionPme < 0) ? "less"     : "more",
                (lossFractionPme < 0) ? "decrease" : "increase",
                (lossFractionPme < 0) ? "decrease" : "increase");
        fprintf(fplog, "%s\n", buf);
        fprintf(stderr, "%s\n", buf);
    }
}

static float dd_vol_min(gmx_domdec_t *dd)
{
    return dd->comm->load[0].cvol_min*dd->nnodes;
}

static gmx_bool dd_load_flags(gmx_domdec_t *dd)
{
    return dd->comm->load[0].flags;
}

static float dd_f_imbal(gmx_domdec_t *dd)
{
    if (dd->comm->load[0].sum > 0)
    {
        return dd->comm->load[0].max*dd->nnodes/dd->comm->load[0].sum - 1.0f;
    }
    else
    {
        /* Something is wrong in the cycle counting, report no load imbalance */
        return 0.0f;
    }
}

float dd_pme_f_ratio(gmx_domdec_t *dd)
{
    /* Should only be called on the DD master rank */
    assert(DDMASTER(dd));

    if (dd->comm->load[0].mdf > 0 && dd->comm->cycl_n[ddCyclPME] > 0)
    {
        return dd->comm->load[0].pme/dd->comm->load[0].mdf;
    }
    else
    {
        return -1.0;
    }
}

static void dd_print_load(FILE *fplog, gmx_domdec_t *dd, int64_t step)
{
    int  flags, d;
    char buf[22];

    flags = dd_load_flags(dd);
    if (flags)
    {
        fprintf(fplog,
                "DD  load balancing is limited by minimum cell size in dimension");
        for (d = 0; d < dd->ndim; d++)
        {
            if (flags & (1<<d))
            {
                fprintf(fplog, " %c", dim2char(dd->dim[d]));
            }
        }
        fprintf(fplog, "\n");
    }
    fprintf(fplog, "DD  step %s", gmx_step_str(step, buf));
    if (isDlbOn(dd->comm))
    {
        fprintf(fplog, "  vol min/aver %5.3f%c",
                dd_vol_min(dd), flags ? '!' : ' ');
    }
    if (dd->nnodes > 1)
    {
        fprintf(fplog, " load imb.: force %4.1f%%", dd_f_imbal(dd)*100);
    }
    if (dd->comm->cycl_n[ddCyclPME])
    {
        fprintf(fplog, "  pme mesh/force %5.3f", dd_pme_f_ratio(dd));
    }
    fprintf(fplog, "\n\n");
}

static void dd_print_load_verbose(gmx_domdec_t *dd)
{
    if (isDlbOn(dd->comm))
    {
        fprintf(stderr, "vol %4.2f%c ",
                dd_vol_min(dd), dd_load_flags(dd) ? '!' : ' ');
    }
    if (dd->nnodes > 1)
    {
        fprintf(stderr, "imb F %2d%% ", (int)(dd_f_imbal(dd)*100+0.5));
    }
    if (dd->comm->cycl_n[ddCyclPME])
    {
        fprintf(stderr, "pme/F %4.2f ", dd_pme_f_ratio(dd));
    }
}

#if GMX_MPI
static void make_load_communicator(gmx_domdec_t *dd, int dim_ind, ivec loc)
{
    MPI_Comm           c_row;
    int                dim, i, rank;
    ivec               loc_c;
    gmx_bool           bPartOfGroup = FALSE;

    dim = dd->dim[dim_ind];
    copy_ivec(loc, loc_c);
    for (i = 0; i < dd->nc[dim]; i++)
    {
        loc_c[dim] = i;
        rank       = dd_index(dd->nc, loc_c);
        if (rank == dd->rank)
        {
            /* This process is part of the group */
            bPartOfGroup = TRUE;
        }
    }
    MPI_Comm_split(dd->mpi_comm_all, bPartOfGroup ? 0 : MPI_UNDEFINED, dd->rank,
                   &c_row);
    if (bPartOfGroup)
    {
        dd->comm->mpi_comm_load[dim_ind] = c_row;
        if (!isDlbDisabled(dd->comm))
        {
            DDCellsizesWithDlb &cellsizes = dd->comm->cellsizesWithDlb[dim_ind];

            if (dd->ci[dim] == dd->master_ci[dim])
            {
                /* This is the root process of this row */
                cellsizes.rowMaster  = gmx::compat::make_unique<RowMaster>();

                RowMaster &rowMaster = *cellsizes.rowMaster;
                rowMaster.cellFrac.resize(ddCellFractionBufferSize(dd, dim_ind));
                rowMaster.oldCellFrac.resize(dd->nc[dim] + 1);
                rowMaster.isCellMin.resize(dd->nc[dim]);
                if (dim_ind > 0)
                {
                    rowMaster.bounds.resize(dd->nc[dim]);
                }
                rowMaster.buf_ncd.resize(dd->nc[dim]);
            }
            else
            {
                /* This is not a root process, we only need to receive cell_f */
                cellsizes.fracRow.resize(ddCellFractionBufferSize(dd, dim_ind));
            }
        }
        if (dd->ci[dim] == dd->master_ci[dim])
        {
            snew(dd->comm->load[dim_ind].load, dd->nc[dim]*DD_NLOAD_MAX);
        }
    }
}
#endif

void dd_setup_dlb_resource_sharing(t_commrec            *cr,
                                   int                   gpu_id)
{
#if GMX_MPI
    int           physicalnode_id_hash;
    gmx_domdec_t *dd;
    MPI_Comm      mpi_comm_pp_physicalnode;

    if (!thisRankHasDuty(cr, DUTY_PP) || gpu_id < 0)
    {
        /* Only ranks with short-ranged tasks (currently) use GPUs.
         * If we don't have GPUs assigned, there are no resources to share.
         */
        return;
    }

    physicalnode_id_hash = gmx_physicalnode_id_hash();

    dd = cr->dd;

    if (debug)
    {
        fprintf(debug, "dd_setup_dd_dlb_gpu_sharing:\n");
        fprintf(debug, "DD PP rank %d physical node hash %d gpu_id %d\n",
                dd->rank, physicalnode_id_hash, gpu_id);
    }
    /* Split the PP communicator over the physical nodes */
    /* TODO: See if we should store this (before), as it's also used for
     * for the nodecomm summation.
     */
    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.
    MPI_Comm_split(dd->mpi_comm_all, physicalnode_id_hash, dd->rank,
                   &mpi_comm_pp_physicalnode);
    MPI_Comm_split(mpi_comm_pp_physicalnode, gpu_id, dd->rank,
                   &dd->comm->mpi_comm_gpu_shared);
    MPI_Comm_free(&mpi_comm_pp_physicalnode);
    MPI_Comm_size(dd->comm->mpi_comm_gpu_shared, &dd->comm->nrank_gpu_shared);

    if (debug)
    {
        fprintf(debug, "nrank_gpu_shared %d\n", dd->comm->nrank_gpu_shared);
    }

    /* Note that some ranks could share a GPU, while others don't */

    if (dd->comm->nrank_gpu_shared == 1)
    {
        MPI_Comm_free(&dd->comm->mpi_comm_gpu_shared);
    }
#else
    GMX_UNUSED_VALUE(cr);
    GMX_UNUSED_VALUE(gpu_id);
#endif
}

static void make_load_communicators(gmx_domdec_t gmx_unused *dd)
{
#if GMX_MPI
    int  dim0, dim1, i, j;
    ivec loc;

    if (debug)
    {
        fprintf(debug, "Making load communicators\n");
    }

    snew(dd->comm->load,          std::max(dd->ndim, 1));
    snew(dd->comm->mpi_comm_load, std::max(dd->ndim, 1));

    if (dd->ndim == 0)
    {
        return;
    }

    clear_ivec(loc);
    make_load_communicator(dd, 0, loc);
    if (dd->ndim > 1)
    {
        dim0 = dd->dim[0];
        for (i = 0; i < dd->nc[dim0]; i++)
        {
            loc[dim0] = i;
            make_load_communicator(dd, 1, loc);
        }
    }
    if (dd->ndim > 2)
    {
        dim0 = dd->dim[0];
        for (i = 0; i < dd->nc[dim0]; i++)
        {
            loc[dim0] = i;
            dim1      = dd->dim[1];
            for (j = 0; j < dd->nc[dim1]; j++)
            {
                loc[dim1] = j;
                make_load_communicator(dd, 2, loc);
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Finished making load communicators\n");
    }
#endif
}

/*! \brief Sets up the relation between neighboring domains and zones */
static void setup_neighbor_relations(gmx_domdec_t *dd)
{
    int                     d, dim, i, j, m;
    ivec                    tmp, s;
    gmx_domdec_zones_t     *zones;
    gmx_domdec_ns_ranges_t *izone;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        copy_ivec(dd->ci, tmp);
        tmp[dim]           = (tmp[dim] + 1) % dd->nc[dim];
        dd->neighbor[d][0] = ddcoord2ddnodeid(dd, tmp);
        copy_ivec(dd->ci, tmp);
        tmp[dim]           = (tmp[dim] - 1 + dd->nc[dim]) % dd->nc[dim];
        dd->neighbor[d][1] = ddcoord2ddnodeid(dd, tmp);
        if (debug)
        {
            fprintf(debug, "DD rank %d neighbor ranks in dir %d are + %d - %d\n",
                    dd->rank, dim,
                    dd->neighbor[d][0],
                    dd->neighbor[d][1]);
        }
    }

    int nzone  = (1 << dd->ndim);
    int nizone = (1 << std::max(dd->ndim - 1, 0));
    assert(nizone >= 1 && nizone <= DD_MAXIZONE);

    zones = &dd->comm->zones;

    for (i = 0; i < nzone; i++)
    {
        m = 0;
        clear_ivec(zones->shift[i]);
        for (d = 0; d < dd->ndim; d++)
        {
            zones->shift[i][dd->dim[d]] = dd_zo[i][m++];
        }
    }

    zones->n = nzone;
    for (i = 0; i < nzone; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            s[d] = dd->ci[d] - zones->shift[i][d];
            if (s[d] < 0)
            {
                s[d] += dd->nc[d];
            }
            else if (s[d] >= dd->nc[d])
            {
                s[d] -= dd->nc[d];
            }
        }
    }
    zones->nizone = nizone;
    for (i = 0; i < zones->nizone; i++)
    {
        assert(ddNonbondedZonePairRanges[i][0] == i);

        izone     = &zones->izone[i];
        /* dd_zp3 is for 3D decomposition, for fewer dimensions use only
         * j-zones up to nzone.
         */
        izone->j0 = std::min(ddNonbondedZonePairRanges[i][1], nzone);
        izone->j1 = std::min(ddNonbondedZonePairRanges[i][2], nzone);
        for (dim = 0; dim < DIM; dim++)
        {
            if (dd->nc[dim] == 1)
            {
                /* All shifts should be allowed */
                izone->shift0[dim] = -1;
                izone->shift1[dim] = 1;
            }
            else
            {
                /* Determine the min/max j-zone shift wrt the i-zone */
                izone->shift0[dim] = 1;
                izone->shift1[dim] = -1;
                for (j = izone->j0; j < izone->j1; j++)
                {
                    int shift_diff = zones->shift[j][dim] - zones->shift[i][dim];
                    if (shift_diff < izone->shift0[dim])
                    {
                        izone->shift0[dim] = shift_diff;
                    }
                    if (shift_diff > izone->shift1[dim])
                    {
                        izone->shift1[dim] = shift_diff;
                    }
                }
            }
        }
    }

    if (!isDlbDisabled(dd->comm))
    {
        dd->comm->cellsizesWithDlb.resize(dd->ndim);
    }

    if (dd->comm->bRecordLoad)
    {
        make_load_communicators(dd);
    }
}

static void make_pp_communicator(FILE                 *fplog,
                                 gmx_domdec_t         *dd,
                                 t_commrec gmx_unused *cr,
                                 int gmx_unused        reorder)
{
#if GMX_MPI
    gmx_domdec_comm_t *comm;
    int                rank, *buf;
    ivec               periods;
    MPI_Comm           comm_cart;

    comm = dd->comm;

    if (comm->bCartesianPP)
    {
        /* Set up cartesian communication for the particle-particle part */
        if (fplog)
        {
            fprintf(fplog, "Will use a Cartesian communicator: %d x %d x %d\n",
                    dd->nc[XX], dd->nc[YY], dd->nc[ZZ]);
        }

        for (int i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Cart_create(cr->mpi_comm_mygroup, DIM, dd->nc, periods, reorder,
                        &comm_cart);
        /* We overwrite the old communicator with the new cartesian one */
        cr->mpi_comm_mygroup = comm_cart;
    }

    dd->mpi_comm_all = cr->mpi_comm_mygroup;
    MPI_Comm_rank(dd->mpi_comm_all, &dd->rank);

    if (comm->bCartesianPP_PME)
    {
        /* Since we want to use the original cartesian setup for sim,
         * and not the one after split, we need to make an index.
         */
        snew(comm->ddindex2ddnodeid, dd->nnodes);
        comm->ddindex2ddnodeid[dd_index(dd->nc, dd->ci)] = dd->rank;
        gmx_sumi(dd->nnodes, comm->ddindex2ddnodeid, cr);
        /* Get the rank of the DD master,
         * above we made sure that the master node is a PP node.
         */
        if (MASTER(cr))
        {
            rank = dd->rank;
        }
        else
        {
            rank = 0;
        }
        MPI_Allreduce(&rank, &dd->masterrank, 1, MPI_INT, MPI_SUM, dd->mpi_comm_all);
    }
    else if (comm->bCartesianPP)
    {
        if (cr->npmenodes == 0)
        {
            /* The PP communicator is also
             * the communicator for this simulation
             */
            cr->mpi_comm_mysim = cr->mpi_comm_mygroup;
        }
        cr->nodeid = dd->rank;

        MPI_Cart_coords(dd->mpi_comm_all, dd->rank, DIM, dd->ci);

        /* We need to make an index to go from the coordinates
         * to the nodeid of this simulation.
         */
        snew(comm->ddindex2simnodeid, dd->nnodes);
        snew(buf, dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->nc, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf, comm->ddindex2simnodeid, dd->nnodes, MPI_INT, MPI_SUM,
                      cr->mpi_comm_mysim);
        sfree(buf);

        /* Determine the master coordinates and rank.
         * The DD master should be the same node as the master of this sim.
         */
        for (int i = 0; i < dd->nnodes; i++)
        {
            if (comm->ddindex2simnodeid[i] == 0)
            {
                ddindex2xyz(dd->nc, i, dd->master_ci);
                MPI_Cart_rank(dd->mpi_comm_all, dd->master_ci, &dd->masterrank);
            }
        }
        if (debug)
        {
            fprintf(debug, "The master rank is %d\n", dd->masterrank);
        }
    }
    else
    {
        /* No Cartesian communicators */
        /* We use the rank in dd->comm->all as DD index */
        ddindex2xyz(dd->nc, dd->rank, dd->ci);
        /* The simulation master nodeid is 0, so the DD master rank is also 0 */
        dd->masterrank = 0;
        clear_ivec(dd->master_ci);
    }
#endif

    if (fplog)
    {
        fprintf(fplog,
                "Domain decomposition rank %d, coordinates %d %d %d\n\n",
                dd->rank, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
    }
    if (debug)
    {
        fprintf(debug,
                "Domain decomposition rank %d, coordinates %d %d %d\n\n",
                dd->rank, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
    }
}

static void receive_ddindex2simnodeid(gmx_domdec_t         *dd,
                                      t_commrec            *cr)
{
#if GMX_MPI
    gmx_domdec_comm_t *comm = dd->comm;

    if (!comm->bCartesianPP_PME && comm->bCartesianPP)
    {
        int *buf;
        snew(comm->ddindex2simnodeid, dd->nnodes);
        snew(buf, dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->nc, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf, comm->ddindex2simnodeid, dd->nnodes, MPI_INT, MPI_SUM,
                      cr->mpi_comm_mysim);
        sfree(buf);
    }
#else
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(cr);
#endif
}

static void split_communicator(FILE *fplog, t_commrec *cr, gmx_domdec_t *dd,
                               DdRankOrder gmx_unused rankOrder,
                               int gmx_unused reorder)
{
    gmx_domdec_comm_t *comm;
    int                i;
    gmx_bool           bDiv[DIM];
#if GMX_MPI
    MPI_Comm           comm_cart;
#endif

    comm = dd->comm;

    if (comm->bCartesianPP)
    {
        for (i = 1; i < DIM; i++)
        {
            bDiv[i] = ((cr->npmenodes*dd->nc[i]) % (dd->nnodes) == 0);
        }
        if (bDiv[YY] || bDiv[ZZ])
        {
            comm->bCartesianPP_PME = TRUE;
            /* If we have 2D PME decomposition, which is always in x+y,
             * we stack the PME only nodes in z.
             * Otherwise we choose the direction that provides the thinnest slab
             * of PME only nodes as this will have the least effect
             * on the PP communication.
             * But for the PME communication the opposite might be better.
             */
            if (bDiv[ZZ] && (comm->npmenodes_y > 1 ||
                             !bDiv[YY] ||
                             dd->nc[YY] > dd->nc[ZZ]))
            {
                comm->cartpmedim = ZZ;
            }
            else
            {
                comm->cartpmedim = YY;
            }
            comm->ntot[comm->cartpmedim]
                += (cr->npmenodes*dd->nc[comm->cartpmedim])/dd->nnodes;
        }
        else if (fplog)
        {
            fprintf(fplog, "Number of PME-only ranks (%d) is not a multiple of nx*ny (%d*%d) or nx*nz (%d*%d)\n", cr->npmenodes, dd->nc[XX], dd->nc[YY], dd->nc[XX], dd->nc[ZZ]);
            fprintf(fplog,
                    "Will not use a Cartesian communicator for PP <-> PME\n\n");
        }
    }

    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        int  rank;
        ivec periods;

        if (fplog)
        {
            fprintf(fplog, "Will use a Cartesian communicator for PP <-> PME: %d x %d x %d\n", comm->ntot[XX], comm->ntot[YY], comm->ntot[ZZ]);
        }

        for (i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Cart_create(cr->mpi_comm_mysim, DIM, comm->ntot, periods, reorder,
                        &comm_cart);
        MPI_Comm_rank(comm_cart, &rank);
        if (MASTER(cr) && rank != 0)
        {
            gmx_fatal(FARGS, "MPI rank 0 was renumbered by MPI_Cart_create, we do not allow this");
        }

        /* With this assigment we loose the link to the original communicator
         * which will usually be MPI_COMM_WORLD, unless have multisim.
         */
        cr->mpi_comm_mysim = comm_cart;
        cr->sim_nodeid     = rank;

        MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, dd->ci);

        if (fplog)
        {
            fprintf(fplog, "Cartesian rank %d, coordinates %d %d %d\n\n",
                    cr->sim_nodeid, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
        }

        if (dd->ci[comm->cartpmedim] < dd->nc[comm->cartpmedim])
        {
            cr->duty = DUTY_PP;
        }
        if (cr->npmenodes == 0 ||
            dd->ci[comm->cartpmedim] >= dd->nc[comm->cartpmedim])
        {
            cr->duty = DUTY_PME;
        }

        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       getThisRankDuties(cr),
                       dd_index(comm->ntot, dd->ci),
                       &cr->mpi_comm_mygroup);
#endif
    }
    else
    {
        switch (rankOrder)
        {
            case DdRankOrder::pp_pme:
                if (fplog)
                {
                    fprintf(fplog, "Order of the ranks: PP first, PME last\n");
                }
                break;
            case DdRankOrder::interleave:
                /* Interleave the PP-only and PME-only ranks */
                if (fplog)
                {
                    fprintf(fplog, "Interleaving PP and PME ranks\n");
                }
                comm->pmenodes = dd_interleaved_pme_ranks(dd);
                break;
            case DdRankOrder::cartesian:
                break;
            default:
                gmx_fatal(FARGS, "Invalid ddRankOrder=%d", static_cast<int>(rankOrder));
        }

        if (dd_simnode2pmenode(dd, cr, cr->sim_nodeid) == -1)
        {
            cr->duty = DUTY_PME;
        }
        else
        {
            cr->duty = DUTY_PP;
        }
#if GMX_MPI
        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       getThisRankDuties(cr),
                       cr->nodeid,
                       &cr->mpi_comm_mygroup);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
#endif
    }

    if (fplog)
    {
        fprintf(fplog, "This rank does only %s work.\n\n",
                thisRankHasDuty(cr, DUTY_PP) ? "particle-particle" : "PME-mesh");
    }
}

/*! \brief Generates the MPI communicators for domain decomposition */
static void make_dd_communicators(FILE *fplog, t_commrec *cr,
                                  gmx_domdec_t *dd, DdRankOrder ddRankOrder)
{
    gmx_domdec_comm_t *comm;
    int                CartReorder;

    comm = dd->comm;

    copy_ivec(dd->nc, comm->ntot);

    comm->bCartesianPP     = (ddRankOrder == DdRankOrder::cartesian);
    comm->bCartesianPP_PME = FALSE;

    /* Reorder the nodes by default. This might change the MPI ranks.
     * Real reordering is only supported on very few architectures,
     * Blue Gene is one of them.
     */
    CartReorder = (getenv("GMX_NO_CART_REORDER") == nullptr);

    if (cr->npmenodes > 0)
    {
        /* Split the communicator into a PP and PME part */
        split_communicator(fplog, cr, dd, ddRankOrder, CartReorder);
        if (comm->bCartesianPP_PME)
        {
            /* We (possibly) reordered the nodes in split_communicator,
             * so it is no longer required in make_pp_communicator.
             */
            CartReorder = FALSE;
        }
    }
    else
    {
        /* All nodes do PP and PME */
#if GMX_MPI
        /* We do not require separate communicators */
        cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
#endif
    }

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Copy or make a new PP communicator */
        make_pp_communicator(fplog, dd, cr, CartReorder);
    }
    else
    {
        receive_ddindex2simnodeid(dd, cr);
    }

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Set up the commnuication to our PME node */
        dd->pme_nodeid           = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
        dd->pme_receive_vir_ener = receive_vir_ener(dd, cr);
        if (debug)
        {
            fprintf(debug, "My pme_nodeid %d receive ener %d\n",
                    dd->pme_nodeid, dd->pme_receive_vir_ener);
        }
    }
    else
    {
        dd->pme_nodeid = -1;
    }

    if (DDMASTER(dd))
    {
        dd->ma = gmx::compat::make_unique<AtomDistribution>(dd->nc,
                                                            comm->cgs_gl.nr,
                                                            comm->cgs_gl.index[comm->cgs_gl.nr]);
    }
}

static real *get_slb_frac(FILE *fplog, const char *dir, int nc, const char *size_string)
{
    real  *slb_frac, tot;
    int    i, n;
    double dbl;

    slb_frac = nullptr;
    if (nc > 1 && size_string != nullptr)
    {
        if (fplog)
        {
            fprintf(fplog, "Using static load balancing for the %s direction\n",
                    dir);
        }
        snew(slb_frac, nc);
        tot = 0;
        for (i = 0; i < nc; i++)
        {
            dbl = 0;
            sscanf(size_string, "%20lf%n", &dbl, &n);
            if (dbl == 0)
            {
                gmx_fatal(FARGS, "Incorrect or not enough DD cell size entries for direction %s: '%s'", dir, size_string);
            }
            slb_frac[i]  = dbl;
            size_string += n;
            tot         += slb_frac[i];
        }
        /* Normalize */
        if (fplog)
        {
            fprintf(fplog, "Relative cell sizes:");
        }
        for (i = 0; i < nc; i++)
        {
            slb_frac[i] /= tot;
            if (fplog)
            {
                fprintf(fplog, " %5.3f", slb_frac[i]);
            }
        }
        if (fplog)
        {
            fprintf(fplog, "\n");
        }
    }

    return slb_frac;
}

static int multi_body_bondeds_count(const gmx_mtop_t *mtop)
{
    int                  n, nmol, ftype;
    gmx_mtop_ilistloop_t iloop;
    const t_ilist       *il;

    n     = 0;
    iloop = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & IF_BOND) &&
                NRAL(ftype) >  2)
            {
                n += nmol*il[ftype].nr/(1 + NRAL(ftype));
            }
        }
    }

    return n;
}

static int dd_getenv(FILE *fplog, const char *env_var, int def)
{
    char *val;
    int   nst;

    nst = def;
    val = getenv(env_var);
    if (val)
    {
        if (sscanf(val, "%20d", &nst) <= 0)
        {
            nst = 1;
        }
        if (fplog)
        {
            fprintf(fplog, "Found env.var. %s = %s, using value %d\n",
                    env_var, val, nst);
        }
    }

    return nst;
}

static void dd_warning(const t_commrec *cr, FILE *fplog, const char *warn_string)
{
    if (MASTER(cr))
    {
        fprintf(stderr, "\n%s\n", warn_string);
    }
    if (fplog)
    {
        fprintf(fplog, "\n%s\n", warn_string);
    }
}

static void check_dd_restrictions(t_commrec *cr, const gmx_domdec_t *dd,
                                  const t_inputrec *ir, FILE *fplog)
{
    if (ir->ePBC == epbcSCREW &&
        (dd->nc[XX] == 1 || dd->nc[YY] > 1 || dd->nc[ZZ] > 1))
    {
        gmx_fatal(FARGS, "With pbc=%s can only do domain decomposition in the x-direction", epbc_names[ir->ePBC]);
    }

    if (ir->ns_type == ensSIMPLE)
    {
        gmx_fatal(FARGS, "Domain decomposition does not support simple neighbor searching, use grid searching or run with one MPI rank");
    }

    if (ir->nstlist == 0)
    {
        gmx_fatal(FARGS, "Domain decomposition does not work with nstlist=0");
    }

    if (ir->comm_mode == ecmANGULAR && ir->ePBC != epbcNONE)
    {
        dd_warning(cr, fplog, "comm-mode angular will give incorrect results when the comm group partially crosses a periodic boundary");
    }
}

static real average_cellsize_min(gmx_domdec_t *dd, gmx_ddbox_t *ddbox)
{
    int  di, d;
    real r;

    r = ddbox->box_size[XX];
    for (di = 0; di < dd->ndim; di++)
    {
        d = dd->dim[di];
        /* Check using the initial average cell size */
        r = std::min(r, ddbox->box_size[d]*ddbox->skew_fac[d]/dd->nc[d]);
    }

    return r;
}

/*! \brief Depending on the DLB initial value return the DLB switched off state or issue an error.
 */
static int forceDlbOffOrBail(int                cmdlineDlbState,
                             const std::string &reasonStr,
                             t_commrec         *cr,
                             FILE              *fplog)
{
    std::string dlbNotSupportedErr  = "Dynamic load balancing requested, but ";
    std::string dlbDisableNote      = "NOTE: disabling dynamic load balancing as ";

    if (cmdlineDlbState == edlbsOnUser)
    {
        gmx_fatal(FARGS, "%s", (dlbNotSupportedErr + reasonStr).c_str());
    }
    else if (cmdlineDlbState == edlbsOffCanTurnOn)
    {
        dd_warning(cr, fplog, (dlbDisableNote + reasonStr + "\n").c_str());
    }
    return edlbsOffForever;
}

/*! \brief Return the dynamic load balancer's initial state based on initial conditions and user inputs.
 *
 * This function parses the parameters of "-dlb" command line option setting
 * corresponding state values. Then it checks the consistency of the determined
 * state with other run parameters and settings. As a result, the initial state
 * may be altered or an error may be thrown if incompatibility of options is detected.
 *
 * \param [in] fplog       Pointer to mdrun log file.
 * \param [in] cr          Pointer to MPI communication object.
 * \param [in] dlbOption   Enum value for the DLB option.
 * \param [in] bRecordLoad True if the load balancer is recording load information.
 * \param [in] mdrunOptions  Options for mdrun.
 * \param [in] ir          Pointer mdrun to input parameters.
 * \returns                DLB initial/startup state.
 */
static int determineInitialDlbState(FILE *fplog, t_commrec *cr,
                                    DlbOption dlbOption, gmx_bool bRecordLoad,
                                    const MdrunOptions &mdrunOptions,
                                    const t_inputrec *ir)
{
    int dlbState = edlbsOffCanTurnOn;

    switch (dlbOption)
    {
        case DlbOption::turnOnWhenUseful: dlbState = edlbsOffCanTurnOn; break;
        case DlbOption::no:               dlbState = edlbsOffUser;      break;
        case DlbOption::yes:              dlbState = edlbsOnUser;       break;
        default: gmx_incons("Invalid dlbOption enum value");
    }

    /* Reruns don't support DLB: bail or override auto mode */
    if (mdrunOptions.rerun)
    {
        std::string reasonStr = "it is not supported in reruns.";
        return forceDlbOffOrBail(dlbState, reasonStr, cr, fplog);
    }

    /* Unsupported integrators */
    if (!EI_DYNAMICS(ir->eI))
    {
        auto reasonStr = gmx::formatString("it is only supported with dynamics, not with integrator '%s'.", EI(ir->eI));
        return forceDlbOffOrBail(dlbState, reasonStr, cr, fplog);
    }

    /* Without cycle counters we can't time work to balance on */
    if (!bRecordLoad)
    {
        std::string reasonStr = "cycle counters unsupported or not enabled in the operating system kernel.";
        return forceDlbOffOrBail(dlbState, reasonStr, cr, fplog);
    }

    if (mdrunOptions.reproducible)
    {
        std::string reasonStr = "you started a reproducible run.";
        switch (dlbState)
        {
            case edlbsOffUser:
                break;
            case edlbsOffForever:
                GMX_RELEASE_ASSERT(false, "edlbsOffForever is not a valid initial state");
                break;
            case edlbsOffCanTurnOn:
                return forceDlbOffOrBail(dlbState, reasonStr, cr, fplog);
            case edlbsOnCanTurnOff:
                GMX_RELEASE_ASSERT(false, "edlbsOffCanTurnOff is not a valid initial state");
                break;
            case edlbsOnUser:
                return forceDlbOffOrBail(dlbState, reasonStr + " In load balanced runs binary reproducibility cannot be ensured.", cr, fplog);
            default:
                gmx_fatal(FARGS, "Death horror: undefined case (%d) for load balancing choice", dlbState);
        }
    }

    return dlbState;
}

static void set_dd_dim(FILE *fplog, gmx_domdec_t *dd)
{
    dd->ndim = 0;
    if (getenv("GMX_DD_ORDER_ZYX") != nullptr)
    {
        /* Decomposition order z,y,x */
        if (fplog)
        {
            fprintf(fplog, "Using domain decomposition order z, y, x\n");
        }
        for (int dim = DIM-1; dim >= 0; dim--)
        {
            if (dd->nc[dim] > 1)
            {
                dd->dim[dd->ndim++] = dim;
            }
        }
    }
    else
    {
        /* Decomposition order x,y,z */
        for (int dim = 0; dim < DIM; dim++)
        {
            if (dd->nc[dim] > 1)
            {
                dd->dim[dd->ndim++] = dim;
            }
        }
    }

    if (dd->ndim == 0)
    {
        /* Set dim[0] to avoid extra checks on ndim in several places */
        dd->dim[0] = XX;
    }
}

static gmx_domdec_comm_t *init_dd_comm()
{
    gmx_domdec_comm_t *comm = new gmx_domdec_comm_t;

    comm->n_load_have      = 0;
    comm->n_load_collect   = 0;

    comm->haveTurnedOffDlb = false;

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
    comm->load_mdf  = 0;
    comm->load_pme  = 0;

    /* This should be replaced by a unique pointer */
    comm->balanceRegion = ddBalanceRegionAllocate();

    return comm;
}

/*! \brief Set the cell size and interaction limits, as well as the DD grid */
static void set_dd_limits_and_grid(FILE *fplog, t_commrec *cr, gmx_domdec_t *dd,
                                   const DomdecOptions &options,
                                   const MdrunOptions &mdrunOptions,
                                   const gmx_mtop_t *mtop,
                                   const t_inputrec *ir,
                                   const matrix box,
                                   gmx::ArrayRef<const gmx::RVec> xGlobal,
                                   gmx_ddbox_t *ddbox)
{
    real               r_bonded         = -1;
    real               r_bonded_limit   = -1;
    const real         tenPercentMargin = 1.1;
    gmx_domdec_comm_t *comm             = dd->comm;

    dd->npbcdim   = ePBC2npbcdim(ir->ePBC);
    dd->bScrewPBC = (ir->ePBC == epbcSCREW);

    dd->pme_recv_f_alloc = 0;
    dd->pme_recv_f_buf   = nullptr;

    /* Initialize to GPU share count to 0, might change later */
    comm->nrank_gpu_shared = 0;

    comm->dlbState         = determineInitialDlbState(fplog, cr, options.dlbOption, comm->bRecordLoad, mdrunOptions, ir);
    dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, TRUE);
    /* To consider turning DLB on after 2*nstlist steps we need to check
     * at partitioning count 3. Thus we need to increase the first count by 2.
     */
    comm->ddPartioningCountFirstDlbOff += 2;

    if (fplog)
    {
        fprintf(fplog, "Dynamic load balancing: %s\n",
                edlbs_names[comm->dlbState]);
    }
    comm->bPMELoadBalDLBLimits = FALSE;

    /* Allocate the charge group/atom sorting struct */
    comm->sort = gmx::compat::make_unique<gmx_domdec_sort_t>();

    comm->bCGs = (ncg_mtop(mtop) < mtop->natoms);

    comm->bInterCGBondeds = ((ncg_mtop(mtop) > gmx_mtop_num_molecules(*mtop)) ||
                             mtop->bIntermolecularInteractions);
    if (comm->bInterCGBondeds)
    {
        comm->bInterCGMultiBody = (multi_body_bondeds_count(mtop) > 0);
    }
    else
    {
        comm->bInterCGMultiBody = FALSE;
    }

    dd->bInterCGcons    = gmx::inter_charge_group_constraints(*mtop);
    dd->bInterCGsettles = gmx::inter_charge_group_settles(*mtop);

    if (ir->rlist == 0)
    {
        /* Set the cut-off to some very large value,
         * so we don't need if statements everywhere in the code.
         * We use sqrt, since the cut-off is squared in some places.
         */
        comm->cutoff   = GMX_CUTOFF_INF;
    }
    else
    {
        comm->cutoff   = ir->rlist;
    }
    comm->cutoff_mbody = 0;

    comm->cellsize_limit = 0;
    comm->bBondComm      = FALSE;

    /* Atoms should be able to move by up to half the list buffer size (if > 0)
     * within nstlist steps. Since boundaries are allowed to displace by half
     * a cell size, DD cells should be at least the size of the list buffer.
     */
    comm->cellsize_limit = std::max(comm->cellsize_limit,
                                    ir->rlist - std::max(ir->rvdw, ir->rcoulomb));

    if (comm->bInterCGBondeds)
    {
        if (options.minimumCommunicationRange > 0)
        {
            comm->cutoff_mbody = options.minimumCommunicationRange;
            if (options.useBondedCommunication)
            {
                comm->bBondComm = (comm->cutoff_mbody > comm->cutoff);
            }
            else
            {
                comm->cutoff = std::max(comm->cutoff, comm->cutoff_mbody);
            }
            r_bonded_limit = comm->cutoff_mbody;
        }
        else if (ir->bPeriodicMols)
        {
            /* Can not easily determine the required cut-off */
            dd_warning(cr, fplog, "NOTE: Periodic molecules are present in this system. Because of this, the domain decomposition algorithm cannot easily determine the minimum cell size that it requires for treating bonded interactions. Instead, domain decomposition will assume that half the non-bonded cut-off will be a suitable lower bound.\n");
            comm->cutoff_mbody = comm->cutoff/2;
            r_bonded_limit     = comm->cutoff_mbody;
        }
        else
        {
            real r_2b, r_mb;

            if (MASTER(cr))
            {
                dd_bonded_cg_distance(fplog, mtop, ir, as_rvec_array(xGlobal.data()), box,
                                      options.checkBondedInteractions,
                                      &r_2b, &r_mb);
            }
            gmx_bcast(sizeof(r_2b), &r_2b, cr);
            gmx_bcast(sizeof(r_mb), &r_mb, cr);

            /* We use an initial margin of 10% for the minimum cell size,
             * except when we are just below the non-bonded cut-off.
             */
            if (options.useBondedCommunication)
            {
                if (std::max(r_2b, r_mb) > comm->cutoff)
                {
                    r_bonded        = std::max(r_2b, r_mb);
                    r_bonded_limit  = tenPercentMargin*r_bonded;
                    comm->bBondComm = TRUE;
                }
                else
                {
                    r_bonded       = r_mb;
                    r_bonded_limit = std::min(tenPercentMargin*r_bonded, comm->cutoff);
                }
                /* We determine cutoff_mbody later */
            }
            else
            {
                /* No special bonded communication,
                 * simply increase the DD cut-off.
                 */
                r_bonded_limit     = tenPercentMargin*std::max(r_2b, r_mb);
                comm->cutoff_mbody = r_bonded_limit;
                comm->cutoff       = std::max(comm->cutoff, comm->cutoff_mbody);
            }
        }
        if (fplog)
        {
            fprintf(fplog,
                    "Minimum cell size due to bonded interactions: %.3f nm\n",
                    r_bonded_limit);
        }
        comm->cellsize_limit = std::max(comm->cellsize_limit, r_bonded_limit);
    }

    real rconstr = 0;
    if (dd->bInterCGcons && options.constraintCommunicationRange <= 0)
    {
        /* There is a cell size limit due to the constraints (P-LINCS) */
        rconstr = gmx::constr_r_max(fplog, mtop, ir);
        if (fplog)
        {
            fprintf(fplog,
                    "Estimated maximum distance required for P-LINCS: %.3f nm\n",
                    rconstr);
            if (rconstr > comm->cellsize_limit)
            {
                fprintf(fplog, "This distance will limit the DD cell size, you can override this with -rcon\n");
            }
        }
    }
    else if (options.constraintCommunicationRange > 0 && fplog)
    {
        /* Here we do not check for dd->bInterCGcons,
         * because one can also set a cell size limit for virtual sites only
         * and at this point we don't know yet if there are intercg v-sites.
         */
        fprintf(fplog,
                "User supplied maximum distance required for P-LINCS: %.3f nm\n",
                options.constraintCommunicationRange);
        rconstr = options.constraintCommunicationRange;
    }
    comm->cellsize_limit = std::max(comm->cellsize_limit, rconstr);

    comm->cgs_gl = gmx_mtop_global_cgs(mtop);

    if (options.numCells[XX] > 0)
    {
        copy_ivec(options.numCells, dd->nc);
        set_dd_dim(fplog, dd);
        set_ddbox_cr(cr, &dd->nc, ir, box, xGlobal, ddbox);

        if (options.numPmeRanks >= 0)
        {
            cr->npmenodes = options.numPmeRanks;
        }
        else
        {
            /* When the DD grid is set explicitly and -npme is set to auto,
             * don't use PME ranks. We check later if the DD grid is
             * compatible with the total number of ranks.
             */
            cr->npmenodes = 0;
        }

        real acs = average_cellsize_min(dd, ddbox);
        if (acs < comm->cellsize_limit)
        {
            if (fplog)
            {
                fprintf(fplog, "ERROR: The initial cell size (%f) is smaller than the cell size limit (%f)\n", acs, comm->cellsize_limit);
            }
            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "The initial cell size (%f) is smaller than the cell size limit (%f), change options -dd, -rdd or -rcon, see the log file for details",
                                 acs, comm->cellsize_limit);
        }
    }
    else
    {
        set_ddbox_cr(cr, nullptr, ir, box, xGlobal, ddbox);

        /* We need to choose the optimal DD grid and possibly PME nodes */
        real limit =
            dd_choose_grid(fplog, cr, dd, ir, mtop, box, ddbox,
                           options.numPmeRanks,
                           !isDlbDisabled(comm),
                           options.dlbScaling,
                           comm->cellsize_limit, comm->cutoff,
                           comm->bInterCGBondeds);

        if (dd->nc[XX] == 0)
        {
            char     buf[STRLEN];
            gmx_bool bC = (dd->bInterCGcons && rconstr > r_bonded_limit);
            sprintf(buf, "Change the number of ranks or mdrun option %s%s%s",
                    !bC ? "-rdd" : "-rcon",
                    comm->dlbState != edlbsOffUser ? " or -dds" : "",
                    bC ? " or your LINCS settings" : "");

            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "There is no domain decomposition for %d ranks that is compatible with the given box and a minimum cell size of %g nm\n"
                                 "%s\n"
                                 "Look in the log file for details on the domain decomposition",
                                 cr->nnodes-cr->npmenodes, limit, buf);
        }
        set_dd_dim(fplog, dd);
    }

    if (fplog)
    {
        fprintf(fplog,
                "Domain decomposition grid %d x %d x %d, separate PME ranks %d\n",
                dd->nc[XX], dd->nc[YY], dd->nc[ZZ], cr->npmenodes);
    }

    dd->nnodes = dd->nc[XX]*dd->nc[YY]*dd->nc[ZZ];
    if (cr->nnodes - dd->nnodes != cr->npmenodes)
    {
        gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                             "The size of the domain decomposition grid (%d) does not match the number of ranks (%d). The total number of ranks is %d",
                             dd->nnodes, cr->nnodes - cr->npmenodes, cr->nnodes);
    }
    if (cr->npmenodes > dd->nnodes)
    {
        gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                             "The number of separate PME ranks (%d) is larger than the number of PP ranks (%d), this is not supported.", cr->npmenodes, dd->nnodes);
    }
    if (cr->npmenodes > 0)
    {
        comm->npmenodes = cr->npmenodes;
    }
    else
    {
        comm->npmenodes = dd->nnodes;
    }

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        /* The following choices should match those
         * in comm_cost_est in domdec_setup.c.
         * Note that here the checks have to take into account
         * that the decomposition might occur in a different order than xyz
         * (for instance through the env.var. GMX_DD_ORDER_ZYX),
         * in which case they will not match those in comm_cost_est,
         * but since that is mainly for testing purposes that's fine.
         */
        if (dd->ndim >= 2 && dd->dim[0] == XX && dd->dim[1] == YY &&
            comm->npmenodes > dd->nc[XX] && comm->npmenodes % dd->nc[XX] == 0 &&
            getenv("GMX_PMEONEDD") == nullptr)
        {
            comm->npmedecompdim = 2;
            comm->npmenodes_x   = dd->nc[XX];
            comm->npmenodes_y   = comm->npmenodes/comm->npmenodes_x;
        }
        else
        {
            /* In case nc is 1 in both x and y we could still choose to
             * decompose pme in y instead of x, but we use x for simplicity.
             */
            comm->npmedecompdim = 1;
            if (dd->dim[0] == YY)
            {
                comm->npmenodes_x = 1;
                comm->npmenodes_y = comm->npmenodes;
            }
            else
            {
                comm->npmenodes_x = comm->npmenodes;
                comm->npmenodes_y = 1;
            }
        }
        if (fplog)
        {
            fprintf(fplog, "PME domain decomposition: %d x %d x %d\n",
                    comm->npmenodes_x, comm->npmenodes_y, 1);
        }
    }
    else
    {
        comm->npmedecompdim = 0;
        comm->npmenodes_x   = 0;
        comm->npmenodes_y   = 0;
    }

    snew(comm->slb_frac, DIM);
    if (isDlbDisabled(comm))
    {
        comm->slb_frac[XX] = get_slb_frac(fplog, "x", dd->nc[XX], options.cellSizeX);
        comm->slb_frac[YY] = get_slb_frac(fplog, "y", dd->nc[YY], options.cellSizeY);
        comm->slb_frac[ZZ] = get_slb_frac(fplog, "z", dd->nc[ZZ], options.cellSizeZ);
    }

    if (comm->bInterCGBondeds && comm->cutoff_mbody == 0)
    {
        if (comm->bBondComm || !isDlbDisabled(comm))
        {
            /* Set the bonded communication distance to halfway
             * the minimum and the maximum,
             * since the extra communication cost is nearly zero.
             */
            real acs           = average_cellsize_min(dd, ddbox);
            comm->cutoff_mbody = 0.5*(r_bonded + acs);
            if (!isDlbDisabled(comm))
            {
                /* Check if this does not limit the scaling */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody,
                                              options.dlbScaling*acs);
            }
            if (!comm->bBondComm)
            {
                /* Without bBondComm do not go beyond the n.b. cut-off */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody, comm->cutoff);
                if (comm->cellsize_limit >= comm->cutoff)
                {
                    /* We don't loose a lot of efficieny
                     * when increasing it to the n.b. cut-off.
                     * It can even be slightly faster, because we need
                     * less checks for the communication setup.
                     */
                    comm->cutoff_mbody = comm->cutoff;
                }
            }
            /* Check if we did not end up below our original limit */
            comm->cutoff_mbody = std::max(comm->cutoff_mbody, r_bonded_limit);

            if (comm->cutoff_mbody > comm->cellsize_limit)
            {
                comm->cellsize_limit = comm->cutoff_mbody;
            }
        }
        /* Without DLB and cutoff_mbody<cutoff, cutoff_mbody is dynamic */
    }

    if (debug)
    {
        fprintf(debug, "Bonded atom communication beyond the cut-off: %d\n"
                "cellsize limit %f\n",
                comm->bBondComm, comm->cellsize_limit);
    }

    if (MASTER(cr))
    {
        check_dd_restrictions(cr, dd, ir, fplog);
    }
}

static void set_dlb_limits(gmx_domdec_t *dd)

{
    for (int d = 0; d < dd->ndim; d++)
    {
        /* Set the number of pulses to the value for DLB */
        dd->comm->cd[d].ind.resize(dd->comm->cd[d].np_dlb);

        dd->comm->cellsize_min[dd->dim[d]] =
            dd->comm->cellsize_min_dlb[dd->dim[d]];
    }
}


static void turn_on_dlb(FILE *fplog, const t_commrec *cr, int64_t step)
{
    gmx_domdec_t      *dd           = cr->dd;
    gmx_domdec_comm_t *comm         = dd->comm;

    real               cellsize_min = comm->cellsize_min[dd->dim[0]];
    for (int d = 1; d < dd->ndim; d++)
    {
        cellsize_min = std::min(cellsize_min, comm->cellsize_min[dd->dim[d]]);
    }

    /* Turn off DLB if we're too close to the cell size limit. */
    if (cellsize_min < comm->cellsize_limit*1.05)
    {
        auto str = gmx::formatString("step %" PRId64 " Measured %.1f %% performance loss due to load imbalance, "
                                     "but the minimum cell size is smaller than 1.05 times the cell size limit."
                                     "Will no longer try dynamic load balancing.\n", step, dd_force_imb_perf_loss(dd)*100);
        dd_warning(cr, fplog, str.c_str());

        comm->dlbState = edlbsOffForever;
        return;
    }

    char buf[STRLEN];
    sprintf(buf, "step %" PRId64 " Turning on dynamic load balancing, because the performance loss due to load imbalance is %.1f %%.\n", step, dd_force_imb_perf_loss(dd)*100);
    dd_warning(cr, fplog, buf);
    comm->dlbState = edlbsOnCanTurnOff;

    /* Store the non-DLB performance, so we can check if DLB actually
     * improves performance.
     */
    GMX_RELEASE_ASSERT(comm->cycl_n[ddCyclStep] > 0, "When we turned on DLB, we should have measured cycles");
    comm->cyclesPerStepBeforeDLB = comm->cycl[ddCyclStep]/comm->cycl_n[ddCyclStep];

    set_dlb_limits(dd);

    /* We can set the required cell size info here,
     * so we do not need to communicate this.
     * The grid is completely uniform.
     */
    for (int d = 0; d < dd->ndim; d++)
    {
        RowMaster *rowMaster = comm->cellsizesWithDlb[d].rowMaster.get();

        if (rowMaster)
        {
            comm->load[d].sum_m = comm->load[d].sum;

            int nc = dd->nc[dd->dim[d]];
            for (int i = 0; i < nc; i++)
            {
                rowMaster->cellFrac[i] = i/static_cast<real>(nc);
                if (d > 0)
                {
                    rowMaster->bounds[i].cellFracLowerMax =  i     /static_cast<real>(nc);
                    rowMaster->bounds[i].cellFracUpperMin = (i + 1)/static_cast<real>(nc);
                }
            }
            rowMaster->cellFrac[nc] = 1.0;
        }
    }
}

static void turn_off_dlb(FILE *fplog, const t_commrec *cr, int64_t step)
{
    gmx_domdec_t *dd = cr->dd;

    char          buf[STRLEN];
    sprintf(buf, "step %" PRId64 " Turning off dynamic load balancing, because it is degrading performance.\n", step);
    dd_warning(cr, fplog, buf);
    dd->comm->dlbState                     = edlbsOffCanTurnOn;
    dd->comm->haveTurnedOffDlb             = true;
    dd->comm->ddPartioningCountFirstDlbOff = dd->ddp_count;
}

static void turn_off_dlb_forever(FILE *fplog, const t_commrec *cr, int64_t step)
{
    GMX_RELEASE_ASSERT(cr->dd->comm->dlbState == edlbsOffCanTurnOn, "Can only turn off DLB forever when it was in the can-turn-on state");
    char buf[STRLEN];
    sprintf(buf, "step %" PRId64 " Will no longer try dynamic load balancing, as it degraded performance.\n", step);
    dd_warning(cr, fplog, buf);
    cr->dd->comm->dlbState = edlbsOffForever;
}

static char *init_bLocalCG(const gmx_mtop_t *mtop)
{
    int   ncg, cg;
    char *bLocalCG;

    ncg = ncg_mtop(mtop);
    snew(bLocalCG, ncg);
    for (cg = 0; cg < ncg; cg++)
    {
        bLocalCG[cg] = FALSE;
    }

    return bLocalCG;
}

void dd_init_bondeds(FILE *fplog,
                     gmx_domdec_t *dd,
                     const gmx_mtop_t *mtop,
                     const gmx_vsite_t *vsite,
                     const t_inputrec *ir,
                     gmx_bool bBCheck, cginfo_mb_t *cginfo_mb)
{
    gmx_domdec_comm_t *comm;

    dd_make_reverse_top(fplog, dd, mtop, vsite, ir, bBCheck);

    comm = dd->comm;

    if (comm->bBondComm)
    {
        /* Communicate atoms beyond the cut-off for bonded interactions */
        comm = dd->comm;

        comm->cglink = make_charge_group_links(mtop, dd, cginfo_mb);

        comm->bLocalCG = init_bLocalCG(mtop);
    }
    else
    {
        /* Only communicate atoms based on cut-off */
        comm->cglink   = nullptr;
        comm->bLocalCG = nullptr;
    }
}

static void print_dd_settings(FILE *fplog, gmx_domdec_t *dd,
                              const gmx_mtop_t *mtop, const t_inputrec *ir,
                              gmx_bool bDynLoadBal, real dlb_scale,
                              const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                d;
    ivec               np;
    real               limit, shrink;
    char               buf[64];

    if (fplog == nullptr)
    {
        return;
    }

    comm = dd->comm;

    if (bDynLoadBal)
    {
        fprintf(fplog, "The maximum number of communication pulses is:");
        for (d = 0; d < dd->ndim; d++)
        {
            fprintf(fplog, " %c %d", dim2char(dd->dim[d]), comm->cd[d].np_dlb);
        }
        fprintf(fplog, "\n");
        fprintf(fplog, "The minimum size for domain decomposition cells is %.3f nm\n", comm->cellsize_limit);
        fprintf(fplog, "The requested allowed shrink of DD cells (option -dds) is: %.2f\n", dlb_scale);
        fprintf(fplog, "The allowed shrink of domain decomposition cells is:");
        for (d = 0; d < DIM; d++)
        {
            if (dd->nc[d] > 1)
            {
                if (d >= ddbox->npbcdim && dd->nc[d] == 2)
                {
                    shrink = 0;
                }
                else
                {
                    shrink =
                        comm->cellsize_min_dlb[d]/
                        (ddbox->box_size[d]*ddbox->skew_fac[d]/dd->nc[d]);
                }
                fprintf(fplog, " %c %.2f", dim2char(d), shrink);
            }
        }
        fprintf(fplog, "\n");
    }
    else
    {
        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbPULSE_ONLY, np);
        fprintf(fplog, "The initial number of communication pulses is:");
        for (d = 0; d < dd->ndim; d++)
        {
            fprintf(fplog, " %c %d", dim2char(dd->dim[d]), np[dd->dim[d]]);
        }
        fprintf(fplog, "\n");
        fprintf(fplog, "The initial domain decomposition cell size is:");
        for (d = 0; d < DIM; d++)
        {
            if (dd->nc[d] > 1)
            {
                fprintf(fplog, " %c %.2f nm",
                        dim2char(d), dd->comm->cellsize_min[d]);
            }
        }
        fprintf(fplog, "\n\n");
    }

    gmx_bool bInterCGVsites = count_intercg_vsites(mtop);

    if (comm->bInterCGBondeds ||
        bInterCGVsites ||
        dd->bInterCGcons || dd->bInterCGsettles)
    {
        fprintf(fplog, "The maximum allowed distance for charge groups involved in interactions is:\n");
        fprintf(fplog, "%40s  %-7s %6.3f nm\n",
                "non-bonded interactions", "", comm->cutoff);

        if (bDynLoadBal)
        {
            limit = dd->comm->cellsize_limit;
        }
        else
        {
            if (dynamic_dd_box(ddbox, ir))
            {
                fprintf(fplog, "(the following are initial values, they could change due to box deformation)\n");
            }
            limit = dd->comm->cellsize_min[XX];
            for (d = 1; d < DIM; d++)
            {
                limit = std::min(limit, dd->comm->cellsize_min[d]);
            }
        }

        if (comm->bInterCGBondeds)
        {
            fprintf(fplog, "%40s  %-7s %6.3f nm\n",
                    "two-body bonded interactions", "(-rdd)",
                    std::max(comm->cutoff, comm->cutoff_mbody));
            fprintf(fplog, "%40s  %-7s %6.3f nm\n",
                    "multi-body bonded interactions", "(-rdd)",
                    (comm->bBondComm || isDlbOn(dd->comm)) ? comm->cutoff_mbody : std::min(comm->cutoff, limit));
        }
        if (bInterCGVsites)
        {
            fprintf(fplog, "%40s  %-7s %6.3f nm\n",
                    "virtual site constructions", "(-rcon)", limit);
        }
        if (dd->bInterCGcons || dd->bInterCGsettles)
        {
            sprintf(buf, "atoms separated by up to %d constraints",
                    1+ir->nProjOrder);
            fprintf(fplog, "%40s  %-7s %6.3f nm\n",
                    buf, "(-rcon)", limit);
        }
        fprintf(fplog, "\n");
    }

    fflush(fplog);
}

static void set_cell_limits_dlb(gmx_domdec_t      *dd,
                                real               dlb_scale,
                                const t_inputrec  *ir,
                                const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                d, dim, npulse, npulse_d_max, npulse_d;
    gmx_bool           bNoCutOff;

    comm = dd->comm;

    bNoCutOff = (ir->rvdw == 0 || ir->rcoulomb == 0);

    /* Determine the maximum number of comm. pulses in one dimension */

    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);

    /* Determine the maximum required number of grid pulses */
    if (comm->cellsize_limit >= comm->cutoff)
    {
        /* Only a single pulse is required */
        npulse = 1;
    }
    else if (!bNoCutOff && comm->cellsize_limit > 0)
    {
        /* We round down slightly here to avoid overhead due to the latency
         * of extra communication calls when the cut-off
         * would be only slightly longer than the cell size.
         * Later cellsize_limit is redetermined,
         * so we can not miss interactions due to this rounding.
         */
        npulse = (int)(0.96 + comm->cutoff/comm->cellsize_limit);
    }
    else
    {
        /* There is no cell size limit */
        npulse = std::max(dd->nc[XX]-1, std::max(dd->nc[YY]-1, dd->nc[ZZ]-1));
    }

    if (!bNoCutOff && npulse > 1)
    {
        /* See if we can do with less pulses, based on dlb_scale */
        npulse_d_max = 0;
        for (d = 0; d < dd->ndim; d++)
        {
            dim      = dd->dim[d];
            npulse_d = (int)(1 + dd->nc[dim]*comm->cutoff
                             /(ddbox->box_size[dim]*ddbox->skew_fac[dim]*dlb_scale));
            npulse_d_max = std::max(npulse_d_max, npulse_d);
        }
        npulse = std::min(npulse, npulse_d_max);
    }

    /* This env var can override npulse */
    d = dd_getenv(debug, "GMX_DD_NPULSE", 0);
    if (d > 0)
    {
        npulse = d;
    }

    comm->maxpulse       = 1;
    comm->bVacDLBNoLimit = (ir->ePBC == epbcNONE);
    for (d = 0; d < dd->ndim; d++)
    {
        comm->cd[d].np_dlb    = std::min(npulse, dd->nc[dd->dim[d]]-1);
        comm->maxpulse        = std::max(comm->maxpulse, comm->cd[d].np_dlb);
        if (comm->cd[d].np_dlb < dd->nc[dd->dim[d]]-1)
        {
            comm->bVacDLBNoLimit = FALSE;
        }
    }

    /* cellsize_limit is set for LINCS in init_domain_decomposition */
    if (!comm->bVacDLBNoLimit)
    {
        comm->cellsize_limit = std::max(comm->cellsize_limit,
                                        comm->cutoff/comm->maxpulse);
    }
    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);
    /* Set the minimum cell size for each DD dimension */
    for (d = 0; d < dd->ndim; d++)
    {
        if (comm->bVacDLBNoLimit ||
            comm->cd[d].np_dlb*comm->cellsize_limit >= comm->cutoff)
        {
            comm->cellsize_min_dlb[dd->dim[d]] = comm->cellsize_limit;
        }
        else
        {
            comm->cellsize_min_dlb[dd->dim[d]] =
                comm->cutoff/comm->cd[d].np_dlb;
        }
    }
    if (comm->cutoff_mbody <= 0)
    {
        comm->cutoff_mbody = std::min(comm->cutoff, comm->cellsize_limit);
    }
    if (isDlbOn(comm))
    {
        set_dlb_limits(dd);
    }
}

gmx_bool dd_bonded_molpbc(const gmx_domdec_t *dd, int ePBC)
{
    /* If each molecule is a single charge group
     * or we use domain decomposition for each periodic dimension,
     * we do not need to take pbc into account for the bonded interactions.
     */
    return (ePBC != epbcNONE && dd->comm->bInterCGBondeds &&
            !(dd->nc[XX] > 1 &&
              dd->nc[YY] > 1 &&
              (dd->nc[ZZ] > 1 || ePBC == epbcXY)));
}

/*! \brief Sets grid size limits and PP-PME setup, prints settings to log */
static void set_ddgrid_parameters(FILE *fplog, gmx_domdec_t *dd, real dlb_scale,
                                  const gmx_mtop_t *mtop, const t_inputrec *ir,
                                  const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                natoms_tot;
    real               vol_frac;

    comm = dd->comm;

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        init_ddpme(dd, &comm->ddpme[0], 0);
        if (comm->npmedecompdim >= 2)
        {
            init_ddpme(dd, &comm->ddpme[1], 1);
        }
    }
    else
    {
        comm->npmenodes = 0;
        if (dd->pme_nodeid >= 0)
        {
            gmx_fatal_collective(FARGS, dd->mpi_comm_all, DDMASTER(dd),
                                 "Can not have separate PME ranks without PME electrostatics");
        }
    }

    if (debug)
    {
        fprintf(debug, "The DD cut-off is %f\n", comm->cutoff);
    }
    if (!isDlbDisabled(comm))
    {
        set_cell_limits_dlb(dd, dlb_scale, ir, ddbox);
    }

    print_dd_settings(fplog, dd, mtop, ir, isDlbOn(comm), dlb_scale, ddbox);
    if (comm->dlbState == edlbsOffCanTurnOn)
    {
        if (fplog)
        {
            fprintf(fplog, "When dynamic load balancing gets turned on, these settings will change to:\n");
        }
        print_dd_settings(fplog, dd, mtop, ir, TRUE, dlb_scale, ddbox);
    }

    if (ir->ePBC == epbcNONE)
    {
        vol_frac = 1 - 1/(double)dd->nnodes;
    }
    else
    {
        vol_frac =
            (1 + comm_box_frac(dd->nc, comm->cutoff, ddbox))/(double)dd->nnodes;
    }
    if (debug)
    {
        fprintf(debug, "Volume fraction for all DD zones: %f\n", vol_frac);
    }
    natoms_tot = comm->cgs_gl.index[comm->cgs_gl.nr];

    dd->ga2la = ga2la_init(natoms_tot, static_cast<int>(vol_frac*natoms_tot));
}

/*! \brief Set some important DD parameters that can be modified by env.vars */
static void set_dd_envvar_options(FILE *fplog, gmx_domdec_t *dd, int rank_mysim)
{
    gmx_domdec_comm_t *comm = dd->comm;

    dd->bSendRecv2      = dd_getenv(fplog, "GMX_DD_USE_SENDRECV2", 0);
    comm->dlb_scale_lim = dd_getenv(fplog, "GMX_DLB_MAX_BOX_SCALING", 10);
    comm->eFlop         = dd_getenv(fplog, "GMX_DLB_BASED_ON_FLOPS", 0);
    int recload         = dd_getenv(fplog, "GMX_DD_RECORD_LOAD", 1);
    comm->nstDDDump     = dd_getenv(fplog, "GMX_DD_NST_DUMP", 0);
    comm->nstDDDumpGrid = dd_getenv(fplog, "GMX_DD_NST_DUMP_GRID", 0);
    comm->DD_debug      = dd_getenv(fplog, "GMX_DD_DEBUG", 0);

    if (dd->bSendRecv2 && fplog)
    {
        fprintf(fplog, "Will use two sequential MPI_Sendrecv calls instead of two simultaneous non-blocking MPI_Irecv and MPI_Isend pairs for constraint and vsite communication\n");
    }

    if (comm->eFlop)
    {
        if (fplog)
        {
            fprintf(fplog, "Will load balance based on FLOP count\n");
        }
        if (comm->eFlop > 1)
        {
            srand(1 + rank_mysim);
        }
        comm->bRecordLoad = TRUE;
    }
    else
    {
        comm->bRecordLoad = (wallcycle_have_counter() && recload > 0);
    }
}

DomdecOptions::DomdecOptions() :
    checkBondedInteractions(TRUE),
    useBondedCommunication(TRUE),
    numPmeRanks(-1),
    rankOrder(DdRankOrder::pp_pme),
    minimumCommunicationRange(0),
    constraintCommunicationRange(0),
    dlbOption(DlbOption::turnOnWhenUseful),
    dlbScaling(0.8),
    cellSizeX(nullptr),
    cellSizeY(nullptr),
    cellSizeZ(nullptr)
{
    clear_ivec(numCells);
}

gmx_domdec_t *init_domain_decomposition(FILE *fplog, t_commrec *cr,
                                        const DomdecOptions &options,
                                        const MdrunOptions &mdrunOptions,
                                        const gmx_mtop_t *mtop,
                                        const t_inputrec *ir,
                                        const matrix box,
                                        gmx::ArrayRef<const gmx::RVec> xGlobal)
{
    gmx_domdec_t      *dd;

    if (fplog)
    {
        fprintf(fplog,
                "\nInitializing Domain Decomposition on %d ranks\n", cr->nnodes);
    }

    dd = new gmx_domdec_t;

    dd->comm = init_dd_comm();

    /* Initialize DD paritioning counters */
    dd->comm->partition_step = INT_MIN;
    dd->ddp_count            = 0;

    set_dd_envvar_options(fplog, dd, cr->nodeid);

    gmx_ddbox_t ddbox = {0};
    set_dd_limits_and_grid(fplog, cr, dd, options, mdrunOptions,
                           mtop, ir,
                           box, xGlobal,
                           &ddbox);

    make_dd_communicators(fplog, cr, dd, options.rankOrder);

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        set_ddgrid_parameters(fplog, dd, options.dlbScaling, mtop, ir, &ddbox);

        setup_neighbor_relations(dd);
    }

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);

    clear_dd_cycle_counts(dd);

    return dd;
}

static gmx_bool test_dd_cutoff(t_commrec *cr,
                               t_state *state, const t_inputrec *ir,
                               real cutoff_req)
{
    gmx_domdec_t *dd;
    gmx_ddbox_t   ddbox;
    int           d, dim, np;
    real          inv_cell_size;
    int           LocallyLimited;

    dd = cr->dd;

    set_ddbox(dd, false, ir, state->box, true, state->x, &ddbox);

    LocallyLimited = 0;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];

        inv_cell_size = DD_CELL_MARGIN*dd->nc[dim]/ddbox.box_size[dim];
        if (dynamic_dd_box(&ddbox, ir))
        {
            inv_cell_size *= DD_PRES_SCALE_MARGIN;
        }

        np = 1 + (int)(cutoff_req*inv_cell_size*ddbox.skew_fac[dim]);

        if (!isDlbDisabled(dd->comm) && (dim < ddbox.npbcdim) && (dd->comm->cd[d].np_dlb > 0))
        {
            if (np > dd->comm->cd[d].np_dlb)
            {
                return FALSE;
            }

            /* If a current local cell size is smaller than the requested
             * cut-off, we could still fix it, but this gets very complicated.
             * Without fixing here, we might actually need more checks.
             */
            if ((dd->comm->cell_x1[dim] - dd->comm->cell_x0[dim])*ddbox.skew_fac[dim]*dd->comm->cd[d].np_dlb < cutoff_req)
            {
                LocallyLimited = 1;
            }
        }
    }

    if (!isDlbDisabled(dd->comm))
    {
        /* If DLB is not active yet, we don't need to check the grid jumps.
         * Actually we shouldn't, because then the grid jump data is not set.
         */
        if (isDlbOn(dd->comm) &&
            check_grid_jump(0, dd, cutoff_req, &ddbox, FALSE))
        {
            LocallyLimited = 1;
        }

        gmx_sumi(1, &LocallyLimited, cr);

        if (LocallyLimited > 0)
        {
            return FALSE;
        }
    }

    return TRUE;
}

gmx_bool change_dd_cutoff(t_commrec *cr, t_state *state, const t_inputrec *ir,
                          real cutoff_req)
{
    gmx_bool bCutoffAllowed;

    bCutoffAllowed = test_dd_cutoff(cr, state, ir, cutoff_req);

    if (bCutoffAllowed)
    {
        cr->dd->comm->cutoff = cutoff_req;
    }

    return bCutoffAllowed;
}

void set_dd_dlb_max_cutoff(t_commrec *cr, real cutoff)
{
    gmx_domdec_comm_t *comm;

    comm = cr->dd->comm;

    /* Turn on the DLB limiting (might have been on already) */
    comm->bPMELoadBalDLBLimits = TRUE;

    /* Change the cut-off limit */
    comm->PMELoadBal_max_cutoff = cutoff;

    if (debug)
    {
        fprintf(debug, "PME load balancing set a limit to the DLB staggering such that a %f cut-off will continue to fit\n",
                comm->PMELoadBal_max_cutoff);
    }
}

/* Sets whether we should later check the load imbalance data, so that
 * we can trigger dynamic load balancing if enough imbalance has
 * arisen.
 *
 * Used after PME load balancing unlocks DLB, so that the check
 * whether DLB will be useful can happen immediately.
 */
static void dd_dlb_set_should_check_whether_to_turn_dlb_on(gmx_domdec_t *dd, gmx_bool bValue)
{
    if (dd->comm->dlbState == edlbsOffCanTurnOn)
    {
        dd->comm->bCheckWhetherToTurnDlbOn = bValue;

        if (bValue == TRUE)
        {
            /* Store the DD partitioning count, so we can ignore cycle counts
             * over the next nstlist steps, which are often slower.
             */
            dd->comm->ddPartioningCountFirstDlbOff = dd->ddp_count;
        }
    }
}

/* Returns if we should check whether there has been enough load
 * imbalance to trigger dynamic load balancing.
 */
static gmx_bool dd_dlb_get_should_check_whether_to_turn_dlb_on(gmx_domdec_t *dd)
{
    if (dd->comm->dlbState != edlbsOffCanTurnOn)
    {
        return FALSE;
    }

    if (dd->ddp_count <= dd->comm->ddPartioningCountFirstDlbOff)
    {
        /* We ignore the first nstlist steps at the start of the run
         * or after PME load balancing or after turning DLB off, since
         * these often have extra allocation or cache miss overhead.
         */
        return FALSE;
    }

    if (dd->comm->cycl_n[ddCyclStep] == 0)
    {
        /* We can have zero timed steps when dd_partition_system is called
         * more than once at the same step, e.g. with replica exchange.
         * Turning on DLB would trigger an assertion failure later, but is
         * also useless right after exchanging replicas.
         */
        return FALSE;
    }

    /* We should check whether we should use DLB directly after
     * unlocking DLB. */
    if (dd->comm->bCheckWhetherToTurnDlbOn)
    {
        /* This flag was set when the PME load-balancing routines
           unlocked DLB, and should now be cleared. */
        dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, FALSE);
        return TRUE;
    }
    /* We check whether we should use DLB every c_checkTurnDlbOnInterval
     * partitionings (we do not do this every partioning, so that we
     * avoid excessive communication). */
    if (dd->comm->n_load_have % c_checkTurnDlbOnInterval == c_checkTurnDlbOnInterval - 1)
    {
        return TRUE;
    }

    return FALSE;
}

gmx_bool dd_dlb_is_on(const gmx_domdec_t *dd)
{
    return isDlbOn(dd->comm);
}

gmx_bool dd_dlb_is_locked(const gmx_domdec_t *dd)
{
    return (dd->comm->dlbState == edlbsOffTemporarilyLocked);
}

void dd_dlb_lock(gmx_domdec_t *dd)
{
    /* We can only lock the DLB when it is set to auto, otherwise don't do anything */
    if (dd->comm->dlbState == edlbsOffCanTurnOn)
    {
        dd->comm->dlbState = edlbsOffTemporarilyLocked;
    }
}

void dd_dlb_unlock(gmx_domdec_t *dd)
{
    /* We can only lock the DLB when it is set to auto, otherwise don't do anything */
    if (dd->comm->dlbState == edlbsOffTemporarilyLocked)
    {
        dd->comm->dlbState = edlbsOffCanTurnOn;
        dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, TRUE);
    }
}

static void merge_cg_buffers(int ncell,
                             gmx_domdec_comm_dim_t *cd, int pulse,
                             int  *ncg_cell,
                             gmx::ArrayRef<int> index_gl,
                             const int  *recv_i,
                             rvec *cg_cm,    rvec *recv_vr,
                             gmx::ArrayRef<int> cgindex,
                             cginfo_mb_t *cginfo_mb, int *cginfo)
{
    gmx_domdec_ind_t *ind, *ind_p;
    int               p, cell, c, cg, cg0, cg1, cg_gl, nat;
    int               shift, shift_at;

    ind = &cd->ind[pulse];

    /* First correct the already stored data */
    shift = ind->nrecv[ncell];
    for (cell = ncell-1; cell >= 0; cell--)
    {
        shift -= ind->nrecv[cell];
        if (shift > 0)
        {
            /* Move the cg's present from previous grid pulses */
            cg0                = ncg_cell[ncell+cell];
            cg1                = ncg_cell[ncell+cell+1];
            cgindex[cg1+shift] = cgindex[cg1];
            for (cg = cg1-1; cg >= cg0; cg--)
            {
                index_gl[cg+shift] = index_gl[cg];
                copy_rvec(cg_cm[cg], cg_cm[cg+shift]);
                cgindex[cg+shift] = cgindex[cg];
                cginfo[cg+shift]  = cginfo[cg];
            }
            /* Correct the already stored send indices for the shift */
            for (p = 1; p <= pulse; p++)
            {
                ind_p = &cd->ind[p];
                cg0   = 0;
                for (c = 0; c < cell; c++)
                {
                    cg0 += ind_p->nsend[c];
                }
                cg1 = cg0 + ind_p->nsend[cell];
                for (cg = cg0; cg < cg1; cg++)
                {
                    ind_p->index[cg] += shift;
                }
            }
        }
    }

    /* Merge in the communicated buffers */
    shift    = 0;
    shift_at = 0;
    cg0      = 0;
    for (cell = 0; cell < ncell; cell++)
    {
        cg1 = ncg_cell[ncell+cell+1] + shift;
        if (shift_at > 0)
        {
            /* Correct the old cg indices */
            for (cg = ncg_cell[ncell+cell]; cg < cg1; cg++)
            {
                cgindex[cg+1] += shift_at;
            }
        }
        for (cg = 0; cg < ind->nrecv[cell]; cg++)
        {
            /* Copy this charge group from the buffer */
            index_gl[cg1] = recv_i[cg0];
            copy_rvec(recv_vr[cg0], cg_cm[cg1]);
            /* Add it to the cgindex */
            cg_gl          = index_gl[cg1];
            cginfo[cg1]    = ddcginfo(cginfo_mb, cg_gl);
            nat            = GET_CGINFO_NATOMS(cginfo[cg1]);
            cgindex[cg1+1] = cgindex[cg1] + nat;
            cg0++;
            cg1++;
            shift_at += nat;
        }
        shift                 += ind->nrecv[cell];
        ncg_cell[ncell+cell+1] = cg1;
    }
}

static void make_cell2at_index(gmx_domdec_comm_dim_t        *cd,
                               int                           nzone,
                               int                           atomGroupStart,
                               const gmx::RangePartitioning &atomGroups)
{
    /* Store the atom block boundaries for easy copying of communication buffers
     */
    int g = atomGroupStart;
    for (int zone = 0; zone < nzone; zone++)
    {
        for (gmx_domdec_ind_t &ind : cd->ind)
        {
            const auto range    = atomGroups.subRange(g, g + ind.nrecv[zone]);
            ind.cell2at0[zone]  = range.begin();
            ind.cell2at1[zone]  = range.end();
            g                  += ind.nrecv[zone];
        }
    }
}

static gmx_bool missing_link(t_blocka *link, int cg_gl, const char *bLocalCG)
{
    int      i;
    gmx_bool bMiss;

    bMiss = FALSE;
    for (i = link->index[cg_gl]; i < link->index[cg_gl+1]; i++)
    {
        if (!bLocalCG[link->a[i]])
        {
            bMiss = TRUE;
        }
    }

    return bMiss;
}

/* Domain corners for communication, a maximum of 4 i-zones see a j domain */
typedef struct {
    real c[DIM][4]; /* the corners for the non-bonded communication */
    real cr0;       /* corner for rounding */
    real cr1[4];    /* corners for rounding */
    real bc[DIM];   /* corners for bounded communication */
    real bcr1;      /* corner for rounding for bonded communication */
} dd_corners_t;

/* Determine the corners of the domain(s) we are communicating with */
static void
set_dd_corners(const gmx_domdec_t *dd,
               int dim0, int dim1, int dim2,
               gmx_bool bDistMB,
               dd_corners_t *c)
{
    const gmx_domdec_comm_t  *comm;
    const gmx_domdec_zones_t *zones;
    int i, j;

    comm = dd->comm;

    zones = &comm->zones;

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
        if (isDlbOn(dd->comm))
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
            for (j = 0; j < 4; j++)
            {
                c->c[2][j] = comm->cell_x0[dim2];
            }
            if (isDlbOn(dd->comm))
            {
                /* Use the maximum of the i-cells that see a j-cell */
                for (i = 0; i < zones->nizone; i++)
                {
                    for (j = zones->izone[i].j0; j < zones->izone[i].j1; j++)
                    {
                        if (j >= 4)
                        {
                            c->c[2][j-4] =
                                std::max(c->c[2][j-4],
                                         comm->zone_d2[zones->shift[i][dim0]][zones->shift[i][dim1]].mch0);
                        }
                    }
                }
                if (bDistMB)
                {
                    /* For the multi-body distance we need the maximum */
                    c->bc[2] = comm->cell_x0[dim2];
                    for (i = 0; i < 2; i++)
                    {
                        for (j = 0; j < 2; j++)
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
            if (isDlbOn(dd->comm))
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

/* Add the atom groups we need to send in this pulse from this zone to
 * \p localAtomGroups and \p work
 */
static void
get_zone_pulse_cgs(gmx_domdec_t *dd,
                   int zonei, int zone,
                   int cg0, int cg1,
                   gmx::ArrayRef<const int> globalAtomGroupIndices,
                   const gmx::RangePartitioning &atomGroups,
                   int dim, int dim_ind,
                   int dim0, int dim1, int dim2,
                   real r_comm2, real r_bcomm2,
                   matrix box,
                   bool distanceIsTriclinic,
                   rvec *normal,
                   real skew_fac2_d, real skew_fac_01,
                   rvec *v_d, rvec *v_0, rvec *v_1,
                   const dd_corners_t *c,
                   const rvec sf2_round,
                   gmx_bool bDistBonded,
                   gmx_bool bBondComm,
                   gmx_bool bDist2B,
                   gmx_bool bDistMB,
                   rvec *cg_cm,
                   const int *cginfo,
                   std::vector<int>     *localAtomGroups,
                   dd_comm_setup_work_t *work)
{
    gmx_domdec_comm_t *comm;
    gmx_bool           bScrew;
    gmx_bool           bDistMB_pulse;
    int                cg, i;
    real               r2, rb2, r, tric_sh;
    rvec               rn, rb;
    int                dimd;
    int                nsend_z, nat;

    comm = dd->comm;

    bScrew = (dd->bScrewPBC && dim == XX);

    bDistMB_pulse = (bDistMB && bDistBonded);

    /* Unpack the work data */
    std::vector<int>       &ibuf = work->atomGroupBuffer;
    std::vector<gmx::RVec> &vbuf = work->positionBuffer;
    nsend_z                      = 0;
    nat                          = work->nat;

    for (cg = cg0; cg < cg1; cg++)
    {
        r2  = 0;
        rb2 = 0;
        if (!distanceIsTriclinic)
        {
            /* Rectangular direction, easy */
            r = cg_cm[cg][dim] - c->c[dim_ind][zone];
            if (r > 0)
            {
                r2 += r*r;
            }
            if (bDistMB_pulse)
            {
                r = cg_cm[cg][dim] - c->bc[dim_ind];
                if (r > 0)
                {
                    rb2 += r*r;
                }
            }
            /* Rounding gives at most a 16% reduction
             * in communicated atoms
             */
            if (dim_ind >= 1 && (zonei == 1 || zonei == 2))
            {
                r = cg_cm[cg][dim0] - c->cr0;
                /* This is the first dimension, so always r >= 0 */
                r2 += r*r;
                if (bDistMB_pulse)
                {
                    rb2 += r*r;
                }
            }
            if (dim_ind == 2 && (zonei == 2 || zonei == 3))
            {
                r = cg_cm[cg][dim1] - c->cr1[zone];
                if (r > 0)
                {
                    r2 += r*r;
                }
                if (bDistMB_pulse)
                {
                    r = cg_cm[cg][dim1] - c->bcr1;
                    if (r > 0)
                    {
                        rb2 += r*r;
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
                rn[dim0] = cg_cm[cg][dim0] - c->cr0;
                for (i = dim0+1; i < DIM; i++)
                {
                    rn[dim0] -= cg_cm[cg][i]*v_0[i][dim0];
                }
                r2 = rn[dim0]*rn[dim0]*sf2_round[dim0];
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
                        rn[dimd] -= rn[dim0]*normal[dim0][dimd];
                        if (bDistMB_pulse)
                        {
                            rb[dimd] -= rb[dim0]*normal[dim0][dimd];
                        }
                    }
                }
            }
            if (dim_ind == 2 && (zonei == 2 || zonei == 3))
            {
                rn[dim1] += cg_cm[cg][dim1] - c->cr1[zone];
                tric_sh   = 0;
                for (i = dim1+1; i < DIM; i++)
                {
                    tric_sh -= cg_cm[cg][i]*v_1[i][dim1];
                }
                rn[dim1] += tric_sh;
                if (rn[dim1] > 0)
                {
                    r2 += rn[dim1]*rn[dim1]*sf2_round[dim1];
                    /* Take care of coupling of the distances
                     * to the planes along dim0 and dim1 through dim2.
                     */
                    r2 -= rn[dim0]*rn[dim1]*skew_fac_01;
                    /* Take care that the cell planes along dim1
                     * might not be orthogonal to that along dim2.
                     */
                    if (normal[dim1][dim2] > 0)
                    {
                        rn[dim2] -= rn[dim1]*normal[dim1][dim2];
                    }
                }
                if (bDistMB_pulse)
                {
                    rb[dim1] +=
                        cg_cm[cg][dim1] - c->bcr1 + tric_sh;
                    if (rb[dim1] > 0)
                    {
                        rb2 += rb[dim1]*rb[dim1]*sf2_round[dim1];
                        /* Take care of coupling of the distances
                         * to the planes along dim0 and dim1 through dim2.
                         */
                        rb2 -= rb[dim0]*rb[dim1]*skew_fac_01;
                        /* Take care that the cell planes along dim1
                         * might not be orthogonal to that along dim2.
                         */
                        if (normal[dim1][dim2] > 0)
                        {
                            rb[dim2] -= rb[dim1]*normal[dim1][dim2];
                        }
                    }
                }
            }
            /* The distance along the communication direction */
            rn[dim] += cg_cm[cg][dim] - c->c[dim_ind][zone];
            tric_sh  = 0;
            for (i = dim+1; i < DIM; i++)
            {
                tric_sh -= cg_cm[cg][i]*v_d[i][dim];
            }
            rn[dim] += tric_sh;
            if (rn[dim] > 0)
            {
                r2 += rn[dim]*rn[dim]*skew_fac2_d;
                /* Take care of coupling of the distances
                 * to the planes along dim0 and dim1 through dim2.
                 */
                if (dim_ind == 1 && zonei == 1)
                {
                    r2 -= rn[dim0]*rn[dim]*skew_fac_01;
                }
            }
            if (bDistMB_pulse)
            {
                clear_rvec(rb);
                rb[dim] += cg_cm[cg][dim] - c->bc[dim_ind] + tric_sh;
                if (rb[dim] > 0)
                {
                    rb2 += rb[dim]*rb[dim]*skew_fac2_d;
                    /* Take care of coupling of the distances
                     * to the planes along dim0 and dim1 through dim2.
                     */
                    if (dim_ind == 1 && zonei == 1)
                    {
                        rb2 -= rb[dim0]*rb[dim]*skew_fac_01;
                    }
                }
            }
        }

        if (r2 < r_comm2 ||
            (bDistBonded &&
             ((bDistMB && rb2 < r_bcomm2) ||
              (bDist2B && r2  < r_bcomm2)) &&
             (!bBondComm ||
              (GET_CGINFO_BOND_INTER(cginfo[cg]) &&
               missing_link(comm->cglink, globalAtomGroupIndices[cg],
                            comm->bLocalCG)))))
        {
            /* Store the local and global atom group indices and position */
            localAtomGroups->push_back(cg);
            ibuf.push_back(globalAtomGroupIndices[cg]);
            nsend_z++;

            rvec posPbc;
            if (dd->ci[dim] == 0)
            {
                /* Correct cg_cm for pbc */
                rvec_add(cg_cm[cg], box[dim], posPbc);
                if (bScrew)
                {
                    posPbc[YY] = box[YY][YY] - posPbc[YY];
                    posPbc[ZZ] = box[ZZ][ZZ] - posPbc[ZZ];
                }
            }
            else
            {
                copy_rvec(cg_cm[cg], posPbc);
            }
            vbuf.emplace_back(posPbc[XX], posPbc[YY], posPbc[ZZ]);

            nat += atomGroups.block(cg).size();
        }
    }

    work->nat        = nat;
    work->nsend_zone = nsend_z;
}

static void clearCommSetupData(dd_comm_setup_work_t *work)
{
    work->localAtomGroupBuffer.clear();
    work->atomGroupBuffer.clear();
    work->positionBuffer.clear();
    work->nat        = 0;
    work->nsend_zone = 0;
}

static void setup_dd_communication(gmx_domdec_t *dd,
                                   matrix box, gmx_ddbox_t *ddbox,
                                   t_forcerec *fr,
                                   t_state *state, PaddedRVecVector *f)
{
    int                    dim_ind, dim, dim0, dim1, dim2, dimd, nat_tot;
    int                    nzone, nzone_send, zone, zonei, cg0, cg1;
    int                    c, i, cg, cg_gl, nrcg;
    int                   *zone_cg_range, pos_cg;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_zones_t    *zones;
    gmx_domdec_comm_dim_t *cd;
    cginfo_mb_t           *cginfo_mb;
    gmx_bool               bBondComm, bDist2B, bDistMB, bDistBonded;
    real                   r_comm2, r_bcomm2;
    dd_corners_t           corners;
    rvec                  *cg_cm, *normal, *v_d, *v_0 = nullptr, *v_1 = nullptr;
    real                   skew_fac2_d, skew_fac_01;
    rvec                   sf2_round;

    if (debug)
    {
        fprintf(debug, "Setting up DD communication\n");
    }

    comm  = dd->comm;

    if (comm->dth.empty())
    {
        /* Initialize the thread data.
         * This can not be done in init_domain_decomposition,
         * as the numbers of threads is determined later.
         */
        int numThreads = gmx_omp_nthreads_get(emntDomdec);
        comm->dth.resize(numThreads);
    }

    switch (fr->cutoff_scheme)
    {
        case ecutsGROUP:
            cg_cm = fr->cg_cm;
            break;
        case ecutsVERLET:
            cg_cm = as_rvec_array(state->x.data());
            break;
        default:
            gmx_incons("unimplemented");
    }

    bBondComm = comm->bBondComm;

    /* Do we need to determine extra distances for multi-body bondeds? */
    bDistMB = (comm->bInterCGMultiBody && isDlbOn(dd->comm) && dd->ndim > 1);

    /* Do we need to determine extra distances for only two-body bondeds? */
    bDist2B = (bBondComm && !bDistMB);

    r_comm2  = gmx::square(comm->cutoff);
    r_bcomm2 = gmx::square(comm->cutoff_mbody);

    if (debug)
    {
        fprintf(debug, "bBondComm %d, r_bc %f\n", bBondComm, std::sqrt(r_bcomm2));
    }

    zones = &comm->zones;

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
            skew_fac_01 =
                ddbox->v[dim0][dim1+1][dim0]*ddbox->v[dim1][dim1+1][dim1];
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

    zone_cg_range = zones->cg_range;
    cginfo_mb     = fr->cginfo_mb;

    zone_cg_range[0]   = 0;
    zone_cg_range[1]   = dd->ncg_home;
    comm->zone_ncg1[0] = dd->ncg_home;
    pos_cg             = dd->ncg_home;

    nat_tot = comm->atomRanges.numHomeAtoms();
    nzone   = 1;
    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        dim = dd->dim[dim_ind];
        cd  = &comm->cd[dim_ind];

        /* Check if we need to compute triclinic distances along this dim */
        bool distanceIsTriclinic = false;
        for (i = 0; i <= dim_ind; i++)
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

        v_d                = ddbox->v[dim];
        skew_fac2_d        = gmx::square(ddbox->skew_fac[dim]);

        cd->receiveInPlace = true;
        for (int p = 0; p < cd->numPulses(); p++)
        {
            /* Only atoms communicated in the first pulse are used
             * for multi-body bonded interactions or for bBondComm.
             */
            bDistBonded = ((bDistMB || bDist2B) && p == 0);

            gmx_domdec_ind_t *ind = &cd->ind[p];

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
                            for (i = dd->dim[dimd]+1; i < DIM; i++)
                            {
                                /* If we are shifted in dimension i
                                 * and the cell plane is tilted forward
                                 * in dimension i, skip this coupling.
                                 */
                                if (!(zones->shift[nzone+zone][i] &&
                                      ddbox->v[dimd][i][dimd] >= 0))
                                {
                                    sf2_round[dimd] +=
                                        gmx::square(ddbox->v[dimd][i][dimd]);
                                }
                            }
                            sf2_round[dimd] = 1/sf2_round[dimd];
                        }
                    }
                }

                zonei = zone_perm[dim_ind][zone];
                if (p == 0)
                {
                    /* Here we permutate the zones to obtain a convenient order
                     * for neighbor searching
                     */
                    cg0 = zone_cg_range[zonei];
                    cg1 = zone_cg_range[zonei+1];
                }
                else
                {
                    /* Look only at the cg's received in the previous grid pulse
                     */
                    cg1 = zone_cg_range[nzone+zone+1];
                    cg0 = cg1 - cd->ind[p-1].nrecv[zone];
                }

                const int numThreads = static_cast<int>(comm->dth.size());
#pragma omp parallel for num_threads(numThreads) schedule(static)
                for (int th = 0; th < numThreads; th++)
                {
                    try
                    {
                        dd_comm_setup_work_t &work = comm->dth[th];

                        /* Retain data accumulated into buffers of thread 0 */
                        if (th > 0)
                        {
                            clearCommSetupData(&work);
                        }

                        int cg0_th = cg0 + ((cg1 - cg0)* th   )/numThreads;
                        int cg1_th = cg0 + ((cg1 - cg0)*(th+1))/numThreads;

                        /* Get the cg's for this pulse in this zone */
                        get_zone_pulse_cgs(dd, zonei, zone, cg0_th, cg1_th,
                                           dd->globalAtomGroupIndices,
                                           dd->atomGrouping(),
                                           dim, dim_ind, dim0, dim1, dim2,
                                           r_comm2, r_bcomm2,
                                           box, distanceIsTriclinic,
                                           normal, skew_fac2_d, skew_fac_01,
                                           v_d, v_0, v_1, &corners, sf2_round,
                                           bDistBonded, bBondComm,
                                           bDist2B, bDistMB,
                                           cg_cm, fr->cginfo,
                                           th == 0 ? &ind->index : &work.localAtomGroupBuffer,
                                           &work);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                } // END

                std::vector<int>       &atomGroups = comm->dth[0].atomGroupBuffer;
                std::vector<gmx::RVec> &positions  = comm->dth[0].positionBuffer;
                ind->nsend[zone]  = comm->dth[0].nsend_zone;
                /* Append data of threads>=1 to the communication buffers */
                for (int th = 1; th < numThreads; th++)
                {
                    const dd_comm_setup_work_t &dth = comm->dth[th];

                    ind->index.insert(ind->index.end(), dth.localAtomGroupBuffer.begin(), dth.localAtomGroupBuffer.end());
                    atomGroups.insert(atomGroups.end(), dth.atomGroupBuffer.begin(), dth.atomGroupBuffer.end());
                    positions.insert(positions.end(), dth.positionBuffer.begin(), dth.positionBuffer.end());
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
            /* Communicate the number of cg's and atoms to receive */
            ddSendrecv(dd, dim_ind, dddirBackward,
                       ind->nsend, nzone+2,
                       ind->nrecv, nzone+2);

            if (p > 0)
            {
                /* We can receive in place if only the last zone is not empty */
                for (zone = 0; zone < nzone-1; zone++)
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
            DDBufferAccess<int>       globalAtomGroupBuffer(comm->intBuffer, receiveBufferSize);
            DDBufferAccess<gmx::RVec> rvecBuffer(comm->rvecBuffer, receiveBufferSize);

            dd_comm_setup_work_t     &work = comm->dth[0];

            /* Make space for the global cg indices */
            int numAtomGroupsNew = pos_cg + ind->nrecv[nzone];
            dd->globalAtomGroupIndices.resize(numAtomGroupsNew);
            /* Communicate the global cg indices */
            gmx::ArrayRef<int> integerBufferRef;
            if (cd->receiveInPlace)
            {
                integerBufferRef = gmx::arrayRefFromArray(dd->globalAtomGroupIndices.data() + pos_cg, ind->nrecv[nzone]);
            }
            else
            {
                integerBufferRef = globalAtomGroupBuffer.buffer;
            }
            ddSendrecv<int>(dd, dim_ind, dddirBackward,
                            work.atomGroupBuffer, integerBufferRef);

            /* Make space for cg_cm */
            dd_check_alloc_ncg(fr, state, f, pos_cg + ind->nrecv[nzone]);
            if (fr->cutoff_scheme == ecutsGROUP)
            {
                cg_cm = fr->cg_cm;
            }
            else
            {
                cg_cm = as_rvec_array(state->x.data());
            }
            /* Communicate cg_cm */
            gmx::ArrayRef<gmx::RVec> rvecBufferRef;
            if (cd->receiveInPlace)
            {
                rvecBufferRef = gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(cg_cm + pos_cg), ind->nrecv[nzone]);
            }
            else
            {
                rvecBufferRef = rvecBuffer.buffer;
            }
            ddSendrecv<gmx::RVec>(dd, dim_ind, dddirBackward,
                                  work.positionBuffer, rvecBufferRef);

            /* Make the charge group index */
            if (cd->receiveInPlace)
            {
                zone = (p == 0 ? 0 : nzone - 1);
                while (zone < nzone)
                {
                    for (cg = 0; cg < ind->nrecv[zone]; cg++)
                    {
                        cg_gl                              = dd->globalAtomGroupIndices[pos_cg];
                        fr->cginfo[pos_cg]                 = ddcginfo(cginfo_mb, cg_gl);
                        nrcg                               = GET_CGINFO_NATOMS(fr->cginfo[pos_cg]);
                        dd->atomGrouping_.appendBlock(nrcg);
                        if (bBondComm)
                        {
                            /* Update the charge group presence,
                             * so we can use it in the next pass of the loop.
                             */
                            comm->bLocalCG[cg_gl] = TRUE;
                        }
                        pos_cg++;
                    }
                    if (p == 0)
                    {
                        comm->zone_ncg1[nzone+zone] = ind->nrecv[zone];
                    }
                    zone++;
                    zone_cg_range[nzone+zone] = pos_cg;
                }
            }
            else
            {
                /* This part of the code is never executed with bBondComm. */
                std::vector<int> &atomGroupsIndex = dd->atomGrouping_.rawIndex();
                atomGroupsIndex.resize(numAtomGroupsNew + 1);

                merge_cg_buffers(nzone, cd, p, zone_cg_range,
                                 dd->globalAtomGroupIndices, integerBufferRef.data(),
                                 cg_cm, as_rvec_array(rvecBufferRef.data()),
                                 atomGroupsIndex,
                                 fr->cginfo_mb, fr->cginfo);
                pos_cg += ind->nrecv[nzone];
            }
            nat_tot += ind->nrecv[nzone+1];
        }
        if (!cd->receiveInPlace)
        {
            /* Store the atom block for easy copying of communication buffers */
            make_cell2at_index(cd, nzone, zone_cg_range[nzone], dd->atomGrouping());
        }
        nzone += nzone;
    }

    comm->atomRanges.setEnd(DDAtomRanges::Type::Zones, nat_tot);

    if (!bBondComm)
    {
        /* We don't need to update cginfo, since that was alrady done above.
         * So we pass NULL for the forcerec.
         */
        dd_set_cginfo(dd->globalAtomGroupIndices,
                      dd->ncg_home, dd->globalAtomGroupIndices.size(),
                      nullptr, comm->bLocalCG);
    }

    if (debug)
    {
        fprintf(debug, "Finished setting up DD communication, zones:");
        for (c = 0; c < zones->n; c++)
        {
            fprintf(debug, " %d", zones->cg_range[c+1]-zones->cg_range[c]);
        }
        fprintf(debug, "\n");
    }
}

static void set_cg_boundaries(gmx_domdec_zones_t *zones)
{
    int c;

    for (c = 0; c < zones->nizone; c++)
    {
        zones->izone[c].cg1  = zones->cg_range[c+1];
        zones->izone[c].jcg0 = zones->cg_range[zones->izone[c].j0];
        zones->izone[c].jcg1 = zones->cg_range[zones->izone[c].j1];
    }
}

/* \brief Set zone dimensions for zones \p zone_start to \p zone_end-1
 *
 * Also sets the atom density for the home zone when \p zone_start=0.
 * For this \p numMovedChargeGroupsInHomeZone needs to be passed to tell
 * how many charge groups will move but are still part of the current range.
 * \todo When converting domdec to use proper classes, all these variables
 *       should be private and a method should return the correct count
 *       depending on an internal state.
 *
 * \param[in,out] dd          The domain decomposition struct
 * \param[in]     box         The box
 * \param[in]     ddbox       The domain decomposition box struct
 * \param[in]     zone_start  The start of the zone range to set sizes for
 * \param[in]     zone_end    The end of the zone range to set sizes for
 * \param[in]     numMovedChargeGroupsInHomeZone  The number of charge groups in the home zone that should moved but are still present in dd->comm->zones.cg_range
 */
static void set_zones_size(gmx_domdec_t *dd,
                           matrix box, const gmx_ddbox_t *ddbox,
                           int zone_start, int zone_end,
                           int numMovedChargeGroupsInHomeZone)
{
    gmx_domdec_comm_t  *comm;
    gmx_domdec_zones_t *zones;
    gmx_bool            bDistMB;
    int                 z, zi, d, dim;
    real                rcs, rcmbs;
    int                 i, j;
    real                vol;

    comm = dd->comm;

    zones = &comm->zones;

    /* Do we need to determine extra distances for multi-body bondeds? */
    bDistMB = (comm->bInterCGMultiBody && isDlbOn(dd->comm) && dd->ndim > 1);

    for (z = zone_start; z < zone_end; z++)
    {
        /* Copy cell limits to zone limits.
         * Valid for non-DD dims and non-shifted dims.
         */
        copy_rvec(comm->cell_x0, zones->size[z].x0);
        copy_rvec(comm->cell_x1, zones->size[z].x1);
    }

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];

        for (z = 0; z < zones->n; z++)
        {
            /* With a staggered grid we have different sizes
             * for non-shifted dimensions.
             */
            if (isDlbOn(dd->comm) && zones->shift[z][dim] == 0)
            {
                if (d == 1)
                {
                    zones->size[z].x0[dim] = comm->zone_d1[zones->shift[z][dd->dim[d-1]]].min0;
                    zones->size[z].x1[dim] = comm->zone_d1[zones->shift[z][dd->dim[d-1]]].max1;
                }
                else if (d == 2)
                {
                    zones->size[z].x0[dim] = comm->zone_d2[zones->shift[z][dd->dim[d-2]]][zones->shift[z][dd->dim[d-1]]].min0;
                    zones->size[z].x1[dim] = comm->zone_d2[zones->shift[z][dd->dim[d-2]]][zones->shift[z][dd->dim[d-1]]].max1;
                }
            }
        }

        rcs   = comm->cutoff;
        rcmbs = comm->cutoff_mbody;
        if (ddbox->tric_dir[dim])
        {
            rcs   /= ddbox->skew_fac[dim];
            rcmbs /= ddbox->skew_fac[dim];
        }

        /* Set the lower limit for the shifted zone dimensions */
        for (z = zone_start; z < zone_end; z++)
        {
            if (zones->shift[z][dim] > 0)
            {
                dim = dd->dim[d];
                if (!isDlbOn(dd->comm) || d == 0)
                {
                    zones->size[z].x0[dim] = comm->cell_x1[dim];
                    zones->size[z].x1[dim] = comm->cell_x1[dim] + rcs;
                }
                else
                {
                    /* Here we take the lower limit of the zone from
                     * the lowest domain of the zone below.
                     */
                    if (z < 4)
                    {
                        zones->size[z].x0[dim] =
                            comm->zone_d1[zones->shift[z][dd->dim[d-1]]].min1;
                    }
                    else
                    {
                        if (d == 1)
                        {
                            zones->size[z].x0[dim] =
                                zones->size[zone_perm[2][z-4]].x0[dim];
                        }
                        else
                        {
                            zones->size[z].x0[dim] =
                                comm->zone_d2[zones->shift[z][dd->dim[d-2]]][zones->shift[z][dd->dim[d-1]]].min1;
                        }
                    }
                    /* A temporary limit, is updated below */
                    zones->size[z].x1[dim] = zones->size[z].x0[dim];

                    if (bDistMB)
                    {
                        for (zi = 0; zi < zones->nizone; zi++)
                        {
                            if (zones->shift[zi][dim] == 0)
                            {
                                /* This takes the whole zone into account.
                                 * With multiple pulses this will lead
                                 * to a larger zone then strictly necessary.
                                 */
                                zones->size[z].x1[dim] = std::max(zones->size[z].x1[dim],
                                                                  zones->size[zi].x1[dim]+rcmbs);
                            }
                        }
                    }
                }
            }
        }

        /* Loop over the i-zones to set the upper limit of each
         * j-zone they see.
         */
        for (zi = 0; zi < zones->nizone; zi++)
        {
            if (zones->shift[zi][dim] == 0)
            {
                /* We should only use zones up to zone_end */
                int jZoneEnd = std::min(zones->izone[zi].j1, zone_end);
                for (z = zones->izone[zi].j0; z < jZoneEnd; z++)
                {
                    if (zones->shift[z][dim] > 0)
                    {
                        zones->size[z].x1[dim] = std::max(zones->size[z].x1[dim],
                                                          zones->size[zi].x1[dim]+rcs);
                    }
                }
            }
        }
    }

    for (z = zone_start; z < zone_end; z++)
    {
        /* Initialization only required to keep the compiler happy */
        rvec corner_min = {0, 0, 0}, corner_max = {0, 0, 0}, corner;
        int  nc, c;

        /* To determine the bounding box for a zone we need to find
         * the extreme corners of 4, 2 or 1 corners.
         */
        nc = 1 << (ddbox->nboundeddim - 1);

        for (c = 0; c < nc; c++)
        {
            /* Set up a zone corner at x=0, ignoring trilinic couplings */
            corner[XX] = 0;
            if ((c & 1) == 0)
            {
                corner[YY] = zones->size[z].x0[YY];
            }
            else
            {
                corner[YY] = zones->size[z].x1[YY];
            }
            if ((c & 2) == 0)
            {
                corner[ZZ] = zones->size[z].x0[ZZ];
            }
            else
            {
                corner[ZZ] = zones->size[z].x1[ZZ];
            }
            if (dd->ndim == 1 && dd->dim[0] < ZZ && ZZ < dd->npbcdim &&
                box[ZZ][1 - dd->dim[0]] != 0)
            {
                /* With 1D domain decomposition the cg's are not in
                 * the triclinic box, but triclinic x-y and rectangular y/x-z.
                 * Shift the corner of the z-vector back to along the box
                 * vector of dimension d, so it will later end up at 0 along d.
                 * This can affect the location of this corner along dd->dim[0]
                 * through the matrix operation below if box[d][dd->dim[0]]!=0.
                 */
                int d = 1 - dd->dim[0];

                corner[d] -= corner[ZZ]*box[ZZ][d]/box[ZZ][ZZ];
            }
            /* Apply the triclinic couplings */
            assert(ddbox->npbcdim <= DIM);
            for (i = YY; i < ddbox->npbcdim; i++)
            {
                for (j = XX; j < i; j++)
                {
                    corner[j] += corner[i]*box[i][j]/box[i][i];
                }
            }
            if (c == 0)
            {
                copy_rvec(corner, corner_min);
                copy_rvec(corner, corner_max);
            }
            else
            {
                for (i = 0; i < DIM; i++)
                {
                    corner_min[i] = std::min(corner_min[i], corner[i]);
                    corner_max[i] = std::max(corner_max[i], corner[i]);
                }
            }
        }
        /* Copy the extreme cornes without offset along x */
        for (i = 0; i < DIM; i++)
        {
            zones->size[z].bb_x0[i] = corner_min[i];
            zones->size[z].bb_x1[i] = corner_max[i];
        }
        /* Add the offset along x */
        zones->size[z].bb_x0[XX] += zones->size[z].x0[XX];
        zones->size[z].bb_x1[XX] += zones->size[z].x1[XX];
    }

    if (zone_start == 0)
    {
        vol = 1;
        for (dim = 0; dim < DIM; dim++)
        {
            vol *= zones->size[0].x1[dim] - zones->size[0].x0[dim];
        }
        zones->dens_zone0 = (zones->cg_range[1] - zones->cg_range[0] - numMovedChargeGroupsInHomeZone)/vol;
    }

    if (debug)
    {
        for (z = zone_start; z < zone_end; z++)
        {
            fprintf(debug, "zone %d    %6.3f - %6.3f  %6.3f - %6.3f  %6.3f - %6.3f\n",
                    z,
                    zones->size[z].x0[XX], zones->size[z].x1[XX],
                    zones->size[z].x0[YY], zones->size[z].x1[YY],
                    zones->size[z].x0[ZZ], zones->size[z].x1[ZZ]);
            fprintf(debug, "zone %d bb %6.3f - %6.3f  %6.3f - %6.3f  %6.3f - %6.3f\n",
                    z,
                    zones->size[z].bb_x0[XX], zones->size[z].bb_x1[XX],
                    zones->size[z].bb_x0[YY], zones->size[z].bb_x1[YY],
                    zones->size[z].bb_x0[ZZ], zones->size[z].bb_x1[ZZ]);
        }
    }
}

static int comp_cgsort(const void *a, const void *b)
{
    int           comp;

    gmx_cgsort_t *cga, *cgb;
    cga = (gmx_cgsort_t *)a;
    cgb = (gmx_cgsort_t *)b;

    comp = cga->nsc - cgb->nsc;
    if (comp == 0)
    {
        comp = cga->ind_gl - cgb->ind_gl;
    }

    return comp;
}

/* Order data in \p dataToSort according to \p sort
 *
 * Note: both buffers should have at least \p sort.size() elements.
 */
template <typename T>
static void
orderVector(gmx::ArrayRef<const gmx_cgsort_t>  sort,
            gmx::ArrayRef<T>                   dataToSort,
            gmx::ArrayRef<T>                   sortBuffer)
{
    GMX_ASSERT(dataToSort.size() >= sort.size(), "The vector needs to be sufficiently large");
    GMX_ASSERT(sortBuffer.size() >= sort.size(), "The sorting buffer needs to be sufficiently large");

    /* Order the data into the temporary buffer */
    size_t i = 0;
    for (const gmx_cgsort_t &entry : sort)
    {
        sortBuffer[i++] = dataToSort[entry.ind];
    }

    /* Copy back to the original array */
    std::copy(sortBuffer.begin(), sortBuffer.begin() + sort.size(),
              dataToSort.begin());
}

/* Order data in \p dataToSort according to \p sort
 *
 * Note: \p vectorToSort should have at least \p sort.size() elements,
 *       \p workVector is resized when it is too small.
 */
template <typename T>
static void
orderVector(gmx::ArrayRef<const gmx_cgsort_t>  sort,
            gmx::ArrayRef<T>                   vectorToSort,
            std::vector<T>                    *workVector)
{
    if (gmx::index(workVector->size()) < sort.size())
    {
        workVector->resize(sort.size());
    }
    orderVector<T>(sort, vectorToSort, *workVector);
}

static void order_vec_atom(const gmx::RangePartitioning      *atomGroups,
                           gmx::ArrayRef<const gmx_cgsort_t>  sort,
                           gmx::ArrayRef<gmx::RVec>           v,
                           gmx::ArrayRef<gmx::RVec>           buf)
{
    if (atomGroups == nullptr)
    {
        /* Avoid the useless loop of the atoms within a cg */
        orderVector(sort, v, buf);

        return;
    }

    /* Order the data */
    int a = 0;
    for (const gmx_cgsort_t &entry : sort)
    {
        for (int i : atomGroups->block(entry.ind))
        {
            copy_rvec(v[i], buf[a]);
            a++;
        }
    }
    int atot = a;

    /* Copy back to the original array */
    for (int a = 0; a < atot; a++)
    {
        copy_rvec(buf[a], v[a]);
    }
}

/* Returns whether a < b */
static bool compareCgsort(const gmx_cgsort_t &a,
                          const gmx_cgsort_t &b)
{
    return (a.nsc < b.nsc ||
            (a.nsc == b.nsc && a.ind_gl < b.ind_gl));
}

static void orderedSort(gmx::ArrayRef<const gmx_cgsort_t>  stationary,
                        gmx::ArrayRef<gmx_cgsort_t>        moved,
                        std::vector<gmx_cgsort_t>         *sort1)
{
    /* The new indices are not very ordered, so we qsort them */
    gmx_qsort_threadsafe(moved.data(), moved.size(), sizeof(moved[0]), comp_cgsort);

    /* stationary is already ordered, so now we can merge the two arrays */
    sort1->resize(stationary.size() + moved.size());
    std::merge(stationary.begin(), stationary.end(),
               moved.begin(), moved.end(),
               sort1->begin(),
               compareCgsort);
}

/* Set the sorting order for systems with charge groups, returned in sort->sort.
 * The order is according to the global charge group index.
 * This adds and removes charge groups that moved between domains.
 */
static void dd_sort_order(const gmx_domdec_t *dd,
                          const t_forcerec   *fr,
                          int                 ncg_home_old,
                          gmx_domdec_sort_t  *sort)
{
    const int *a          = fr->ns->grid->cell_index;

    const int  movedValue = NSGRID_SIGNAL_MOVED_FAC*fr->ns->grid->ncells;

    if (ncg_home_old >= 0)
    {
        std::vector<gmx_cgsort_t> &stationary = sort->stationary;
        std::vector<gmx_cgsort_t> &moved      = sort->moved;

        /* The charge groups that remained in the same ns grid cell
         * are completely ordered. So we can sort efficiently by sorting
         * the charge groups that did move into the stationary list.
         * Note: push_back() seems to be slightly slower than direct access.
         */
        stationary.clear();
        moved.clear();
        for (int i = 0; i < dd->ncg_home; i++)
        {
            /* Check if this cg did not move to another node */
            if (a[i] < movedValue)
            {
                gmx_cgsort_t entry;
                entry.nsc    = a[i];
                entry.ind_gl = dd->globalAtomGroupIndices[i];
                entry.ind    = i;

                if (i >= ncg_home_old || a[i] != sort->sorted[i].nsc)
                {
                    /* This cg is new on this node or moved ns grid cell */
                    moved.push_back(entry);
                }
                else
                {
                    /* This cg did not move */
                    stationary.push_back(entry);
                }
            }
        }

        if (debug)
        {
            fprintf(debug, "ordered sort cgs: stationary %zu moved %zu\n",
                    stationary.size(), moved.size());
        }
        /* Sort efficiently */
        orderedSort(stationary, moved, &sort->sorted);
    }
    else
    {
        std::vector<gmx_cgsort_t> &cgsort   = sort->sorted;
        cgsort.clear();
        cgsort.reserve(dd->ncg_home);
        int                        numCGNew = 0;
        for (int i = 0; i < dd->ncg_home; i++)
        {
            /* Sort on the ns grid cell indices
             * and the global topology index
             */
            gmx_cgsort_t entry;
            entry.nsc    = a[i];
            entry.ind_gl = dd->globalAtomGroupIndices[i];
            entry.ind    = i;
            cgsort.push_back(entry);
            if (cgsort[i].nsc < movedValue)
            {
                numCGNew++;
            }
        }
        if (debug)
        {
            fprintf(debug, "qsort cgs: %d new home %d\n", dd->ncg_home, numCGNew);
        }
        /* Determine the order of the charge groups using qsort */
        gmx_qsort_threadsafe(cgsort.data(), dd->ncg_home, sizeof(cgsort[0]), comp_cgsort);

        /* Remove the charge groups which are no longer at home here */
        cgsort.resize(numCGNew);
    }
}

/* Returns the sorting order for atoms based on the nbnxn grid order in sort */
static void dd_sort_order_nbnxn(const t_forcerec          *fr,
                                std::vector<gmx_cgsort_t> *sort)
{
    gmx::ArrayRef<const int> atomOrder = nbnxn_get_atomorder(fr->nbv->nbs.get());

    /* Using push_back() instead of this resize results in much slower code */
    sort->resize(atomOrder.size());
    gmx::ArrayRef<gmx_cgsort_t> buffer    = *sort;
    size_t                      numSorted = 0;
    for (int i : atomOrder)
    {
        if (i >= 0)
        {
            /* The values of nsc and ind_gl are not used in this case */
            buffer[numSorted++].ind = i;
        }
    }
    sort->resize(numSorted);
}

static void dd_sort_state(gmx_domdec_t *dd, rvec *cgcm, t_forcerec *fr, t_state *state,
                          int ncg_home_old)
{
    gmx_domdec_sort_t *sort = dd->comm->sort.get();

    switch (fr->cutoff_scheme)
    {
        case ecutsGROUP:
            dd_sort_order(dd, fr, ncg_home_old, sort);
            break;
        case ecutsVERLET:
            dd_sort_order_nbnxn(fr, &sort->sorted);
            break;
        default:
            gmx_incons("unimplemented");
    }

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    /* We alloc with the old size, since cgindex is still old */
    GMX_ASSERT(atomGrouping.numBlocks() == dd->ncg_home, "atomGroups and dd should be consistent");
    DDBufferAccess<gmx::RVec>     rvecBuffer(dd->comm->rvecBuffer, atomGrouping.fullRange().end());

    const gmx::RangePartitioning *atomGroupsPtr = (dd->comm->bCGs ? &atomGrouping : nullptr);

    /* Set the new home atom/charge group count */
    dd->ncg_home = sort->sorted.size();
    if (debug)
    {
        fprintf(debug, "Set the new home charge group count to %d\n",
                dd->ncg_home);
    }

    /* Reorder the state */
    gmx::ArrayRef<const gmx_cgsort_t> cgsort = sort->sorted;
    GMX_RELEASE_ASSERT(cgsort.size() == dd->ncg_home, "We should sort all the home atom groups");

    if (state->flags & (1 << estX))
    {
        order_vec_atom(atomGroupsPtr, cgsort, state->x, rvecBuffer.buffer);
    }
    if (state->flags & (1 << estV))
    {
        order_vec_atom(atomGroupsPtr, cgsort, state->v, rvecBuffer.buffer);
    }
    if (state->flags & (1 << estCGP))
    {
        order_vec_atom(atomGroupsPtr, cgsort, state->cg_p, rvecBuffer.buffer);
    }

    if (fr->cutoff_scheme == ecutsGROUP)
    {
        /* Reorder cgcm */
        gmx::ArrayRef<gmx::RVec> cgcmRef = gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(cgcm[0]), cgsort.size());
        orderVector(cgsort, cgcmRef, rvecBuffer.buffer);
    }

    /* Reorder the global cg index */
    orderVector<int>(cgsort, dd->globalAtomGroupIndices, &sort->intBuffer);
    /* Reorder the cginfo */
    orderVector<int>(cgsort, gmx::arrayRefFromArray(fr->cginfo, cgsort.size()), &sort->intBuffer);
    /* Rebuild the local cg index */
    if (dd->comm->bCGs)
    {
        /* We make a new, ordered atomGroups object and assign it to
         * the old one. This causes some allocation overhead, but saves
         * a copy back of the whole index.
         */
        gmx::RangePartitioning ordered;
        for (const gmx_cgsort_t &entry : cgsort)
        {
            ordered.appendBlock(atomGrouping.block(entry.ind).size());
        }
        dd->atomGrouping_ = ordered;
    }
    else
    {
        dd->atomGrouping_.setAllBlocksSizeOne(dd->ncg_home);
    }
    /* Set the home atom number */
    dd->comm->atomRanges.setEnd(DDAtomRanges::Type::Home, dd->atomGrouping().fullRange().end());

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        /* The atoms are now exactly in grid order, update the grid order */
        nbnxn_set_atomorder(fr->nbv->nbs.get());
    }
    else
    {
        /* Copy the sorted ns cell indices back to the ns grid struct */
        for (gmx::index i = 0; i < cgsort.size(); i++)
        {
            fr->ns->grid->cell_index[i] = cgsort[i].nsc;
        }
        fr->ns->grid->nr = cgsort.size();
    }
}

static void add_dd_statistics(gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm = dd->comm;

    for (int i = 0; i < static_cast<int>(DDAtomRanges::Type::Number); i++)
    {
        auto range = static_cast<DDAtomRanges::Type>(i);
        comm->sum_nat[i] +=
            comm->atomRanges.end(range) - comm->atomRanges.start(range);
    }
    comm->ndecomp++;
}

void reset_dd_statistics_counters(gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm = dd->comm;

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

void print_dd_statistics(const t_commrec *cr, const t_inputrec *ir, FILE *fplog)
{
    gmx_domdec_comm_t *comm      = cr->dd->comm;

    const int          numRanges = static_cast<int>(DDAtomRanges::Type::Number);
    gmx_sumd(numRanges, comm->sum_nat, cr);

    if (fplog == nullptr)
    {
        return;
    }

    fprintf(fplog, "\n    D O M A I N   D E C O M P O S I T I O N   S T A T I S T I C S\n\n");

    for (int i = static_cast<int>(DDAtomRanges::Type::Zones); i < numRanges; i++)
    {
        auto   range = static_cast<DDAtomRanges::Type>(i);
        double av    = comm->sum_nat[i]/comm->ndecomp;
        switch (range)
        {
            case DDAtomRanges::Type::Zones:
                fprintf(fplog,
                        " av. #atoms communicated per step for force:  %d x %.1f\n",
                        2, av);
                break;
            case DDAtomRanges::Type::Vsites:
                if (cr->dd->vsite_comm)
                {
                    fprintf(fplog,
                            " av. #atoms communicated per step for vsites: %d x %.1f\n",
                            (EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD) ? 3 : 2,
                            av);
                }
                break;
            case DDAtomRanges::Type::Constraints:
                if (cr->dd->constraint_comm)
                {
                    fprintf(fplog,
                            " av. #atoms communicated per step for LINCS:  %d x %.1f\n",
                            1 + ir->nLincsIter, av);
                }
                break;
            default:
                gmx_incons(" Unknown type for DD statistics");
        }
    }
    fprintf(fplog, "\n");

    if (comm->bRecordLoad && EI_DYNAMICS(ir->eI))
    {
        print_dd_load_av(fplog, cr->dd);
    }
}

void dd_partition_system(FILE                *fplog,
                         int64_t              step,
                         const t_commrec     *cr,
                         gmx_bool             bMasterState,
                         int                  nstglobalcomm,
                         t_state             *state_global,
                         const gmx_mtop_t    *top_global,
                         const t_inputrec    *ir,
                         t_state             *state_local,
                         PaddedRVecVector    *f,
                         gmx::MDAtoms        *mdAtoms,
                         gmx_localtop_t      *top_local,
                         t_forcerec          *fr,
                         gmx_vsite_t         *vsite,
                         gmx::Constraints    *constr,
                         t_nrnb              *nrnb,
                         gmx_wallcycle       *wcycle,
                         gmx_bool             bVerbose)
{
    gmx_domdec_t      *dd;
    gmx_domdec_comm_t *comm;
    gmx_ddbox_t        ddbox = {0};
    t_block           *cgs_gl;
    int64_t            step_pcoupl;
    rvec               cell_ns_x0, cell_ns_x1;
    int                ncgindex_set, ncg_home_old = -1, ncg_moved, nat_f_novirsum;
    gmx_bool           bBoxChanged, bNStGlobalComm, bDoDLB, bCheckWhetherToTurnDlbOn, bLogLoad;
    gmx_bool           bRedist, bSortCG, bResortAll;
    ivec               ncells_old = {0, 0, 0}, ncells_new = {0, 0, 0}, np;
    real               grid_density;
    char               sbuf[22];

    wallcycle_start(wcycle, ewcDOMDEC);

    dd   = cr->dd;
    comm = dd->comm;

    // TODO if the update code becomes accessible here, use
    // upd->deform for this logic.
    bBoxChanged = (bMasterState || inputrecDeform(ir));
    if (ir->epc != epcNO)
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
        int n = ir->nstpcouple;
        if (n == 1)
        {
            step_pcoupl = step - 1;
        }
        else
        {
            step_pcoupl = ((step - 1)/n)*n + 1;
        }
        if (step_pcoupl >= comm->partition_step)
        {
            bBoxChanged = TRUE;
        }
    }

    bNStGlobalComm = (step % nstglobalcomm == 0);

    if (!isDlbOn(comm))
    {
        bDoDLB = FALSE;
    }
    else
    {
        /* Should we do dynamic load balacing this step?
         * Since it requires (possibly expensive) global communication,
         * we might want to do DLB less frequently.
         */
        if (bBoxChanged || ir->epc != epcNO)
        {
            bDoDLB = bBoxChanged;
        }
        else
        {
            bDoDLB = bNStGlobalComm;
        }
    }

    /* Check if we have recorded loads on the nodes */
    if (comm->bRecordLoad && dd_load_count(comm) > 0)
    {
        bCheckWhetherToTurnDlbOn = dd_dlb_get_should_check_whether_to_turn_dlb_on(dd);

        /* Print load every nstlog, first and last step to the log file */
        bLogLoad = ((ir->nstlog > 0 && step % ir->nstlog == 0) ||
                    comm->n_load_collect == 0 ||
                    (ir->nsteps >= 0 &&
                     (step + ir->nstlist > ir->init_step + ir->nsteps)));

        /* Avoid extra communication due to verbose screen output
         * when nstglobalcomm is set.
         */
        if (bDoDLB || bLogLoad || bCheckWhetherToTurnDlbOn ||
            (bVerbose && (ir->nstlist == 0 || nstglobalcomm <= ir->nstlist)))
        {
            get_load_distribution(dd, wcycle);
            if (DDMASTER(dd))
            {
                if (bLogLoad)
                {
                    dd_print_load(fplog, dd, step-1);
                }
                if (bVerbose)
                {
                    dd_print_load_verbose(dd);
                }
            }
            comm->n_load_collect++;

            if (isDlbOn(comm))
            {
                if (DDMASTER(dd))
                {
                    /* Add the measured cycles to the running average */
                    const float averageFactor        = 0.1f;
                    comm->cyclesPerStepDlbExpAverage =
                        (1 - averageFactor)*comm->cyclesPerStepDlbExpAverage +
                        averageFactor*comm->cycl[ddCyclStep]/comm->cycl_n[ddCyclStep];
                }
                if (comm->dlbState == edlbsOnCanTurnOff &&
                    dd->comm->n_load_have % c_checkTurnDlbOffInterval == c_checkTurnDlbOffInterval - 1)
                {
                    gmx_bool turnOffDlb;
                    if (DDMASTER(dd))
                    {
                        /* If the running averaged cycles with DLB are more
                         * than before we turned on DLB, turn off DLB.
                         * We will again run and check the cycles without DLB
                         * and we can then decide if to turn off DLB forever.
                         */
                        turnOffDlb = (comm->cyclesPerStepDlbExpAverage >
                                      comm->cyclesPerStepBeforeDLB);
                    }
                    dd_bcast(dd, sizeof(turnOffDlb), &turnOffDlb);
                    if (turnOffDlb)
                    {
                        /* To turn off DLB, we need to redistribute the atoms */
                        dd_collect_state(dd, state_local, state_global);
                        bMasterState = TRUE;
                        turn_off_dlb(fplog, cr, step);
                    }
                }
            }
            else if (bCheckWhetherToTurnDlbOn)
            {
                gmx_bool turnOffDlbForever = FALSE;
                gmx_bool turnOnDlb         = FALSE;

                /* Since the timings are node dependent, the master decides */
                if (DDMASTER(dd))
                {
                    /* If we recently turned off DLB, we want to check if
                     * performance is better without DLB. We want to do this
                     * ASAP to minimize the chance that external factors
                     * slowed down the DLB step are gone here and we
                     * incorrectly conclude that DLB was causing the slowdown.
                     * So we measure one nstlist block, no running average.
                     */
                    if (comm->haveTurnedOffDlb &&
                        comm->cycl[ddCyclStep]/comm->cycl_n[ddCyclStep] <
                        comm->cyclesPerStepDlbExpAverage)
                    {
                        /* After turning off DLB we ran nstlist steps in fewer
                         * cycles than with DLB. This likely means that DLB
                         * in not benefical, but this could be due to a one
                         * time unlucky fluctuation, so we require two such
                         * observations in close succession to turn off DLB
                         * forever.
                         */
                        if (comm->dlbSlowerPartitioningCount > 0 &&
                            dd->ddp_count < comm->dlbSlowerPartitioningCount + 10*c_checkTurnDlbOnInterval)
                        {
                            turnOffDlbForever = TRUE;
                        }
                        comm->haveTurnedOffDlb           = false;
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
                        if (cr->npmenodes > 0 &&
                            dd_pme_f_ratio(dd) > 1 - DD_PERF_LOSS_DLB_ON)
                        {
                            turnOnDlb = FALSE;
                        }
                        else
                        {
                            turnOnDlb = (dd_force_imb_perf_loss(dd) >= DD_PERF_LOSS_DLB_ON);
                        }
                    }
                }
                struct
                {
                    gmx_bool turnOffDlbForever;
                    gmx_bool turnOnDlb;
                }
                bools {
                    turnOffDlbForever, turnOnDlb
                };
                dd_bcast(dd, sizeof(bools), &bools);
                if (bools.turnOffDlbForever)
                {
                    turn_off_dlb_forever(fplog, cr, step);
                }
                else if (bools.turnOnDlb)
                {
                    turn_on_dlb(fplog, cr, step);
                    bDoDLB = TRUE;
                }
            }
        }
        comm->n_load_have++;
    }

    cgs_gl = &comm->cgs_gl;

    bRedist = FALSE;
    if (bMasterState)
    {
        /* Clear the old state */
        clearDDStateIndices(dd, 0, 0);
        ncgindex_set = 0;

        auto xGlobal = positionsFromStatePointer(state_global);

        set_ddbox(dd, true, ir,
                  DDMASTER(dd) ? state_global->box : nullptr,
                  true, xGlobal,
                  &ddbox);

        distributeState(fplog, dd, state_global, ddbox, state_local, f);

        dd_make_local_cgs(dd, &top_local->cgs);

        /* Ensure that we have space for the new distribution */
        dd_check_alloc_ncg(fr, state_local, f, dd->ncg_home);

        if (fr->cutoff_scheme == ecutsGROUP)
        {
            calc_cgcm(fplog, 0, dd->ncg_home,
                      &top_local->cgs, as_rvec_array(state_local->x.data()), fr->cg_cm);
        }

        inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());

        dd_set_cginfo(dd->globalAtomGroupIndices, 0, dd->ncg_home, fr, comm->bLocalCG);
    }
    else if (state_local->ddp_count != dd->ddp_count)
    {
        if (state_local->ddp_count > dd->ddp_count)
        {
            gmx_fatal(FARGS, "Internal inconsistency state_local->ddp_count (%d) > dd->ddp_count (%ld)", state_local->ddp_count, dd->ddp_count);
        }

        if (state_local->ddp_count_cg_gl != state_local->ddp_count)
        {
            gmx_fatal(FARGS, "Internal inconsistency state_local->ddp_count_cg_gl (%d) != state_local->ddp_count (%d)", state_local->ddp_count_cg_gl, state_local->ddp_count);
        }

        /* Clear the old state */
        clearDDStateIndices(dd, 0, 0);

        /* Restore the atom group indices from state_local */
        restoreAtomGroups(dd, cgs_gl->index, state_local);
        make_dd_indices(dd, cgs_gl->index, 0);
        ncgindex_set = dd->ncg_home;

        if (fr->cutoff_scheme == ecutsGROUP)
        {
            /* Redetermine the cg COMs */
            calc_cgcm(fplog, 0, dd->ncg_home,
                      &top_local->cgs, as_rvec_array(state_local->x.data()), fr->cg_cm);
        }

        inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());

        dd_set_cginfo(dd->globalAtomGroupIndices, 0, dd->ncg_home, fr, comm->bLocalCG);

        set_ddbox(dd, bMasterState, ir, state_local->box,
                  true, state_local->x, &ddbox);

        bRedist = isDlbOn(comm);
    }
    else
    {
        /* We have the full state, only redistribute the cgs */

        /* Clear the non-home indices */
        clearDDStateIndices(dd, dd->ncg_home, comm->atomRanges.numHomeAtoms());
        ncgindex_set = 0;

        /* Avoid global communication for dim's without pbc and -gcom */
        if (!bNStGlobalComm)
        {
            copy_rvec(comm->box0, ddbox.box0    );
            copy_rvec(comm->box_size, ddbox.box_size);
        }
        set_ddbox(dd, bMasterState, ir, state_local->box,
                  bNStGlobalComm, state_local->x, &ddbox);

        bBoxChanged = TRUE;
        bRedist     = TRUE;
    }
    /* For dim's without pbc and -gcom */
    copy_rvec(ddbox.box0, comm->box0    );
    copy_rvec(ddbox.box_size, comm->box_size);

    set_dd_cell_sizes(dd, &ddbox, dynamic_dd_box(&ddbox, ir), bMasterState, bDoDLB,
                      step, wcycle);

    if (comm->nstDDDumpGrid > 0 && step % comm->nstDDDumpGrid == 0)
    {
        write_dd_grid_pdb("dd_grid", step, dd, state_local->box, &ddbox);
    }

    /* Check if we should sort the charge groups */
    bSortCG = (bMasterState || bRedist);

    ncg_home_old = dd->ncg_home;

    /* When repartitioning we mark charge groups that will move to neighboring
     * DD cells, but we do not move them right away for performance reasons.
     * Thus we need to keep track of how many charge groups will move for
     * obtaining correct local charge group / atom counts.
     */
    ncg_moved = 0;
    if (bRedist)
    {
        wallcycle_sub_start(wcycle, ewcsDD_REDIST);

        dd_redistribute_cg(fplog, step, dd, ddbox.tric_dir,
                           state_local, f, fr,
                           !bSortCG, nrnb, &ncgindex_set, &ncg_moved);

        wallcycle_sub_stop(wcycle, ewcsDD_REDIST);
    }

    get_nsgrid_boundaries(ddbox.nboundeddim, state_local->box,
                          dd, &ddbox,
                          &comm->cell_x0, &comm->cell_x1,
                          dd->ncg_home, fr->cg_cm,
                          cell_ns_x0, cell_ns_x1, &grid_density);

    if (bBoxChanged)
    {
        comm_dd_ns_cell_sizes(dd, &ddbox, cell_ns_x0, cell_ns_x1, step);
    }

    switch (fr->cutoff_scheme)
    {
        case ecutsGROUP:
            copy_ivec(fr->ns->grid->n, ncells_old);
            grid_first(fplog, fr->ns->grid, dd, &ddbox,
                       state_local->box, cell_ns_x0, cell_ns_x1,
                       fr->rlist, grid_density);
            break;
        case ecutsVERLET:
            nbnxn_get_ncells(fr->nbv->nbs.get(), &ncells_old[XX], &ncells_old[YY]);
            break;
        default:
            gmx_incons("unimplemented");
    }
    /* We need to store tric_dir for dd_get_ns_ranges called from ns.c */
    copy_ivec(ddbox.tric_dir, comm->tric_dir);

    if (bSortCG)
    {
        wallcycle_sub_start(wcycle, ewcsDD_GRID);

        /* Sort the state on charge group position.
         * This enables exact restarts from this step.
         * It also improves performance by about 15% with larger numbers
         * of atoms per node.
         */

        /* Fill the ns grid with the home cell,
         * so we can sort with the indices.
         */
        set_zones_ncg_home(dd);

        switch (fr->cutoff_scheme)
        {
            case ecutsVERLET:
                set_zones_size(dd, state_local->box, &ddbox, 0, 1, ncg_moved);

                nbnxn_put_on_grid(fr->nbv->nbs.get(), fr->ePBC, state_local->box,
                                  0,
                                  comm->zones.size[0].bb_x0,
                                  comm->zones.size[0].bb_x1,
                                  0, dd->ncg_home,
                                  comm->zones.dens_zone0,
                                  fr->cginfo,
                                  as_rvec_array(state_local->x.data()),
                                  ncg_moved, bRedist ? comm->movedBuffer.data() : nullptr,
                                  fr->nbv->grp[eintLocal].kernel_type,
                                  fr->nbv->nbat);

                nbnxn_get_ncells(fr->nbv->nbs.get(), &ncells_new[XX], &ncells_new[YY]);
                break;
            case ecutsGROUP:
                fill_grid(&comm->zones, fr->ns->grid, dd->ncg_home,
                          0, dd->ncg_home, fr->cg_cm);

                copy_ivec(fr->ns->grid->n, ncells_new);
                break;
            default:
                gmx_incons("unimplemented");
        }

        bResortAll = bMasterState;

        /* Check if we can user the old order and ns grid cell indices
         * of the charge groups to sort the charge groups efficiently.
         */
        if (ncells_new[XX] != ncells_old[XX] ||
            ncells_new[YY] != ncells_old[YY] ||
            ncells_new[ZZ] != ncells_old[ZZ])
        {
            bResortAll = TRUE;
        }

        if (debug)
        {
            fprintf(debug, "Step %s, sorting the %d home charge groups\n",
                    gmx_step_str(step, sbuf), dd->ncg_home);
        }
        dd_sort_state(dd, fr->cg_cm, fr, state_local,
                      bResortAll ? -1 : ncg_home_old);

        /* After sorting and compacting we set the correct size */
        dd_resize_state(state_local, f, comm->atomRanges.numHomeAtoms());

        /* Rebuild all the indices */
        ga2la_clear(dd->ga2la);
        ncgindex_set = 0;

        wallcycle_sub_stop(wcycle, ewcsDD_GRID);
    }

    wallcycle_sub_start(wcycle, ewcsDD_SETUPCOMM);

    /* Setup up the communication and communicate the coordinates */
    setup_dd_communication(dd, state_local->box, &ddbox, fr, state_local, f);

    /* Set the indices */
    make_dd_indices(dd, cgs_gl->index, ncgindex_set);

    /* Set the charge group boundaries for neighbor searching */
    set_cg_boundaries(&comm->zones);

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        /* When bSortCG=true, we have already set the size for zone 0 */
        set_zones_size(dd, state_local->box, &ddbox,
                       bSortCG ? 1 : 0, comm->zones.n,
                       0);
    }

    wallcycle_sub_stop(wcycle, ewcsDD_SETUPCOMM);

    /*
       write_dd_pdb("dd_home",step,"dump",top_global,cr,
                 -1,as_rvec_array(state_local->x.data()),state_local->box);
     */

    wallcycle_sub_start(wcycle, ewcsDD_MAKETOP);

    /* Extract a local topology from the global topology */
    for (int i = 0; i < dd->ndim; i++)
    {
        np[dd->dim[i]] = comm->cd[i].numPulses();
    }
    dd_make_local_top(dd, &comm->zones, dd->npbcdim, state_local->box,
                      comm->cellsize_min, np,
                      fr,
                      fr->cutoff_scheme == ecutsGROUP ? fr->cg_cm : as_rvec_array(state_local->x.data()),
                      vsite, top_global, top_local);

    wallcycle_sub_stop(wcycle, ewcsDD_MAKETOP);

    wallcycle_sub_start(wcycle, ewcsDD_MAKECONSTR);

    /* Set up the special atom communication */
    int n = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int i = static_cast<int>(DDAtomRanges::Type::Zones) + 1; i < static_cast<int>(DDAtomRanges::Type::Number); i++)
    {
        auto range = static_cast<DDAtomRanges::Type>(i);
        switch (range)
        {
            case DDAtomRanges::Type::Vsites:
                if (vsite && vsite->n_intercg_vsite)
                {
                    n = dd_make_local_vsites(dd, n, top_local->idef.il);
                }
                break;
            case DDAtomRanges::Type::Constraints:
                if (dd->bInterCGcons || dd->bInterCGsettles)
                {
                    /* Only for inter-cg constraints we need special code */
                    n = dd_make_local_constraints(dd, n, top_global, fr->cginfo,
                                                  constr, ir->nProjOrder,
                                                  top_local->idef.il);
                }
                break;
            default:
                gmx_incons("Unknown special atom type setup");
        }
        comm->atomRanges.setEnd(range, n);
    }

    wallcycle_sub_stop(wcycle, ewcsDD_MAKECONSTR);

    wallcycle_sub_start(wcycle, ewcsDD_TOPOTHER);

    /* Make space for the extra coordinates for virtual site
     * or constraint communication.
     */
    state_local->natoms = comm->atomRanges.numAtomsTotal();

    dd_resize_state(state_local, f, state_local->natoms);

    if (fr->haveDirectVirialContributions)
    {
        if (vsite && vsite->n_intercg_vsite)
        {
            nat_f_novirsum = comm->atomRanges.end(DDAtomRanges::Type::Vsites);
        }
        else
        {
            if (EEL_FULL(ir->coulombtype) && dd->n_intercg_excl > 0)
            {
                nat_f_novirsum = comm->atomRanges.end(DDAtomRanges::Type::Zones);
            }
            else
            {
                nat_f_novirsum = comm->atomRanges.numHomeAtoms();
            }
        }
    }
    else
    {
        nat_f_novirsum = 0;
    }

    /* Set the number of atoms required for the force calculation.
     * Forces need to be constrained when doing energy
     * minimization. For simple simulations we could avoid some
     * allocation, zeroing and copying, but this is probably not worth
     * the complications and checking.
     */
    forcerec_set_ranges(fr, dd->ncg_home, dd->globalAtomGroupIndices.size(),
                        comm->atomRanges.end(DDAtomRanges::Type::Zones),
                        comm->atomRanges.end(DDAtomRanges::Type::Constraints),
                        nat_f_novirsum);

    /* Update atom data for mdatoms and several algorithms */
    mdAlgorithmsSetupAtomData(cr, ir, top_global, top_local, fr,
                              nullptr, mdAtoms, constr, vsite, nullptr);

    auto mdatoms = mdAtoms->mdatoms();
    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Send the charges and/or c6/sigmas to our PME only node */
        gmx_pme_send_parameters(cr,
                                fr->ic,
                                mdatoms->nChargePerturbed, mdatoms->nTypePerturbed,
                                mdatoms->chargeA, mdatoms->chargeB,
                                mdatoms->sqrt_c6A, mdatoms->sqrt_c6B,
                                mdatoms->sigmaA, mdatoms->sigmaB,
                                dd_pme_maxshift_x(dd), dd_pme_maxshift_y(dd));
    }

    if (ir->bPull)
    {
        /* Update the local pull groups */
        dd_make_local_pull_groups(cr, ir->pull_work, mdatoms);
    }

    if (ir->bRot)
    {
        /* Update the local rotation groups */
        dd_make_local_rotation_groups(dd, ir->rot);
    }

    if (ir->eSwapCoords != eswapNO)
    {
        /* Update the local groups needed for ion swapping */
        dd_make_local_swap_groups(dd, ir->swap);
    }

    /* Update the local atoms to be communicated via the IMD protocol if bIMD is TRUE. */
    dd_make_local_IMD_atoms(ir->bIMD, dd, ir->imd);

    add_dd_statistics(dd);

    /* Make sure we only count the cycles for this DD partitioning */
    clear_dd_cycle_counts(dd);

    /* Because the order of the atoms might have changed since
     * the last vsite construction, we need to communicate the constructing
     * atom coordinates again (for spreading the forces this MD step).
     */
    dd_move_x_vsites(dd, state_local->box, as_rvec_array(state_local->x.data()));

    wallcycle_sub_stop(wcycle, ewcsDD_TOPOTHER);

    if (comm->nstDDDump > 0 && step % comm->nstDDDump == 0)
    {
        dd_move_x(dd, state_local->box, state_local->x, nullWallcycle);
        write_dd_pdb("dd_dump", step, "dump", top_global, cr,
                     -1, as_rvec_array(state_local->x.data()), state_local->box);
    }

    /* Store the partitioning step */
    comm->partition_step = step;

    /* Increase the DD partitioning counter */
    dd->ddp_count++;
    /* The state currently matches this DD partitioning count, store it */
    state_local->ddp_count = dd->ddp_count;
    if (bMasterState)
    {
        /* The DD master node knows the complete cg distribution,
         * store the count so we can possibly skip the cg info communication.
         */
        comm->master_cg_ddp_count = (bSortCG ? 0 : dd->ddp_count);
    }

    if (comm->DD_debug > 0)
    {
        /* Set the env var GMX_DD_DEBUG if you suspect corrupted indices */
        check_index_consistency(dd, top_global->natoms, ncg_mtop(top_global),
                                "after partitioning");
    }

    wallcycle_stop(wcycle, ewcDOMDEC);
}

/*! \brief Check whether bonded interactions are missing, if appropriate */
void checkNumberOfBondedInteractions(FILE                 *fplog,
                                     t_commrec            *cr,
                                     int                   totalNumberOfBondedInteractions,
                                     const gmx_mtop_t     *top_global,
                                     const gmx_localtop_t *top_local,
                                     const t_state        *state,
                                     bool                 *shouldCheckNumberOfBondedInteractions)
{
    if (*shouldCheckNumberOfBondedInteractions)
    {
        if (totalNumberOfBondedInteractions != cr->dd->nbonded_global)
        {
            dd_print_missing_interactions(fplog, cr, totalNumberOfBondedInteractions, top_global, top_local, state); // Does not return
        }
        *shouldCheckNumberOfBondedInteractions = false;
    }
}
