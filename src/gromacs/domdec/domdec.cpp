/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

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
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/genborn.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
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

#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_vsite.h"

#define DDRANK(dd, rank)    (rank)
#define DDMASTERRANK(dd)   (dd->masterrank)

struct gmx_domdec_master_t
{
    /* The cell boundaries */
    real **cell_x;
    /* The global charge group division */
    int   *ncg;    /* Number of home charge groups for each node */
    int   *index;  /* Index of nnodes+1 into cg */
    int   *cg;     /* Global charge group index */
    int   *nat;    /* Number of home atoms for each node. */
    int   *ibuf;   /* Buffer for communication */
    rvec  *vbuf;   /* Buffer for state scattering and gathering */
};

#define DD_NLOAD_MAX 9

const char *edlbs_names[edlbsNR] = { "off", "auto", "locked", "on", "on" };

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

/* Factors used to avoid problems due to rounding issues */
#define DD_CELL_MARGIN       1.0001
#define DD_CELL_MARGIN2      1.00005
/* Factor to account for pressure scaling during nstlist steps */
#define DD_PRES_SCALE_MARGIN 1.02

/* Turn on DLB when the load imbalance causes this amount of total loss.
 * There is a bit of overhead with DLB and it's difficult to achieve
 * a load imbalance of less than 2% with DLB.
 */
#define DD_PERF_LOSS_DLB_ON  0.02

/* Warn about imbalance due to PP or PP/PME load imbalance at this loss */
#define DD_PERF_LOSS_WARN    0.05

#define DD_CELL_F_SIZE(dd, di) ((dd)->nc[(dd)->dim[(di)]]+1+(di)*2+1+(di))

/* Use separate MPI send and receive commands
 * when nnodes <= GMX_DD_NNODES_SENDRECV.
 * This saves memory (and some copying for small nnodes).
 * For high parallelization scatter and gather calls are used.
 */
#define GMX_DD_NNODES_SENDRECV 4


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

/* This order is required to minimize the coordinate communication in PME
 * which uses decomposition in the x direction.
 */
#define dd_index(n, i) ((((i)[XX]*(n)[YY] + (i)[YY])*(n)[ZZ]) + (i)[ZZ])

static void ddindex2xyz(ivec nc, int ind, ivec xyz)
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
        if (i >= dd->comm->nat[ddnatNR-1])
        {
            gmx_fatal(FARGS, "glatnr called with %d, which is larger than the local number of atoms (%d)", i, dd->comm->nat[ddnatNR-1]);
        }
        atnr = dd->gatindex[i] + 1;
    }

    return atnr;
}

t_block *dd_charge_groups_global(gmx_domdec_t *dd)
{
    return &dd->comm->cgs_gl;
}

/*! \brief Returns true if the DLB state indicates that the balancer is on. */
static bool isDlbOn(const gmx_domdec_comm_t *comm)
{
    return (comm->dlbState == edlbsOnCanTurnOff ||
            comm->dlbState == edlbsOnUser);
}
/*! \brief Returns true if the DLB state indicates that the balancer is off/disabled.
 */
static bool isDlbDisabled(const gmx_domdec_comm_t *comm)
{
    return (comm->dlbState == edlbsOffUser ||
            comm->dlbState == edlbsOffForever);
}

static void vec_rvec_init(vec_rvec_t *v)
{
    v->nalloc = 0;
    v->v      = nullptr;
}

static void vec_rvec_check_alloc(vec_rvec_t *v, int n)
{
    if (n > v->nalloc)
    {
        v->nalloc = over_alloc_dd(n);
        srenew(v->v, v->nalloc);
    }
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
        state->cg_gl[i] = dd->index_gl[i];
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

int dd_natoms_mdatoms(const gmx_domdec_t *dd)
{
    /* We currently set mdatoms entries for all atoms:
     * local + non-local + communicated for vsite + constraints
     */

    return dd->comm->nat[ddnatNR - 1];
}

int dd_natoms_vsite(const gmx_domdec_t *dd)
{
    return dd->comm->nat[ddnatVSITE];
}

void dd_get_constraint_range(const gmx_domdec_t *dd, int *at_start, int *at_end)
{
    *at_start = dd->comm->nat[ddnatCON-1];
    *at_end   = dd->comm->nat[ddnatCON];
}

void dd_move_x(gmx_domdec_t *dd, matrix box, rvec x[])
{
    int                    nzone, nat_tot, n, d, p, i, j, at0, at1, zone;
    int                   *index, *cgindex;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    gmx_domdec_ind_t      *ind;
    rvec                   shift = {0, 0, 0}, *buf, *rbuf;
    gmx_bool               bPBC, bScrew;

    comm = dd->comm;

    cgindex = dd->cgindex;

    buf = comm->vbuf.v;

    nzone   = 1;
    nat_tot = dd->nat_home;
    for (d = 0; d < dd->ndim; d++)
    {
        bPBC   = (dd->ci[dd->dim[d]] == 0);
        bScrew = (bPBC && dd->bScrewPBC && dd->dim[d] == XX);
        if (bPBC)
        {
            copy_rvec(box[dd->dim[d]], shift);
        }
        cd = &comm->cd[d];
        for (p = 0; p < cd->np; p++)
        {
            ind   = &cd->ind[p];
            index = ind->index;
            n     = 0;
            if (!bPBC)
            {
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        copy_rvec(x[j], buf[n]);
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        /* We need to shift the coordinates */
                        rvec_add(x[j], shift, buf[n]);
                        n++;
                    }
                }
            }
            else
            {
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        /* Shift x */
                        buf[n][XX] = x[j][XX] + shift[XX];
                        /* Rotate y and z.
                         * This operation requires a special shift force
                         * treatment, which is performed in calc_vir.
                         */
                        buf[n][YY] = box[YY][YY] - x[j][YY];
                        buf[n][ZZ] = box[ZZ][ZZ] - x[j][ZZ];
                        n++;
                    }
                }
            }

            if (cd->bInPlace)
            {
                rbuf = x + nat_tot;
            }
            else
            {
                rbuf = comm->vbuf2.v;
            }
            /* Send and receive the coordinates */
            dd_sendrecv_rvec(dd, d, dddirBackward,
                             buf,  ind->nsend[nzone+1],
                             rbuf, ind->nrecv[nzone+1]);
            if (!cd->bInPlace)
            {
                j = 0;
                for (zone = 0; zone < nzone; zone++)
                {
                    for (i = ind->cell2at0[zone]; i < ind->cell2at1[zone]; i++)
                    {
                        copy_rvec(rbuf[j], x[i]);
                        j++;
                    }
                }
            }
            nat_tot += ind->nrecv[nzone+1];
        }
        nzone += nzone;
    }
}

void dd_move_f(gmx_domdec_t *dd, rvec f[], rvec *fshift)
{
    int                    nzone, nat_tot, n, d, p, i, j, at0, at1, zone;
    int                   *index, *cgindex;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    gmx_domdec_ind_t      *ind;
    rvec                  *buf, *sbuf;
    ivec                   vis;
    int                    is;
    gmx_bool               bShiftForcesNeedPbc, bScrew;

    comm = dd->comm;

    cgindex = dd->cgindex;

    buf = comm->vbuf.v;

    nzone   = comm->zones.n/2;
    nat_tot = dd->nat_tot;
    for (d = dd->ndim-1; d >= 0; d--)
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
        for (p = cd->np-1; p >= 0; p--)
        {
            ind      = &cd->ind[p];
            nat_tot -= ind->nrecv[nzone+1];
            if (cd->bInPlace)
            {
                sbuf = f + nat_tot;
            }
            else
            {
                sbuf = comm->vbuf2.v;
                j    = 0;
                for (zone = 0; zone < nzone; zone++)
                {
                    for (i = ind->cell2at0[zone]; i < ind->cell2at1[zone]; i++)
                    {
                        copy_rvec(f[i], sbuf[j]);
                        j++;
                    }
                }
            }
            /* Communicate the forces */
            dd_sendrecv_rvec(dd, d, dddirForward,
                             sbuf, ind->nrecv[nzone+1],
                             buf,  ind->nsend[nzone+1]);
            index = ind->index;
            /* Add the received forces */
            n = 0;
            if (!bShiftForcesNeedPbc)
            {
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        rvec_inc(f[j], buf[n]);
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                /* fshift should always be defined if this function is
                 * called when bShiftForcesNeedPbc is true */
                assert(NULL != fshift);
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        rvec_inc(f[j], buf[n]);
                        /* Add this force to the shift force */
                        rvec_inc(fshift[is], buf[n]);
                        n++;
                    }
                }
            }
            else
            {
                for (i = 0; i < ind->nsend[nzone]; i++)
                {
                    at0 = cgindex[index[i]];
                    at1 = cgindex[index[i]+1];
                    for (j = at0; j < at1; j++)
                    {
                        /* Rotate the force */
                        f[j][XX] += buf[n][XX];
                        f[j][YY] -= buf[n][YY];
                        f[j][ZZ] -= buf[n][ZZ];
                        if (fshift)
                        {
                            /* Add this force to the shift force */
                            rvec_inc(fshift[is], buf[n]);
                        }
                        n++;
                    }
                }
            }
        }
        nzone /= 2;
    }
}

void dd_atom_spread_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot, n, d, p, i, j, at0, at1, zone;
    int                   *index, *cgindex;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    gmx_domdec_ind_t      *ind;
    real                  *buf, *rbuf;

    comm = dd->comm;

    cgindex = dd->cgindex;

    buf = &comm->vbuf.v[0][0];

    nzone   = 1;
    nat_tot = dd->nat_home;
    for (d = 0; d < dd->ndim; d++)
    {
        cd = &comm->cd[d];
        for (p = 0; p < cd->np; p++)
        {
            ind   = &cd->ind[p];
            index = ind->index;
            n     = 0;
            for (i = 0; i < ind->nsend[nzone]; i++)
            {
                at0 = cgindex[index[i]];
                at1 = cgindex[index[i]+1];
                for (j = at0; j < at1; j++)
                {
                    buf[n] = v[j];
                    n++;
                }
            }

            if (cd->bInPlace)
            {
                rbuf = v + nat_tot;
            }
            else
            {
                rbuf = &comm->vbuf2.v[0][0];
            }
            /* Send and receive the coordinates */
            dd_sendrecv_real(dd, d, dddirBackward,
                             buf,  ind->nsend[nzone+1],
                             rbuf, ind->nrecv[nzone+1]);
            if (!cd->bInPlace)
            {
                j = 0;
                for (zone = 0; zone < nzone; zone++)
                {
                    for (i = ind->cell2at0[zone]; i < ind->cell2at1[zone]; i++)
                    {
                        v[i] = rbuf[j];
                        j++;
                    }
                }
            }
            nat_tot += ind->nrecv[nzone+1];
        }
        nzone += nzone;
    }
}

void dd_atom_sum_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot, n, d, p, i, j, at0, at1, zone;
    int                   *index, *cgindex;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    gmx_domdec_ind_t      *ind;
    real                  *buf, *sbuf;

    comm = dd->comm;

    cgindex = dd->cgindex;

    buf = &comm->vbuf.v[0][0];

    nzone   = comm->zones.n/2;
    nat_tot = dd->nat_tot;
    for (d = dd->ndim-1; d >= 0; d--)
    {
        cd = &comm->cd[d];
        for (p = cd->np-1; p >= 0; p--)
        {
            ind      = &cd->ind[p];
            nat_tot -= ind->nrecv[nzone+1];
            if (cd->bInPlace)
            {
                sbuf = v + nat_tot;
            }
            else
            {
                sbuf = &comm->vbuf2.v[0][0];
                j    = 0;
                for (zone = 0; zone < nzone; zone++)
                {
                    for (i = ind->cell2at0[zone]; i < ind->cell2at1[zone]; i++)
                    {
                        sbuf[j] = v[i];
                        j++;
                    }
                }
            }
            /* Communicate the forces */
            dd_sendrecv_real(dd, d, dddirForward,
                             sbuf, ind->nrecv[nzone+1],
                             buf,  ind->nsend[nzone+1]);
            index = ind->index;
            /* Add the received forces */
            n = 0;
            for (i = 0; i < ind->nsend[nzone]; i++)
            {
                at0 = cgindex[index[i]];
                at1 = cgindex[index[i]+1];
                for (j = at0; j < at1; j++)
                {
                    v[j] += buf[n];
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


#define DDZONECOMM_MAXZONE  5
#define DDZONECOMM_BUFSIZE  3

static void dd_sendrecv_ddzone(const gmx_domdec_t *dd,
                               int ddimind, int direction,
                               gmx_ddzone_t *buf_s, int n_s,
                               gmx_ddzone_t *buf_r, int n_r)
{
#define ZBS  DDZONECOMM_BUFSIZE
    rvec vbuf_s[DDZONECOMM_MAXZONE*ZBS];
    rvec vbuf_r[DDZONECOMM_MAXZONE*ZBS];
    int  i;

    for (i = 0; i < n_s; i++)
    {
        vbuf_s[i*ZBS  ][0] = buf_s[i].min0;
        vbuf_s[i*ZBS  ][1] = buf_s[i].max1;
        vbuf_s[i*ZBS  ][2] = buf_s[i].min1;
        vbuf_s[i*ZBS+1][0] = buf_s[i].mch0;
        vbuf_s[i*ZBS+1][1] = buf_s[i].mch1;
        vbuf_s[i*ZBS+1][2] = 0;
        vbuf_s[i*ZBS+2][0] = buf_s[i].p1_0;
        vbuf_s[i*ZBS+2][1] = buf_s[i].p1_1;
        vbuf_s[i*ZBS+2][2] = 0;
    }

    dd_sendrecv_rvec(dd, ddimind, direction,
                     vbuf_s, n_s*ZBS,
                     vbuf_r, n_r*ZBS);

    for (i = 0; i < n_r; i++)
    {
        buf_r[i].min0 = vbuf_r[i*ZBS  ][0];
        buf_r[i].max1 = vbuf_r[i*ZBS  ][1];
        buf_r[i].min1 = vbuf_r[i*ZBS  ][2];
        buf_r[i].mch0 = vbuf_r[i*ZBS+1][0];
        buf_r[i].mch1 = vbuf_r[i*ZBS+1][1];
        buf_r[i].p1_0 = vbuf_r[i*ZBS+2][0];
        buf_r[i].p1_1 = vbuf_r[i*ZBS+2][1];
    }

#undef ZBS
}

static void dd_move_cellx(gmx_domdec_t *dd, gmx_ddbox_t *ddbox,
                          rvec cell_ns_x0, rvec cell_ns_x1)
{
    int                d, d1, dim, pos, buf_size, i, j, p, npulse, npulse_min;
    gmx_ddzone_t      *zp;
    gmx_ddzone_t       buf_s[DDZONECOMM_MAXZONE];
    gmx_ddzone_t       buf_r[DDZONECOMM_MAXZONE];
    gmx_ddzone_t       buf_e[DDZONECOMM_MAXZONE];
    rvec               extr_s[2], extr_r[2];
    rvec               dh;
    real               dist_d, c = 0, det;
    gmx_domdec_comm_t *comm;
    gmx_bool           bPBC, bUse;

    comm = dd->comm;

    for (d = 1; d < dd->ndim; d++)
    {
        dim      = dd->dim[d];
        zp       = (d == 1) ? &comm->zone_d1[0] : &comm->zone_d2[0][0];
        zp->min0 = cell_ns_x0[dim];
        zp->max1 = cell_ns_x1[dim];
        zp->min1 = cell_ns_x1[dim];
        zp->mch0 = cell_ns_x0[dim];
        zp->mch1 = cell_ns_x1[dim];
        zp->p1_0 = cell_ns_x0[dim];
        zp->p1_1 = cell_ns_x1[dim];
    }

    for (d = dd->ndim-2; d >= 0; d--)
    {
        dim  = dd->dim[d];
        bPBC = (dim < ddbox->npbcdim);

        /* Use an rvec to store two reals */
        extr_s[d][0] = comm->cell_f0[d+1];
        extr_s[d][1] = comm->cell_f1[d+1];
        extr_s[d][2] = comm->cell_f1[d+1];

        pos = 0;
        /* Store the extremes in the backward sending buffer,
         * so the get updated separately from the forward communication.
         */
        for (d1 = d; d1 < dd->ndim-1; d1++)
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
        npulse = comm->cd[d].np;
        if (bPBC)
        {
            /* Take the minimum to avoid double communication */
            npulse_min = std::min(npulse, dd->nc[dim]-1-npulse);
        }
        else
        {
            /* Without PBC we should really not communicate over
             * the boundaries, but implementing that complicates
             * the communication setup and therefore we simply
             * do all communication, but ignore some data.
             */
            npulse_min = npulse;
        }
        for (p = 0; p < npulse_min; p++)
        {
            /* Communicate the extremes forward */
            bUse = (bPBC || dd->ci[dim] > 0);

            dd_sendrecv_rvec(dd, d, dddirForward,
                             extr_s+d, dd->ndim-d-1,
                             extr_r+d, dd->ndim-d-1);

            if (bUse)
            {
                for (d1 = d; d1 < dd->ndim-1; d1++)
                {
                    extr_s[d1][0] = std::max(extr_s[d1][0], extr_r[d1][0]);
                    extr_s[d1][1] = std::min(extr_s[d1][1], extr_r[d1][1]);
                    extr_s[d1][2] = std::min(extr_s[d1][2], extr_r[d1][2]);
                }
            }
        }

        buf_size = pos;
        for (p = 0; p < npulse; p++)
        {
            /* Communicate all the zone information backward */
            bUse = (bPBC || dd->ci[dim] < dd->nc[dim] - 1);

            dd_sendrecv_ddzone(dd, d, dddirBackward,
                               buf_s, buf_size,
                               buf_r, buf_size);

            clear_rvec(dh);
            if (p > 0)
            {
                for (d1 = d+1; d1 < dd->ndim; d1++)
                {
                    /* Determine the decrease of maximum required
                     * communication height along d1 due to the distance along d,
                     * this avoids a lot of useless atom communication.
                     */
                    dist_d = comm->cell_x1[dim] - buf_r[0].p1_0;

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
                    det = (1 + c*c)*comm->cutoff*comm->cutoff - dist_d*dist_d;
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
            for (i = 0; i < buf_size; i++)
            {
                if (p == 0)
                {
                    buf_e[i] = buf_r[i];
                }
                else
                {
                    if (bUse)
                    {
                        buf_e[i].min0 = std::min(buf_e[i].min0, buf_r[i].min0);
                        buf_e[i].max1 = std::max(buf_e[i].max1, buf_r[i].max1);
                        buf_e[i].min1 = std::min(buf_e[i].min1, buf_r[i].min1);
                    }

                    if (dd->ndim == 3 && d == 0 && i == buf_size - 1)
                    {
                        d1 = 1;
                    }
                    else
                    {
                        d1 = d + 1;
                    }
                    if (bUse && dh[d1] >= 0)
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
            if (((bPBC || dd->ci[dim]+npulse < dd->nc[dim]) && p == npulse-1) ||
                (!bPBC && dd->ci[dim]+1+p == dd->nc[dim]-1))
            {
                /* Store the extremes */
                pos = 0;

                for (d1 = d; d1 < dd->ndim-1; d1++)
                {
                    extr_s[d1][1] = std::min(extr_s[d1][1], buf_e[pos].min0);
                    extr_s[d1][0] = std::max(extr_s[d1][0], buf_e[pos].max1);
                    extr_s[d1][2] = std::min(extr_s[d1][2], buf_e[pos].min1);
                    pos++;
                }

                if (d == 1 || (d == 0 && dd->ndim == 3))
                {
                    for (i = d; i < 2; i++)
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
        dim = dd->dim[1];
        for (i = 0; i < 2; i++)
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
        dim = dd->dim[2];
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
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
    for (d = 1; d < dd->ndim; d++)
    {
        comm->cell_f_max0[d] = extr_s[d-1][0];
        comm->cell_f_min1[d] = extr_s[d-1][1];
        if (debug)
        {
            fprintf(debug, "Cell fraction d %d, max0 %f, min1 %f\n",
                    d, comm->cell_f_max0[d], comm->cell_f_min1[d]);
        }
    }
}

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

static void dd_collect_vec_sendrecv(gmx_domdec_t *dd,
                                    const rvec *lv, rvec *v)
{
    gmx_domdec_master_t *ma;
    int                  n, i, c, a, nalloc = 0;
    rvec                *buf = nullptr;
    t_block             *cgs_gl;

    ma = dd->ma;

    if (!DDMASTER(dd))
    {
#if GMX_MPI
        MPI_Send(const_cast<void *>(static_cast<const void *>(lv)), dd->nat_home*sizeof(rvec), MPI_BYTE,
                 DDMASTERRANK(dd), dd->rank, dd->mpi_comm_all);
#endif
    }
    else
    {
        /* Copy the master coordinates to the global array */
        cgs_gl = &dd->comm->cgs_gl;

        n = DDMASTERRANK(dd);
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
                MPI_Recv(buf, ma->nat[n]*sizeof(rvec), MPI_BYTE, DDRANK(dd, n),
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

static void get_commbuffer_counts(gmx_domdec_t *dd,
                                  int **counts, int **disps)
{
    gmx_domdec_master_t *ma;
    int                  n;

    ma = dd->ma;

    /* Make the rvec count and displacment arrays */
    *counts  = ma->ibuf;
    *disps   = ma->ibuf + dd->nnodes;
    for (n = 0; n < dd->nnodes; n++)
    {
        (*counts)[n] = ma->nat[n]*sizeof(rvec);
        (*disps)[n]  = (n == 0 ? 0 : (*disps)[n-1] + (*counts)[n-1]);
    }
}

static void dd_collect_vec_gatherv(gmx_domdec_t *dd,
                                   const rvec *lv, rvec *v)
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

    dd_gatherv(dd, dd->nat_home*sizeof(rvec), lv, rcounts, disps, buf);

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

void dd_collect_vec(gmx_domdec_t           *dd,
                    const t_state          *state_local,
                    const PaddedRVecVector *localVector,
                    rvec                   *v)
{
    dd_collect_cg(dd, state_local);

    const rvec *lv = as_rvec_array(localVector->data());

    if (dd->nnodes <= GMX_DD_NNODES_SENDRECV)
    {
        dd_collect_vec_sendrecv(dd, lv, v);
    }
    else
    {
        dd_collect_vec_gatherv(dd, lv, v);
    }
}

void dd_collect_vec(gmx_domdec_t           *dd,
                    const t_state          *state_local,
                    const PaddedRVecVector *localVector,
                    PaddedRVecVector       *vector)
{
    dd_collect_vec(dd, state_local, localVector,
                   DDMASTER(dd) ? as_rvec_array(vector->data()) : nullptr);
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
        dd_collect_vec(dd, state_local, &state_local->x, &state->x);
    }
    if (state_local->flags & (1 << estV))
    {
        dd_collect_vec(dd, state_local, &state_local->v, &state->v);
    }
    if (state_local->flags & (1 << estCGP))
    {
        dd_collect_vec(dd, state_local, &state_local->cg_p, &state->cg_p);
    }
}

static void dd_resize_state(t_state *state, PaddedRVecVector *f, int natoms)
{
    if (debug)
    {
        fprintf(debug, "Resizing state: currently %d, required %d\n", state->natoms, natoms);
    }

    state_change_natoms(state, natoms);

    if (f != nullptr)
    {
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        f->resize(natoms + 1);
    }
}

static void dd_check_alloc_ncg(t_forcerec       *fr,
                               t_state          *state,
                               PaddedRVecVector *f,
                               int               numChargeGroups)
{
    if (numChargeGroups > fr->cg_nalloc)
    {
        if (debug)
        {
            fprintf(debug, "Reallocating forcerec: currently %d, required %d, allocating %d\n", fr->cg_nalloc, numChargeGroups, over_alloc_dd(numChargeGroups));
        }
        fr->cg_nalloc = over_alloc_dd(numChargeGroups);
        srenew(fr->cginfo, fr->cg_nalloc);
        if (fr->cutoff_scheme == ecutsGROUP)
        {
            srenew(fr->cg_cm, fr->cg_nalloc);
        }
    }
    if (fr->cutoff_scheme == ecutsVERLET)
    {
        /* We don't use charge groups, we use x in state to set up
         * the atom communication.
         */
        dd_resize_state(state, f, numChargeGroups);
    }
}

static void dd_distribute_vec_sendrecv(gmx_domdec_t *dd, t_block *cgs,
                                       const rvec *v, rvec *lv)
{
    gmx_domdec_master_t *ma;
    int                  n, i, c, a, nalloc = 0;
    rvec                *buf = nullptr;

    if (DDMASTER(dd))
    {
        ma  = dd->ma;

        for (n = 0; n < dd->nnodes; n++)
        {
            if (n != dd->rank)
            {
                if (ma->nat[n] > nalloc)
                {
                    nalloc = over_alloc_dd(ma->nat[n]);
                    srenew(buf, nalloc);
                }
                /* Use lv as a temporary buffer */
                a = 0;
                for (i = ma->index[n]; i < ma->index[n+1]; i++)
                {
                    for (c = cgs->index[ma->cg[i]]; c < cgs->index[ma->cg[i]+1]; c++)
                    {
                        copy_rvec(v[c], buf[a++]);
                    }
                }
                if (a != ma->nat[n])
                {
                    gmx_fatal(FARGS, "Internal error a (%d) != nat (%d)",
                              a, ma->nat[n]);
                }

#if GMX_MPI
                MPI_Send(buf, ma->nat[n]*sizeof(rvec), MPI_BYTE,
                         DDRANK(dd, n), n, dd->mpi_comm_all);
#endif
            }
        }
        sfree(buf);
        n = DDMASTERRANK(dd);
        a = 0;
        for (i = ma->index[n]; i < ma->index[n+1]; i++)
        {
            for (c = cgs->index[ma->cg[i]]; c < cgs->index[ma->cg[i]+1]; c++)
            {
                copy_rvec(v[c], lv[a++]);
            }
        }
    }
    else
    {
#if GMX_MPI
        MPI_Recv(lv, dd->nat_home*sizeof(rvec), MPI_BYTE, DDMASTERRANK(dd),
                 MPI_ANY_TAG, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
    }
}

static void dd_distribute_vec_scatterv(gmx_domdec_t *dd, t_block *cgs,
                                       const rvec *v, rvec *lv)
{
    gmx_domdec_master_t *ma;
    int                 *scounts = nullptr, *disps = nullptr;
    int                  n, i, c, a;
    rvec                *buf = nullptr;

    if (DDMASTER(dd))
    {
        ma  = dd->ma;

        get_commbuffer_counts(dd, &scounts, &disps);

        buf = ma->vbuf;
        a   = 0;
        for (n = 0; n < dd->nnodes; n++)
        {
            for (i = ma->index[n]; i < ma->index[n+1]; i++)
            {
                for (c = cgs->index[ma->cg[i]]; c < cgs->index[ma->cg[i]+1]; c++)
                {
                    copy_rvec(v[c], buf[a++]);
                }
            }
        }
    }

    dd_scatterv(dd, scounts, disps, buf, dd->nat_home*sizeof(rvec), lv);
}

static void dd_distribute_vec(gmx_domdec_t *dd, t_block *cgs,
                              const rvec *v, rvec *lv)
{
    if (dd->nnodes <= GMX_DD_NNODES_SENDRECV)
    {
        dd_distribute_vec_sendrecv(dd, cgs, v, lv);
    }
    else
    {
        dd_distribute_vec_scatterv(dd, cgs, v, lv);
    }
}

static void dd_distribute_dfhist(gmx_domdec_t *dd, df_history_t *dfhist)
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
        dd_bcast(dd, sizeof(int)*nlam, dfhist->n_at_lam);
        dd_bcast(dd, sizeof(real)*nlam, dfhist->wl_histo);
        dd_bcast(dd, sizeof(real)*nlam, dfhist->sum_weights);
        dd_bcast(dd, sizeof(real)*nlam, dfhist->sum_dg);
        dd_bcast(dd, sizeof(real)*nlam, dfhist->sum_minvar);
        dd_bcast(dd, sizeof(real)*nlam, dfhist->sum_variance);

        for (int i = 0; i < nlam; i++)
        {
            dd_bcast(dd, sizeof(real)*nlam, dfhist->accum_p[i]);
            dd_bcast(dd, sizeof(real)*nlam, dfhist->accum_m[i]);
            dd_bcast(dd, sizeof(real)*nlam, dfhist->accum_p2[i]);
            dd_bcast(dd, sizeof(real)*nlam, dfhist->accum_m2[i]);
            dd_bcast(dd, sizeof(real)*nlam, dfhist->Tij[i]);
            dd_bcast(dd, sizeof(real)*nlam, dfhist->Tij_empirical[i]);
        }
    }
}

static void dd_distribute_state(gmx_domdec_t *dd, t_block *cgs,
                                t_state *state, t_state *state_local,
                                PaddedRVecVector *f)
{
    int nh = state_local->nhchainlength;

    if (DDMASTER(dd))
    {
        GMX_RELEASE_ASSERT(state->nhchainlength == nh, "The global and local Nose-Hoover chain lengths should match");

        for (int i = 0; i < efptNR; i++)
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
                state_local->nosehoover_xi[i*nh+j]        = state->nosehoover_xi[i*nh+j];
                state_local->nosehoover_vxi[i*nh+j]       = state->nosehoover_vxi[i*nh+j];
            }
            state_local->therm_integral[i] = state->therm_integral[i];
        }
        for (int i = 0; i < state_local->nnhpres; i++)
        {
            for (int j = 0; j < nh; j++)
            {
                state_local->nhpres_xi[i*nh+j]        = state->nhpres_xi[i*nh+j];
                state_local->nhpres_vxi[i*nh+j]       = state->nhpres_vxi[i*nh+j];
            }
        }
        state_local->baros_integral = state->baros_integral;
    }
    dd_bcast(dd, ((efptNR)*sizeof(real)), state_local->lambda.data());
    dd_bcast(dd, sizeof(int), &state_local->fep_state);
    dd_bcast(dd, sizeof(real), &state_local->veta);
    dd_bcast(dd, sizeof(real), &state_local->vol0);
    dd_bcast(dd, sizeof(state_local->box), state_local->box);
    dd_bcast(dd, sizeof(state_local->box_rel), state_local->box_rel);
    dd_bcast(dd, sizeof(state_local->boxv), state_local->boxv);
    dd_bcast(dd, sizeof(state_local->svir_prev), state_local->svir_prev);
    dd_bcast(dd, sizeof(state_local->fvir_prev), state_local->fvir_prev);
    dd_bcast(dd, ((state_local->ngtc*nh)*sizeof(double)), state_local->nosehoover_xi.data());
    dd_bcast(dd, ((state_local->ngtc*nh)*sizeof(double)), state_local->nosehoover_vxi.data());
    dd_bcast(dd, state_local->ngtc*sizeof(double), state_local->therm_integral.data());
    dd_bcast(dd, ((state_local->nnhpres*nh)*sizeof(double)), state_local->nhpres_xi.data());
    dd_bcast(dd, ((state_local->nnhpres*nh)*sizeof(double)), state_local->nhpres_vxi.data());

    /* communicate df_history -- required for restarting from checkpoint */
    dd_distribute_dfhist(dd, state_local->dfhist);

    dd_resize_state(state_local, f, dd->nat_home);

    if (state_local->flags & (1 << estX))
    {
        const rvec *xGlobal = (DDMASTER(dd) ? as_rvec_array(state->x.data()) : nullptr);
        dd_distribute_vec(dd, cgs, xGlobal, as_rvec_array(state_local->x.data()));
    }
    if (state_local->flags & (1 << estV))
    {
        const rvec *vGlobal = (DDMASTER(dd) ? as_rvec_array(state->v.data()) : nullptr);
        dd_distribute_vec(dd, cgs, vGlobal, as_rvec_array(state_local->v.data()));
    }
    if (state_local->flags & (1 << estCGP))
    {
        const rvec *cgpGlobal = (DDMASTER(dd) ? as_rvec_array(state->cg_p.data()) : nullptr);
        dd_distribute_vec(dd, cgs, cgpGlobal, as_rvec_array(state_local->cg_p.data()));
    }
}

static char dim2char(int dim)
{
    char c = '?';

    switch (dim)
    {
        case XX: c = 'X'; break;
        case YY: c = 'Y'; break;
        case ZZ: c = 'Z'; break;
        default: gmx_fatal(FARGS, "Unknown dim %d", dim);
    }

    return c;
}

static void write_dd_grid_pdb(const char *fn, gmx_int64_t step,
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

void write_dd_pdb(const char *fn, gmx_int64_t step, const char *title,
                  const gmx_mtop_t *mtop, t_commrec *cr,
                  int natoms, rvec x[], matrix box)
{
    char          fname[STRLEN], buf[22];
    FILE         *out;
    int           i, ii, resnr, c;
    const char   *atomname, *resname;
    real          b;
    gmx_domdec_t *dd;

    dd = cr->dd;
    if (natoms == -1)
    {
        natoms = dd->comm->nat[ddnatVSITE];
    }

    sprintf(fname, "%s_%s_n%d.pdb", fn, gmx_step_str(step, buf), cr->sim_nodeid);

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, dd->bScrewPBC ? epbcSCREW : epbcXYZ, box);
    int molb = 0;
    for (i = 0; i < natoms; i++)
    {
        ii = dd->gatindex[i];
        mtopGetAtomAndResidueName(mtop, ii, &molb, &atomname, &resnr, &resname, nullptr);
        if (i < dd->comm->nat[ddnatZONE])
        {
            c = 0;
            while (i >= dd->cgindex[dd->comm->zones.cg_range[c+1]])
            {
                c++;
            }
            b = c;
        }
        else if (i < dd->comm->nat[ddnatVSITE])
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

static int gmx_ddcoord2pmeindex(t_commrec *cr, int x, int y, int z)
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

static int ddcoord2simnodeid(t_commrec *cr, int x, int y, int z)
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

void get_pme_nnodes(const gmx_domdec_t *dd,
                    int *npmenodes_x, int *npmenodes_y)
{
    if (dd != nullptr)
    {
        *npmenodes_x = dd->comm->npmenodes_x;
        *npmenodes_y = dd->comm->npmenodes_y;
    }
    else
    {
        *npmenodes_x = 1;
        *npmenodes_y = 1;
    }
}

void get_pme_ddnodes(t_commrec *cr, int pmenodeid,
                     int *nmy_ddnodes, int **my_ddnodes, int *node_peer)
{
    gmx_domdec_t *dd;
    int           x, y, z;
    ivec          coord, coord_pme;

    dd = cr->dd;

    snew(*my_ddnodes, (dd->nnodes+cr->npmenodes-1)/cr->npmenodes);

    *nmy_ddnodes = 0;
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
                        (*my_ddnodes)[(*nmy_ddnodes)++] = ddcoord2simnodeid(cr, x, y, z);
                    }
                }
                else
                {
                    /* The slab corresponds to the nodeid in the PME group */
                    if (gmx_ddcoord2pmeindex(cr, x, y, z) == pmenodeid)
                    {
                        (*my_ddnodes)[(*nmy_ddnodes)++] = ddcoord2simnodeid(cr, x, y, z);
                    }
                }
            }
        }
    }

    /* The last PP-only node is the peer node */
    *node_peer = (*my_ddnodes)[*nmy_ddnodes-1];

    if (debug)
    {
        fprintf(debug, "Receive coordinates from PP ranks:");
        for (x = 0; x < *nmy_ddnodes; x++)
        {
            fprintf(debug, " %d", (*my_ddnodes)[x]);
        }
        fprintf(debug, "\n");
    }
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

static void rebuild_cgindex(gmx_domdec_t *dd,
                            const int *gcgs_index, const t_state *state)
{
    int * gmx_restrict dd_cg_gl = dd->index_gl;
    int * gmx_restrict cgindex  = dd->cgindex;
    int                nat      = 0;

    /* Copy back the global charge group indices from state
     * and rebuild the local charge group to atom index.
     */
    cgindex[0] = nat;
    for (unsigned int i = 0; i < state->cg_gl.size(); i++)
    {
        cgindex[i]  = nat;
        int cg_gl   = state->cg_gl[i];
        dd_cg_gl[i] = cg_gl;
        nat        += gcgs_index[cg_gl+1] - gcgs_index[cg_gl];
    }
    cgindex[state->cg_gl.size()] = nat;

    dd->ncg_home = state->cg_gl.size();
    dd->nat_home = nat;

    set_zones_ncg_home(dd);
}

static int ddcginfo(const cginfo_mb_t *cginfo_mb, int cg)
{
    while (cg >= cginfo_mb->cg_end)
    {
        cginfo_mb++;
    }

    return cginfo_mb->cginfo[(cg - cginfo_mb->cg_start) % cginfo_mb->cg_mod];
}

static void dd_set_cginfo(int *index_gl, int cg0, int cg1,
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
    int          nzone, zone, zone1, cg0, cg1, cg1_p1, cg, cg_gl, a, a_gl;
    int         *zone2cg, *zone_ncg1, *index_gl, *gatindex;
    gmx_bool     bCGs;

    if (dd->nat_tot > dd->gatindex_nalloc)
    {
        dd->gatindex_nalloc = over_alloc_dd(dd->nat_tot);
        srenew(dd->gatindex, dd->gatindex_nalloc);
    }

    nzone      = dd->comm->zones.n;
    zone2cg    = dd->comm->zones.cg_range;
    zone_ncg1  = dd->comm->zone_ncg1;
    index_gl   = dd->index_gl;
    gatindex   = dd->gatindex;
    bCGs       = dd->comm->bCGs;

    if (zone2cg[1] != dd->ncg_home)
    {
        gmx_incons("dd->ncg_zone is not up to date");
    }

    /* Make the local to global and global to local atom index */
    a = dd->cgindex[cg_start];
    for (zone = 0; zone < nzone; zone++)
    {
        if (zone == 0)
        {
            cg0 = cg_start;
        }
        else
        {
            cg0 = zone2cg[zone];
        }
        cg1    = zone2cg[zone+1];
        cg1_p1 = cg0 + zone_ncg1[zone];

        for (cg = cg0; cg < cg1; cg++)
        {
            zone1 = zone;
            if (cg >= cg1_p1)
            {
                /* Signal that this cg is from more than one pulse away */
                zone1 += nzone;
            }
            cg_gl = index_gl[cg];
            if (bCGs)
            {
                for (a_gl = gcgs_index[cg_gl]; a_gl < gcgs_index[cg_gl+1]; a_gl++)
                {
                    gatindex[a] = a_gl;
                    ga2la_set(dd->ga2la, a_gl, a, zone1);
                    a++;
                }
            }
            else
            {
                gatindex[a] = cg_gl;
                ga2la_set(dd->ga2la, cg_gl, a, zone1);
                a++;
            }
        }
    }
}

static int check_bLocalCG(gmx_domdec_t *dd, int ncg_sys, const char *bLocalCG,
                          const char *where)
{
    int i, ngl, nerr;

    nerr = 0;
    if (bLocalCG == nullptr)
    {
        return nerr;
    }
    for (i = 0; i < dd->ncg_tot; i++)
    {
        if (!bLocalCG[dd->index_gl[i]])
        {
            fprintf(stderr,
                    "DD rank %d, %s: cg %d, global cg %d is not marked in bLocalCG (ncg_home %d)\n", dd->rank, where, i+1, dd->index_gl[i]+1, dd->ncg_home);
            nerr++;
        }
    }
    ngl = 0;
    for (i = 0; i < ncg_sys; i++)
    {
        if (bLocalCG[i])
        {
            ngl++;
        }
    }
    if (ngl != dd->ncg_tot)
    {
        fprintf(stderr, "DD rank %d, %s: In bLocalCG %d cgs are marked as local, whereas there are %d\n", dd->rank, where, ngl, dd->ncg_tot);
        nerr++;
    }

    return nerr;
}

static void check_index_consistency(gmx_domdec_t *dd,
                                    int natoms_sys, int ncg_sys,
                                    const char *where)
{
    int   nerr, ngl, i, a, cell;
    int  *have;

    nerr = 0;

    if (dd->comm->DD_debug > 1)
    {
        snew(have, natoms_sys);
        for (a = 0; a < dd->nat_tot; a++)
        {
            if (have[dd->gatindex[a]] > 0)
            {
                fprintf(stderr, "DD rank %d: global atom %d occurs twice: index %d and %d\n", dd->rank, dd->gatindex[a]+1, have[dd->gatindex[a]], a+1);
            }
            else
            {
                have[dd->gatindex[a]] = a + 1;
            }
        }
        sfree(have);
    }

    snew(have, dd->nat_tot);

    ngl  = 0;
    for (i = 0; i < natoms_sys; i++)
    {
        if (ga2la_get(dd->ga2la, i, &a, &cell))
        {
            if (a >= dd->nat_tot)
            {
                fprintf(stderr, "DD rank %d: global atom %d marked as local atom %d, which is larger than nat_tot (%d)\n", dd->rank, i+1, a+1, dd->nat_tot);
                nerr++;
            }
            else
            {
                have[a] = 1;
                if (dd->gatindex[a] != i)
                {
                    fprintf(stderr, "DD rank %d: global atom %d marked as local atom %d, which has global atom index %d\n", dd->rank, i+1, a+1, dd->gatindex[a]+1);
                    nerr++;
                }
            }
            ngl++;
        }
    }
    if (ngl != dd->nat_tot)
    {
        fprintf(stderr,
                "DD rank %d, %s: %d global atom indices, %d local atoms\n",
                dd->rank, where, ngl, dd->nat_tot);
    }
    for (a = 0; a < dd->nat_tot; a++)
    {
        if (have[a] == 0)
        {
            fprintf(stderr,
                    "DD rank %d, %s: local atom %d, global %d has no global index\n",
                    dd->rank, where, a+1, dd->gatindex[a]+1);
        }
    }
    sfree(have);

    nerr += check_bLocalCG(dd, ncg_sys, dd->comm->bLocalCG, where);

    if (nerr > 0)
    {
        gmx_fatal(FARGS, "DD rank %d, %s: %d atom/cg index inconsistencies",
                  dd->rank, where, nerr);
    }
}

static void clear_dd_indices(gmx_domdec_t *dd, int cg_start, int a_start)
{
    int   i;
    char *bLocalCG;

    if (a_start == 0)
    {
        /* Clear the whole list without searching */
        ga2la_clear(dd->ga2la);
    }
    else
    {
        for (i = a_start; i < dd->nat_tot; i++)
        {
            ga2la_del(dd->ga2la, dd->gatindex[i]);
        }
    }

    bLocalCG = dd->comm->bLocalCG;
    if (bLocalCG)
    {
        for (i = cg_start; i < dd->ncg_tot; i++)
        {
            bLocalCG[dd->index_gl[i]] = FALSE;
        }
    }

    dd_clear_local_vsite_indices(dd);

    if (dd->constraints)
    {
        dd_clear_local_constraint_indices(dd);
    }
}

/* This function should be used for moving the domain boudaries during DLB,
 * for obtaining the minimum cell size. It checks the initially set limit
 * comm->cellsize_min, for bonded and initial non-bonded cut-offs,
 * and, possibly, a longer cut-off limit set for PME load balancing.
 */
static real cellsize_min_dlb(gmx_domdec_comm_t *comm, int dim_ind, int dim)
{
    real cellsize_min;

    cellsize_min = comm->cellsize_min[dim];

    if (!comm->bVacDLBNoLimit)
    {
        /* The cut-off might have changed, e.g. by PME load balacning,
         * from the value used to set comm->cellsize_min, so check it.
         */
        cellsize_min = std::max(cellsize_min, comm->cutoff/comm->cd[dim_ind].np_dlb);

        if (comm->bPMELoadBalDLBLimits)
        {
            /* Check for the cut-off limit set by the PME load balancing */
            cellsize_min = std::max(cellsize_min, comm->PMELoadBal_max_cutoff/comm->cd[dim_ind].np_dlb);
        }
    }

    return cellsize_min;
}

static real grid_jump_limit(gmx_domdec_comm_t *comm, real cutoff,
                            int dim_ind)
{
    real grid_jump_limit;

    /* The distance between the boundaries of cells at distance
     * x+-1,y+-1 or y+-1,z+-1 is limited by the cut-off restrictions
     * and by the fact that cells should not be shifted by more than
     * half their size, such that cg's only shift by one cell
     * at redecomposition.
     */
    grid_jump_limit = comm->cellsize_limit;
    if (!comm->bVacDLBNoLimit)
    {
        if (comm->bPMELoadBalDLBLimits)
        {
            cutoff = std::max(cutoff, comm->PMELoadBal_max_cutoff);
        }
        grid_jump_limit = std::max(grid_jump_limit,
                                   cutoff/comm->cd[dim_ind].np);
    }

    return grid_jump_limit;
}

static gmx_bool check_grid_jump(gmx_int64_t     step,
                                gmx_domdec_t   *dd,
                                real            cutoff,
                                gmx_ddbox_t    *ddbox,
                                gmx_bool        bFatal)
{
    gmx_domdec_comm_t *comm;
    int                d, dim;
    real               limit, bfac;
    gmx_bool           bInvalid;

    bInvalid = FALSE;

    comm = dd->comm;

    for (d = 1; d < dd->ndim; d++)
    {
        dim   = dd->dim[d];
        limit = grid_jump_limit(comm, cutoff, d);
        bfac  = ddbox->box_size[dim];
        if (ddbox->tric_dir[dim])
        {
            bfac *= ddbox->skew_fac[dim];
        }
        if ((comm->cell_f1[d] - comm->cell_f_max0[d])*bfac <  limit ||
                                                              (comm->cell_f0[d] - comm->cell_f_min1[d])*bfac > -limit)
        {
            bInvalid = TRUE;

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

    return bInvalid;
}

static int dd_load_count(gmx_domdec_comm_t *comm)
{
    return (comm->eFlop ? comm->flop_n : comm->cycl_n[ddCyclF]);
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

static void set_pme_maxshift(gmx_domdec_t *dd, gmx_ddpme_t *ddpme,
                             gmx_bool bUniform, const gmx_ddbox_t *ddbox,
                             const real *cell_f)
{
    gmx_domdec_comm_t *comm;
    int                nc, ns, s;
    int               *xmin, *xmax;
    real               range, pme_boundary;
    int                sh;

    comm = dd->comm;
    nc   = dd->nc[ddpme->dim];
    ns   = ddpme->nslab;

    if (!ddpme->dim_match)
    {
        /* PP decomposition is not along dim: the worst situation */
        sh = ns/2;
    }
    else if (ns <= 3 || (bUniform && ns == nc))
    {
        /* The optimal situation */
        sh = 1;
    }
    else
    {
        /* We need to check for all pme nodes which nodes they
         * could possibly need to communicate with.
         */
        xmin = ddpme->pp_min;
        xmax = ddpme->pp_max;
        /* Allow for atoms to be maximally 2/3 times the cut-off
         * out of their DD cell. This is a reasonable balance between
         * between performance and support for most charge-group/cut-off
         * combinations.
         */
        range  = 2.0/3.0*comm->cutoff/ddbox->box_size[ddpme->dim];
        /* Avoid extra communication when we are exactly at a boundary */
        range *= 0.999;

        sh = 1;
        for (s = 0; s < ns; s++)
        {
            /* PME slab s spreads atoms between box frac. s/ns and (s+1)/ns */
            pme_boundary = (real)s/ns;
            while (sh+1 < ns &&
                   ((s-(sh+1) >= 0 &&
                     cell_f[xmax[s-(sh+1)   ]+1]     + range > pme_boundary) ||
                    (s-(sh+1) <  0 &&
                     cell_f[xmax[s-(sh+1)+ns]+1] - 1 + range > pme_boundary)))
            {
                sh++;
            }
            pme_boundary = (real)(s+1)/ns;
            while (sh+1 < ns &&
                   ((s+(sh+1) <  ns &&
                     cell_f[xmin[s+(sh+1)   ]  ]     - range < pme_boundary) ||
                    (s+(sh+1) >= ns &&
                     cell_f[xmin[s+(sh+1)-ns]  ] + 1 - range < pme_boundary)))
            {
                sh++;
            }
        }
    }

    ddpme->maxshift = sh;

    if (debug)
    {
        fprintf(debug, "PME slab communication range for dim %d is %d\n",
                ddpme->dim, ddpme->maxshift);
    }
}

static void check_box_size(gmx_domdec_t *dd, gmx_ddbox_t *ddbox)
{
    int d, dim;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        if (dim < ddbox->nboundeddim &&
            ddbox->box_size[dim]*ddbox->skew_fac[dim] <
            dd->nc[dim]*dd->comm->cellsize_limit*DD_CELL_MARGIN)
        {
            gmx_fatal(FARGS, "The %c-size of the box (%f) times the triclinic skew factor (%f) is smaller than the number of DD cells (%d) times the smallest allowed cell size (%f)\n",
                      dim2char(dim), ddbox->box_size[dim], ddbox->skew_fac[dim],
                      dd->nc[dim], dd->comm->cellsize_limit);
        }
    }
}

enum {
    setcellsizeslbLOCAL, setcellsizeslbMASTER, setcellsizeslbPULSE_ONLY
};

/* Set the domain boundaries. Use for static (or no) load balancing,
 * and also for the starting state for dynamic load balancing.
 * setmode determine if and where the boundaries are stored, use enum above.
 * Returns the number communication pulses in npulse.
 */
static void set_dd_cell_sizes_slb(gmx_domdec_t *dd, const gmx_ddbox_t *ddbox,
                                  int setmode, ivec npulse)
{
    gmx_domdec_comm_t *comm;
    int                d, j;
    rvec               cellsize_min;
    real              *cell_x, cell_dx, cellsize;

    comm = dd->comm;

    for (d = 0; d < DIM; d++)
    {
        cellsize_min[d] = ddbox->box_size[d]*ddbox->skew_fac[d];
        npulse[d]       = 1;
        if (dd->nc[d] == 1 || comm->slb_frac[d] == nullptr)
        {
            /* Uniform grid */
            cell_dx = ddbox->box_size[d]/dd->nc[d];
            switch (setmode)
            {
                case setcellsizeslbMASTER:
                    for (j = 0; j < dd->nc[d]+1; j++)
                    {
                        dd->ma->cell_x[d][j] = ddbox->box0[d] + j*cell_dx;
                    }
                    break;
                case setcellsizeslbLOCAL:
                    comm->cell_x0[d] = ddbox->box0[d] + (dd->ci[d]  )*cell_dx;
                    comm->cell_x1[d] = ddbox->box0[d] + (dd->ci[d]+1)*cell_dx;
                    break;
                default:
                    break;
            }
            cellsize = cell_dx*ddbox->skew_fac[d];
            while (cellsize*npulse[d] < comm->cutoff)
            {
                npulse[d]++;
            }
            cellsize_min[d] = cellsize;
        }
        else
        {
            /* Statically load balanced grid */
            /* Also when we are not doing a master distribution we determine
             * all cell borders in a loop to obtain identical values
             * to the master distribution case and to determine npulse.
             */
            if (setmode == setcellsizeslbMASTER)
            {
                cell_x = dd->ma->cell_x[d];
            }
            else
            {
                snew(cell_x, dd->nc[d]+1);
            }
            cell_x[0] = ddbox->box0[d];
            for (j = 0; j < dd->nc[d]; j++)
            {
                cell_dx     = ddbox->box_size[d]*comm->slb_frac[d][j];
                cell_x[j+1] = cell_x[j] + cell_dx;
                cellsize    = cell_dx*ddbox->skew_fac[d];
                while (cellsize*npulse[d] < comm->cutoff &&
                       npulse[d] < dd->nc[d]-1)
                {
                    npulse[d]++;
                }
                cellsize_min[d] = std::min(cellsize_min[d], cellsize);
            }
            if (setmode == setcellsizeslbLOCAL)
            {
                comm->cell_x0[d] = cell_x[dd->ci[d]];
                comm->cell_x1[d] = cell_x[dd->ci[d]+1];
            }
            if (setmode != setcellsizeslbMASTER)
            {
                sfree(cell_x);
            }
        }
        /* The following limitation is to avoid that a cell would receive
         * some of its own home charge groups back over the periodic boundary.
         * Double charge groups cause trouble with the global indices.
         */
        if (d < ddbox->npbcdim &&
            dd->nc[d] > 1 && npulse[d] >= dd->nc[d])
        {
            char error_string[STRLEN];

            sprintf(error_string,
                    "The box size in direction %c (%f) times the triclinic skew factor (%f) is too small for a cut-off of %f with %d domain decomposition cells, use 1 or more than %d %s or increase the box size in this direction",
                    dim2char(d), ddbox->box_size[d], ddbox->skew_fac[d],
                    comm->cutoff,
                    dd->nc[d], dd->nc[d],
                    dd->nnodes > dd->nc[d] ? "cells" : "ranks");

            if (setmode == setcellsizeslbLOCAL)
            {
                gmx_fatal_collective(FARGS, dd->mpi_comm_all, DDMASTER(dd),
                                     error_string);
            }
            else
            {
                gmx_fatal(FARGS, error_string);
            }
        }
    }

    if (!isDlbOn(comm))
    {
        copy_rvec(cellsize_min, comm->cellsize_min);
    }

    for (d = 0; d < comm->npmedecompdim; d++)
    {
        set_pme_maxshift(dd, &comm->ddpme[d],
                         comm->slb_frac[dd->dim[d]] == nullptr, ddbox,
                         comm->ddpme[d].slb_dim_f);
    }
}


static void dd_cell_sizes_dlb_root_enforce_limits(gmx_domdec_t *dd,
                                                  int d, int dim, domdec_root_t *root,
                                                  const gmx_ddbox_t *ddbox,
                                                  gmx_bool bUniform, gmx_int64_t step, real cellsize_limit_f, int range[])
{
    gmx_domdec_comm_t *comm;
    int                ncd, i, j, nmin, nmin_old;
    gmx_bool           bLimLo, bLimHi;
    real              *cell_size;
    real               fac, halfway, cellsize_limit_f_i, region_size;
    gmx_bool           bPBC, bLastHi = FALSE;
    int                nrange[] = {range[0], range[1]};

    region_size = root->cell_f[range[1]]-root->cell_f[range[0]];

    comm = dd->comm;

    ncd = dd->nc[dim];

    bPBC = (dim < ddbox->npbcdim);

    cell_size = root->buf_ncd;

    if (debug)
    {
        fprintf(debug, "enforce_limits: %d %d\n", range[0], range[1]);
    }

    /* First we need to check if the scaling does not make cells
     * smaller than the smallest allowed size.
     * We need to do this iteratively, since if a cell is too small,
     * it needs to be enlarged, which makes all the other cells smaller,
     * which could in turn make another cell smaller than allowed.
     */
    for (i = range[0]; i < range[1]; i++)
    {
        root->bCellMin[i] = FALSE;
    }
    nmin = 0;
    do
    {
        nmin_old = nmin;
        /* We need the total for normalization */
        fac = 0;
        for (i = range[0]; i < range[1]; i++)
        {
            if (root->bCellMin[i] == FALSE)
            {
                fac += cell_size[i];
            }
        }
        fac = ( region_size - nmin*cellsize_limit_f)/fac; /* substracting cells already set to cellsize_limit_f */
        /* Determine the cell boundaries */
        for (i = range[0]; i < range[1]; i++)
        {
            if (root->bCellMin[i] == FALSE)
            {
                cell_size[i] *= fac;
                if (!bPBC && (i == 0 || i == dd->nc[dim] -1))
                {
                    cellsize_limit_f_i = 0;
                }
                else
                {
                    cellsize_limit_f_i = cellsize_limit_f;
                }
                if (cell_size[i] < cellsize_limit_f_i)
                {
                    root->bCellMin[i] = TRUE;
                    cell_size[i]      = cellsize_limit_f_i;
                    nmin++;
                }
            }
            root->cell_f[i+1] = root->cell_f[i] + cell_size[i];
        }
    }
    while (nmin > nmin_old);

    i            = range[1]-1;
    cell_size[i] = root->cell_f[i+1] - root->cell_f[i];
    /* For this check we should not use DD_CELL_MARGIN,
     * but a slightly smaller factor,
     * since rounding could get use below the limit.
     */
    if (bPBC && cell_size[i] < cellsize_limit_f*DD_CELL_MARGIN2/DD_CELL_MARGIN)
    {
        char buf[22];
        gmx_fatal(FARGS, "step %s: the dynamic load balancing could not balance dimension %c: box size %f, triclinic skew factor %f, #cells %d, minimum cell size %f\n",
                  gmx_step_str(step, buf),
                  dim2char(dim), ddbox->box_size[dim], ddbox->skew_fac[dim],
                  ncd, comm->cellsize_min[dim]);
    }

    root->bLimited = (nmin > 0) || (range[0] > 0) || (range[1] < ncd);

    if (!bUniform)
    {
        /* Check if the boundary did not displace more than halfway
         * each of the cells it bounds, as this could cause problems,
         * especially when the differences between cell sizes are large.
         * If changes are applied, they will not make cells smaller
         * than the cut-off, as we check all the boundaries which
         * might be affected by a change and if the old state was ok,
         * the cells will at most be shrunk back to their old size.
         */
        for (i = range[0]+1; i < range[1]; i++)
        {
            halfway = 0.5*(root->old_cell_f[i] + root->old_cell_f[i-1]);
            if (root->cell_f[i] < halfway)
            {
                root->cell_f[i] = halfway;
                /* Check if the change also causes shifts of the next boundaries */
                for (j = i+1; j < range[1]; j++)
                {
                    if (root->cell_f[j] < root->cell_f[j-1] + cellsize_limit_f)
                    {
                        root->cell_f[j] =  root->cell_f[j-1] + cellsize_limit_f;
                    }
                }
            }
            halfway = 0.5*(root->old_cell_f[i] + root->old_cell_f[i+1]);
            if (root->cell_f[i] > halfway)
            {
                root->cell_f[i] = halfway;
                /* Check if the change also causes shifts of the next boundaries */
                for (j = i-1; j >= range[0]+1; j--)
                {
                    if (root->cell_f[j] > root->cell_f[j+1] - cellsize_limit_f)
                    {
                        root->cell_f[j] = root->cell_f[j+1] - cellsize_limit_f;
                    }
                }
            }
        }
    }

    /* nrange is defined as [lower, upper) range for new call to enforce_limits */
    /* find highest violation of LimLo (a) and the following violation of LimHi (thus the lowest following) (b)
     * then call enforce_limits for (oldb,a), (a,b). In the next step: (b,nexta). oldb and nexta can be the boundaries.
     * for a and b nrange is used */
    if (d > 0)
    {
        /* Take care of the staggering of the cell boundaries */
        if (bUniform)
        {
            for (i = range[0]; i < range[1]; i++)
            {
                root->cell_f_max0[i] = root->cell_f[i];
                root->cell_f_min1[i] = root->cell_f[i+1];
            }
        }
        else
        {
            for (i = range[0]+1; i < range[1]; i++)
            {
                bLimLo = (root->cell_f[i] < root->bound_min[i]);
                bLimHi = (root->cell_f[i] > root->bound_max[i]);
                if (bLimLo && bLimHi)
                {
                    /* Both limits violated, try the best we can */
                    /* For this case we split the original range (range) in two parts and care about the other limitiations in the next iteration. */
                    root->cell_f[i] = 0.5*(root->bound_min[i] + root->bound_max[i]);
                    nrange[0]       = range[0];
                    nrange[1]       = i;
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);

                    nrange[0] = i;
                    nrange[1] = range[1];
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);

                    return;
                }
                else if (bLimLo)
                {
                    /* root->cell_f[i] = root->bound_min[i]; */
                    nrange[1] = i;  /* only store violation location. There could be a LimLo violation following with an higher index */
                    bLastHi   = FALSE;
                }
                else if (bLimHi && !bLastHi)
                {
                    bLastHi = TRUE;
                    if (nrange[1] < range[1])   /* found a LimLo before */
                    {
                        root->cell_f[nrange[1]] = root->bound_min[nrange[1]];
                        dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);
                        nrange[0] = nrange[1];
                    }
                    root->cell_f[i] = root->bound_max[i];
                    nrange[1]       = i;
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);
                    nrange[0] = i;
                    nrange[1] = range[1];
                }
            }
            if (nrange[1] < range[1])   /* found last a LimLo */
            {
                root->cell_f[nrange[1]] = root->bound_min[nrange[1]];
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);
                nrange[0] = nrange[1];
                nrange[1] = range[1];
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);
            }
            else if (nrange[0] > range[0]) /* found at least one LimHi */
            {
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, nrange);
            }
        }
    }
}


static void set_dd_cell_sizes_dlb_root(gmx_domdec_t *dd,
                                       int d, int dim, domdec_root_t *root,
                                       const gmx_ddbox_t *ddbox,
                                       gmx_bool bDynamicBox,
                                       gmx_bool bUniform, gmx_int64_t step)
{
    gmx_domdec_comm_t *comm;
    int                ncd, d1, i, pos;
    real              *cell_size;
    real               load_aver, load_i, imbalance, change, change_max, sc;
    real               cellsize_limit_f, dist_min_f, dist_min_f_hard, space;
    real               change_limit;
    real               relax = 0.5;
    gmx_bool           bPBC;
    int                range[] = { 0, 0 };

    comm = dd->comm;

    /* Convert the maximum change from the input percentage to a fraction */
    change_limit = comm->dlb_scale_lim*0.01;

    ncd = dd->nc[dim];

    bPBC = (dim < ddbox->npbcdim);

    cell_size = root->buf_ncd;

    /* Store the original boundaries */
    for (i = 0; i < ncd+1; i++)
    {
        root->old_cell_f[i] = root->cell_f[i];
    }
    if (bUniform)
    {
        for (i = 0; i < ncd; i++)
        {
            cell_size[i] = 1.0/ncd;
        }
    }
    else if (dd_load_count(comm) > 0)
    {
        load_aver  = comm->load[d].sum_m/ncd;
        change_max = 0;
        for (i = 0; i < ncd; i++)
        {
            /* Determine the relative imbalance of cell i */
            load_i    = comm->load[d].load[i*comm->load[d].nload+2];
            imbalance = (load_i - load_aver)/(load_aver > 0 ? load_aver : 1);
            /* Determine the change of the cell size using underrelaxation */
            change     = -relax*imbalance;
            change_max = std::max(change_max, std::max(change, -change));
        }
        /* Limit the amount of scaling.
         * We need to use the same rescaling for all cells in one row,
         * otherwise the load balancing might not converge.
         */
        sc = relax;
        if (change_max > change_limit)
        {
            sc *= change_limit/change_max;
        }
        for (i = 0; i < ncd; i++)
        {
            /* Determine the relative imbalance of cell i */
            load_i    = comm->load[d].load[i*comm->load[d].nload+2];
            imbalance = (load_i - load_aver)/(load_aver > 0 ? load_aver : 1);
            /* Determine the change of the cell size using underrelaxation */
            change       = -sc*imbalance;
            cell_size[i] = (root->cell_f[i+1]-root->cell_f[i])*(1 + change);
        }
    }

    cellsize_limit_f  = cellsize_min_dlb(comm, d, dim)/ddbox->box_size[dim];
    cellsize_limit_f *= DD_CELL_MARGIN;
    dist_min_f_hard   = grid_jump_limit(comm, comm->cutoff, d)/ddbox->box_size[dim];
    dist_min_f        = dist_min_f_hard * DD_CELL_MARGIN;
    if (ddbox->tric_dir[dim])
    {
        cellsize_limit_f /= ddbox->skew_fac[dim];
        dist_min_f       /= ddbox->skew_fac[dim];
    }
    if (bDynamicBox && d > 0)
    {
        dist_min_f *= DD_PRES_SCALE_MARGIN;
    }
    if (d > 0 && !bUniform)
    {
        /* Make sure that the grid is not shifted too much */
        for (i = 1; i < ncd; i++)
        {
            if (root->cell_f_min1[i] - root->cell_f_max0[i-1] < 2 * dist_min_f_hard)
            {
                gmx_incons("Inconsistent DD boundary staggering limits!");
            }
            root->bound_min[i] = root->cell_f_max0[i-1] + dist_min_f;
            space              = root->cell_f[i] - (root->cell_f_max0[i-1] + dist_min_f);
            if (space > 0)
            {
                root->bound_min[i] += 0.5*space;
            }
            root->bound_max[i] = root->cell_f_min1[i] - dist_min_f;
            space              = root->cell_f[i] - (root->cell_f_min1[i] - dist_min_f);
            if (space < 0)
            {
                root->bound_max[i] += 0.5*space;
            }
            if (debug)
            {
                fprintf(debug,
                        "dim %d boundary %d %.3f < %.3f < %.3f < %.3f < %.3f\n",
                        d, i,
                        root->cell_f_max0[i-1] + dist_min_f,
                        root->bound_min[i], root->cell_f[i], root->bound_max[i],
                        root->cell_f_min1[i] - dist_min_f);
            }
        }
    }
    range[1]          = ncd;
    root->cell_f[0]   = 0;
    root->cell_f[ncd] = 1;
    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, root, ddbox, bUniform, step, cellsize_limit_f, range);


    /* After the checks above, the cells should obey the cut-off
     * restrictions, but it does not hurt to check.
     */
    for (i = 0; i < ncd; i++)
    {
        if (debug)
        {
            fprintf(debug, "Relative bounds dim %d  cell %d: %f %f\n",
                    dim, i, root->cell_f[i], root->cell_f[i+1]);
        }

        if ((bPBC || (i != 0 && i != dd->nc[dim]-1)) &&
            root->cell_f[i+1] - root->cell_f[i] <
            cellsize_limit_f/DD_CELL_MARGIN)
        {
            char buf[22];
            fprintf(stderr,
                    "\nWARNING step %s: direction %c, cell %d too small: %f\n",
                    gmx_step_str(step, buf), dim2char(dim), i,
                    (root->cell_f[i+1] - root->cell_f[i])
                    *ddbox->box_size[dim]*ddbox->skew_fac[dim]);
        }
    }

    pos = ncd + 1;
    /* Store the cell boundaries of the lower dimensions at the end */
    for (d1 = 0; d1 < d; d1++)
    {
        root->cell_f[pos++] = comm->cell_f0[d1];
        root->cell_f[pos++] = comm->cell_f1[d1];
    }

    if (d < comm->npmedecompdim)
    {
        /* The master determines the maximum shift for
         * the coordinate communication between separate PME nodes.
         */
        set_pme_maxshift(dd, &comm->ddpme[d], bUniform, ddbox, root->cell_f);
    }
    root->cell_f[pos++] = comm->ddpme[0].maxshift;
    if (d >= 1)
    {
        root->cell_f[pos++] = comm->ddpme[1].maxshift;
    }
}

static void relative_to_absolute_cell_bounds(gmx_domdec_t      *dd,
                                             const gmx_ddbox_t *ddbox,
                                             int                dimind)
{
    gmx_domdec_comm_t *comm;
    int                dim;

    comm = dd->comm;

    /* Set the cell dimensions */
    dim                = dd->dim[dimind];
    comm->cell_x0[dim] = comm->cell_f0[dimind]*ddbox->box_size[dim];
    comm->cell_x1[dim] = comm->cell_f1[dimind]*ddbox->box_size[dim];
    if (dim >= ddbox->nboundeddim)
    {
        comm->cell_x0[dim] += ddbox->box0[dim];
        comm->cell_x1[dim] += ddbox->box0[dim];
    }
}

static void distribute_dd_cell_sizes_dlb(gmx_domdec_t *dd,
                                         int d, int dim, real *cell_f_row,
                                         const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                d1, pos;

    comm = dd->comm;

#if GMX_MPI
    /* Each node would only need to know two fractions,
     * but it is probably cheaper to broadcast the whole array.
     */
    MPI_Bcast(cell_f_row, DD_CELL_F_SIZE(dd, d)*sizeof(real), MPI_BYTE,
              0, comm->mpi_comm_load[d]);
#endif
    /* Copy the fractions for this dimension from the buffer */
    comm->cell_f0[d] = cell_f_row[dd->ci[dim]  ];
    comm->cell_f1[d] = cell_f_row[dd->ci[dim]+1];
    /* The whole array was communicated, so set the buffer position */
    pos = dd->nc[dim] + 1;
    for (d1 = 0; d1 <= d; d1++)
    {
        if (d1 < d)
        {
            /* Copy the cell fractions of the lower dimensions */
            comm->cell_f0[d1] = cell_f_row[pos++];
            comm->cell_f1[d1] = cell_f_row[pos++];
        }
        relative_to_absolute_cell_bounds(dd, ddbox, d1);
    }
    /* Convert the communicated shift from float to int */
    comm->ddpme[0].maxshift = (int)(cell_f_row[pos++] + 0.5);
    if (d >= 1)
    {
        comm->ddpme[1].maxshift = (int)(cell_f_row[pos++] + 0.5);
    }
}

static void set_dd_cell_sizes_dlb_change(gmx_domdec_t *dd,
                                         const gmx_ddbox_t *ddbox,
                                         gmx_bool bDynamicBox,
                                         gmx_bool bUniform, gmx_int64_t step)
{
    gmx_domdec_comm_t *comm;
    int                d, dim, d1;
    gmx_bool           bRowMember, bRowRoot;
    real              *cell_f_row;

    comm = dd->comm;

    for (d = 0; d < dd->ndim; d++)
    {
        dim        = dd->dim[d];
        bRowMember = TRUE;
        bRowRoot   = TRUE;
        for (d1 = d; d1 < dd->ndim; d1++)
        {
            if (dd->ci[dd->dim[d1]] > 0)
            {
                if (d1 != d)
                {
                    bRowMember = FALSE;
                }
                bRowRoot = FALSE;
            }
        }
        if (bRowMember)
        {
            if (bRowRoot)
            {
                set_dd_cell_sizes_dlb_root(dd, d, dim, comm->root[d],
                                           ddbox, bDynamicBox, bUniform, step);
                cell_f_row = comm->root[d]->cell_f;
            }
            else
            {
                cell_f_row = comm->cell_f_row;
            }
            distribute_dd_cell_sizes_dlb(dd, d, dim, cell_f_row, ddbox);
        }
    }
}

static void set_dd_cell_sizes_dlb_nochange(gmx_domdec_t      *dd,
                                           const gmx_ddbox_t *ddbox)
{
    int d;

    /* This function assumes the box is static and should therefore
     * not be called when the box has changed since the last
     * call to dd_partition_system.
     */
    for (d = 0; d < dd->ndim; d++)
    {
        relative_to_absolute_cell_bounds(dd, ddbox, d);
    }
}



static void set_dd_cell_sizes_dlb(gmx_domdec_t *dd,
                                  const gmx_ddbox_t *ddbox, gmx_bool bDynamicBox,
                                  gmx_bool bUniform, gmx_bool bDoDLB, gmx_int64_t step,
                                  gmx_wallcycle_t wcycle)
{
    gmx_domdec_comm_t *comm;
    int                dim;

    comm = dd->comm;

    if (bDoDLB)
    {
        wallcycle_start(wcycle, ewcDDCOMMBOUND);
        set_dd_cell_sizes_dlb_change(dd, ddbox, bDynamicBox, bUniform, step);
        wallcycle_stop(wcycle, ewcDDCOMMBOUND);
    }
    else if (bDynamicBox)
    {
        set_dd_cell_sizes_dlb_nochange(dd, ddbox);
    }

    /* Set the dimensions for which no DD is used */
    for (dim = 0; dim < DIM; dim++)
    {
        if (dd->nc[dim] == 1)
        {
            comm->cell_x0[dim] = 0;
            comm->cell_x1[dim] = ddbox->box_size[dim];
            if (dim >= ddbox->nboundeddim)
            {
                comm->cell_x0[dim] += ddbox->box0[dim];
                comm->cell_x1[dim] += ddbox->box0[dim];
            }
        }
    }
}

static void realloc_comm_ind(gmx_domdec_t *dd, ivec npulse)
{
    int                    d, np, i;
    gmx_domdec_comm_dim_t *cd;

    for (d = 0; d < dd->ndim; d++)
    {
        cd = &dd->comm->cd[d];
        np = npulse[dd->dim[d]];
        if (np > cd->np_nalloc)
        {
            if (debug)
            {
                fprintf(debug, "(Re)allocing cd for %c to %d pulses\n",
                        dim2char(dd->dim[d]), np);
            }
            if (DDMASTER(dd) && cd->np_nalloc > 0)
            {
                fprintf(stderr, "\nIncreasing the number of cell to communicate in dimension %c to %d for the first time\n", dim2char(dd->dim[d]), np);
            }
            srenew(cd->ind, np);
            for (i = cd->np_nalloc; i < np; i++)
            {
                cd->ind[i].index  = nullptr;
                cd->ind[i].nalloc = 0;
            }
            cd->np_nalloc = np;
        }
        cd->np = np;
    }
}


static void set_dd_cell_sizes(gmx_domdec_t *dd,
                              gmx_ddbox_t *ddbox, gmx_bool bDynamicBox,
                              gmx_bool bUniform, gmx_bool bDoDLB, gmx_int64_t step,
                              gmx_wallcycle_t wcycle)
{
    gmx_domdec_comm_t *comm;
    int                d;
    ivec               npulse;

    comm = dd->comm;

    /* Copy the old cell boundaries for the cg displacement check */
    copy_rvec(comm->cell_x0, comm->old_cell_x0);
    copy_rvec(comm->cell_x1, comm->old_cell_x1);

    if (isDlbOn(comm))
    {
        if (DDMASTER(dd))
        {
            check_box_size(dd, ddbox);
        }
        set_dd_cell_sizes_dlb(dd, ddbox, bDynamicBox, bUniform, bDoDLB, step, wcycle);
    }
    else
    {
        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbLOCAL, npulse);
        realloc_comm_ind(dd, npulse);
    }

    if (debug)
    {
        for (d = 0; d < DIM; d++)
        {
            fprintf(debug, "cell_x[%d] %f - %f skew_fac %f\n",
                    d, comm->cell_x0[d], comm->cell_x1[d], ddbox->skew_fac[d]);
        }
    }
}

static void comm_dd_ns_cell_sizes(gmx_domdec_t *dd,
                                  gmx_ddbox_t *ddbox,
                                  rvec cell_ns_x0, rvec cell_ns_x1,
                                  gmx_int64_t step)
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

static void make_tric_corr_matrix(int npbcdim, const matrix box, matrix tcm)
{
    if (YY < npbcdim)
    {
        tcm[YY][XX] = -box[YY][XX]/box[YY][YY];
    }
    else
    {
        tcm[YY][XX] = 0;
    }
    if (ZZ < npbcdim)
    {
        tcm[ZZ][XX] = -(box[ZZ][YY]*tcm[YY][XX] + box[ZZ][XX])/box[ZZ][ZZ];
        tcm[ZZ][YY] = -box[ZZ][YY]/box[ZZ][ZZ];
    }
    else
    {
        tcm[ZZ][XX] = 0;
        tcm[ZZ][YY] = 0;
    }
}

static void check_screw_box(const matrix box)
{
    /* Mathematical limitation */
    if (box[YY][XX] != 0 || box[ZZ][XX] != 0)
    {
        gmx_fatal(FARGS, "With screw pbc the unit cell can not have non-zero off-diagonal x-components");
    }

    /* Limitation due to the asymmetry of the eighth shell method */
    if (box[ZZ][YY] != 0)
    {
        gmx_fatal(FARGS, "pbc=screw with non-zero box_zy is not supported");
    }
}

static void distribute_cg(FILE *fplog,
                          const matrix box, ivec tric_dir, t_block *cgs, rvec pos[],
                          gmx_domdec_t *dd)
{
    gmx_domdec_master_t *ma;
    int                **tmp_ind = nullptr, *tmp_nalloc = nullptr;
    int                  i, icg, j, k, k0, k1, d;
    matrix               tcm;
    rvec                 cg_cm;
    ivec                 ind;
    real                 nrcg, inv_ncg, pos_d;
    int                 *cgindex;
    gmx_bool             bScrew;

    ma = dd->ma;

    snew(tmp_nalloc, dd->nnodes);
    snew(tmp_ind, dd->nnodes);
    for (i = 0; i < dd->nnodes; i++)
    {
        tmp_nalloc[i] = over_alloc_large(cgs->nr/dd->nnodes+1);
        snew(tmp_ind[i], tmp_nalloc[i]);
    }

    /* Clear the count */
    for (i = 0; i < dd->nnodes; i++)
    {
        ma->ncg[i] = 0;
        ma->nat[i] = 0;
    }

    make_tric_corr_matrix(dd->npbcdim, box, tcm);

    cgindex = cgs->index;

    /* Compute the center of geometry for all charge groups */
    for (icg = 0; icg < cgs->nr; icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1 - k0;
        if (nrcg == 1)
        {
            copy_rvec(pos[k0], cg_cm);
        }
        else
        {
            inv_ncg = 1.0/nrcg;

            clear_rvec(cg_cm);
            for (k = k0; (k < k1); k++)
            {
                rvec_inc(cg_cm, pos[k]);
            }
            for (d = 0; (d < DIM); d++)
            {
                cg_cm[d] *= inv_ncg;
            }
        }
        /* Put the charge group in the box and determine the cell index */
        for (d = DIM-1; d >= 0; d--)
        {
            pos_d = cg_cm[d];
            if (d < dd->npbcdim)
            {
                bScrew = (dd->bScrewPBC && d == XX);
                if (tric_dir[d] && dd->nc[d] > 1)
                {
                    /* Use triclinic coordintates for this dimension */
                    for (j = d+1; j < DIM; j++)
                    {
                        pos_d += cg_cm[j]*tcm[j][d];
                    }
                }
                while (pos_d >= box[d][d])
                {
                    pos_d -= box[d][d];
                    rvec_dec(cg_cm, box[d]);
                    if (bScrew)
                    {
                        cg_cm[YY] = box[YY][YY] - cg_cm[YY];
                        cg_cm[ZZ] = box[ZZ][ZZ] - cg_cm[ZZ];
                    }
                    for (k = k0; (k < k1); k++)
                    {
                        rvec_dec(pos[k], box[d]);
                        if (bScrew)
                        {
                            pos[k][YY] = box[YY][YY] - pos[k][YY];
                            pos[k][ZZ] = box[ZZ][ZZ] - pos[k][ZZ];
                        }
                    }
                }
                while (pos_d < 0)
                {
                    pos_d += box[d][d];
                    rvec_inc(cg_cm, box[d]);
                    if (bScrew)
                    {
                        cg_cm[YY] = box[YY][YY] - cg_cm[YY];
                        cg_cm[ZZ] = box[ZZ][ZZ] - cg_cm[ZZ];
                    }
                    for (k = k0; (k < k1); k++)
                    {
                        rvec_inc(pos[k], box[d]);
                        if (bScrew)
                        {
                            pos[k][YY] = box[YY][YY] - pos[k][YY];
                            pos[k][ZZ] = box[ZZ][ZZ] - pos[k][ZZ];
                        }
                    }
                }
            }
            /* This could be done more efficiently */
            ind[d] = 0;
            while (ind[d]+1 < dd->nc[d] && pos_d >= ma->cell_x[d][ind[d]+1])
            {
                ind[d]++;
            }
        }
        i = dd_index(dd->nc, ind);
        if (ma->ncg[i] == tmp_nalloc[i])
        {
            tmp_nalloc[i] = over_alloc_large(ma->ncg[i]+1);
            srenew(tmp_ind[i], tmp_nalloc[i]);
        }
        tmp_ind[i][ma->ncg[i]] = icg;
        ma->ncg[i]++;
        ma->nat[i] += cgindex[icg+1] - cgindex[icg];
    }

    k1 = 0;
    for (i = 0; i < dd->nnodes; i++)
    {
        ma->index[i] = k1;
        for (k = 0; k < ma->ncg[i]; k++)
        {
            ma->cg[k1++] = tmp_ind[i][k];
        }
    }
    ma->index[dd->nnodes] = k1;

    for (i = 0; i < dd->nnodes; i++)
    {
        sfree(tmp_ind[i]);
    }
    sfree(tmp_ind);
    sfree(tmp_nalloc);

    if (fplog)
    {
        // Use double for the sums to avoid natoms^2 overflowing
        // (65537^2 > 2^32)
        int    nat_sum, nat_min, nat_max;
        double nat2_sum;

        nat_sum  = 0;
        nat2_sum = 0;
        nat_min  = ma->nat[0];
        nat_max  = ma->nat[0];
        for (i = 0; i < dd->nnodes; i++)
        {
            nat_sum  += ma->nat[i];
            // cast to double to avoid integer overflows when squaring
            nat2_sum += gmx::square(static_cast<double>(ma->nat[i]));
            nat_min   = std::min(nat_min, ma->nat[i]);
            nat_max   = std::max(nat_max, ma->nat[i]);
        }
        nat_sum  /= dd->nnodes;
        nat2_sum /= dd->nnodes;

        fprintf(fplog, "Atom distribution over %d domains: av %d stddev %d min %d max %d\n",
                dd->nnodes,
                nat_sum,
                static_cast<int>(std::sqrt(nat2_sum - gmx::square(static_cast<double>(nat_sum)) + 0.5)),
                nat_min, nat_max);
    }
}

static void get_cg_distribution(FILE *fplog, gmx_domdec_t *dd,
                                t_block *cgs, const matrix box, gmx_ddbox_t *ddbox,
                                rvec pos[])
{
    gmx_domdec_master_t *ma = nullptr;
    ivec                 npulse;
    int                  i, cg_gl;
    int                 *ibuf, buf2[2] = { 0, 0 };
    gmx_bool             bMaster = DDMASTER(dd);

    if (bMaster)
    {
        ma = dd->ma;

        if (dd->bScrewPBC)
        {
            check_screw_box(box);
        }

        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbMASTER, npulse);

        distribute_cg(fplog, box, ddbox->tric_dir, cgs, pos, dd);
        for (i = 0; i < dd->nnodes; i++)
        {
            ma->ibuf[2*i]   = ma->ncg[i];
            ma->ibuf[2*i+1] = ma->nat[i];
        }
        ibuf = ma->ibuf;
    }
    else
    {
        ibuf = nullptr;
    }
    dd_scatter(dd, 2*sizeof(int), ibuf, buf2);

    dd->ncg_home = buf2[0];
    dd->nat_home = buf2[1];
    dd->ncg_tot  = dd->ncg_home;
    dd->nat_tot  = dd->nat_home;
    if (dd->ncg_home > dd->cg_nalloc || dd->cg_nalloc == 0)
    {
        dd->cg_nalloc = over_alloc_dd(dd->ncg_home);
        srenew(dd->index_gl, dd->cg_nalloc);
        srenew(dd->cgindex, dd->cg_nalloc+1);
    }
    if (bMaster)
    {
        for (i = 0; i < dd->nnodes; i++)
        {
            ma->ibuf[i]            = ma->ncg[i]*sizeof(int);
            ma->ibuf[dd->nnodes+i] = ma->index[i]*sizeof(int);
        }
    }

    dd_scatterv(dd,
                bMaster ? ma->ibuf : nullptr,
                bMaster ? ma->ibuf+dd->nnodes : nullptr,
                bMaster ? ma->cg : nullptr,
                dd->ncg_home*sizeof(int), dd->index_gl);

    /* Determine the home charge group sizes */
    dd->cgindex[0] = 0;
    for (i = 0; i < dd->ncg_home; i++)
    {
        cg_gl            = dd->index_gl[i];
        dd->cgindex[i+1] =
            dd->cgindex[i] + cgs->index[cg_gl+1] - cgs->index[cg_gl];
    }

    if (debug)
    {
        fprintf(debug, "Home charge groups:\n");
        for (i = 0; i < dd->ncg_home; i++)
        {
            fprintf(debug, " %d", dd->index_gl[i]);
            if (i % 10 == 9)
            {
                fprintf(debug, "\n");
            }
        }
        fprintf(debug, "\n");
    }
}

static int compact_and_copy_vec_at(int ncg, int *move,
                                   int *cgindex,
                                   int nvec, int vec,
                                   rvec *src, gmx_domdec_comm_t *comm,
                                   gmx_bool bCompact)
{
    int m, icg, i, i0, i1, nrcg;
    int home_pos;
    int pos_vec[DIM*2];

    home_pos = 0;

    for (m = 0; m < DIM*2; m++)
    {
        pos_vec[m] = 0;
    }

    i0 = 0;
    for (icg = 0; icg < ncg; icg++)
    {
        i1 = cgindex[icg+1];
        m  = move[icg];
        if (m == -1)
        {
            if (bCompact)
            {
                /* Compact the home array in place */
                for (i = i0; i < i1; i++)
                {
                    copy_rvec(src[i], src[home_pos++]);
                }
            }
        }
        else
        {
            /* Copy to the communication buffer */
            nrcg        = i1 - i0;
            pos_vec[m] += 1 + vec*nrcg;
            for (i = i0; i < i1; i++)
            {
                copy_rvec(src[i], comm->cgcm_state[m][pos_vec[m]++]);
            }
            pos_vec[m] += (nvec - vec - 1)*nrcg;
        }
        if (!bCompact)
        {
            home_pos += i1 - i0;
        }
        i0 = i1;
    }

    return home_pos;
}

static int compact_and_copy_vec_cg(int ncg, int *move,
                                   int *cgindex,
                                   int nvec, rvec *src, gmx_domdec_comm_t *comm,
                                   gmx_bool bCompact)
{
    int m, icg, i0, i1, nrcg;
    int home_pos;
    int pos_vec[DIM*2];

    home_pos = 0;

    for (m = 0; m < DIM*2; m++)
    {
        pos_vec[m] = 0;
    }

    i0 = 0;
    for (icg = 0; icg < ncg; icg++)
    {
        i1 = cgindex[icg+1];
        m  = move[icg];
        if (m == -1)
        {
            if (bCompact)
            {
                /* Compact the home array in place */
                copy_rvec(src[icg], src[home_pos++]);
            }
        }
        else
        {
            nrcg = i1 - i0;
            /* Copy to the communication buffer */
            copy_rvec(src[icg], comm->cgcm_state[m][pos_vec[m]]);
            pos_vec[m] += 1 + nrcg*nvec;
        }
        i0 = i1;
    }
    if (!bCompact)
    {
        home_pos = ncg;
    }

    return home_pos;
}

static int compact_ind(int ncg, int *move,
                       int *index_gl, int *cgindex,
                       int *gatindex,
                       gmx_ga2la_t *ga2la, char *bLocalCG,
                       int *cginfo)
{
    int cg, nat, a0, a1, a, a_gl;
    int home_pos;

    home_pos = 0;
    nat      = 0;
    for (cg = 0; cg < ncg; cg++)
    {
        a0 = cgindex[cg];
        a1 = cgindex[cg+1];
        if (move[cg] == -1)
        {
            /* Compact the home arrays in place.
             * Anything that can be done here avoids access to global arrays.
             */
            cgindex[home_pos] = nat;
            for (a = a0; a < a1; a++)
            {
                a_gl          = gatindex[a];
                gatindex[nat] = a_gl;
                /* The cell number stays 0, so we don't need to set it */
                ga2la_change_la(ga2la, a_gl, nat);
                nat++;
            }
            index_gl[home_pos] = index_gl[cg];
            cginfo[home_pos]   = cginfo[cg];
            /* The charge group remains local, so bLocalCG does not change */
            home_pos++;
        }
        else
        {
            /* Clear the global indices */
            for (a = a0; a < a1; a++)
            {
                ga2la_del(ga2la, gatindex[a]);
            }
            if (bLocalCG)
            {
                bLocalCG[index_gl[cg]] = FALSE;
            }
        }
    }
    cgindex[home_pos] = nat;

    return home_pos;
}

static void clear_and_mark_ind(int ncg, int *move,
                               int *index_gl, int *cgindex, int *gatindex,
                               gmx_ga2la_t *ga2la, char *bLocalCG,
                               int *cell_index)
{
    int cg, a0, a1, a;

    for (cg = 0; cg < ncg; cg++)
    {
        if (move[cg] >= 0)
        {
            a0 = cgindex[cg];
            a1 = cgindex[cg+1];
            /* Clear the global indices */
            for (a = a0; a < a1; a++)
            {
                ga2la_del(ga2la, gatindex[a]);
            }
            if (bLocalCG)
            {
                bLocalCG[index_gl[cg]] = FALSE;
            }
            /* Signal that this cg has moved using the ns cell index.
             * Here we set it to -1. fill_grid will change it
             * from -1 to NSGRID_SIGNAL_MOVED_FAC*grid->ncells.
             */
            cell_index[cg] = -1;
        }
    }
}

static void print_cg_move(FILE *fplog,
                          gmx_domdec_t *dd,
                          gmx_int64_t step, int cg, int dim, int dir,
                          gmx_bool bHaveCgcmOld, real limitd,
                          rvec cm_old, rvec cm_new, real pos_d)
{
    gmx_domdec_comm_t *comm;
    char               buf[22];

    comm = dd->comm;

    fprintf(fplog, "\nStep %s:\n", gmx_step_str(step, buf));
    if (limitd > 0)
    {
        fprintf(fplog, "%s %d moved more than the distance allowed by the domain decomposition (%f) in direction %c\n",
                dd->comm->bCGs ? "The charge group starting at atom" : "Atom",
                ddglatnr(dd, dd->cgindex[cg]), limitd, dim2char(dim));
    }
    else
    {
        /* We don't have a limiting distance available: don't print it */
        fprintf(fplog, "%s %d moved more than the distance allowed by the domain decomposition in direction %c\n",
                dd->comm->bCGs ? "The charge group starting at atom" : "Atom",
                ddglatnr(dd, dd->cgindex[cg]), dim2char(dim));
    }
    fprintf(fplog, "distance out of cell %f\n",
            dir == 1 ? pos_d - comm->cell_x1[dim] : pos_d - comm->cell_x0[dim]);
    if (bHaveCgcmOld)
    {
        fprintf(fplog, "Old coordinates: %8.3f %8.3f %8.3f\n",
                cm_old[XX], cm_old[YY], cm_old[ZZ]);
    }
    fprintf(fplog, "New coordinates: %8.3f %8.3f %8.3f\n",
            cm_new[XX], cm_new[YY], cm_new[ZZ]);
    fprintf(fplog, "Old cell boundaries in direction %c: %8.3f %8.3f\n",
            dim2char(dim),
            comm->old_cell_x0[dim], comm->old_cell_x1[dim]);
    fprintf(fplog, "New cell boundaries in direction %c: %8.3f %8.3f\n",
            dim2char(dim),
            comm->cell_x0[dim], comm->cell_x1[dim]);
}

static void cg_move_error(FILE *fplog,
                          gmx_domdec_t *dd,
                          gmx_int64_t step, int cg, int dim, int dir,
                          gmx_bool bHaveCgcmOld, real limitd,
                          rvec cm_old, rvec cm_new, real pos_d)
{
    if (fplog)
    {
        print_cg_move(fplog, dd, step, cg, dim, dir,
                      bHaveCgcmOld, limitd, cm_old, cm_new, pos_d);
    }
    print_cg_move(stderr, dd, step, cg, dim, dir,
                  bHaveCgcmOld, limitd, cm_old, cm_new, pos_d);
    gmx_fatal(FARGS,
              "%s moved too far between two domain decomposition steps\n"
              "This usually means that your system is not well equilibrated",
              dd->comm->bCGs ? "A charge group" : "An atom");
}

static void rotate_state_atom(t_state *state, int a)
{
    if (state->flags & (1 << estX))
    {
        /* Rotate the complete state; for a rectangular box only */
        state->x[a][YY] = state->box[YY][YY] - state->x[a][YY];
        state->x[a][ZZ] = state->box[ZZ][ZZ] - state->x[a][ZZ];
    }
    if (state->flags & (1 << estV))
    {
        state->v[a][YY] = -state->v[a][YY];
        state->v[a][ZZ] = -state->v[a][ZZ];
    }
    if (state->flags & (1 << estCGP))
    {
        state->cg_p[a][YY] = -state->cg_p[a][YY];
        state->cg_p[a][ZZ] = -state->cg_p[a][ZZ];
    }
}

static int *get_moved(gmx_domdec_comm_t *comm, int natoms)
{
    if (natoms > comm->moved_nalloc)
    {
        /* Contents should be preserved here */
        comm->moved_nalloc = over_alloc_dd(natoms);
        srenew(comm->moved, comm->moved_nalloc);
    }

    return comm->moved;
}

static void calc_cg_move(FILE *fplog, gmx_int64_t step,
                         gmx_domdec_t *dd,
                         t_state *state,
                         ivec tric_dir, matrix tcm,
                         rvec cell_x0, rvec cell_x1,
                         rvec limitd, rvec limit0, rvec limit1,
                         const int *cgindex,
                         int cg_start, int cg_end,
                         rvec *cg_cm,
                         int *move)
{
    int      npbcdim;
    int      cg, k, k0, k1, d, dim, d2;
    int      mc, nrcg;
    int      flag;
    gmx_bool bScrew;
    ivec     dev;
    real     inv_ncg, pos_d;
    rvec     cm_new;

    npbcdim = dd->npbcdim;

    for (cg = cg_start; cg < cg_end; cg++)
    {
        k0   = cgindex[cg];
        k1   = cgindex[cg+1];
        nrcg = k1 - k0;
        if (nrcg == 1)
        {
            copy_rvec(state->x[k0], cm_new);
        }
        else
        {
            inv_ncg = 1.0/nrcg;

            clear_rvec(cm_new);
            for (k = k0; (k < k1); k++)
            {
                rvec_inc(cm_new, state->x[k]);
            }
            for (d = 0; (d < DIM); d++)
            {
                cm_new[d] = inv_ncg*cm_new[d];
            }
        }

        clear_ivec(dev);
        /* Do pbc and check DD cell boundary crossings */
        for (d = DIM-1; d >= 0; d--)
        {
            if (dd->nc[d] > 1)
            {
                bScrew = (dd->bScrewPBC && d == XX);
                /* Determine the location of this cg in lattice coordinates */
                pos_d = cm_new[d];
                if (tric_dir[d])
                {
                    for (d2 = d+1; d2 < DIM; d2++)
                    {
                        pos_d += cm_new[d2]*tcm[d2][d];
                    }
                }
                /* Put the charge group in the triclinic unit-cell */
                if (pos_d >= cell_x1[d])
                {
                    if (pos_d >= limit1[d])
                    {
                        cg_move_error(fplog, dd, step, cg, d, 1,
                                      cg_cm != as_rvec_array(state->x.data()), limitd[d],
                                      cg_cm[cg], cm_new, pos_d);
                    }
                    dev[d] = 1;
                    if (dd->ci[d] == dd->nc[d] - 1)
                    {
                        rvec_dec(cm_new, state->box[d]);
                        if (bScrew)
                        {
                            cm_new[YY] = state->box[YY][YY] - cm_new[YY];
                            cm_new[ZZ] = state->box[ZZ][ZZ] - cm_new[ZZ];
                        }
                        for (k = k0; (k < k1); k++)
                        {
                            rvec_dec(state->x[k], state->box[d]);
                            if (bScrew)
                            {
                                rotate_state_atom(state, k);
                            }
                        }
                    }
                }
                else if (pos_d < cell_x0[d])
                {
                    if (pos_d < limit0[d])
                    {
                        cg_move_error(fplog, dd, step, cg, d, -1,
                                      cg_cm != as_rvec_array(state->x.data()), limitd[d],
                                      cg_cm[cg], cm_new, pos_d);
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
                        for (k = k0; (k < k1); k++)
                        {
                            rvec_inc(state->x[k], state->box[d]);
                            if (bScrew)
                            {
                                rotate_state_atom(state, k);
                            }
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
                    for (k = k0; (k < k1); k++)
                    {
                        rvec_dec(state->x[k], state->box[d]);
                    }
                }
                while (cm_new[d] < 0)
                {
                    rvec_inc(cm_new, state->box[d]);
                    for (k = k0; (k < k1); k++)
                    {
                        rvec_inc(state->x[k], state->box[d]);
                    }
                }
            }
        }

        copy_rvec(cm_new, cg_cm[cg]);

        /* Determine where this cg should go */
        flag = 0;
        mc   = -1;
        for (d = 0; d < dd->ndim; d++)
        {
            dim = dd->dim[d];
            if (dev[dim] == 1)
            {
                flag |= DD_FLAG_FW(d);
                if (mc == -1)
                {
                    mc = d*2;
                }
            }
            else if (dev[dim] == -1)
            {
                flag |= DD_FLAG_BW(d);
                if (mc == -1)
                {
                    if (dd->nc[dim] > 2)
                    {
                        mc = d*2 + 1;
                    }
                    else
                    {
                        mc = d*2;
                    }
                }
            }
        }
        /* Temporarily store the flag in move */
        move[cg] = mc + flag;
    }
}

static void dd_redistribute_cg(FILE *fplog, gmx_int64_t step,
                               gmx_domdec_t *dd, ivec tric_dir,
                               t_state *state, PaddedRVecVector *f,
                               t_forcerec *fr,
                               gmx_bool bCompact,
                               t_nrnb *nrnb,
                               int *ncg_stay_home,
                               int *ncg_moved)
{
    int               *move;
    int                npbcdim;
    int                ncg[DIM*2] = { 0 }, nat[DIM*2] = { 0 };
    int                i, cg, k, d, dim, dim2, dir, d2, d3;
    int                mc, cdd, nrcg, ncg_recv, nvs, nvr, nvec, vec;
    int                sbuf[2], rbuf[2];
    int                home_pos_cg, home_pos_at, buf_pos;
    int                flag;
    real               pos_d;
    matrix             tcm;
    rvec              *cg_cm = nullptr, cell_x0, cell_x1, limitd, limit0, limit1;
    int               *cgindex;
    cginfo_mb_t       *cginfo_mb;
    gmx_domdec_comm_t *comm;
    int               *moved;
    int                nthread, thread;

    if (dd->bScrewPBC)
    {
        check_screw_box(state->box);
    }

    comm  = dd->comm;
    if (fr->cutoff_scheme == ecutsGROUP)
    {
        cg_cm = fr->cg_cm;
    }

    // Positions are always present, so there's nothing to flag
    bool bV   = state->flags & (1<<estV);
    bool bCGP = state->flags & (1<<estCGP);

    if (dd->ncg_tot > comm->nalloc_int)
    {
        comm->nalloc_int = over_alloc_dd(dd->ncg_tot);
        srenew(comm->buf_int, comm->nalloc_int);
    }
    move = comm->buf_int;

    npbcdim = dd->npbcdim;

    for (d = 0; (d < DIM); d++)
    {
        limitd[d] = dd->comm->cellsize_min[d];
        if (d >= npbcdim && dd->ci[d] == 0)
        {
            cell_x0[d] = -GMX_FLOAT_MAX;
        }
        else
        {
            cell_x0[d] = comm->cell_x0[d];
        }
        if (d >= npbcdim && dd->ci[d] == dd->nc[d] - 1)
        {
            cell_x1[d] = GMX_FLOAT_MAX;
        }
        else
        {
            cell_x1[d] = comm->cell_x1[d];
        }
        if (d < npbcdim)
        {
            limit0[d] = comm->old_cell_x0[d] - limitd[d];
            limit1[d] = comm->old_cell_x1[d] + limitd[d];
        }
        else
        {
            /* We check after communication if a charge group moved
             * more than one cell. Set the pre-comm check limit to float_max.
             */
            limit0[d] = -GMX_FLOAT_MAX;
            limit1[d] =  GMX_FLOAT_MAX;
        }
    }

    make_tric_corr_matrix(npbcdim, state->box, tcm);

    cgindex = dd->cgindex;

    nthread = gmx_omp_nthreads_get(emntDomdec);

    /* Compute the center of geometry for all home charge groups
     * and put them in the box and determine where they should go.
     */
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        try
        {
            calc_cg_move(fplog, step, dd, state, tric_dir, tcm,
                         cell_x0, cell_x1, limitd, limit0, limit1,
                         cgindex,
                         ( thread   *dd->ncg_home)/nthread,
                         ((thread+1)*dd->ncg_home)/nthread,
                         fr->cutoff_scheme == ecutsGROUP ? cg_cm : as_rvec_array(state->x.data()),
                         move);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    for (cg = 0; cg < dd->ncg_home; cg++)
    {
        if (move[cg] >= 0)
        {
            mc       = move[cg];
            flag     = mc & ~DD_FLAG_NRCG;
            mc       = mc & DD_FLAG_NRCG;
            move[cg] = mc;

            if (ncg[mc]+1 > comm->cggl_flag_nalloc[mc])
            {
                comm->cggl_flag_nalloc[mc] = over_alloc_dd(ncg[mc]+1);
                srenew(comm->cggl_flag[mc], comm->cggl_flag_nalloc[mc]*DD_CGIBS);
            }
            comm->cggl_flag[mc][ncg[mc]*DD_CGIBS  ] = dd->index_gl[cg];
            /* We store the cg size in the lower 16 bits
             * and the place where the charge group should go
             * in the next 6 bits. This saves some communication volume.
             */
            nrcg = cgindex[cg+1] - cgindex[cg];
            comm->cggl_flag[mc][ncg[mc]*DD_CGIBS+1] = nrcg | flag;
            ncg[mc] += 1;
            nat[mc] += nrcg;
        }
    }

    inc_nrnb(nrnb, eNR_CGCM, dd->nat_home);
    inc_nrnb(nrnb, eNR_RESETX, dd->ncg_home);

    *ncg_moved = 0;
    for (i = 0; i < dd->ndim*2; i++)
    {
        *ncg_moved += ncg[i];
    }

    nvec = 1;
    if (bV)
    {
        nvec++;
    }
    if (bCGP)
    {
        nvec++;
    }

    /* Make sure the communication buffers are large enough */
    for (mc = 0; mc < dd->ndim*2; mc++)
    {
        nvr = ncg[mc] + nat[mc]*nvec;
        if (nvr > comm->cgcm_state_nalloc[mc])
        {
            comm->cgcm_state_nalloc[mc] = over_alloc_dd(nvr);
            srenew(comm->cgcm_state[mc], comm->cgcm_state_nalloc[mc]);
        }
    }

    switch (fr->cutoff_scheme)
    {
        case ecutsGROUP:
            /* Recalculating cg_cm might be cheaper than communicating,
             * but that could give rise to rounding issues.
             */
            home_pos_cg =
                compact_and_copy_vec_cg(dd->ncg_home, move, cgindex,
                                        nvec, cg_cm, comm, bCompact);
            break;
        case ecutsVERLET:
            /* Without charge groups we send the moved atom coordinates
             * over twice. This is so the code below can be used without
             * many conditionals for both for with and without charge groups.
             */
            home_pos_cg =
                compact_and_copy_vec_cg(dd->ncg_home, move, cgindex,
                                        nvec, as_rvec_array(state->x.data()), comm, FALSE);
            if (bCompact)
            {
                home_pos_cg -= *ncg_moved;
            }
            break;
        default:
            gmx_incons("unimplemented");
            home_pos_cg = 0;
    }

    vec         = 0;
    home_pos_at =
        compact_and_copy_vec_at(dd->ncg_home, move, cgindex,
                                nvec, vec++, as_rvec_array(state->x.data()),
                                comm, bCompact);
    if (bV)
    {
        compact_and_copy_vec_at(dd->ncg_home, move, cgindex,
                                nvec, vec++, as_rvec_array(state->v.data()),
                                comm, bCompact);
    }
    if (bCGP)
    {
        compact_and_copy_vec_at(dd->ncg_home, move, cgindex,
                                nvec, vec++, as_rvec_array(state->cg_p.data()),
                                comm, bCompact);
    }

    if (bCompact)
    {
        compact_ind(dd->ncg_home, move,
                    dd->index_gl, dd->cgindex, dd->gatindex,
                    dd->ga2la, comm->bLocalCG,
                    fr->cginfo);
    }
    else
    {
        if (fr->cutoff_scheme == ecutsVERLET)
        {
            moved = get_moved(comm, dd->ncg_home);

            for (k = 0; k < dd->ncg_home; k++)
            {
                moved[k] = 0;
            }
        }
        else
        {
            moved = fr->ns->grid->cell_index;
        }

        clear_and_mark_ind(dd->ncg_home, move,
                           dd->index_gl, dd->cgindex, dd->gatindex,
                           dd->ga2la, comm->bLocalCG,
                           moved);
    }

    cginfo_mb = fr->cginfo_mb;

    *ncg_stay_home = home_pos_cg;
    for (d = 0; d < dd->ndim; d++)
    {
        dim      = dd->dim[d];
        ncg_recv = 0;
        nvr      = 0;
        for (dir = 0; dir < (dd->nc[dim] == 2 ? 1 : 2); dir++)
        {
            cdd = d*2 + dir;
            /* Communicate the cg and atom counts */
            sbuf[0] = ncg[cdd];
            sbuf[1] = nat[cdd];
            if (debug)
            {
                fprintf(debug, "Sending ddim %d dir %d: ncg %d nat %d\n",
                        d, dir, sbuf[0], sbuf[1]);
            }
            dd_sendrecv_int(dd, d, dir, sbuf, 2, rbuf, 2);

            if ((ncg_recv+rbuf[0])*DD_CGIBS > comm->nalloc_int)
            {
                comm->nalloc_int = over_alloc_dd((ncg_recv+rbuf[0])*DD_CGIBS);
                srenew(comm->buf_int, comm->nalloc_int);
            }

            /* Communicate the charge group indices, sizes and flags */
            dd_sendrecv_int(dd, d, dir,
                            comm->cggl_flag[cdd], sbuf[0]*DD_CGIBS,
                            comm->buf_int+ncg_recv*DD_CGIBS, rbuf[0]*DD_CGIBS);

            nvs = ncg[cdd] + nat[cdd]*nvec;
            i   = rbuf[0]  + rbuf[1] *nvec;
            vec_rvec_check_alloc(&comm->vbuf, nvr+i);

            /* Communicate cgcm and state */
            dd_sendrecv_rvec(dd, d, dir,
                             comm->cgcm_state[cdd], nvs,
                             comm->vbuf.v+nvr, i);
            ncg_recv += rbuf[0];
            nvr      += i;
        }

        dd_check_alloc_ncg(fr, state, f, home_pos_cg + ncg_recv);
        if (fr->cutoff_scheme == ecutsGROUP)
        {
            /* Here we resize to more than necessary and shrink later */
            dd_resize_state(state, f, home_pos_at + ncg_recv*MAX_CGCGSIZE);
        }

        /* Process the received charge groups */
        buf_pos = 0;
        for (cg = 0; cg < ncg_recv; cg++)
        {
            flag = comm->buf_int[cg*DD_CGIBS+1];

            if (dim >= npbcdim && dd->nc[dim] > 2)
            {
                /* No pbc in this dim and more than one domain boundary.
                 * We do a separate check if a charge group didn't move too far.
                 */
                if (((flag & DD_FLAG_FW(d)) &&
                     comm->vbuf.v[buf_pos][dim] > cell_x1[dim]) ||
                    ((flag & DD_FLAG_BW(d)) &&
                     comm->vbuf.v[buf_pos][dim] < cell_x0[dim]))
                {
                    cg_move_error(fplog, dd, step, cg, dim,
                                  (flag & DD_FLAG_FW(d)) ? 1 : 0,
                                  fr->cutoff_scheme == ecutsGROUP, 0,
                                  comm->vbuf.v[buf_pos],
                                  comm->vbuf.v[buf_pos],
                                  comm->vbuf.v[buf_pos][dim]);
                }
            }

            mc = -1;
            if (d < dd->ndim-1)
            {
                /* Check which direction this cg should go */
                for (d2 = d+1; (d2 < dd->ndim && mc == -1); d2++)
                {
                    if (isDlbOn(dd->comm))
                    {
                        /* The cell boundaries for dimension d2 are not equal
                         * for each cell row of the lower dimension(s),
                         * therefore we might need to redetermine where
                         * this cg should go.
                         */
                        dim2 = dd->dim[d2];
                        /* If this cg crosses the box boundary in dimension d2
                         * we can use the communicated flag, so we do not
                         * have to worry about pbc.
                         */
                        if (!((dd->ci[dim2] == dd->nc[dim2]-1 &&
                               (flag & DD_FLAG_FW(d2))) ||
                              (dd->ci[dim2] == 0 &&
                               (flag & DD_FLAG_BW(d2)))))
                        {
                            /* Clear the two flags for this dimension */
                            flag &= ~(DD_FLAG_FW(d2) | DD_FLAG_BW(d2));
                            /* Determine the location of this cg
                             * in lattice coordinates
                             */
                            pos_d = comm->vbuf.v[buf_pos][dim2];
                            if (tric_dir[dim2])
                            {
                                for (d3 = dim2+1; d3 < DIM; d3++)
                                {
                                    pos_d +=
                                        comm->vbuf.v[buf_pos][d3]*tcm[d3][dim2];
                                }
                            }
                            /* Check of we are not at the box edge.
                             * pbc is only handled in the first step above,
                             * but this check could move over pbc while
                             * the first step did not due to different rounding.
                             */
                            if (pos_d >= cell_x1[dim2] &&
                                dd->ci[dim2] != dd->nc[dim2]-1)
                            {
                                flag |= DD_FLAG_FW(d2);
                            }
                            else if (pos_d < cell_x0[dim2] &&
                                     dd->ci[dim2] != 0)
                            {
                                flag |= DD_FLAG_BW(d2);
                            }
                            comm->buf_int[cg*DD_CGIBS+1] = flag;
                        }
                    }
                    /* Set to which neighboring cell this cg should go */
                    if (flag & DD_FLAG_FW(d2))
                    {
                        mc = d2*2;
                    }
                    else if (flag & DD_FLAG_BW(d2))
                    {
                        if (dd->nc[dd->dim[d2]] > 2)
                        {
                            mc = d2*2+1;
                        }
                        else
                        {
                            mc = d2*2;
                        }
                    }
                }
            }

            nrcg = flag & DD_FLAG_NRCG;
            if (mc == -1)
            {
                if (home_pos_cg+1 > dd->cg_nalloc)
                {
                    dd->cg_nalloc = over_alloc_dd(home_pos_cg+1);
                    srenew(dd->index_gl, dd->cg_nalloc);
                    srenew(dd->cgindex, dd->cg_nalloc+1);
                }
                /* Set the global charge group index and size */
                dd->index_gl[home_pos_cg]  = comm->buf_int[cg*DD_CGIBS];
                dd->cgindex[home_pos_cg+1] = dd->cgindex[home_pos_cg] + nrcg;
                /* Copy the state from the buffer */
                if (fr->cutoff_scheme == ecutsGROUP)
                {
                    cg_cm = fr->cg_cm;
                    copy_rvec(comm->vbuf.v[buf_pos], cg_cm[home_pos_cg]);
                }
                buf_pos++;

                /* Set the cginfo */
                fr->cginfo[home_pos_cg] = ddcginfo(cginfo_mb,
                                                   dd->index_gl[home_pos_cg]);
                if (comm->bLocalCG)
                {
                    comm->bLocalCG[dd->index_gl[home_pos_cg]] = TRUE;
                }

                for (i = 0; i < nrcg; i++)
                {
                    copy_rvec(comm->vbuf.v[buf_pos++],
                              state->x[home_pos_at+i]);
                }
                if (bV)
                {
                    for (i = 0; i < nrcg; i++)
                    {
                        copy_rvec(comm->vbuf.v[buf_pos++],
                                  state->v[home_pos_at+i]);
                    }
                }
                if (bCGP)
                {
                    for (i = 0; i < nrcg; i++)
                    {
                        copy_rvec(comm->vbuf.v[buf_pos++],
                                  state->cg_p[home_pos_at+i]);
                    }
                }
                home_pos_cg += 1;
                home_pos_at += nrcg;
            }
            else
            {
                /* Reallocate the buffers if necessary  */
                if (ncg[mc]+1 > comm->cggl_flag_nalloc[mc])
                {
                    comm->cggl_flag_nalloc[mc] = over_alloc_dd(ncg[mc]+1);
                    srenew(comm->cggl_flag[mc], comm->cggl_flag_nalloc[mc]*DD_CGIBS);
                }
                nvr = ncg[mc] + nat[mc]*nvec;
                if (nvr + 1 + nrcg*nvec > comm->cgcm_state_nalloc[mc])
                {
                    comm->cgcm_state_nalloc[mc] = over_alloc_dd(nvr + 1 + nrcg*nvec);
                    srenew(comm->cgcm_state[mc], comm->cgcm_state_nalloc[mc]);
                }
                /* Copy from the receive to the send buffers */
                memcpy(comm->cggl_flag[mc] + ncg[mc]*DD_CGIBS,
                       comm->buf_int + cg*DD_CGIBS,
                       DD_CGIBS*sizeof(int));
                memcpy(comm->cgcm_state[mc][nvr],
                       comm->vbuf.v[buf_pos],
                       (1+nrcg*nvec)*sizeof(rvec));
                buf_pos += 1 + nrcg*nvec;
                ncg[mc] += 1;
                nat[mc] += nrcg;
            }
        }
    }

    /* With sorting (!bCompact) the indices are now only partially up to date
     * and ncg_home and nat_home are not the real count, since there are
     * "holes" in the arrays for the charge groups that moved to neighbors.
     */
    if (fr->cutoff_scheme == ecutsVERLET)
    {
        moved = get_moved(comm, home_pos_cg);

        for (i = dd->ncg_home; i < home_pos_cg; i++)
        {
            moved[i] = 0;
        }
    }
    dd->ncg_home = home_pos_cg;
    dd->nat_home = home_pos_at;

    if (fr->cutoff_scheme == ecutsGROUP && !bCompact)
    {
        /* We overallocated before, we need to set the right size here */
        dd_resize_state(state, f, dd->nat_home);
    }

    if (debug)
    {
        fprintf(debug,
                "Finished repartitioning: cgs moved out %d, new home %d\n",
                *ncg_moved, dd->ncg_home-*ncg_moved);

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
    domdec_root_t     *root = nullptr;
    int                d, dim, i, pos;
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

    for (d = dd->ndim-1; d >= 0; d--)
    {
        dim = dd->dim[d];
        /* Check if we participate in the communication in this dimension */
        if (d == dd->ndim-1 ||
            (dd->ci[dd->dim[d+1]] == 0 && dd->ci[dd->dim[dd->ndim-1]] == 0))
        {
            load = &comm->load[d];
            if (isDlbOn(dd->comm))
            {
                cell_frac = comm->cell_f1[d] - comm->cell_f0[d];
            }
            pos = 0;
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
                        sbuf[pos++] = comm->cell_f_max0[d];
                        sbuf[pos++] = comm->cell_f_min1[d];
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
                        sbuf[pos++] = comm->cell_f_max0[d];
                        sbuf[pos++] = comm->cell_f_min1[d];
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
                /* We are the root, process this row */
                if (isDlbOn(comm))
                {
                    root = comm->root[d];
                }
                load->sum      = 0;
                load->max      = 0;
                load->sum_m    = 0;
                load->cvol_min = 1;
                load->flags    = 0;
                load->mdf      = 0;
                load->pme      = 0;
                pos            = 0;
                for (i = 0; i < dd->nc[dim]; i++)
                {
                    load->sum += load->load[pos++];
                    load->max  = std::max(load->max, load->load[pos]);
                    pos++;
                    if (isDlbOn(dd->comm))
                    {
                        if (root->bLimited)
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
                            root->cell_f_max0[i] = load->load[pos++];
                            root->cell_f_min1[i] = load->load[pos++];
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
                if (isDlbOn(comm) && root->bLimited)
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
            for (d = 0; d < dd->ndim; d++)
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
        std::string dlbStateStr = "";

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
        sprintf(buf, " Part of the total run time spent waiting due to PP/PME imbalance: %.1f %%\n", fabs(lossFractionPme)*100);
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
    if (numPmeRanks > 0 && fabs(lossFractionPme) >= DD_PERF_LOSS_WARN)
    {
        sprintf(buf,
                "NOTE: %.1f %% performance was lost because the PME ranks\n"
                "      had %s work to do than the PP ranks.\n"
                "      You might want to %s the number of PME ranks\n"
                "      or %s the cut-off and the grid spacing.\n",
                fabs(lossFractionPme*100),
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

static void dd_print_load(FILE *fplog, gmx_domdec_t *dd, gmx_int64_t step)
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
    domdec_root_t     *root;
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
            if (dd->ci[dim] == dd->master_ci[dim])
            {
                /* This is the root process of this row */
                snew(dd->comm->root[dim_ind], 1);
                root = dd->comm->root[dim_ind];
                snew(root->cell_f, DD_CELL_F_SIZE(dd, dim_ind));
                snew(root->old_cell_f, dd->nc[dim]+1);
                snew(root->bCellMin, dd->nc[dim]);
                if (dim_ind > 0)
                {
                    snew(root->cell_f_max0, dd->nc[dim]);
                    snew(root->cell_f_min1, dd->nc[dim]);
                    snew(root->bound_min, dd->nc[dim]);
                    snew(root->bound_max, dd->nc[dim]);
                }
                snew(root->buf_ncd, dd->nc[dim]);
            }
            else
            {
                /* This is not a root process, we only need to receive cell_f */
                snew(dd->comm->cell_f_row, DD_CELL_F_SIZE(dd, dim_ind));
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

    if (!(cr->duty & DUTY_PP) || gpu_id < 0)
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
     * for the nodecomm summution.
     */
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
        snew(dd->comm->root, dd->ndim);
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
        if (cr->duty & DUTY_PP)
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
        if (cr->duty & DUTY_PP)
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

static gmx_domdec_master_t *init_gmx_domdec_master_t(gmx_domdec_t *dd,
                                                     int ncg, int natoms)
{
    gmx_domdec_master_t *ma;
    int                  i;

    snew(ma, 1);

    snew(ma->ncg, dd->nnodes);
    snew(ma->index, dd->nnodes+1);
    snew(ma->cg, ncg);
    snew(ma->nat, dd->nnodes);
    snew(ma->ibuf, dd->nnodes*2);
    snew(ma->cell_x, DIM);
    for (i = 0; i < DIM; i++)
    {
        snew(ma->cell_x[i], dd->nc[i]+1);
    }

    if (dd->nnodes <= GMX_DD_NNODES_SENDRECV)
    {
        ma->vbuf = nullptr;
    }
    else
    {
        snew(ma->vbuf, natoms);
    }

    return ma;
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

#if GMX_MPI
    if (comm->bCartesianPP_PME)
    {
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
                       cr->duty,
                       dd_index(comm->ntot, dd->ci),
                       &cr->mpi_comm_mygroup);
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

        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       cr->duty,
                       cr->nodeid,
                       &cr->mpi_comm_mygroup);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
    }
#endif

    if (fplog)
    {
        fprintf(fplog, "This rank does only %s work.\n\n",
                (cr->duty & DUTY_PP) ? "particle-particle" : "PME-mesh");
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

    if (cr->duty & DUTY_PP)
    {
        /* Copy or make a new PP communicator */
        make_pp_communicator(fplog, dd, cr, CartReorder);
    }
    else
    {
        receive_ddindex2simnodeid(dd, cr);
    }

    if (!(cr->duty & DUTY_PME))
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
        dd->ma = init_gmx_domdec_master_t(dd,
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
    t_ilist             *il;

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

static void dd_warning(t_commrec *cr, FILE *fplog, const char *warn_string)
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
        gmx_fatal(FARGS, (dlbNotSupportedErr + reasonStr).c_str());
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
                break;
            case edlbsOnCanTurnOff:
                GMX_RELEASE_ASSERT(false, "edlbsOffCanTurnOff is not a valid initial state");
                break;
            case edlbsOnUser:
                return forceDlbOffOrBail(dlbState, reasonStr + " In load balanced runs binary reproducibility cannot be ensured.", cr, fplog);
                break;
            default:
                gmx_fatal(FARGS, "Death horror: undefined case (%d) for load balancing choice", dlbState);
                break;
        }
    }

    return dlbState;
}

static void set_dd_dim(FILE *fplog, gmx_domdec_t *dd)
{
    int dim;

    dd->ndim = 0;
    if (getenv("GMX_DD_ORDER_ZYX") != nullptr)
    {
        /* Decomposition order z,y,x */
        if (fplog)
        {
            fprintf(fplog, "Using domain decomposition order z, y, x\n");
        }
        for (dim = DIM-1; dim >= 0; dim--)
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
        for (dim = 0; dim < DIM; dim++)
        {
            if (dd->nc[dim] > 1)
            {
                dd->dim[dd->ndim++] = dim;
            }
        }
    }
}

static gmx_domdec_comm_t *init_dd_comm()
{
    gmx_domdec_comm_t *comm;
    int                i;

    snew(comm, 1);
    snew(comm->cggl_flag, DIM*2);
    snew(comm->cgcm_state, DIM*2);
    for (i = 0; i < DIM*2; i++)
    {
        comm->cggl_flag_nalloc[i]  = 0;
        comm->cgcm_state_nalloc[i] = 0;
    }

    comm->nalloc_int = 0;
    comm->buf_int    = nullptr;

    vec_rvec_init(&comm->vbuf);

    comm->n_load_have    = 0;
    comm->n_load_collect = 0;

    for (i = 0; i < ddnatNR-ddnatZONE; i++)
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
                                   const matrix box, const rvec *xGlobal,
                                   gmx_ddbox_t *ddbox,
                                   int *npme_x, int *npme_y)
{
    real               r_bonded         = -1;
    real               r_bonded_limit   = -1;
    const real         tenPercentMargin = 1.1;
    gmx_domdec_comm_t *comm             = dd->comm;

    snew(comm->cggl_flag, DIM*2);
    snew(comm->cgcm_state, DIM*2);

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
    snew(comm->sort, 1);

    comm->bCGs = (ncg_mtop(mtop) < mtop->natoms);

    comm->bInterCGBondeds = ((ncg_mtop(mtop) > mtop->mols.nr) ||
                             mtop->bIntermolecularInteractions);
    if (comm->bInterCGBondeds)
    {
        comm->bInterCGMultiBody = (multi_body_bondeds_count(mtop) > 0);
    }
    else
    {
        comm->bInterCGMultiBody = FALSE;
    }

    dd->bInterCGcons    = inter_charge_group_constraints(mtop);
    dd->bInterCGsettles = inter_charge_group_settles(mtop);

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
                dd_bonded_cg_distance(fplog, mtop, ir, xGlobal, box,
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
        rconstr = constr_r_max(fplog, mtop, ir);
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
        set_ddbox_cr(cr, &dd->nc, ir, box, &comm->cgs_gl, xGlobal, ddbox);

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
        set_ddbox_cr(cr, nullptr, ir, box, &comm->cgs_gl, xGlobal, ddbox);

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

    /* Technically we don't need both of these,
     * but it simplifies code not having to recalculate it.
     */
    *npme_x = comm->npmenodes_x;
    *npme_y = comm->npmenodes_y;

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
    int d;

    for (d = 0; d < dd->ndim; d++)
    {
        dd->comm->cd[d].np                 = dd->comm->cd[d].np_dlb;
        dd->comm->cellsize_min[dd->dim[d]] =
            dd->comm->cellsize_min_dlb[dd->dim[d]];
    }
}


static void turn_on_dlb(FILE *fplog, t_commrec *cr, gmx_int64_t step)
{
    gmx_domdec_t      *dd;
    gmx_domdec_comm_t *comm;
    real               cellsize_min;
    int                d, nc, i;

    dd   = cr->dd;
    comm = dd->comm;

    cellsize_min = comm->cellsize_min[dd->dim[0]];
    for (d = 1; d < dd->ndim; d++)
    {
        cellsize_min = std::min(cellsize_min, comm->cellsize_min[dd->dim[d]]);
    }

    /* Turn off DLB if we're too close to the cell size limit. */
    if (cellsize_min < comm->cellsize_limit*1.05)
    {
        auto str = gmx::formatString("step %" GMX_PRId64 " Measured %.1f %% performance loss due to load imbalance, "
                                     "but the minimum cell size is smaller than 1.05 times the cell size limit."
                                     "Will no longer try dynamic load balancing.\n", step, dd_force_imb_perf_loss(dd)*100);
        dd_warning(cr, fplog, str.c_str());

        comm->dlbState = edlbsOffForever;
        return;
    }

    char buf[STRLEN];
    sprintf(buf, "step %" GMX_PRId64 " Turning on dynamic load balancing, because the performance loss due to load imbalance is %.1f %%.\n", step, dd_force_imb_perf_loss(dd)*100);
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
    for (d = 0; d < dd->ndim; d++)
    {
        if (comm->root[d])
        {
            comm->load[d].sum_m = comm->load[d].sum;

            nc = dd->nc[dd->dim[d]];
            for (i = 0; i < nc; i++)
            {
                comm->root[d]->cell_f[i]    = i/(real)nc;
                if (d > 0)
                {
                    comm->root[d]->cell_f_max0[i] =  i   /(real)nc;
                    comm->root[d]->cell_f_min1[i] = (i+1)/(real)nc;
                }
            }
            comm->root[d]->cell_f[nc] = 1.0;
        }
    }
}

static void turn_off_dlb(FILE *fplog, t_commrec *cr, gmx_int64_t step)
{
    gmx_domdec_t *dd = cr->dd;

    char          buf[STRLEN];
    sprintf(buf, "step %" GMX_PRId64 " Turning off dynamic load balancing, because it is degrading performance.\n", step);
    dd_warning(cr, fplog, buf);
    dd->comm->dlbState                     = edlbsOffCanTurnOn;
    dd->comm->haveTurnedOffDlb             = true;
    dd->comm->ddPartioningCountFirstDlbOff = dd->ddp_count;
}

static void turn_off_dlb_forever(FILE *fplog, t_commrec *cr, gmx_int64_t step)
{
    GMX_RELEASE_ASSERT(cr->dd->comm->dlbState == edlbsOffCanTurnOn, "Can only turn off DLB forever when it was in the can-turn-on state");
    char buf[STRLEN];
    sprintf(buf, "step %" GMX_PRId64 " Will no longer try dynamic load balancing, as it degraded performance.\n", step);
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
        comm->cd[d].np_nalloc = comm->cd[d].np_dlb;
        snew(comm->cd[d].ind, comm->cd[d].np_nalloc);
        comm->maxpulse = std::max(comm->maxpulse, comm->cd[d].np_dlb);
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
                                        const rvec *xGlobal,
                                        gmx_ddbox_t *ddbox,
                                        int *npme_x, int *npme_y)
{
    gmx_domdec_t      *dd;

    if (fplog)
    {
        fprintf(fplog,
                "\nInitializing Domain Decomposition on %d ranks\n", cr->nnodes);
    }

    snew(dd, 1);

    dd->comm = init_dd_comm();

    set_dd_envvar_options(fplog, dd, cr->nodeid);

    set_dd_limits_and_grid(fplog, cr, dd, options, mdrunOptions,
                           mtop, ir,
                           box, xGlobal,
                           ddbox,
                           npme_x, npme_y);

    make_dd_communicators(fplog, cr, dd, options.rankOrder);

    if (cr->duty & DUTY_PP)
    {
        set_ddgrid_parameters(fplog, dd, options.dlbScaling, mtop, ir, ddbox);

        setup_neighbor_relations(dd);
    }

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);

    /* Initialize DD paritioning counters */
    dd->comm->partition_step = INT_MIN;
    dd->ddp_count            = 0;

    /* We currently don't know the number of threads yet, we set this later */
    dd->comm->nth = 0;

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

    set_ddbox(dd, FALSE, cr, ir, state->box,
              TRUE, &dd->comm->cgs_gl, as_rvec_array(state->x.data()), &ddbox);

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
                             int  *index_gl, int  *recv_i,
                             rvec *cg_cm,    rvec *recv_vr,
                             int *cgindex,
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

static void make_cell2at_index(gmx_domdec_comm_dim_t *cd,
                               int nzone, int cg0, const int *cgindex)
{
    int cg, zone, p;

    /* Store the atom block boundaries for easy copying of communication buffers
     */
    cg = cg0;
    for (zone = 0; zone < nzone; zone++)
    {
        for (p = 0; p < cd->np; p++)
        {
            cd->ind[p].cell2at0[zone] = cgindex[cg];
            cg += cd->ind[p].nrecv[zone];
            cd->ind[p].cell2at1[zone] = cgindex[cg];
        }
    }
}

static gmx_bool missing_link(t_blocka *link, int cg_gl, char *bLocalCG)
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

/* Determine which cg's we need to send in this pulse from this zone */
static void
get_zone_pulse_cgs(gmx_domdec_t *dd,
                   int zonei, int zone,
                   int cg0, int cg1,
                   const int *index_gl,
                   const int *cgindex,
                   int dim, int dim_ind,
                   int dim0, int dim1, int dim2,
                   real r_comm2, real r_bcomm2,
                   matrix box,
                   ivec tric_dist,
                   rvec *normal,
                   real skew_fac2_d, real skew_fac_01,
                   rvec *v_d, rvec *v_0, rvec *v_1,
                   const dd_corners_t *c,
                   rvec sf2_round,
                   gmx_bool bDistBonded,
                   gmx_bool bBondComm,
                   gmx_bool bDist2B,
                   gmx_bool bDistMB,
                   rvec *cg_cm,
                   int *cginfo,
                   gmx_domdec_ind_t *ind,
                   int **ibuf, int *ibuf_nalloc,
                   vec_rvec_t *vbuf,
                   int *nsend_ptr,
                   int *nat_ptr,
                   int *nsend_z_ptr)
{
    gmx_domdec_comm_t *comm;
    gmx_bool           bScrew;
    gmx_bool           bDistMB_pulse;
    int                cg, i;
    real               r2, rb2, r, tric_sh;
    rvec               rn, rb;
    int                dimd;
    int                nsend_z, nsend, nat;

    comm = dd->comm;

    bScrew = (dd->bScrewPBC && dim == XX);

    bDistMB_pulse = (bDistMB && bDistBonded);

    nsend_z = 0;
    nsend   = *nsend_ptr;
    nat     = *nat_ptr;

    for (cg = cg0; cg < cg1; cg++)
    {
        r2  = 0;
        rb2 = 0;
        if (tric_dist[dim_ind] == 0)
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
               missing_link(comm->cglink, index_gl[cg],
                            comm->bLocalCG)))))
        {
            /* Make an index to the local charge groups */
            if (nsend+1 > ind->nalloc)
            {
                ind->nalloc = over_alloc_large(nsend+1);
                srenew(ind->index, ind->nalloc);
            }
            if (nsend+1 > *ibuf_nalloc)
            {
                *ibuf_nalloc = over_alloc_large(nsend+1);
                srenew(*ibuf, *ibuf_nalloc);
            }
            ind->index[nsend] = cg;
            (*ibuf)[nsend]    = index_gl[cg];
            nsend_z++;
            vec_rvec_check_alloc(vbuf, nsend+1);

            if (dd->ci[dim] == 0)
            {
                /* Correct cg_cm for pbc */
                rvec_add(cg_cm[cg], box[dim], vbuf->v[nsend]);
                if (bScrew)
                {
                    vbuf->v[nsend][YY] = box[YY][YY] - vbuf->v[nsend][YY];
                    vbuf->v[nsend][ZZ] = box[ZZ][ZZ] - vbuf->v[nsend][ZZ];
                }
            }
            else
            {
                copy_rvec(cg_cm[cg], vbuf->v[nsend]);
            }
            nsend++;
            nat += cgindex[cg+1] - cgindex[cg];
        }
    }

    *nsend_ptr   = nsend;
    *nat_ptr     = nat;
    *nsend_z_ptr = nsend_z;
}

static void setup_dd_communication(gmx_domdec_t *dd,
                                   matrix box, gmx_ddbox_t *ddbox,
                                   t_forcerec *fr,
                                   t_state *state, PaddedRVecVector *f)
{
    int                    dim_ind, dim, dim0, dim1, dim2, dimd, p, nat_tot;
    int                    nzone, nzone_send, zone, zonei, cg0, cg1;
    int                    c, i, cg, cg_gl, nrcg;
    int                   *zone_cg_range, pos_cg, *index_gl, *cgindex, *recv_i;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_zones_t    *zones;
    gmx_domdec_comm_dim_t *cd;
    gmx_domdec_ind_t      *ind;
    cginfo_mb_t           *cginfo_mb;
    gmx_bool               bBondComm, bDist2B, bDistMB, bDistBonded;
    real                   r_comm2, r_bcomm2;
    dd_corners_t           corners;
    ivec                   tric_dist;
    rvec                  *cg_cm, *normal, *v_d, *v_0 = nullptr, *v_1 = nullptr, *recv_vr;
    real                   skew_fac2_d, skew_fac_01;
    rvec                   sf2_round;
    int                    nsend, nat;
    int                    th;

    if (debug)
    {
        fprintf(debug, "Setting up DD communication\n");
    }

    comm  = dd->comm;

    if (comm->nth == 0)
    {
        /* Initialize the thread data.
         * This can not be done in init_domain_decomposition,
         * as the numbers of threads is determined later.
         */
        comm->nth = gmx_omp_nthreads_get(emntDomdec);
        if (comm->nth > 1)
        {
            snew(comm->dth, comm->nth);
        }
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
            cg_cm = nullptr;
    }

    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        /* Check if we need to use triclinic distances */
        tric_dist[dim_ind] = 0;
        for (i = 0; i <= dim_ind; i++)
        {
            if (ddbox->tric_dir[dd->dim[i]])
            {
                tric_dist[dim_ind] = 1;
            }
        }
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
    index_gl      = dd->index_gl;
    cgindex       = dd->cgindex;
    cginfo_mb     = fr->cginfo_mb;

    zone_cg_range[0]   = 0;
    zone_cg_range[1]   = dd->ncg_home;
    comm->zone_ncg1[0] = dd->ncg_home;
    pos_cg             = dd->ncg_home;

    nat_tot = dd->nat_home;
    nzone   = 1;
    for (dim_ind = 0; dim_ind < dd->ndim; dim_ind++)
    {
        dim = dd->dim[dim_ind];
        cd  = &comm->cd[dim_ind];

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

        cd->bInPlace = TRUE;
        for (p = 0; p < cd->np; p++)
        {
            /* Only atoms communicated in the first pulse are used
             * for multi-body bonded interactions or for bBondComm.
             */
            bDistBonded = ((bDistMB || bDist2B) && p == 0);

            ind   = &cd->ind[p];
            nsend = 0;
            nat   = 0;
            for (zone = 0; zone < nzone_send; zone++)
            {
                if (tric_dist[dim_ind] && dim_ind > 0)
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

#pragma omp parallel for num_threads(comm->nth) schedule(static)
                for (th = 0; th < comm->nth; th++)
                {
                    try
                    {
                        gmx_domdec_ind_t *ind_p;
                        int             **ibuf_p, *ibuf_nalloc_p;
                        vec_rvec_t       *vbuf_p;
                        int              *nsend_p, *nat_p;
                        int              *nsend_zone_p;
                        int               cg0_th, cg1_th;

                        if (th == 0)
                        {
                            /* Thread 0 writes in the comm buffers */
                            ind_p         = ind;
                            ibuf_p        = &comm->buf_int;
                            ibuf_nalloc_p = &comm->nalloc_int;
                            vbuf_p        = &comm->vbuf;
                            nsend_p       = &nsend;
                            nat_p         = &nat;
                            nsend_zone_p  = &ind->nsend[zone];
                        }
                        else
                        {
                            /* Other threads write into temp buffers */
                            ind_p         = &comm->dth[th].ind;
                            ibuf_p        = &comm->dth[th].ibuf;
                            ibuf_nalloc_p = &comm->dth[th].ibuf_nalloc;
                            vbuf_p        = &comm->dth[th].vbuf;
                            nsend_p       = &comm->dth[th].nsend;
                            nat_p         = &comm->dth[th].nat;
                            nsend_zone_p  = &comm->dth[th].nsend_zone;

                            comm->dth[th].nsend      = 0;
                            comm->dth[th].nat        = 0;
                            comm->dth[th].nsend_zone = 0;
                        }

                        if (comm->nth == 1)
                        {
                            cg0_th = cg0;
                            cg1_th = cg1;
                        }
                        else
                        {
                            cg0_th = cg0 + ((cg1 - cg0)* th   )/comm->nth;
                            cg1_th = cg0 + ((cg1 - cg0)*(th+1))/comm->nth;
                        }

                        /* Get the cg's for this pulse in this zone */
                        get_zone_pulse_cgs(dd, zonei, zone, cg0_th, cg1_th,
                                           index_gl, cgindex,
                                           dim, dim_ind, dim0, dim1, dim2,
                                           r_comm2, r_bcomm2,
                                           box, tric_dist,
                                           normal, skew_fac2_d, skew_fac_01,
                                           v_d, v_0, v_1, &corners, sf2_round,
                                           bDistBonded, bBondComm,
                                           bDist2B, bDistMB,
                                           cg_cm, fr->cginfo,
                                           ind_p,
                                           ibuf_p, ibuf_nalloc_p,
                                           vbuf_p,
                                           nsend_p, nat_p,
                                           nsend_zone_p);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                } // END

                /* Append data of threads>=1 to the communication buffers */
                for (th = 1; th < comm->nth; th++)
                {
                    dd_comm_setup_work_t *dth;
                    int                   i, ns1;

                    dth = &comm->dth[th];

                    ns1 = nsend + dth->nsend_zone;
                    if (ns1 > ind->nalloc)
                    {
                        ind->nalloc = over_alloc_dd(ns1);
                        srenew(ind->index, ind->nalloc);
                    }
                    if (ns1 > comm->nalloc_int)
                    {
                        comm->nalloc_int = over_alloc_dd(ns1);
                        srenew(comm->buf_int, comm->nalloc_int);
                    }
                    if (ns1 > comm->vbuf.nalloc)
                    {
                        comm->vbuf.nalloc = over_alloc_dd(ns1);
                        srenew(comm->vbuf.v, comm->vbuf.nalloc);
                    }

                    for (i = 0; i < dth->nsend_zone; i++)
                    {
                        ind->index[nsend]    = dth->ind.index[i];
                        comm->buf_int[nsend] = dth->ibuf[i];
                        copy_rvec(dth->vbuf.v[i],
                                  comm->vbuf.v[nsend]);
                        nsend++;
                    }
                    nat              += dth->nat;
                    ind->nsend[zone] += dth->nsend_zone;
                }
            }
            /* Clear the counts in case we do not have pbc */
            for (zone = nzone_send; zone < nzone; zone++)
            {
                ind->nsend[zone] = 0;
            }
            ind->nsend[nzone]   = nsend;
            ind->nsend[nzone+1] = nat;
            /* Communicate the number of cg's and atoms to receive */
            dd_sendrecv_int(dd, dim_ind, dddirBackward,
                            ind->nsend, nzone+2,
                            ind->nrecv, nzone+2);

            /* The rvec buffer is also required for atom buffers of size nsend
             * in dd_move_x and dd_move_f.
             */
            vec_rvec_check_alloc(&comm->vbuf, ind->nsend[nzone+1]);

            if (p > 0)
            {
                /* We can receive in place if only the last zone is not empty */
                for (zone = 0; zone < nzone-1; zone++)
                {
                    if (ind->nrecv[zone] > 0)
                    {
                        cd->bInPlace = FALSE;
                    }
                }
                if (!cd->bInPlace)
                {
                    /* The int buffer is only required here for the cg indices */
                    if (ind->nrecv[nzone] > comm->nalloc_int2)
                    {
                        comm->nalloc_int2 = over_alloc_dd(ind->nrecv[nzone]);
                        srenew(comm->buf_int2, comm->nalloc_int2);
                    }
                    /* The rvec buffer is also required for atom buffers
                     * of size nrecv in dd_move_x and dd_move_f.
                     */
                    i = std::max(cd->ind[0].nrecv[nzone+1], ind->nrecv[nzone+1]);
                    vec_rvec_check_alloc(&comm->vbuf2, i);
                }
            }

            /* Make space for the global cg indices */
            if (pos_cg + ind->nrecv[nzone] > dd->cg_nalloc
                || dd->cg_nalloc == 0)
            {
                dd->cg_nalloc = over_alloc_dd(pos_cg + ind->nrecv[nzone]);
                srenew(index_gl, dd->cg_nalloc);
                srenew(cgindex, dd->cg_nalloc+1);
            }
            /* Communicate the global cg indices */
            if (cd->bInPlace)
            {
                recv_i = index_gl + pos_cg;
            }
            else
            {
                recv_i = comm->buf_int2;
            }
            dd_sendrecv_int(dd, dim_ind, dddirBackward,
                            comm->buf_int, nsend,
                            recv_i,        ind->nrecv[nzone]);

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
            if (cd->bInPlace)
            {
                recv_vr = cg_cm + pos_cg;
            }
            else
            {
                recv_vr = comm->vbuf2.v;
            }
            dd_sendrecv_rvec(dd, dim_ind, dddirBackward,
                             comm->vbuf.v, nsend,
                             recv_vr,      ind->nrecv[nzone]);

            /* Make the charge group index */
            if (cd->bInPlace)
            {
                zone = (p == 0 ? 0 : nzone - 1);
                while (zone < nzone)
                {
                    for (cg = 0; cg < ind->nrecv[zone]; cg++)
                    {
                        cg_gl              = index_gl[pos_cg];
                        fr->cginfo[pos_cg] = ddcginfo(cginfo_mb, cg_gl);
                        nrcg               = GET_CGINFO_NATOMS(fr->cginfo[pos_cg]);
                        cgindex[pos_cg+1]  = cgindex[pos_cg] + nrcg;
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
                merge_cg_buffers(nzone, cd, p, zone_cg_range,
                                 index_gl, recv_i, cg_cm, recv_vr,
                                 cgindex, fr->cginfo_mb, fr->cginfo);
                pos_cg += ind->nrecv[nzone];
            }
            nat_tot += ind->nrecv[nzone+1];
        }
        if (!cd->bInPlace)
        {
            /* Store the atom block for easy copying of communication buffers */
            make_cell2at_index(cd, nzone, zone_cg_range[nzone], cgindex);
        }
        nzone += nzone;
    }
    dd->index_gl = index_gl;
    dd->cgindex  = cgindex;

    dd->ncg_tot          = zone_cg_range[zones->n];
    dd->nat_tot          = nat_tot;
    comm->nat[ddnatHOME] = dd->nat_home;
    for (i = ddnatZONE; i < ddnatNR; i++)
    {
        comm->nat[i] = dd->nat_tot;
    }

    if (!bBondComm)
    {
        /* We don't need to update cginfo, since that was alrady done above.
         * So we pass NULL for the forcerec.
         */
        dd_set_cginfo(dd->index_gl, dd->ncg_home, dd->ncg_tot,
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

static void set_zones_size(gmx_domdec_t *dd,
                           matrix box, const gmx_ddbox_t *ddbox,
                           int zone_start, int zone_end)
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
                for (z = zones->izone[zi].j0; z < zones->izone[zi].j1; z++)
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
        zones->dens_zone0 = (zones->cg_range[1] - zones->cg_range[0])/vol;
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

static void order_int_cg(int n, const gmx_cgsort_t *sort,
                         int *a, int *buf)
{
    int i;

    /* Order the data */
    for (i = 0; i < n; i++)
    {
        buf[i] = a[sort[i].ind];
    }

    /* Copy back to the original array */
    for (i = 0; i < n; i++)
    {
        a[i] = buf[i];
    }
}

static void order_vec_cg(int n, const gmx_cgsort_t *sort,
                         rvec *v, rvec *buf)
{
    int i;

    /* Order the data */
    for (i = 0; i < n; i++)
    {
        copy_rvec(v[sort[i].ind], buf[i]);
    }

    /* Copy back to the original array */
    for (i = 0; i < n; i++)
    {
        copy_rvec(buf[i], v[i]);
    }
}

static void order_vec_atom(int ncg, const int *cgindex, const gmx_cgsort_t *sort,
                           rvec *v, rvec *buf)
{
    int a, atot, cg, cg0, cg1, i;

    if (cgindex == nullptr)
    {
        /* Avoid the useless loop of the atoms within a cg */
        order_vec_cg(ncg, sort, v, buf);

        return;
    }

    /* Order the data */
    a = 0;
    for (cg = 0; cg < ncg; cg++)
    {
        cg0 = cgindex[sort[cg].ind];
        cg1 = cgindex[sort[cg].ind+1];
        for (i = cg0; i < cg1; i++)
        {
            copy_rvec(v[i], buf[a]);
            a++;
        }
    }
    atot = a;

    /* Copy back to the original array */
    for (a = 0; a < atot; a++)
    {
        copy_rvec(buf[a], v[a]);
    }
}

static void ordered_sort(int nsort2, gmx_cgsort_t *sort2,
                         int nsort_new, gmx_cgsort_t *sort_new,
                         gmx_cgsort_t *sort1)
{
    int i1, i2, i_new;

    /* The new indices are not very ordered, so we qsort them */
    gmx_qsort_threadsafe(sort_new, nsort_new, sizeof(sort_new[0]), comp_cgsort);

    /* sort2 is already ordered, so now we can merge the two arrays */
    i1    = 0;
    i2    = 0;
    i_new = 0;
    while (i2 < nsort2 || i_new < nsort_new)
    {
        if (i2 == nsort2)
        {
            sort1[i1++] = sort_new[i_new++];
        }
        else if (i_new == nsort_new)
        {
            sort1[i1++] = sort2[i2++];
        }
        else if (sort2[i2].nsc < sort_new[i_new].nsc ||
                 (sort2[i2].nsc == sort_new[i_new].nsc &&
                  sort2[i2].ind_gl < sort_new[i_new].ind_gl))
        {
            sort1[i1++] = sort2[i2++];
        }
        else
        {
            sort1[i1++] = sort_new[i_new++];
        }
    }
}

static int dd_sort_order(gmx_domdec_t *dd, t_forcerec *fr, int ncg_home_old)
{
    gmx_domdec_sort_t *sort;
    gmx_cgsort_t      *cgsort, *sort_i;
    int                ncg_new, nsort2, nsort_new, i, *a, moved;

    sort = dd->comm->sort;

    a = fr->ns->grid->cell_index;

    moved = NSGRID_SIGNAL_MOVED_FAC*fr->ns->grid->ncells;

    if (ncg_home_old >= 0)
    {
        /* The charge groups that remained in the same ns grid cell
         * are completely ordered. So we can sort efficiently by sorting
         * the charge groups that did move into the stationary list.
         */
        ncg_new   = 0;
        nsort2    = 0;
        nsort_new = 0;
        for (i = 0; i < dd->ncg_home; i++)
        {
            /* Check if this cg did not move to another node */
            if (a[i] < moved)
            {
                if (i >= ncg_home_old || a[i] != sort->sort[i].nsc)
                {
                    /* This cg is new on this node or moved ns grid cell */
                    if (nsort_new >= sort->sort_new_nalloc)
                    {
                        sort->sort_new_nalloc = over_alloc_dd(nsort_new+1);
                        srenew(sort->sort_new, sort->sort_new_nalloc);
                    }
                    sort_i = &(sort->sort_new[nsort_new++]);
                }
                else
                {
                    /* This cg did not move */
                    sort_i = &(sort->sort2[nsort2++]);
                }
                /* Sort on the ns grid cell indices
                 * and the global topology index.
                 * index_gl is irrelevant with cell ns,
                 * but we set it here anyhow to avoid a conditional.
                 */
                sort_i->nsc    = a[i];
                sort_i->ind_gl = dd->index_gl[i];
                sort_i->ind    = i;
                ncg_new++;
            }
        }
        if (debug)
        {
            fprintf(debug, "ordered sort cgs: stationary %d moved %d\n",
                    nsort2, nsort_new);
        }
        /* Sort efficiently */
        ordered_sort(nsort2, sort->sort2, nsort_new, sort->sort_new,
                     sort->sort);
    }
    else
    {
        cgsort  = sort->sort;
        ncg_new = 0;
        for (i = 0; i < dd->ncg_home; i++)
        {
            /* Sort on the ns grid cell indices
             * and the global topology index
             */
            cgsort[i].nsc    = a[i];
            cgsort[i].ind_gl = dd->index_gl[i];
            cgsort[i].ind    = i;
            if (cgsort[i].nsc < moved)
            {
                ncg_new++;
            }
        }
        if (debug)
        {
            fprintf(debug, "qsort cgs: %d new home %d\n", dd->ncg_home, ncg_new);
        }
        /* Determine the order of the charge groups using qsort */
        gmx_qsort_threadsafe(cgsort, dd->ncg_home, sizeof(cgsort[0]), comp_cgsort);
    }

    return ncg_new;
}

static int dd_sort_order_nbnxn(gmx_domdec_t *dd, t_forcerec *fr)
{
    gmx_cgsort_t *sort;
    int           ncg_new, i, na;
    const int    *a;

    sort = dd->comm->sort->sort;

    nbnxn_get_atomorder(fr->nbv->nbs, &a, &na);

    ncg_new = 0;
    for (i = 0; i < na; i++)
    {
        if (a[i] >= 0)
        {
            sort[ncg_new].ind = a[i];
            ncg_new++;
        }
    }

    return ncg_new;
}

static void dd_sort_state(gmx_domdec_t *dd, rvec *cgcm, t_forcerec *fr, t_state *state,
                          int ncg_home_old)
{
    gmx_domdec_sort_t *sort;
    gmx_cgsort_t      *cgsort;
    int               *cgindex;
    int                ncg_new, i, *ibuf, cgsize;
    rvec              *vbuf;

    sort = dd->comm->sort;

    if (dd->ncg_home > sort->sort_nalloc)
    {
        sort->sort_nalloc = over_alloc_dd(dd->ncg_home);
        srenew(sort->sort, sort->sort_nalloc);
        srenew(sort->sort2, sort->sort_nalloc);
    }
    cgsort = sort->sort;

    switch (fr->cutoff_scheme)
    {
        case ecutsGROUP:
            ncg_new = dd_sort_order(dd, fr, ncg_home_old);
            break;
        case ecutsVERLET:
            ncg_new = dd_sort_order_nbnxn(dd, fr);
            break;
        default:
            gmx_incons("unimplemented");
            ncg_new = 0;
    }

    /* We alloc with the old size, since cgindex is still old */
    vec_rvec_check_alloc(&dd->comm->vbuf, dd->cgindex[dd->ncg_home]);
    vbuf = dd->comm->vbuf.v;

    if (dd->comm->bCGs)
    {
        cgindex = dd->cgindex;
    }
    else
    {
        cgindex = nullptr;
    }

    /* Remove the charge groups which are no longer at home here */
    dd->ncg_home = ncg_new;
    if (debug)
    {
        fprintf(debug, "Set the new home charge group count to %d\n",
                dd->ncg_home);
    }

    /* Reorder the state */
    if (state->flags & (1 << estX))
    {
        order_vec_atom(dd->ncg_home, cgindex, cgsort, as_rvec_array(state->x.data()), vbuf);
    }
    if (state->flags & (1 << estV))
    {
        order_vec_atom(dd->ncg_home, cgindex, cgsort, as_rvec_array(state->v.data()), vbuf);
    }
    if (state->flags & (1 << estCGP))
    {
        order_vec_atom(dd->ncg_home, cgindex, cgsort, as_rvec_array(state->cg_p.data()), vbuf);
    }

    if (fr->cutoff_scheme == ecutsGROUP)
    {
        /* Reorder cgcm */
        order_vec_cg(dd->ncg_home, cgsort, cgcm, vbuf);
    }

    if (dd->ncg_home+1 > sort->ibuf_nalloc)
    {
        sort->ibuf_nalloc = over_alloc_dd(dd->ncg_home+1);
        srenew(sort->ibuf, sort->ibuf_nalloc);
    }
    ibuf = sort->ibuf;
    /* Reorder the global cg index */
    order_int_cg(dd->ncg_home, cgsort, dd->index_gl, ibuf);
    /* Reorder the cginfo */
    order_int_cg(dd->ncg_home, cgsort, fr->cginfo, ibuf);
    /* Rebuild the local cg index */
    if (dd->comm->bCGs)
    {
        ibuf[0] = 0;
        for (i = 0; i < dd->ncg_home; i++)
        {
            cgsize    = dd->cgindex[cgsort[i].ind+1] - dd->cgindex[cgsort[i].ind];
            ibuf[i+1] = ibuf[i] + cgsize;
        }
        for (i = 0; i < dd->ncg_home+1; i++)
        {
            dd->cgindex[i] = ibuf[i];
        }
    }
    else
    {
        for (i = 0; i < dd->ncg_home+1; i++)
        {
            dd->cgindex[i] = i;
        }
    }
    /* Set the home atom number */
    dd->nat_home = dd->cgindex[dd->ncg_home];

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        /* The atoms are now exactly in grid order, update the grid order */
        nbnxn_set_atomorder(fr->nbv->nbs);
    }
    else
    {
        /* Copy the sorted ns cell indices back to the ns grid struct */
        for (i = 0; i < dd->ncg_home; i++)
        {
            fr->ns->grid->cell_index[i] = cgsort[i].nsc;
        }
        fr->ns->grid->nr = dd->ncg_home;
    }
}

static void add_dd_statistics(gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm;
    int                ddnat;

    comm = dd->comm;

    for (ddnat = ddnatZONE; ddnat < ddnatNR; ddnat++)
    {
        comm->sum_nat[ddnat-ddnatZONE] +=
            comm->nat[ddnat] - comm->nat[ddnat-1];
    }
    comm->ndecomp++;
}

void reset_dd_statistics_counters(gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm;
    int                ddnat;

    comm = dd->comm;

    /* Reset all the statistics and counters for total run counting */
    for (ddnat = ddnatZONE; ddnat < ddnatNR; ddnat++)
    {
        comm->sum_nat[ddnat-ddnatZONE] = 0;
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

void print_dd_statistics(t_commrec *cr, const t_inputrec *ir, FILE *fplog)
{
    gmx_domdec_comm_t *comm;
    int                ddnat;
    double             av;

    comm = cr->dd->comm;

    gmx_sumd(ddnatNR-ddnatZONE, comm->sum_nat, cr);

    if (fplog == nullptr)
    {
        return;
    }

    fprintf(fplog, "\n    D O M A I N   D E C O M P O S I T I O N   S T A T I S T I C S\n\n");

    for (ddnat = ddnatZONE; ddnat < ddnatNR; ddnat++)
    {
        av = comm->sum_nat[ddnat-ddnatZONE]/comm->ndecomp;
        switch (ddnat)
        {
            case ddnatZONE:
                fprintf(fplog,
                        " av. #atoms communicated per step for force:  %d x %.1f\n",
                        2, av);
                break;
            case ddnatVSITE:
                if (cr->dd->vsite_comm)
                {
                    fprintf(fplog,
                            " av. #atoms communicated per step for vsites: %d x %.1f\n",
                            (EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD) ? 3 : 2,
                            av);
                }
                break;
            case ddnatCON:
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
                         gmx_int64_t          step,
                         t_commrec           *cr,
                         gmx_bool             bMasterState,
                         int                  nstglobalcomm,
                         t_state             *state_global,
                         const gmx_mtop_t    *top_global,
                         const t_inputrec    *ir,
                         t_state             *state_local,
                         PaddedRVecVector    *f,
                         t_mdatoms           *mdatoms,
                         gmx_localtop_t      *top_local,
                         t_forcerec          *fr,
                         gmx_vsite_t         *vsite,
                         gmx_constr_t         constr,
                         t_nrnb              *nrnb,
                         gmx_wallcycle_t      wcycle,
                         gmx_bool             bVerbose)
{
    gmx_domdec_t      *dd;
    gmx_domdec_comm_t *comm;
    gmx_ddbox_t        ddbox = {0};
    t_block           *cgs_gl;
    gmx_int64_t        step_pcoupl;
    rvec               cell_ns_x0, cell_ns_x1;
    int                i, n, ncgindex_set, ncg_home_old = -1, ncg_moved, nat_f_novirsum;
    gmx_bool           bBoxChanged, bNStGlobalComm, bDoDLB, bCheckWhetherToTurnDlbOn, bLogLoad;
    gmx_bool           bRedist, bSortCG, bResortAll;
    ivec               ncells_old = {0, 0, 0}, ncells_new = {0, 0, 0}, np;
    real               grid_density;
    char               sbuf[22];

    wallcycle_start(wcycle, ewcDOMDEC);

    dd   = cr->dd;
    comm = dd->comm;

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
        n = ir->nstpcouple;
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
        clear_dd_indices(dd, 0, 0);
        ncgindex_set = 0;

        rvec *xGlobal = (SIMMASTER(cr) ? as_rvec_array(state_global->x.data()) : nullptr);

        set_ddbox(dd, bMasterState, cr, ir,
                  SIMMASTER(cr) ? state_global->box : nullptr,
                  TRUE, cgs_gl, xGlobal,
                  &ddbox);

        get_cg_distribution(fplog, dd, cgs_gl,
                            SIMMASTER(cr) ? state_global->box : nullptr,
                            &ddbox, xGlobal);

        dd_distribute_state(dd, cgs_gl,
                            state_global, state_local, f);

        dd_make_local_cgs(dd, &top_local->cgs);

        /* Ensure that we have space for the new distribution */
        dd_check_alloc_ncg(fr, state_local, f, dd->ncg_home);

        if (fr->cutoff_scheme == ecutsGROUP)
        {
            calc_cgcm(fplog, 0, dd->ncg_home,
                      &top_local->cgs, as_rvec_array(state_local->x.data()), fr->cg_cm);
        }

        inc_nrnb(nrnb, eNR_CGCM, dd->nat_home);

        dd_set_cginfo(dd->index_gl, 0, dd->ncg_home, fr, comm->bLocalCG);
    }
    else if (state_local->ddp_count != dd->ddp_count)
    {
        if (state_local->ddp_count > dd->ddp_count)
        {
            gmx_fatal(FARGS, "Internal inconsistency state_local->ddp_count (%d) > dd->ddp_count (%d)", state_local->ddp_count, dd->ddp_count);
        }

        if (state_local->ddp_count_cg_gl != state_local->ddp_count)
        {
            gmx_fatal(FARGS, "Internal inconsistency state_local->ddp_count_cg_gl (%d) != state_local->ddp_count (%d)", state_local->ddp_count_cg_gl, state_local->ddp_count);
        }

        /* Clear the old state */
        clear_dd_indices(dd, 0, 0);

        /* Build the new indices */
        rebuild_cgindex(dd, cgs_gl->index, state_local);
        make_dd_indices(dd, cgs_gl->index, 0);
        ncgindex_set = dd->ncg_home;

        if (fr->cutoff_scheme == ecutsGROUP)
        {
            /* Redetermine the cg COMs */
            calc_cgcm(fplog, 0, dd->ncg_home,
                      &top_local->cgs, as_rvec_array(state_local->x.data()), fr->cg_cm);
        }

        inc_nrnb(nrnb, eNR_CGCM, dd->nat_home);

        dd_set_cginfo(dd->index_gl, 0, dd->ncg_home, fr, comm->bLocalCG);

        set_ddbox(dd, bMasterState, cr, ir, state_local->box,
                  TRUE, &top_local->cgs, as_rvec_array(state_local->x.data()), &ddbox);

        bRedist = isDlbOn(comm);
    }
    else
    {
        /* We have the full state, only redistribute the cgs */

        /* Clear the non-home indices */
        clear_dd_indices(dd, dd->ncg_home, dd->nat_home);
        ncgindex_set = 0;

        /* Avoid global communication for dim's without pbc and -gcom */
        if (!bNStGlobalComm)
        {
            copy_rvec(comm->box0, ddbox.box0    );
            copy_rvec(comm->box_size, ddbox.box_size);
        }
        set_ddbox(dd, bMasterState, cr, ir, state_local->box,
                  bNStGlobalComm, &top_local->cgs, as_rvec_array(state_local->x.data()), &ddbox);

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
            nbnxn_get_ncells(fr->nbv->nbs, &ncells_old[XX], &ncells_old[YY]);
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
                set_zones_size(dd, state_local->box, &ddbox, 0, 1);

                nbnxn_put_on_grid(fr->nbv->nbs, fr->ePBC, state_local->box,
                                  0,
                                  comm->zones.size[0].bb_x0,
                                  comm->zones.size[0].bb_x1,
                                  0, dd->ncg_home,
                                  comm->zones.dens_zone0,
                                  fr->cginfo,
                                  as_rvec_array(state_local->x.data()),
                                  ncg_moved, bRedist ? comm->moved : nullptr,
                                  fr->nbv->grp[eintLocal].kernel_type,
                                  fr->nbv->grp[eintLocal].nbat);

                nbnxn_get_ncells(fr->nbv->nbs, &ncells_new[XX], &ncells_new[YY]);
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
        dd_resize_state(state_local, f, dd->nat_home);

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
        set_zones_size(dd, state_local->box, &ddbox,
                       bSortCG ? 1 : 0, comm->zones.n);
    }

    wallcycle_sub_stop(wcycle, ewcsDD_SETUPCOMM);

    /*
       write_dd_pdb("dd_home",step,"dump",top_global,cr,
                 -1,as_rvec_array(state_local->x.data()),state_local->box);
     */

    wallcycle_sub_start(wcycle, ewcsDD_MAKETOP);

    /* Extract a local topology from the global topology */
    for (i = 0; i < dd->ndim; i++)
    {
        np[dd->dim[i]] = comm->cd[i].np;
    }
    dd_make_local_top(dd, &comm->zones, dd->npbcdim, state_local->box,
                      comm->cellsize_min, np,
                      fr,
                      fr->cutoff_scheme == ecutsGROUP ? fr->cg_cm : as_rvec_array(state_local->x.data()),
                      vsite, top_global, top_local);

    wallcycle_sub_stop(wcycle, ewcsDD_MAKETOP);

    wallcycle_sub_start(wcycle, ewcsDD_MAKECONSTR);

    /* Set up the special atom communication */
    n = comm->nat[ddnatZONE];
    for (i = ddnatZONE+1; i < ddnatNR; i++)
    {
        switch (i)
        {
            case ddnatVSITE:
                if (vsite && vsite->n_intercg_vsite)
                {
                    n = dd_make_local_vsites(dd, n, top_local->idef.il);
                }
                break;
            case ddnatCON:
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
        comm->nat[i] = n;
    }

    wallcycle_sub_stop(wcycle, ewcsDD_MAKECONSTR);

    wallcycle_sub_start(wcycle, ewcsDD_TOPOTHER);

    /* Make space for the extra coordinates for virtual site
     * or constraint communication.
     */
    state_local->natoms = comm->nat[ddnatNR-1];

    dd_resize_state(state_local, f, state_local->natoms);

    if (fr->bF_NoVirSum)
    {
        if (vsite && vsite->n_intercg_vsite)
        {
            nat_f_novirsum = comm->nat[ddnatVSITE];
        }
        else
        {
            if (EEL_FULL(ir->coulombtype) && dd->n_intercg_excl > 0)
            {
                nat_f_novirsum = dd->nat_tot;
            }
            else
            {
                nat_f_novirsum = dd->nat_home;
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
    forcerec_set_ranges(fr, dd->ncg_home, dd->ncg_tot,
                        dd->nat_tot, comm->nat[ddnatCON], nat_f_novirsum);

    /* Update atom data for mdatoms and several algorithms */
    mdAlgorithmsSetupAtomData(cr, ir, top_global, top_local, fr,
                              nullptr, mdatoms, vsite, nullptr);

    if (ir->implicit_solvent)
    {
        make_local_gb(cr, fr->born, ir->gb_algorithm);
    }

    if (!(cr->duty & DUTY_PME))
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

    if (constr)
    {
        set_constraints(constr, top_local, ir, mdatoms, cr);
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
        dd_move_x(dd, state_local->box, as_rvec_array(state_local->x.data()));
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
