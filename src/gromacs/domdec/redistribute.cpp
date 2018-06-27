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
/* \internal \file
 *
 * \brief Implements atom redistribution functions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "redistribute.h"

#include <cstring>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nsgrid.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_internal.h"
#include "utility.h"

/* The size per charge group of the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_CGIBS 2

/* The flags for the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_FLAG_NRCG  65535
#define DD_FLAG_FW(d) (1<<(16+(d)*2))
#define DD_FLAG_BW(d) (1<<(16+(d)*2+1))

static int compact_and_copy_vec_at(int ncg, int *move,
                                   const gmx::BlockRanges &atomGroups,
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
        i1 = atomGroups.index[icg+1];
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
                                   const gmx::BlockRanges &atomGroups,
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
        i1 = atomGroups.index[icg + 1];
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

static int compact_ind(int                 ncg,
                       const int          *move,
                       gmx::ArrayRef<int>  globalAtomGroupIndices,
                       gmx::BlockRanges   *atomGroups,
                       gmx::ArrayRef<int>  globalAtomIndices,
                       gmx_ga2la_t        *ga2la,
                       char               *bLocalCG,
                       int                *cginfo)
{
    int home_pos = 0;
    int nat      = 0;
    for (int cg = 0; cg < ncg; cg++)
    {
        int a0 = atomGroups->index[cg];
        int a1 = atomGroups->index[cg+1];
        if (move[cg] == -1)
        {
            /* Compact the home arrays in place.
             * Anything that can be done here avoids access to global arrays.
             */
            atomGroups->index[home_pos] = nat;
            for (int a = a0; a < a1; a++)
            {
                const int a_gl         = globalAtomIndices[a];
                globalAtomIndices[nat] = a_gl;
                /* The cell number stays 0, so we don't need to set it */
                ga2la_change_la(ga2la, a_gl, nat);
                nat++;
            }
            globalAtomGroupIndices[home_pos] = globalAtomGroupIndices[cg];
            cginfo[home_pos]                 = cginfo[cg];
            /* The charge group remains local, so bLocalCG does not change */
            home_pos++;
        }
        else
        {
            /* Clear the global indices */
            for (int a = a0; a < a1; a++)
            {
                ga2la_del(ga2la, globalAtomIndices[a]);
            }
            if (bLocalCG)
            {
                bLocalCG[globalAtomGroupIndices[cg]] = FALSE;
            }
        }
    }
    atomGroups->index[home_pos] = nat;

    return home_pos;
}

static void clear_and_mark_ind(int                       ncg,
                               const int                *move,
                               gmx::ArrayRef<const int>  globalAtomGroupIndices,
                               const gmx::BlockRanges   &atomGroups,
                               gmx::ArrayRef<const int>  globalAtomIndices,
                               gmx_ga2la_t              *ga2la,
                               char                     *bLocalCG,
                               int                      *cell_index)
{
    for (int cg = 0; cg < ncg; cg++)
    {
        if (move[cg] >= 0)
        {
            int a0 = atomGroups.index[cg];
            int a1 = atomGroups.index[cg + 1];
            /* Clear the global indices */
            for (int a = a0; a < a1; a++)
            {
                ga2la_del(ga2la, globalAtomIndices[a]);
            }
            if (bLocalCG)
            {
                bLocalCG[globalAtomGroupIndices[cg]] = FALSE;
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
                ddglatnr(dd, dd->atomGroups().index[cg]), limitd, dim2char(dim));
    }
    else
    {
        /* We don't have a limiting distance available: don't print it */
        fprintf(fplog, "%s %d moved more than the distance allowed by the domain decomposition in direction %c\n",
                dd->comm->bCGs ? "The charge group starting at atom" : "Atom",
                ddglatnr(dd, dd->atomGroups().index[cg]), dim2char(dim));
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
                         const gmx::BlockRanges &atomGroups,
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
        k0   = atomGroups.index[cg];
        k1   = atomGroups.index[cg + 1];
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

void dd_redistribute_cg(FILE *fplog, gmx_int64_t step,
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

    if (dd->globalAtomGroupIndices.size() > static_cast<size_t>(comm->nalloc_int))
    {
        comm->nalloc_int = over_alloc_dd(dd->globalAtomGroupIndices.size());
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

    const gmx::BlockRanges &atomGroups = dd->atomGroups();

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
                         atomGroups,
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
            comm->cggl_flag[mc][ncg[mc]*DD_CGIBS  ] = dd->globalAtomGroupIndices[cg];
            /* We store the cg size in the lower 16 bits
             * and the place where the charge group should go
             * in the next 6 bits. This saves some communication volume.
             */
            nrcg = atomGroups.index[cg+1] - atomGroups.index[cg];
            comm->cggl_flag[mc][ncg[mc]*DD_CGIBS+1] = nrcg | flag;
            ncg[mc] += 1;
            nat[mc] += nrcg;
        }
    }

    inc_nrnb(nrnb, eNR_CGCM, comm->atomRanges.numHomeAtoms());
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
                compact_and_copy_vec_cg(dd->ncg_home, move, dd->atomGroups(),
                                        nvec, cg_cm, comm, bCompact);
            break;
        case ecutsVERLET:
            /* Without charge groups we send the moved atom coordinates
             * over twice. This is so the code below can be used without
             * many conditionals for both for with and without charge groups.
             */
            home_pos_cg =
                compact_and_copy_vec_cg(dd->ncg_home, move, dd->atomGroups(),
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
        compact_and_copy_vec_at(dd->ncg_home, move, dd->atomGroups(),
                                nvec, vec++, as_rvec_array(state->x.data()),
                                comm, bCompact);
    if (bV)
    {
        compact_and_copy_vec_at(dd->ncg_home, move, dd->atomGroups(),
                                nvec, vec++, as_rvec_array(state->v.data()),
                                comm, bCompact);
    }
    if (bCGP)
    {
        compact_and_copy_vec_at(dd->ncg_home, move, dd->atomGroups(),
                                nvec, vec++, as_rvec_array(state->cg_p.data()),
                                comm, bCompact);
    }

    if (bCompact)
    {
        compact_ind(dd->ncg_home, move,
                    dd->globalAtomGroupIndices, &dd->atomGroups_, dd->globalAtomIndices,
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
                           dd->globalAtomGroupIndices, dd->atomGroups(), dd->globalAtomIndices,
                           dd->ga2la, comm->bLocalCG,
                           moved);
    }

    /* Now we can remove the excess global atom-group indices from the list */
    dd->globalAtomGroupIndices.resize(home_pos_cg);
    dd->atomGroups_.index.resize(home_pos_cg + 1);

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
            ddSendrecv(dd, d, dir, sbuf, 2, rbuf, 2);

            if ((ncg_recv+rbuf[0])*DD_CGIBS > comm->nalloc_int)
            {
                comm->nalloc_int = over_alloc_dd((ncg_recv+rbuf[0])*DD_CGIBS);
                srenew(comm->buf_int, comm->nalloc_int);
            }

            /* Communicate the charge group indices, sizes and flags */
            ddSendrecv(dd, d, dir,
                       comm->cggl_flag[cdd], sbuf[0]*DD_CGIBS,
                       comm->buf_int+ncg_recv*DD_CGIBS, rbuf[0]*DD_CGIBS);

            nvs = ncg[cdd] + nat[cdd]*nvec;
            i   = rbuf[0]  + rbuf[1] *nvec;
            vec_rvec_check_alloc(&comm->vbuf, nvr+i);

            /* Communicate cgcm and state */
            ddSendrecv(dd, d, dir,
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

                            GMX_ASSERT(dim2 >= 0 && dim2 < DIM, "Keep the static analyzer happy");

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
                /* Set the global charge group index and size */
                const int globalAtomGroupIndex = comm->buf_int[cg*DD_CGIBS];
                dd->globalAtomGroupIndices.push_back(globalAtomGroupIndex);
                dd->atomGroups_.index.push_back(dd->atomGroups_.index[home_pos_cg] + nrcg);
                /* Copy the state from the buffer */
                if (fr->cutoff_scheme == ecutsGROUP)
                {
                    cg_cm = fr->cg_cm;
                    copy_rvec(comm->vbuf.v[buf_pos], cg_cm[home_pos_cg]);
                }
                buf_pos++;

                /* Set the cginfo */
                fr->cginfo[home_pos_cg] = ddcginfo(cginfo_mb,
                                                   globalAtomGroupIndex);
                if (comm->bLocalCG)
                {
                    comm->bLocalCG[globalAtomGroupIndex] = TRUE;
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
    comm->atomRanges.setEnd(DDAtomRanges::Type::Home, home_pos_at);

    if (fr->cutoff_scheme == ecutsGROUP && !bCompact)
    {
        /* We overallocated before, we need to set the right size here */
        dd_resize_state(state, f, comm->atomRanges.numHomeAtoms());
    }

    if (debug)
    {
        fprintf(debug,
                "Finished repartitioning: cgs moved out %d, new home %d\n",
                *ncg_moved, dd->ncg_home-*ncg_moved);

    }
}
