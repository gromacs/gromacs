/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "nbnxn_search.h"

#include "config.h"

#include <assert.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_grid.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/mdlib/nbnxn_util.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace


/* We shift the i-particles backward for PBC.
 * This leads to more conditionals than shifting forward.
 * We do this to get more balanced pair lists.
 */
constexpr bool c_pbcShiftBackward = true;


static void nbs_cycle_clear(nbnxn_cycle_t *cc)
{
    for (int i = 0; i < enbsCCnr; i++)
    {
        cc[i].count = 0;
        cc[i].c     = 0;
    }
}

static double Mcyc_av(const nbnxn_cycle_t *cc)
{
    return (double)cc->c*1e-6/cc->count;
}

static void nbs_cycle_print(FILE *fp, const nbnxn_search_t nbs)
{
    fprintf(fp, "\n");
    fprintf(fp, "ns %4d grid %4.1f search %4.1f red.f %5.3f",
            nbs->cc[enbsCCgrid].count,
            Mcyc_av(&nbs->cc[enbsCCgrid]),
            Mcyc_av(&nbs->cc[enbsCCsearch]),
            Mcyc_av(&nbs->cc[enbsCCreducef]));

    if (nbs->nthread_max > 1)
    {
        if (nbs->cc[enbsCCcombine].count > 0)
        {
            fprintf(fp, " comb %5.2f",
                    Mcyc_av(&nbs->cc[enbsCCcombine]));
        }
        fprintf(fp, " s. th");
        for (int t = 0; t < nbs->nthread_max; t++)
        {
            fprintf(fp, " %4.1f",
                    Mcyc_av(&nbs->work[t].cc[enbsCCsearch]));
        }
    }
    fprintf(fp, "\n");
}

/* Layout for the nonbonded NxN pair lists */
enum class NbnxnLayout
{
    NoSimd4x4, // i-cluster size 4, j-cluster size 4
    Simd4xN,   // i-cluster size 4, j-cluster size SIMD width
    Simd2xNN,  // i-cluster size 4, j-cluster size half SIMD width
    Gpu8x8x8   // i-cluster size 8, j-cluster size 8 + super-clustering
};

/* Returns the j-cluster size */
template <NbnxnLayout layout>
static constexpr int jClusterSize()
{
#if GMX_SIMD
    static_assert(layout == NbnxnLayout::NoSimd4x4 || layout == NbnxnLayout::Simd4xN || layout == NbnxnLayout::Simd2xNN, "Currently jClusterSize only supports CPU layouts");

    return layout == NbnxnLayout::Simd4xN ? GMX_SIMD_REAL_WIDTH : (layout == NbnxnLayout::Simd2xNN ? GMX_SIMD_REAL_WIDTH/2 : NBNXN_CPU_CLUSTER_I_SIZE);
#else
    static_assert(layout == NbnxnLayout::NoSimd4x4, "Currently without SIMD, jClusterSize only supports NoSimd4x4");

    return NBNXN_CPU_CLUSTER_I_SIZE;
#endif
}

/* Returns the j-cluster index given the i-cluster index */
template <int jClusterSize>
static inline int cjFromCi(int ci)
{
    static_assert(jClusterSize == NBNXN_CPU_CLUSTER_I_SIZE/2 || jClusterSize == NBNXN_CPU_CLUSTER_I_SIZE || jClusterSize == NBNXN_CPU_CLUSTER_I_SIZE*2, "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    if (jClusterSize == NBNXN_CPU_CLUSTER_I_SIZE/2)
    {
        return ci << 1;
    }
    else if (jClusterSize == NBNXN_CPU_CLUSTER_I_SIZE)
    {
        return ci;
    }
    else
    {
        return ci >> 1;
    }
}

/* Returns the j-cluster index given the i-cluster index */
template <NbnxnLayout layout>
static inline int cjFromCi(int ci)
{
    constexpr int clusterSize = jClusterSize<layout>();

    return cjFromCi<clusterSize>(ci);
}

/* Returns the nbnxn coordinate data index given the i-cluster index */
template <NbnxnLayout layout>
static inline int xIndexFromCi(int ci)
{
    constexpr int clusterSize = jClusterSize<layout>();

    static_assert(clusterSize == NBNXN_CPU_CLUSTER_I_SIZE/2 || clusterSize == NBNXN_CPU_CLUSTER_I_SIZE || clusterSize == NBNXN_CPU_CLUSTER_I_SIZE*2, "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    if (clusterSize <= NBNXN_CPU_CLUSTER_I_SIZE)
    {
        /* Coordinates are stored packed in groups of 4 */
        return ci*STRIDE_P4;
    }
    else
    {
        /* Coordinates packed in 8, i-cluster size is half the packing width */
        return (ci >> 1)*STRIDE_P8 + (ci & 1)*(c_packX8 >> 1);
    }
}

/* Returns the nbnxn coordinate data index given the j-cluster index */
template <NbnxnLayout layout>
static inline int xIndexFromCj(int cj)
{
    constexpr int clusterSize = jClusterSize<layout>();

    static_assert(clusterSize == NBNXN_CPU_CLUSTER_I_SIZE/2 || clusterSize == NBNXN_CPU_CLUSTER_I_SIZE || clusterSize == NBNXN_CPU_CLUSTER_I_SIZE*2, "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    if (clusterSize == NBNXN_CPU_CLUSTER_I_SIZE/2)
    {
        /* Coordinates are stored packed in groups of 4 */
        return (cj >> 1)*STRIDE_P4 + (cj & 1)*(c_packX4 >> 1);
    }
    else if (clusterSize == NBNXN_CPU_CLUSTER_I_SIZE)
    {
        /* Coordinates are stored packed in groups of 4 */
        return cj*STRIDE_P4;
    }
    else
    {
        /* Coordinates are stored packed in groups of 8 */
        return cj*STRIDE_P8;
    }
}

gmx_bool nbnxn_kernel_pairlist_simple(int nb_kernel_type)
{
    if (nb_kernel_type == nbnxnkNotSet)
    {
        gmx_fatal(FARGS, "Non-bonded kernel type not set for Verlet-style pair-list.");
    }

    switch (nb_kernel_type)
    {
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            return FALSE;

        case nbnxnk4x4_PlainC:
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
            return TRUE;

        default:
            gmx_incons("Invalid nonbonded kernel type passed!");
            return FALSE;
    }
}

/* Initializes a single nbnxn_pairlist_t data structure */
static void nbnxn_init_pairlist_fep(t_nblist *nl)
{
    nl->type        = GMX_NBLIST_INTERACTION_FREE_ENERGY;
    nl->igeometry   = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
    /* The interaction functions are set in the free energy kernel fuction */
    nl->ivdw        = -1;
    nl->ivdwmod     = -1;
    nl->ielec       = -1;
    nl->ielecmod    = -1;

    nl->maxnri      = 0;
    nl->maxnrj      = 0;
    nl->nri         = 0;
    nl->nrj         = 0;
    nl->iinr        = nullptr;
    nl->gid         = nullptr;
    nl->shift       = nullptr;
    nl->jindex      = nullptr;
    nl->jjnr        = nullptr;
    nl->excl_fep    = nullptr;

}

void nbnxn_init_search(nbnxn_search_t           * nbs_ptr,
                       ivec                      *n_dd_cells,
                       struct gmx_domdec_zones_t *zones,
                       gmx_bool                   bFEP,
                       int                        nthread_max)
{
    nbnxn_search_t nbs;
    int            ngrid;

    snew(nbs, 1);
    *nbs_ptr = nbs;

    nbs->bFEP   = bFEP;

    nbs->DomDec = (n_dd_cells != nullptr);

    clear_ivec(nbs->dd_dim);
    ngrid = 1;
    if (nbs->DomDec)
    {
        nbs->zones = zones;

        for (int d = 0; d < DIM; d++)
        {
            if ((*n_dd_cells)[d] > 1)
            {
                nbs->dd_dim[d] = 1;
                /* Each grid matches a DD zone */
                ngrid *= 2;
            }
        }
    }

    nbnxn_grids_init(nbs, ngrid);

    nbs->cell        = nullptr;
    nbs->cell_nalloc = 0;
    nbs->a           = nullptr;
    nbs->a_nalloc    = 0;

    nbs->nthread_max = nthread_max;

    /* Initialize the work data structures for each thread */
    snew(nbs->work, nbs->nthread_max);
    for (int t = 0; t < nbs->nthread_max; t++)
    {
        nbs->work[t].cxy_na           = nullptr;
        nbs->work[t].cxy_na_nalloc    = 0;
        nbs->work[t].sort_work        = nullptr;
        nbs->work[t].sort_work_nalloc = 0;

        snew(nbs->work[t].nbl_fep, 1);
        nbnxn_init_pairlist_fep(nbs->work[t].nbl_fep);
    }

    /* Initialize detailed nbsearch cycle counting */
    nbs->print_cycles = (getenv("GMX_NBNXN_CYCLE") != nullptr);
    nbs->search_count = 0;
    nbs_cycle_clear(nbs->cc);
    for (int t = 0; t < nbs->nthread_max; t++)
    {
        nbs_cycle_clear(nbs->work[t].cc);
    }
}

static void init_buffer_flags(nbnxn_buffer_flags_t *flags,
                              int                   natoms)
{
    flags->nflag = (natoms + NBNXN_BUFFERFLAG_SIZE - 1)/NBNXN_BUFFERFLAG_SIZE;
    if (flags->nflag > flags->flag_nalloc)
    {
        flags->flag_nalloc = over_alloc_large(flags->nflag);
        srenew(flags->flag, flags->flag_nalloc);
    }
    for (int b = 0; b < flags->nflag; b++)
    {
        bitmask_clear(&(flags->flag[b]));
    }
}

/* Determines the cell range along one dimension that
 * the bounding box b0 - b1 sees.
 */
static void get_cell_range(real b0, real b1,
                           int nc, real c0, real s, real invs,
                           real d2, real r2, int *cf, int *cl)
{
    *cf = std::max(static_cast<int>((b0 - c0)*invs), 0);

    while (*cf > 0 && d2 + gmx::square((b0 - c0) - (*cf-1+1)*s) < r2)
    {
        (*cf)--;
    }

    *cl = std::min(static_cast<int>((b1 - c0)*invs), nc-1);
    while (*cl < nc-1 && d2 + gmx::square((*cl+1)*s - (b1 - c0)) < r2)
    {
        (*cl)++;
    }
}

/* Reference code calculating the distance^2 between two bounding boxes */
/*
   static float box_dist2(float bx0, float bx1, float by0,
                       float by1, float bz0, float bz1,
                       const nbnxn_bb_t *bb)
   {
    float d2;
    float dl, dh, dm, dm0;

    d2 = 0;

    dl  = bx0 - bb->upper[BB_X];
    dh  = bb->lower[BB_X] - bx1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = by0 - bb->upper[BB_Y];
    dh  = bb->lower[BB_Y] - by1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = bz0 - bb->upper[BB_Z];
    dh  = bb->lower[BB_Z] - bz1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    return d2;
   }
 */

/* Plain C code calculating the distance^2 between two bounding boxes */
static float subc_bb_dist2(int si, const nbnxn_bb_t *bb_i_ci,
                           int csj, const nbnxn_bb_t *bb_j_all)
{
    const nbnxn_bb_t *bb_i, *bb_j;
    float             d2;
    float             dl, dh, dm, dm0;

    bb_i = bb_i_ci  +  si;
    bb_j = bb_j_all + csj;

    d2 = 0;

    dl  = bb_i->lower[BB_X] - bb_j->upper[BB_X];
    dh  = bb_j->lower[BB_X] - bb_i->upper[BB_X];
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Y] - bb_j->upper[BB_Y];
    dh  = bb_j->lower[BB_Y] - bb_i->upper[BB_Y];
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Z] - bb_j->upper[BB_Z];
    dh  = bb_j->lower[BB_Z] - bb_i->upper[BB_Z];
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    return d2;
}

#if NBNXN_SEARCH_BB_SIMD4

/* 4-wide SIMD code for bb distance for bb format xyz0 */
static float subc_bb_dist2_simd4(int si, const nbnxn_bb_t *bb_i_ci,
                                 int csj, const nbnxn_bb_t *bb_j_all)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    Simd4Float bb_i_S0, bb_i_S1;
    Simd4Float bb_j_S0, bb_j_S1;
    Simd4Float dl_S;
    Simd4Float dh_S;
    Simd4Float dm_S;
    Simd4Float dm0_S;

    bb_i_S0 = load4(&bb_i_ci[si].lower[0]);
    bb_i_S1 = load4(&bb_i_ci[si].upper[0]);
    bb_j_S0 = load4(&bb_j_all[csj].lower[0]);
    bb_j_S1 = load4(&bb_j_all[csj].upper[0]);

    dl_S    = bb_i_S0 - bb_j_S1;
    dh_S    = bb_j_S0 - bb_i_S1;

    dm_S    = max(dl_S, dh_S);
    dm0_S   = max(dm_S, simd4SetZeroF());

    return dotProduct(dm0_S, dm0_S);
}

/* Calculate bb bounding distances of bb_i[si,...,si+3] and store them in d2 */
#define SUBC_BB_DIST2_SIMD4_XXXX_INNER(si, bb_i, d2) \
    {                                                \
        int               shi;                                  \
                                                 \
        Simd4Float        dx_0, dy_0, dz_0;                    \
        Simd4Float        dx_1, dy_1, dz_1;                    \
                                                 \
        Simd4Float        mx, my, mz;                          \
        Simd4Float        m0x, m0y, m0z;                       \
                                                 \
        Simd4Float        d2x, d2y, d2z;                       \
        Simd4Float        d2s, d2t;                            \
                                                 \
        shi = si*NNBSBB_D*DIM;                       \
                                                 \
        xi_l = load4(bb_i+shi+0*STRIDE_PBB);   \
        yi_l = load4(bb_i+shi+1*STRIDE_PBB);   \
        zi_l = load4(bb_i+shi+2*STRIDE_PBB);   \
        xi_h = load4(bb_i+shi+3*STRIDE_PBB);   \
        yi_h = load4(bb_i+shi+4*STRIDE_PBB);   \
        zi_h = load4(bb_i+shi+5*STRIDE_PBB);   \
                                                 \
        dx_0 = xi_l - xj_h;                 \
        dy_0 = yi_l - yj_h;                 \
        dz_0 = zi_l - zj_h;                 \
                                                 \
        dx_1 = xj_l - xi_h;                 \
        dy_1 = yj_l - yi_h;                 \
        dz_1 = zj_l - zi_h;                 \
                                                 \
        mx   = max(dx_0, dx_1);                 \
        my   = max(dy_0, dy_1);                 \
        mz   = max(dz_0, dz_1);                 \
                                                 \
        m0x  = max(mx, zero);                   \
        m0y  = max(my, zero);                   \
        m0z  = max(mz, zero);                   \
                                                 \
        d2x  = m0x * m0x;                   \
        d2y  = m0y * m0y;                   \
        d2z  = m0z * m0z;                   \
                                                 \
        d2s  = d2x + d2y;                   \
        d2t  = d2s + d2z;                   \
                                                 \
        store4(d2+si, d2t);                      \
    }

/* 4-wide SIMD code for nsi bb distances for bb format xxxxyyyyzzzz */
static void subc_bb_dist2_simd4_xxxx(const float *bb_j,
                                     int nsi, const float *bb_i,
                                     float *d2)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    Simd4Float xj_l, yj_l, zj_l;
    Simd4Float xj_h, yj_h, zj_h;
    Simd4Float xi_l, yi_l, zi_l;
    Simd4Float xi_h, yi_h, zi_h;

    Simd4Float zero;

    zero = setZero();

    xj_l = Simd4Float(bb_j[0*STRIDE_PBB]);
    yj_l = Simd4Float(bb_j[1*STRIDE_PBB]);
    zj_l = Simd4Float(bb_j[2*STRIDE_PBB]);
    xj_h = Simd4Float(bb_j[3*STRIDE_PBB]);
    yj_h = Simd4Float(bb_j[4*STRIDE_PBB]);
    zj_h = Simd4Float(bb_j[5*STRIDE_PBB]);

    /* Here we "loop" over si (0,STRIDE_PBB) from 0 to nsi with step STRIDE_PBB.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    SUBC_BB_DIST2_SIMD4_XXXX_INNER(0, bb_i, d2);
    if (STRIDE_PBB < nsi)
    {
        SUBC_BB_DIST2_SIMD4_XXXX_INNER(STRIDE_PBB, bb_i, d2);
    }
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */


/* Returns if any atom pair from two clusters is within distance sqrt(rlist2) */
static gmx_inline gmx_bool
clusterpair_in_range(const nbnxn_list_work_t *work,
                     int si,
                     int csj, int stride, const real *x_j,
                     real rlist2)
{
#if !GMX_SIMD4_HAVE_REAL

    /* Plain C version.
     * All coordinates are stored as xyzxyz...
     */

    const real *x_i = work->x_ci;

    for (int i = 0; i < c_nbnxnGpuClusterSize; i++)
    {
        int i0 = (si*c_nbnxnGpuClusterSize + i)*DIM;
        for (int j = 0; j < c_nbnxnGpuClusterSize; j++)
        {
            int  j0 = (csj*c_nbnxnGpuClusterSize + j)*stride;

            real d2 = gmx::square(x_i[i0  ] - x_j[j0  ]) + gmx::square(x_i[i0+1] - x_j[j0+1]) + gmx::square(x_i[i0+2] - x_j[j0+2]);

            if (d2 < rlist2)
            {
                return TRUE;
            }
        }
    }

    return FALSE;

#else /* !GMX_SIMD4_HAVE_REAL */

    /* 4-wide SIMD version.
     * A cluster is hard-coded to 8 atoms.
     * The coordinates x_i are stored as xxxxyyyy..., x_j is stored xyzxyz...
     * Using 8-wide AVX(2) is not faster on Intel Sandy Bridge and Haswell.
     */
    assert(c_nbnxnGpuClusterSize == 8);

    Simd4Real   rc2_S      = Simd4Real(rlist2);

    const real *x_i        = work->x_ci_simd;

    int         dim_stride = c_nbnxnGpuClusterSize*DIM;
    Simd4Real   ix_S0      = load4(x_i + si*dim_stride + 0*GMX_SIMD4_WIDTH);
    Simd4Real   iy_S0      = load4(x_i + si*dim_stride + 1*GMX_SIMD4_WIDTH);
    Simd4Real   iz_S0      = load4(x_i + si*dim_stride + 2*GMX_SIMD4_WIDTH);
    Simd4Real   ix_S1      = load4(x_i + si*dim_stride + 3*GMX_SIMD4_WIDTH);
    Simd4Real   iy_S1      = load4(x_i + si*dim_stride + 4*GMX_SIMD4_WIDTH);
    Simd4Real   iz_S1      = load4(x_i + si*dim_stride + 5*GMX_SIMD4_WIDTH);

    /* We loop from the outer to the inner particles to maximize
     * the chance that we find a pair in range quickly and return.
     */
    int j0 = csj*c_nbnxnGpuClusterSize;
    int j1 = j0 + c_nbnxnGpuClusterSize - 1;
    while (j0 < j1)
    {
        Simd4Real jx0_S, jy0_S, jz0_S;
        Simd4Real jx1_S, jy1_S, jz1_S;

        Simd4Real dx_S0, dy_S0, dz_S0;
        Simd4Real dx_S1, dy_S1, dz_S1;
        Simd4Real dx_S2, dy_S2, dz_S2;
        Simd4Real dx_S3, dy_S3, dz_S3;

        Simd4Real rsq_S0;
        Simd4Real rsq_S1;
        Simd4Real rsq_S2;
        Simd4Real rsq_S3;

        Simd4Bool wco_S0;
        Simd4Bool wco_S1;
        Simd4Bool wco_S2;
        Simd4Bool wco_S3;
        Simd4Bool wco_any_S01, wco_any_S23, wco_any_S;

        jx0_S = Simd4Real(x_j[j0*stride+0]);
        jy0_S = Simd4Real(x_j[j0*stride+1]);
        jz0_S = Simd4Real(x_j[j0*stride+2]);

        jx1_S = Simd4Real(x_j[j1*stride+0]);
        jy1_S = Simd4Real(x_j[j1*stride+1]);
        jz1_S = Simd4Real(x_j[j1*stride+2]);

        /* Calculate distance */
        dx_S0            = ix_S0 - jx0_S;
        dy_S0            = iy_S0 - jy0_S;
        dz_S0            = iz_S0 - jz0_S;
        dx_S1            = ix_S1 - jx0_S;
        dy_S1            = iy_S1 - jy0_S;
        dz_S1            = iz_S1 - jz0_S;
        dx_S2            = ix_S0 - jx1_S;
        dy_S2            = iy_S0 - jy1_S;
        dz_S2            = iz_S0 - jz1_S;
        dx_S3            = ix_S1 - jx1_S;
        dy_S3            = iy_S1 - jy1_S;
        dz_S3            = iz_S1 - jz1_S;

        /* rsq = dx*dx+dy*dy+dz*dz */
        rsq_S0           = norm2(dx_S0, dy_S0, dz_S0);
        rsq_S1           = norm2(dx_S1, dy_S1, dz_S1);
        rsq_S2           = norm2(dx_S2, dy_S2, dz_S2);
        rsq_S3           = norm2(dx_S3, dy_S3, dz_S3);

        wco_S0           = (rsq_S0 < rc2_S);
        wco_S1           = (rsq_S1 < rc2_S);
        wco_S2           = (rsq_S2 < rc2_S);
        wco_S3           = (rsq_S3 < rc2_S);

        wco_any_S01      = wco_S0 || wco_S1;
        wco_any_S23      = wco_S2 || wco_S3;
        wco_any_S        = wco_any_S01 || wco_any_S23;

        if (anyTrue(wco_any_S))
        {
            return TRUE;
        }

        j0++;
        j1--;
    }

    return FALSE;

#endif /* !GMX_SIMD4_HAVE_REAL */
}

/* Returns the j-cluster index for index cjIndex in a cj list */
static inline int nblCj(const nbnxn_cj_t *cjList, int cjIndex)
{
    return cjList[cjIndex].cj;
}

/* Returns the j-cluster index for index cjIndex in a cj4 list */
static inline int nblCj(const nbnxn_cj4_t *cj4List, int cjIndex)
{
    return cj4List[cjIndex/c_nbnxnGpuJgroupSize].cj[cjIndex & (c_nbnxnGpuJgroupSize - 1)];
}

/* Returns the i-interaction mask of the j sub-cell for index cj_ind */
static unsigned int nbl_imask0(const nbnxn_pairlist_t *nbl, int cj_ind)
{
    return nbl->cj4[cj_ind/c_nbnxnGpuJgroupSize].imei[0].imask;
}

/* Ensures there is enough space for extra extra exclusion masks */
static void check_excl_space(nbnxn_pairlist_t *nbl, int extra)
{
    if (nbl->nexcl+extra > nbl->excl_nalloc)
    {
        nbl->excl_nalloc = over_alloc_small(nbl->nexcl+extra);
        nbnxn_realloc_void((void **)&nbl->excl,
                           nbl->nexcl*sizeof(*nbl->excl),
                           nbl->excl_nalloc*sizeof(*nbl->excl),
                           nbl->alloc, nbl->free);
    }
}

/* Ensures there is enough space for maxNumExtraClusters extra j-clusters in the list */
static void check_cell_list_space_simple(nbnxn_pairlist_t *nbl,
                                         int               maxNumExtraClusters)
{
    int cj_max;

    cj_max = nbl->ncj + maxNumExtraClusters;

    if (cj_max > nbl->cj_nalloc)
    {
        nbl->cj_nalloc = over_alloc_small(cj_max);
        nbnxn_realloc_void((void **)&nbl->cj,
                           nbl->ncj*sizeof(*nbl->cj),
                           nbl->cj_nalloc*sizeof(*nbl->cj),
                           nbl->alloc, nbl->free);

        nbnxn_realloc_void((void **)&nbl->cjOuter,
                           nbl->ncj*sizeof(*nbl->cjOuter),
                           nbl->cj_nalloc*sizeof(*nbl->cjOuter),
                           nbl->alloc, nbl->free);
    }
}

/* Ensures there is enough space for ncell extra j-clusters in the list */
static void check_cell_list_space_supersub(nbnxn_pairlist_t *nbl,
                                           int               ncell)
{
    int ncj4_max, w;

    /* We can have maximally nsupercell*c_gpuNumClusterPerCell sj lists */
    /* We can store 4 j-subcell - i-supercell pairs in one struct.
     * since we round down, we need one extra entry.
     */
    ncj4_max = ((nbl->work->cj_ind + ncell*c_gpuNumClusterPerCell + c_nbnxnGpuJgroupSize - 1)/c_nbnxnGpuJgroupSize);

    if (ncj4_max > nbl->cj4_nalloc)
    {
        nbl->cj4_nalloc = over_alloc_small(ncj4_max);
        nbnxn_realloc_void((void **)&nbl->cj4,
                           nbl->work->cj4_init*sizeof(*nbl->cj4),
                           nbl->cj4_nalloc*sizeof(*nbl->cj4),
                           nbl->alloc, nbl->free);
    }

    if (ncj4_max > nbl->work->cj4_init)
    {
        for (int j4 = nbl->work->cj4_init; j4 < ncj4_max; j4++)
        {
            /* No i-subcells and no excl's in the list initially */
            for (w = 0; w < c_nbnxnGpuClusterpairSplit; w++)
            {
                nbl->cj4[j4].imei[w].imask    = 0U;
                nbl->cj4[j4].imei[w].excl_ind = 0;

            }
        }
        nbl->work->cj4_init = ncj4_max;
    }
}

/* Set all excl masks for one GPU warp no exclusions */
static void set_no_excls(nbnxn_excl_t *excl)
{
    for (int t = 0; t < c_nbnxnGpuExclSize; t++)
    {
        /* Turn all interaction bits on */
        excl->pair[t] = NBNXN_INTERACTION_MASK_ALL;
    }
}

/* Initializes a single nbnxn_pairlist_t data structure */
static void nbnxn_init_pairlist(nbnxn_pairlist_t *nbl,
                                gmx_bool          bSimple,
                                nbnxn_alloc_t    *alloc,
                                nbnxn_free_t     *free)
{
    if (alloc == nullptr)
    {
        nbl->alloc = nbnxn_alloc_aligned;
    }
    else
    {
        nbl->alloc = alloc;
    }
    if (free == nullptr)
    {
        nbl->free = nbnxn_free_aligned;
    }
    else
    {
        nbl->free = free;
    }

    nbl->bSimple     = bSimple;
    nbl->na_sc       = 0;
    nbl->na_ci       = 0;
    nbl->na_cj       = 0;
    nbl->nci         = 0;
    nbl->ci          = nullptr;
    nbl->ci_nalloc   = 0;
    nbl->nsci        = 0;
    nbl->sci         = nullptr;
    nbl->sci_nalloc  = 0;
    nbl->ncj         = 0;
    nbl->ncjInUse    = 0;
    nbl->cj          = nullptr;
    nbl->cj_nalloc   = 0;
    nbl->ncj4        = 0;
    /* We need one element extra in sj, so alloc initially with 1 */
    nbl->cj4_nalloc  = 0;
    nbl->cj4         = nullptr;
    nbl->nci_tot     = 0;

    if (!nbl->bSimple)
    {
        GMX_ASSERT(c_nbnxnGpuNumClusterPerSupercluster == c_gpuNumClusterPerCell, "The search code assumes that the a super-cluster matches a search grid cell");

        GMX_ASSERT(sizeof(nbl->cj4[0].imei[0].imask)*8 >= c_nbnxnGpuJgroupSize*c_gpuNumClusterPerCell, "The i super-cluster cluster interaction mask does not contain a sufficient number of bits");
        GMX_ASSERT(sizeof(nbl->excl[0])*8 >= c_nbnxnGpuJgroupSize*c_gpuNumClusterPerCell, "The GPU exclusion mask does not contain a sufficient number of bits");

        nbl->excl        = nullptr;
        nbl->excl_nalloc = 0;
        nbl->nexcl       = 0;
        check_excl_space(nbl, 1);
        nbl->nexcl       = 1;
        set_no_excls(&nbl->excl[0]);
    }

    snew(nbl->work, 1);
    if (nbl->bSimple)
    {
        snew_aligned(nbl->work->bb_ci, 1, NBNXN_SEARCH_BB_MEM_ALIGN);
    }
    else
    {
#if NBNXN_BBXXXX
        snew_aligned(nbl->work->pbb_ci, c_gpuNumClusterPerCell/STRIDE_PBB*NNBSBB_XXXX, NBNXN_SEARCH_BB_MEM_ALIGN);
#else
        snew_aligned(nbl->work->bb_ci, c_gpuNumClusterPerCell, NBNXN_SEARCH_BB_MEM_ALIGN);
#endif
    }
    int gpu_clusterpair_nc = c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize*DIM;
    snew(nbl->work->x_ci, gpu_clusterpair_nc);
#if GMX_SIMD
    snew_aligned(nbl->work->x_ci_simd,
                 std::max(NBNXN_CPU_CLUSTER_I_SIZE*DIM*GMX_SIMD_REAL_WIDTH,
                          gpu_clusterpair_nc),
                 GMX_SIMD_REAL_WIDTH);
#endif
    snew_aligned(nbl->work->d2, c_gpuNumClusterPerCell, NBNXN_SEARCH_BB_MEM_ALIGN);

    nbl->work->sort            = nullptr;
    nbl->work->sort_nalloc     = 0;
    nbl->work->sci_sort        = nullptr;
    nbl->work->sci_sort_nalloc = 0;
}

void nbnxn_init_pairlist_set(nbnxn_pairlist_set_t *nbl_list,
                             gmx_bool bSimple, gmx_bool bCombined,
                             nbnxn_alloc_t *alloc,
                             nbnxn_free_t  *free)
{
    nbl_list->bSimple   = bSimple;
    nbl_list->bCombined = bCombined;

    nbl_list->nnbl = gmx_omp_nthreads_get(emntNonbonded);

    if (!nbl_list->bCombined &&
        nbl_list->nnbl > NBNXN_BUFFERFLAG_MAX_THREADS)
    {
        gmx_fatal(FARGS, "%d OpenMP threads were requested. Since the non-bonded force buffer reduction is prohibitively slow with more than %d threads, we do not allow this. Use %d or less OpenMP threads.",
                  nbl_list->nnbl, NBNXN_BUFFERFLAG_MAX_THREADS, NBNXN_BUFFERFLAG_MAX_THREADS);
    }

    snew(nbl_list->nbl, nbl_list->nnbl);
    if (bSimple && nbl_list->nnbl > 1)
    {
        snew(nbl_list->nbl_work, nbl_list->nnbl);
    }
    snew(nbl_list->nbl_fep, nbl_list->nnbl);
    /* Execute in order to avoid memory interleaving between threads */
#pragma omp parallel for num_threads(nbl_list->nnbl) schedule(static)
    for (int i = 0; i < nbl_list->nnbl; i++)
    {
        try
        {
            /* Allocate the nblist data structure locally on each thread
             * to optimize memory access for NUMA architectures.
             */
            snew(nbl_list->nbl[i], 1);

            /* Only list 0 is used on the GPU, use normal allocation for i>0 */
            if (!bSimple && i == 0)
            {
                nbnxn_init_pairlist(nbl_list->nbl[i], nbl_list->bSimple, alloc, free);
            }
            else
            {
                nbnxn_init_pairlist(nbl_list->nbl[i], nbl_list->bSimple, nullptr, nullptr);
                if (bSimple && nbl_list->nnbl > 1)
                {
                    snew(nbl_list->nbl_work[i], 1);
                    nbnxn_init_pairlist(nbl_list->nbl_work[i], nbl_list->bSimple, nullptr, nullptr);
                }
            }

            snew(nbl_list->nbl_fep[i], 1);
            nbnxn_init_pairlist_fep(nbl_list->nbl_fep[i]);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

/* Print statistics of a pair list, used for debug output */
static void print_nblist_statistics_simple(FILE *fp, const nbnxn_pairlist_t *nbl,
                                           const nbnxn_search_t nbs, real rl)
{
    const nbnxn_grid_t *grid;
    int                 cs[SHIFTS];
    int                 npexcl;

    grid = &nbs->grid[0];

    fprintf(fp, "nbl nci %d ncj %d\n",
            nbl->nci, nbl->ncjInUse);
    fprintf(fp, "nbl na_sc %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_sc, rl, nbl->ncjInUse, nbl->ncjInUse/(double)grid->nc,
            nbl->ncjInUse/(double)grid->nc*grid->na_sc,
            nbl->ncjInUse/(double)grid->nc*grid->na_sc/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nc*grid->na_sc/(grid->size[XX]*grid->size[YY]*grid->size[ZZ])));

    fprintf(fp, "nbl average j cell list length %.1f\n",
            0.25*nbl->ncjInUse/(double)std::max(nbl->nci, 1));

    for (int s = 0; s < SHIFTS; s++)
    {
        cs[s] = 0;
    }
    npexcl = 0;
    for (int i = 0; i < nbl->nci; i++)
    {
        cs[nbl->ci[i].shift & NBNXN_CI_SHIFT] +=
            nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start;

        int j = nbl->ci[i].cj_ind_start;
        while (j < nbl->ci[i].cj_ind_end &&
               nbl->cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            npexcl++;
            j++;
        }
    }
    fprintf(fp, "nbl cell pairs, total: %d excl: %d %.1f%%\n",
            nbl->ncj, npexcl, 100*npexcl/(double)std::max(nbl->ncj, 1));
    for (int s = 0; s < SHIFTS; s++)
    {
        if (cs[s] > 0)
        {
            fprintf(fp, "nbl shift %2d ncj %3d\n", s, cs[s]);
        }
    }
}

/* Print statistics of a pair lists, used for debug output */
static void print_nblist_statistics_supersub(FILE *fp, const nbnxn_pairlist_t *nbl,
                                             const nbnxn_search_t nbs, real rl)
{
    const nbnxn_grid_t *grid;
    int                 b;
    int                 c[c_gpuNumClusterPerCell + 1];
    double              sum_nsp, sum_nsp2;
    int                 nsp_max;

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp, "nbl nsci %d ncj4 %d nsi %d excl4 %d\n",
            nbl->nsci, nbl->ncj4, nbl->nci_tot, nbl->nexcl);
    fprintf(fp, "nbl na_c %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_ci, rl, nbl->nci_tot, nbl->nci_tot/(double)grid->nsubc_tot,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nsubc_tot*grid->na_c/(grid->size[XX]*grid->size[YY]*grid->size[ZZ])));

    sum_nsp  = 0;
    sum_nsp2 = 0;
    nsp_max  = 0;
    for (int si = 0; si <= c_gpuNumClusterPerCell; si++)
    {
        c[si] = 0;
    }
    for (int i = 0; i < nbl->nsci; i++)
    {
        int nsp;

        nsp = 0;
        for (int j4 = nbl->sci[i].cj4_ind_start; j4 < nbl->sci[i].cj4_ind_end; j4++)
        {
            for (int j = 0; j < c_nbnxnGpuJgroupSize; j++)
            {
                b = 0;
                for (int si = 0; si < c_gpuNumClusterPerCell; si++)
                {
                    if (nbl->cj4[j4].imei[0].imask & (1U << (j*c_gpuNumClusterPerCell + si)))
                    {
                        b++;
                    }
                }
                nsp += b;
                c[b]++;
            }
        }
        sum_nsp  += nsp;
        sum_nsp2 += nsp*nsp;
        nsp_max   = std::max(nsp_max, nsp);
    }
    if (nbl->nsci > 0)
    {
        sum_nsp  /= nbl->nsci;
        sum_nsp2 /= nbl->nsci;
    }
    fprintf(fp, "nbl #cluster-pairs: av %.1f stddev %.1f max %d\n",
            sum_nsp, std::sqrt(sum_nsp2 - sum_nsp*sum_nsp), nsp_max);

    if (nbl->ncj4 > 0)
    {
        for (b = 0; b <= c_gpuNumClusterPerCell; b++)
        {
            fprintf(fp, "nbl j-list #i-subcell %d %7d %4.1f\n",
                    b, c[b],
                    100.0*c[b]/(double)(nbl->ncj4*c_nbnxnGpuJgroupSize));
        }
    }
}

/* Returns a pointer to the exclusion mask for cj4-unit cj4, warp warp */
static void low_get_nbl_exclusions(nbnxn_pairlist_t *nbl, int cj4,
                                   int warp, nbnxn_excl_t **excl)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* No exclusions set, make a new list entry */
        nbl->cj4[cj4].imei[warp].excl_ind = nbl->nexcl;
        nbl->nexcl++;
        *excl = &nbl->excl[nbl->cj4[cj4].imei[warp].excl_ind];
        set_no_excls(*excl);
    }
    else
    {
        /* We already have some exclusions, new ones can be added to the list */
        *excl = &nbl->excl[nbl->cj4[cj4].imei[warp].excl_ind];
    }
}

/* Returns a pointer to the exclusion mask for cj4-unit cj4, warp warp,
 * generates a new element and allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_1(nbnxn_pairlist_t *nbl, int cj4,
                                 int warp, nbnxn_excl_t **excl)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* We need to make a new list entry, check if we have space */
        check_excl_space(nbl, 1);
    }
    low_get_nbl_exclusions(nbl, cj4, warp, excl);
}

/* Returns pointers to the exclusion masks for cj4-unit cj4 for both warps,
 * generates a new element and allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_2(nbnxn_pairlist_t *nbl, int cj4,
                                 nbnxn_excl_t **excl_w0,
                                 nbnxn_excl_t **excl_w1)
{
    /* Check for space we might need */
    check_excl_space(nbl, 2);

    low_get_nbl_exclusions(nbl, cj4, 0, excl_w0);
    low_get_nbl_exclusions(nbl, cj4, 1, excl_w1);
}

/* Sets the self exclusions i=j and pair exclusions i>j */
static void set_self_and_newton_excls_supersub(nbnxn_pairlist_t *nbl,
                                               int cj4_ind, int sj_offset,
                                               int i_cluster_in_cell)
{
    nbnxn_excl_t *excl[c_nbnxnGpuClusterpairSplit];

    /* Here we only set the set self and double pair exclusions */

    assert(c_nbnxnGpuClusterpairSplit == 2);

    get_nbl_exclusions_2(nbl, cj4_ind, &excl[0], &excl[1]);

    /* Only minor < major bits set */
    for (int ej = 0; ej < nbl->na_ci; ej++)
    {
        int w = (ej>>2);
        for (int ei = ej; ei < nbl->na_ci; ei++)
        {
            excl[w]->pair[(ej & (c_nbnxnGpuJgroupSize-1))*nbl->na_ci + ei] &=
                ~(1U << (sj_offset*c_gpuNumClusterPerCell + i_cluster_in_cell));
        }
    }
}

/* Returns a diagonal or off-diagonal interaction mask for plain C lists */
static unsigned int get_imask(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=2 */
gmx_unused static unsigned int get_imask_simd_j2(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci*2 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_0 :
            (rdiag && ci*2+1 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_1 :
             NBNXN_INTERACTION_MASK_ALL));
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=4 */
gmx_unused static unsigned int get_imask_simd_j4(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=8 */
gmx_unused static unsigned int get_imask_simd_j8(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj*2 ? NBNXN_INTERACTION_MASK_DIAG_J8_0 :
            (rdiag && ci == cj*2+1 ? NBNXN_INTERACTION_MASK_DIAG_J8_1 :
             NBNXN_INTERACTION_MASK_ALL));
}

#if GMX_SIMD
#if GMX_SIMD_REAL_WIDTH == 2
#define get_imask_simd_4xn  get_imask_simd_j2
#endif
#if GMX_SIMD_REAL_WIDTH == 4
#define get_imask_simd_4xn  get_imask_simd_j4
#endif
#if GMX_SIMD_REAL_WIDTH == 8
#define get_imask_simd_4xn  get_imask_simd_j8
#define get_imask_simd_2xnn get_imask_simd_j4
#endif
#if GMX_SIMD_REAL_WIDTH == 16
#define get_imask_simd_2xnn get_imask_simd_j8
#endif
#endif

/* Plain C code for checking and adding cluster-pairs to the list.
 *
 * \param[in]     gridj               The j-grid
 * \param[in,out] nbl                 The pair-list to store the cluster pairs in
 * \param[in]     icluster            The index of the i-cluster
 * \param[in]     jclusterFirst       The first cluster in the j-range
 * \param[in]     jclusterLast        The last cluster in the j-range
 * \param[in]     excludeSubDiagonal  Exclude atom pairs with i-index > j-index
 * \param[in]     x_j                 Coordinates for the j-atom, in xyz format
 * \param[in]     rlist2              The squared list cut-off
 * \param[in]     rbb2                The squared cut-off for putting cluster-pairs in the list based on bounding box distance only
 * \param[in,out] numDistanceChecks   The number of distance checks performed
 */
static void
makeClusterListSimple(const nbnxn_grid_t *      gridj,
                      nbnxn_pairlist_t *        nbl,
                      int                       icluster,
                      int                       jclusterFirst,
                      int                       jclusterLast,
                      bool                      excludeSubDiagonal,
                      const real * gmx_restrict x_j,
                      real                      rlist2,
                      float                     rbb2,
                      int * gmx_restrict        numDistanceChecks)
{
    const nbnxn_bb_t * gmx_restrict bb_ci = nbl->work->bb_ci;
    const real * gmx_restrict       x_ci  = nbl->work->x_ci;

    gmx_bool                        InRange;

    InRange = FALSE;
    while (!InRange && jclusterFirst <= jclusterLast)
    {
        real d2  = subc_bb_dist2(0, bb_ci, jclusterFirst, gridj->bb);
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rlist2)
        {
            int cjf_gl = gridj->cell0 + jclusterFirst;
            for (int i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (int j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (gmx::square(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                         gmx::square(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                         gmx::square(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rlist2);
                }
            }
            *numDistanceChecks += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
        }
        if (!InRange)
        {
            jclusterFirst++;
        }
    }
    if (!InRange)
    {
        return;
    }

    InRange = FALSE;
    while (!InRange && jclusterLast > jclusterFirst)
    {
        real d2  = subc_bb_dist2(0, bb_ci, jclusterLast, gridj->bb);
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rlist2)
        {
            int cjl_gl = gridj->cell0 + jclusterLast;
            for (int i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (int j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (gmx::square(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                         gmx::square(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                         gmx::square(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rlist2);
                }
            }
            *numDistanceChecks += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
        }
        if (!InRange)
        {
            jclusterLast--;
        }
    }

    if (jclusterFirst <= jclusterLast)
    {
        for (int jcluster = jclusterFirst; jcluster <= jclusterLast; jcluster++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[nbl->ncj].cj   = gridj->cell0 + jcluster;
            nbl->cj[nbl->ncj].excl = get_imask(excludeSubDiagonal, icluster, jcluster);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_search_simd_4xn.h"
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_search_simd_2xnn.h"
#endif

/* Plain C or SIMD4 code for making a pair list of super-cell sci vs scj.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_supersub(const nbnxn_grid_t *gridi,
                                       const nbnxn_grid_t *gridj,
                                       nbnxn_pairlist_t *nbl,
                                       int sci, int scj,
                                       gmx_bool sci_equals_scj,
                                       int stride, const real *x,
                                       real rlist2, float rbb2,
                                       int *numDistanceChecks)
{
    nbnxn_list_work_t *work   = nbl->work;

#if NBNXN_BBXXXX
    const float       *pbb_ci = work->pbb_ci;
#else
    const nbnxn_bb_t  *bb_ci  = work->bb_ci;
#endif

    assert(c_nbnxnGpuClusterSize == gridi->na_c);
    assert(c_nbnxnGpuClusterSize == gridj->na_c);

    /* We generate the pairlist mainly based on bounding-box distances
     * and do atom pair distance based pruning on the GPU.
     * Only if a j-group contains a single cluster-pair, we try to prune
     * that pair based on atom distances on the CPU to avoid empty j-groups.
     */
#define PRUNE_LIST_CPU_ONE 1
#define PRUNE_LIST_CPU_ALL 0

#if PRUNE_LIST_CPU_ONE
    int  ci_last = -1;
#endif

    float *d2l = work->d2;

    for (int subc = 0; subc < gridj->nsubc[scj]; subc++)
    {
        int          cj4_ind   = nbl->work->cj_ind/c_nbnxnGpuJgroupSize;
        int          cj_offset = nbl->work->cj_ind - cj4_ind*c_nbnxnGpuJgroupSize;
        nbnxn_cj4_t *cj4       = &nbl->cj4[cj4_ind];

        int          cj        = scj*c_gpuNumClusterPerCell + subc;

        int          cj_gl     = gridj->cell0*c_gpuNumClusterPerCell + cj;

        /* Initialize this j-subcell i-subcell list */
        cj4->cj[cj_offset] = cj_gl;

        int ci1;
        if (sci_equals_scj)
        {
            ci1 = subc + 1;
        }
        else
        {
            ci1 = gridi->nsubc[sci];
        }

#if NBNXN_BBXXXX
        /* Determine all ci1 bb distances in one call with SIMD4 */
        subc_bb_dist2_simd4_xxxx(gridj->pbb+(cj>>STRIDE_PBB_2LOG)*NNBSBB_XXXX+(cj & (STRIDE_PBB-1)),
                                 ci1, pbb_ci, d2l);
        *numDistanceChecks += c_nbnxnGpuClusterSize*2;
#endif

        int          npair = 0;
        unsigned int imask = 0;
        /* We use a fixed upper-bound instead of ci1 to help optimization */
        for (int ci = 0; ci < c_gpuNumClusterPerCell; ci++)
        {
            if (ci == ci1)
            {
                break;
            }

#if !NBNXN_BBXXXX
            /* Determine the bb distance between ci and cj */
            d2l[ci]             = subc_bb_dist2(ci, bb_ci, cj, gridj->bb);
            *numDistanceChecks += 2;
#endif
            float d2 = d2l[ci];

#if PRUNE_LIST_CPU_ALL
            /* Check if the distance is within the distance where
             * we use only the bounding box distance rbb,
             * or within the cut-off and there is at least one atom pair
             * within the cut-off. This check is very costly.
             */
            *numDistanceChecks += c_nbnxnGpuClusterSize*c_nbnxnGpuClusterSize;
            if (d2 < rbb2 ||
                (d2 < rlist2 &&
                 clusterpair_in_range(work, ci, cj_gl, stride, x, rlist2)))
#else
            /* Check if the distance between the two bounding boxes
             * in within the pair-list cut-off.
             */
            if (d2 < rlist2)
#endif
            {
                /* Flag this i-subcell to be taken into account */
                imask |= (1U << (cj_offset*c_gpuNumClusterPerCell + ci));

#if PRUNE_LIST_CPU_ONE
                ci_last = ci;
#endif

                npair++;
            }
        }

#if PRUNE_LIST_CPU_ONE
        /* If we only found 1 pair, check if any atoms are actually
         * within the cut-off, so we could get rid of it.
         */
        if (npair == 1 && d2l[ci_last] >= rbb2 &&
            !clusterpair_in_range(work, ci_last, cj_gl, stride, x, rlist2))
        {
            imask &= ~(1U << (cj_offset*c_gpuNumClusterPerCell + ci_last));
            npair--;
        }
#endif

        if (npair > 0)
        {
            /* We have a useful sj entry, close it now */

            /* Set the exclusions for the ci==sj entry.
             * Here we don't bother to check if this entry is actually flagged,
             * as it will nearly always be in the list.
             */
            if (sci_equals_scj)
            {
                set_self_and_newton_excls_supersub(nbl, cj4_ind, cj_offset, subc);
            }

            /* Copy the cluster interaction mask to the list */
            for (int w = 0; w < c_nbnxnGpuClusterpairSplit; w++)
            {
                cj4->imei[w].imask |= imask;
            }

            nbl->work->cj_ind++;

            /* Keep the count */
            nbl->nci_tot += npair;

            /* Increase the closing index in i super-cell list */
            nbl->sci[nbl->nsci].cj4_ind_end =
                (nbl->work->cj_ind + c_nbnxnGpuJgroupSize - 1)/c_nbnxnGpuJgroupSize;
        }
    }
}

/* Returns how many contiguous j-clusters we have starting in the i-list */
template <typename CjListType>
static int numContiguousJClusters(const int         cjIndexStart,
                                  const int         cjIndexEnd,
                                  const CjListType &cjList)
{
    const int firstJCluster = nblCj(cjList, cjIndexStart);

    int       numContiguous = 0;

    while (cjIndexStart + numContiguous < cjIndexEnd &&
           nblCj(cjList, cjIndexStart + numContiguous) == firstJCluster + numContiguous)
    {
        numContiguous++;
    }

    return numContiguous;
}

/* Helper struct for efficient searching for excluded atoms in a j-list */
struct JListRanges
{
    /* Constructor */
    template <typename CjListType>
    JListRanges(int               cjIndexStart,
                int               cjIndexEnd,
                const CjListType &cjList);

    int cjIndexStart; // The start index in the j-list
    int cjIndexEnd;   // The end index in the j-list
    int cjFirst;      // The j-cluster with index cjIndexStart
    int cjLast;       // The j-cluster with index cjIndexEnd-1
    int numDirect;    // Up to cjIndexStart+numDirect the j-clusters are cjFirst + the index offset
};

template <typename CjListType>
JListRanges::JListRanges(int               cjIndexStart,
                         int               cjIndexEnd,
                         const CjListType &cjList) :
    cjIndexStart(cjIndexStart),
    cjIndexEnd(cjIndexEnd)
{
    GMX_ASSERT(cjIndexEnd > cjIndexStart, "JListRanges should only be called with non-empty lists");

    cjFirst   = nblCj(cjList, cjIndexStart);
    cjLast    = nblCj(cjList, cjIndexEnd - 1);

    /* Determine how many contiguous j-cells we have starting
     * from the first i-cell. This number can be used to directly
     * calculate j-cell indices for excluded atoms.
     */
    numDirect = numContiguousJClusters(cjIndexStart, cjIndexEnd, cjList);
}

/* Return the index of \p jCluster in the given range or -1 when not present
 *
 * Note: This code is executed very often and therefore performance is
 *       important. It should be inlined and fully optimized.
 */
template <typename CjListType>
static inline int findJClusterInJList(int                jCluster,
                                      const JListRanges &ranges,
                                      const CjListType  &cjList)
{
    int index;

    if (jCluster < ranges.cjFirst + ranges.numDirect)
    {
        /* We can calculate the index directly using the offset */
        index = ranges.cjIndexStart + jCluster - ranges.cjFirst;
    }
    else
    {
        /* Search for jCluster using bisection */
        index           = -1;
        int rangeStart  = ranges.cjIndexStart + ranges.numDirect;
        int rangeEnd    = ranges.cjIndexEnd;
        int rangeMiddle;
        while (index == -1 && rangeStart < rangeEnd)
        {
            rangeMiddle = (rangeStart + rangeEnd) >> 1;

            const int clusterMiddle = nblCj(cjList, rangeMiddle);

            if (jCluster == clusterMiddle)
            {
                index      = rangeMiddle;
            }
            else if (jCluster < clusterMiddle)
            {
                rangeEnd   = rangeMiddle;
            }
            else
            {
                rangeStart = rangeMiddle + 1;
            }
        }
    }

    return index;
}

/* Set all atom-pair exclusions for a simple type list i-entry
 *
 * Set all atom-pair exclusions from the topology stored in exclusions
 * as masks in the pair-list for simple list entry iEntry.
 */
static void
setExclusionsForSimpleIentry(const nbnxn_search_t  nbs,
                             nbnxn_pairlist_t     *nbl,
                             gmx_bool              diagRemoved,
                             int                   na_cj_2log,
                             const nbnxn_ci_t     &iEntry,
                             const t_blocka       &exclusions)
{
    if (iEntry.cj_ind_end == iEntry.cj_ind_start)
    {
        /* Empty list: no exclusions */
        return;
    }

    const JListRanges  ranges(iEntry.cj_ind_start, iEntry.cj_ind_end, nbl->cj);

    const int          iCluster = iEntry.ci;

    const int         *cell = nbs->cell;

    /* Loop over the atoms in the i-cluster */
    for (int i = 0; i < nbl->na_sc; i++)
    {
        const int iIndex = iCluster*nbl->na_sc + i;
        const int iAtom  = nbs->a[iIndex];
        if (iAtom >= 0)
        {
            /* Loop over the topology-based exclusions for this i-atom */
            for (int exclIndex = exclusions.index[iAtom]; exclIndex < exclusions.index[iAtom + 1]; exclIndex++)
            {
                const int jAtom = exclusions.a[exclIndex];

                if (jAtom == iAtom)
                {
                    /* The self exclusion are already set, save some time */
                    continue;
                }

                /* Get the index of the j-atom in the nbnxn atom data */
                const int jIndex = cell[jAtom];

                /* Without shifts we only calculate interactions j>i
                 * for one-way pair-lists.
                 */
                if (diagRemoved && jIndex <= iIndex)
                {
                    continue;
                }

                const int jCluster = (jIndex >> na_cj_2log);

                /* Could the cluster se be in our list? */
                if (jCluster >= ranges.cjFirst && jCluster <= ranges.cjLast)
                {
                    const int index =
                        findJClusterInJList(jCluster, ranges, nbl->cj);

                    if (index >= 0)
                    {
                        /* We found an exclusion, clear the corresponding
                         * interaction bit.
                         */
                        const int innerJ     = jIndex - (jCluster << na_cj_2log);

                        nbl->cj[index].excl &= ~(1U << ((i << na_cj_2log) + innerJ));
                    }
                }
            }
        }
    }
}

/* Add a new i-entry to the FEP list and copy the i-properties */
static gmx_inline void fep_list_new_nri_copy(t_nblist *nlist)
{
    /* Add a new i-entry */
    nlist->nri++;

    assert(nlist->nri < nlist->maxnri);

    /* Duplicate the last i-entry, except for jindex, which continues */
    nlist->iinr[nlist->nri]   = nlist->iinr[nlist->nri-1];
    nlist->shift[nlist->nri]  = nlist->shift[nlist->nri-1];
    nlist->gid[nlist->nri]    = nlist->gid[nlist->nri-1];
    nlist->jindex[nlist->nri] = nlist->nrj;
}

/* For load balancing of the free-energy lists over threads, we set
 * the maximum nrj size of an i-entry to 40. This leads to good
 * load balancing in the worst case scenario of a single perturbed
 * particle on 16 threads, while not introducing significant overhead.
 * Note that half of the perturbed pairs will anyhow end up in very small lists,
 * since non perturbed i-particles will see few perturbed j-particles).
 */
const int max_nrj_fep = 40;

/* Exclude the perturbed pairs from the Verlet list. This is only done to avoid
 * singularities for overlapping particles (0/0), since the charges and
 * LJ parameters have been zeroed in the nbnxn data structure.
 * Simultaneously make a group pair list for the perturbed pairs.
 */
static void make_fep_list(const nbnxn_search_t    nbs,
                          const nbnxn_atomdata_t *nbat,
                          nbnxn_pairlist_t       *nbl,
                          gmx_bool                bDiagRemoved,
                          nbnxn_ci_t             *nbl_ci,
                          const nbnxn_grid_t     *gridi,
                          const nbnxn_grid_t     *gridj,
                          t_nblist               *nlist)
{
    int      ci, cj_ind_start, cj_ind_end, cja, cjr;
    int      nri_max;
    int      ngid, gid_i = 0, gid_j, gid;
    int      egp_shift, egp_mask;
    int      gid_cj = 0;
    int      ind_i, ind_j, ai, aj;
    int      nri;
    gmx_bool bFEP_i, bFEP_i_all;

    if (nbl_ci->cj_ind_end == nbl_ci->cj_ind_start)
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    cj_ind_start = nbl_ci->cj_ind_start;
    cj_ind_end   = nbl_ci->cj_ind_end;

    /* In worst case we have alternating energy groups
     * and create #atom-pair lists, which means we need the size
     * of a cluster pair (na_ci*na_cj) times the number of cj's.
     */
    nri_max = nbl->na_ci*nbl->na_cj*(cj_ind_end - cj_ind_start);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    ngid = nbat->nenergrp;

    if (static_cast<std::size_t>(ngid*gridj->na_cj) > sizeof(gid_cj)*8)
    {
        gmx_fatal(FARGS, "The Verlet scheme with %dx%d kernels and free-energy only supports up to %d energy groups",
                  gridi->na_c, gridj->na_cj, (sizeof(gid_cj)*8)/gridj->na_cj);
    }

    egp_shift = nbat->neg_2log;
    egp_mask  = (1<<nbat->neg_2log) - 1;

    /* Loop over the atoms in the i sub-cell */
    bFEP_i_all = TRUE;
    for (int i = 0; i < nbl->na_ci; i++)
    {
        ind_i = ci*nbl->na_ci + i;
        ai    = nbs->a[ind_i];
        if (ai >= 0)
        {
            nri                  = nlist->nri;
            nlist->jindex[nri+1] = nlist->jindex[nri];
            nlist->iinr[nri]     = ai;
            /* The actual energy group pair index is set later */
            nlist->gid[nri]      = 0;
            nlist->shift[nri]    = nbl_ci->shift & NBNXN_CI_SHIFT;

            bFEP_i = gridi->fep[ci - gridi->cell0] & (1 << i);

            bFEP_i_all = bFEP_i_all && bFEP_i;

            if (nlist->nrj + (cj_ind_end - cj_ind_start)*nbl->na_cj > nlist->maxnrj)
            {
                nlist->maxnrj = over_alloc_small(nlist->nrj + (cj_ind_end - cj_ind_start)*nbl->na_cj);
                srenew(nlist->jjnr,     nlist->maxnrj);
                srenew(nlist->excl_fep, nlist->maxnrj);
            }

            if (ngid > 1)
            {
                gid_i = (nbat->energrp[ci] >> (egp_shift*i)) & egp_mask;
            }

            for (int cj_ind = cj_ind_start; cj_ind < cj_ind_end; cj_ind++)
            {
                unsigned int fep_cj;

                cja = nbl->cj[cj_ind].cj;

                if (gridj->na_cj == gridj->na_c)
                {
                    cjr    = cja - gridj->cell0;
                    fep_cj = gridj->fep[cjr];
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja];
                    }
                }
                else if (2*gridj->na_cj == gridj->na_c)
                {
                    cjr    = cja - gridj->cell0*2;
                    /* Extract half of the ci fep/energrp mask */
                    fep_cj = (gridj->fep[cjr>>1] >> ((cjr&1)*gridj->na_cj)) & ((1<<gridj->na_cj) - 1);
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja>>1] >> ((cja&1)*gridj->na_cj*egp_shift) & ((1<<(gridj->na_cj*egp_shift)) - 1);
                    }
                }
                else
                {
                    cjr    = cja - (gridj->cell0>>1);
                    /* Combine two ci fep masks/energrp */
                    fep_cj = gridj->fep[cjr*2] + (gridj->fep[cjr*2+1] << gridj->na_c);
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja*2] + (nbat->energrp[cja*2+1] << (gridj->na_c*egp_shift));
                    }
                }

                if (bFEP_i || fep_cj != 0)
                {
                    for (int j = 0; j < nbl->na_cj; j++)
                    {
                        /* Is this interaction perturbed and not excluded? */
                        ind_j = cja*nbl->na_cj + j;
                        aj    = nbs->a[ind_j];
                        if (aj >= 0 &&
                            (bFEP_i || (fep_cj & (1 << j))) &&
                            (!bDiagRemoved || ind_j >= ind_i))
                        {
                            if (ngid > 1)
                            {
                                gid_j = (gid_cj >> (j*egp_shift)) & egp_mask;
                                gid   = GID(gid_i, gid_j, ngid);

                                if (nlist->nrj > nlist->jindex[nri] &&
                                    nlist->gid[nri] != gid)
                                {
                                    /* Energy group pair changed: new list */
                                    fep_list_new_nri_copy(nlist);
                                    nri = nlist->nri;
                                }
                                nlist->gid[nri] = gid;
                            }

                            if (nlist->nrj - nlist->jindex[nri] >= max_nrj_fep)
                            {
                                fep_list_new_nri_copy(nlist);
                                nri = nlist->nri;
                            }

                            /* Add it to the FEP list */
                            nlist->jjnr[nlist->nrj]     = aj;
                            nlist->excl_fep[nlist->nrj] = (nbl->cj[cj_ind].excl >> (i*nbl->na_cj + j)) & 1;
                            nlist->nrj++;

                            /* Exclude it from the normal list.
                             * Note that the charge has been set to zero,
                             * but we need to avoid 0/0, as perturbed atoms
                             * can be on top of each other.
                             */
                            nbl->cj[cj_ind].excl &= ~(1U << (i*nbl->na_cj + j));
                        }
                    }
                }
            }

            if (nlist->nrj > nlist->jindex[nri])
            {
                /* Actually add this new, non-empty, list */
                nlist->nri++;
                nlist->jindex[nlist->nri] = nlist->nrj;
            }
        }
    }

    if (bFEP_i_all)
    {
        /* All interactions are perturbed, we can skip this entry */
        nbl_ci->cj_ind_end = cj_ind_start;
        nbl->ncjInUse     -= cj_ind_end - cj_ind_start;
    }
}

/* Return the index of atom a within a cluster */
static gmx_inline int cj_mod_cj4(int cj)
{
    return cj & (c_nbnxnGpuJgroupSize - 1);
}

/* Convert a j-cluster to a cj4 group */
static gmx_inline int cj_to_cj4(int cj)
{
    return cj/c_nbnxnGpuJgroupSize;
}

/* Return the index of an j-atom within a warp */
static gmx_inline int a_mod_wj(int a)
{
    return a & (c_nbnxnGpuClusterSize/c_nbnxnGpuClusterpairSplit - 1);
}

/* As make_fep_list above, but for super/sub lists. */
static void make_fep_list_supersub(const nbnxn_search_t    nbs,
                                   const nbnxn_atomdata_t *nbat,
                                   nbnxn_pairlist_t       *nbl,
                                   gmx_bool                bDiagRemoved,
                                   const nbnxn_sci_t      *nbl_sci,
                                   real                    shx,
                                   real                    shy,
                                   real                    shz,
                                   real                    rlist_fep2,
                                   const nbnxn_grid_t     *gridi,
                                   const nbnxn_grid_t     *gridj,
                                   t_nblist               *nlist)
{
    int                sci, cj4_ind_start, cj4_ind_end, cjr;
    int                nri_max;
    int                c_abs;
    int                ind_i, ind_j, ai, aj;
    int                nri;
    gmx_bool           bFEP_i;
    real               xi, yi, zi;
    const nbnxn_cj4_t *cj4;

    if (nbl_sci->cj4_ind_end == nbl_sci->cj4_ind_start)
    {
        /* Empty list */
        return;
    }

    sci = nbl_sci->sci;

    cj4_ind_start = nbl_sci->cj4_ind_start;
    cj4_ind_end   = nbl_sci->cj4_ind_end;

    /* Here we process one super-cell, max #atoms na_sc, versus a list
     * cj4 entries, each with max c_nbnxnGpuJgroupSize cj's, each
     * of size na_cj atoms.
     * On the GPU we don't support energy groups (yet).
     * So for each of the na_sc i-atoms, we need max one FEP list
     * for each max_nrj_fep j-atoms.
     */
    nri_max = nbl->na_sc*nbl->na_cj*(1 + ((cj4_ind_end - cj4_ind_start)*c_nbnxnGpuJgroupSize)/max_nrj_fep);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    /* Loop over the atoms in the i super-cluster */
    for (int c = 0; c < c_gpuNumClusterPerCell; c++)
    {
        c_abs = sci*c_gpuNumClusterPerCell + c;

        for (int i = 0; i < nbl->na_ci; i++)
        {
            ind_i = c_abs*nbl->na_ci + i;
            ai    = nbs->a[ind_i];
            if (ai >= 0)
            {
                nri                  = nlist->nri;
                nlist->jindex[nri+1] = nlist->jindex[nri];
                nlist->iinr[nri]     = ai;
                /* With GPUs, energy groups are not supported */
                nlist->gid[nri]      = 0;
                nlist->shift[nri]    = nbl_sci->shift & NBNXN_CI_SHIFT;

                bFEP_i = (gridi->fep[c_abs - gridi->cell0*c_gpuNumClusterPerCell] & (1 << i));

                xi = nbat->x[ind_i*nbat->xstride+XX] + shx;
                yi = nbat->x[ind_i*nbat->xstride+YY] + shy;
                zi = nbat->x[ind_i*nbat->xstride+ZZ] + shz;

                if ((nlist->nrj + cj4_ind_end - cj4_ind_start)*c_nbnxnGpuJgroupSize*nbl->na_cj > nlist->maxnrj)
                {
                    nlist->maxnrj = over_alloc_small((nlist->nrj + cj4_ind_end - cj4_ind_start)*c_nbnxnGpuJgroupSize*nbl->na_cj);
                    srenew(nlist->jjnr,     nlist->maxnrj);
                    srenew(nlist->excl_fep, nlist->maxnrj);
                }

                for (int cj4_ind = cj4_ind_start; cj4_ind < cj4_ind_end; cj4_ind++)
                {
                    cj4 = &nbl->cj4[cj4_ind];

                    for (int gcj = 0; gcj < c_nbnxnGpuJgroupSize; gcj++)
                    {
                        unsigned int fep_cj;

                        if ((cj4->imei[0].imask & (1U << (gcj*c_gpuNumClusterPerCell + c))) == 0)
                        {
                            /* Skip this ci for this cj */
                            continue;
                        }

                        cjr = cj4->cj[gcj] - gridj->cell0*c_gpuNumClusterPerCell;

                        fep_cj = gridj->fep[cjr];

                        if (bFEP_i || fep_cj != 0)
                        {
                            for (int j = 0; j < nbl->na_cj; j++)
                            {
                                /* Is this interaction perturbed and not excluded? */
                                ind_j = (gridj->cell0*c_gpuNumClusterPerCell + cjr)*nbl->na_cj + j;
                                aj    = nbs->a[ind_j];
                                if (aj >= 0 &&
                                    (bFEP_i || (fep_cj & (1 << j))) &&
                                    (!bDiagRemoved || ind_j >= ind_i))
                                {
                                    nbnxn_excl_t *excl;
                                    int           excl_pair;
                                    unsigned int  excl_bit;
                                    real          dx, dy, dz;

                                    get_nbl_exclusions_1(nbl, cj4_ind, j>>2, &excl);

                                    excl_pair = a_mod_wj(j)*nbl->na_ci + i;
                                    excl_bit  = (1U << (gcj*c_gpuNumClusterPerCell + c));

                                    dx = nbat->x[ind_j*nbat->xstride+XX] - xi;
                                    dy = nbat->x[ind_j*nbat->xstride+YY] - yi;
                                    dz = nbat->x[ind_j*nbat->xstride+ZZ] - zi;

                                    /* The unpruned GPU list has more than 2/3
                                     * of the atom pairs beyond rlist. Using
                                     * this list will cause a lot of overhead
                                     * in the CPU FEP kernels, especially
                                     * relative to the fast GPU kernels.
                                     * So we prune the FEP list here.
                                     */
                                    if (dx*dx + dy*dy + dz*dz < rlist_fep2)
                                    {
                                        if (nlist->nrj - nlist->jindex[nri] >= max_nrj_fep)
                                        {
                                            fep_list_new_nri_copy(nlist);
                                            nri = nlist->nri;
                                        }

                                        /* Add it to the FEP list */
                                        nlist->jjnr[nlist->nrj]     = aj;
                                        nlist->excl_fep[nlist->nrj] = (excl->pair[excl_pair] & excl_bit) ? 1 : 0;
                                        nlist->nrj++;
                                    }

                                    /* Exclude it from the normal list.
                                     * Note that the charge and LJ parameters have
                                     * been set to zero, but we need to avoid 0/0,
                                     * as perturbed atoms can be on top of each other.
                                     */
                                    excl->pair[excl_pair] &= ~excl_bit;
                                }
                            }

                            /* Note that we could mask out this pair in imask
                             * if all i- and/or all j-particles are perturbed.
                             * But since the perturbed pairs on the CPU will
                             * take an order of magnitude more time, the GPU
                             * will finish before the CPU and there is no gain.
                             */
                        }
                    }
                }

                if (nlist->nrj > nlist->jindex[nri])
                {
                    /* Actually add this new, non-empty, list */
                    nlist->nri++;
                    nlist->jindex[nlist->nri] = nlist->nrj;
                }
            }
        }
    }
}

/* Set all atom-pair exclusions for a GPU type list i-entry
 *
 * Sets all atom-pair exclusions from the topology stored in exclusions
 * as masks in the pair-list for i-super-cluster list entry iEntry.
 */
static void
setExclusionsForGpuIentry(const nbnxn_search_t  nbs,
                          nbnxn_pairlist_t     *nbl,
                          gmx_bool              diagRemoved,
                          const nbnxn_sci_t    &iEntry,
                          const t_blocka       &exclusions)
{
    if (iEntry.cj4_ind_end == iEntry.cj4_ind_start)
    {
        /* Empty list */
        return;
    }

    /* Set the search ranges using start and end j-cluster indices.
     * Note that here we can not use cj4_ind_end, since the last cj4
     * can be only partially filled, so we use cj_ind.
     */
    const JListRanges ranges(iEntry.cj4_ind_start*c_nbnxnGpuJgroupSize,
                             nbl->work->cj_ind,
                             nbl->cj4);

    GMX_ASSERT(nbl->na_ci == c_nbnxnGpuClusterSize, "na_ci should match the GPU cluster size");
    constexpr int  c_clusterSize      = c_nbnxnGpuClusterSize;
    constexpr int  c_superClusterSize = c_nbnxnGpuNumClusterPerSupercluster*c_nbnxnGpuClusterSize;

    const int      iSuperCluster = iEntry.sci;

    const int     *cell = nbs->cell;

    /* Loop over the atoms in the i super-cluster */
    for (int i = 0; i < c_superClusterSize; i++)
    {
        const int iIndex = iSuperCluster*c_superClusterSize + i;
        const int iAtom  = nbs->a[iIndex];
        if (iAtom >= 0)
        {
            const int iCluster = i/c_clusterSize;

            /* Loop over the topology-based exclusions for this i-atom */
            for (int exclIndex = exclusions.index[iAtom]; exclIndex < exclusions.index[iAtom + 1]; exclIndex++)
            {
                const int jAtom = exclusions.a[exclIndex];

                if (jAtom == iAtom)
                {
                    /* The self exclusions are already set, save some time */
                    continue;
                }

                /* Get the index of the j-atom in the nbnxn atom data */
                const int jIndex = cell[jAtom];

                /* Without shifts we only calculate interactions j>i
                 * for one-way pair-lists.
                 */
                /* NOTE: We would like to use iIndex on the right hand side,
                 * but that makes this routine 25% slower with gcc6/7.
                 * Even using c_superClusterSize makes it slower.
                 * Either of these changes triggers peeling of the exclIndex
                 * loop, which apparently leads to far less efficient code.
                 */
                if (diagRemoved && jIndex <= iSuperCluster*nbl->na_sc + i)
                {
                    continue;
                }

                const int jCluster = jIndex/c_clusterSize;

                /* Check whether the cluster is in our list? */
                if (jCluster >= ranges.cjFirst && jCluster <= ranges.cjLast)
                {
                    const int index =
                        findJClusterInJList(jCluster, ranges, nbl->cj4);

                    if (index >= 0)
                    {
                        /* We found an exclusion, clear the corresponding
                         * interaction bit.
                         */
                        const unsigned int pairMask = (1U << (cj_mod_cj4(index)*c_gpuNumClusterPerCell + iCluster));
                        /* Check if the i-cluster interacts with the j-cluster */
                        if (nbl_imask0(nbl, index) & pairMask)
                        {
                            const int innerI = (i      & (c_clusterSize - 1));
                            const int innerJ = (jIndex & (c_clusterSize - 1));

                            /* Determine which j-half (CUDA warp) we are in */
                            const int     jHalf = innerJ/(c_clusterSize/c_nbnxnGpuClusterpairSplit);

                            nbnxn_excl_t *interactionMask;
                            get_nbl_exclusions_1(nbl, cj_to_cj4(index), jHalf, &interactionMask);

                            interactionMask->pair[a_mod_wj(innerJ)*c_clusterSize + innerI] &= ~pairMask;
                        }
                    }
                }
            }
        }
    }
}

/* Reallocate the simple ci list for at least n entries */
static void nb_realloc_ci(nbnxn_pairlist_t *nbl, int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->ci,
                       nbl->nci*sizeof(*nbl->ci),
                       nbl->ci_nalloc*sizeof(*nbl->ci),
                       nbl->alloc, nbl->free);

    nbnxn_realloc_void((void **)&nbl->ciOuter,
                       nbl->nci*sizeof(*nbl->ciOuter),
                       nbl->ci_nalloc*sizeof(*nbl->ciOuter),
                       nbl->alloc, nbl->free);
}

/* Reallocate the super-cell sci list for at least n entries */
static void nb_realloc_sci(nbnxn_pairlist_t *nbl, int n)
{
    nbl->sci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->sci,
                       nbl->nsci*sizeof(*nbl->sci),
                       nbl->sci_nalloc*sizeof(*nbl->sci),
                       nbl->alloc, nbl->free);
}

/* Make a new ci entry at index nbl->nci */
static void new_ci_entry(nbnxn_pairlist_t *nbl, int ci, int shift, int flags)
{
    if (nbl->nci + 1 > nbl->ci_nalloc)
    {
        nb_realloc_ci(nbl, nbl->nci+1);
    }
    nbl->ci[nbl->nci].ci            = ci;
    nbl->ci[nbl->nci].shift         = shift;
    /* Store the interaction flags along with the shift */
    nbl->ci[nbl->nci].shift        |= flags;
    nbl->ci[nbl->nci].cj_ind_start  = nbl->ncj;
    nbl->ci[nbl->nci].cj_ind_end    = nbl->ncj;
}

/* Make a new sci entry at index nbl->nsci */
static void new_sci_entry(nbnxn_pairlist_t *nbl, int sci, int shift)
{
    if (nbl->nsci + 1 > nbl->sci_nalloc)
    {
        nb_realloc_sci(nbl, nbl->nsci+1);
    }
    nbl->sci[nbl->nsci].sci           = sci;
    nbl->sci[nbl->nsci].shift         = shift;
    nbl->sci[nbl->nsci].cj4_ind_start = nbl->ncj4;
    nbl->sci[nbl->nsci].cj4_ind_end   = nbl->ncj4;
}

/* Sort the simple j-list cj on exclusions.
 * Entries with exclusions will all be sorted to the beginning of the list.
 */
static void sort_cj_excl(nbnxn_cj_t *cj, int ncj,
                         nbnxn_list_work_t *work)
{
    int jnew;

    if (ncj > work->cj_nalloc)
    {
        work->cj_nalloc = over_alloc_large(ncj);
        srenew(work->cj, work->cj_nalloc);
    }

    /* Make a list of the j-cells involving exclusions */
    jnew = 0;
    for (int j = 0; j < ncj; j++)
    {
        if (cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) ||
          (jnew == 1 && cj[0].excl != NBNXN_INTERACTION_MASK_ALL)))
    {
        for (int j = 0; j < ncj; j++)
        {
            if (cj[j].excl == NBNXN_INTERACTION_MASK_ALL)
            {
                work->cj[jnew++] = cj[j];
            }
        }
        for (int j = 0; j < ncj; j++)
        {
            cj[j] = work->cj[j];
        }
    }
}

/* Close this simple list i entry */
static void close_ci_entry_simple(nbnxn_pairlist_t *nbl)
{
    int jlen;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    jlen = nbl->ci[nbl->nci].cj_ind_end - nbl->ci[nbl->nci].cj_ind_start;
    if (jlen > 0)
    {
        sort_cj_excl(nbl->cj+nbl->ci[nbl->nci].cj_ind_start, jlen, nbl->work);

        /* The counts below are used for non-bonded pair/flop counts
         * and should therefore match the available kernel setups.
         */
        if (!(nbl->ci[nbl->nci].shift & NBNXN_CI_DO_COUL(0)))
        {
            nbl->work->ncj_noq += jlen;
        }
        else if ((nbl->ci[nbl->nci].shift & NBNXN_CI_HALF_LJ(0)) ||
                 !(nbl->ci[nbl->nci].shift & NBNXN_CI_DO_LJ(0)))
        {
            nbl->work->ncj_hlj += jlen;
        }

        nbl->nci++;
    }
}

/* Split sci entry for load balancing on the GPU.
 * Splitting ensures we have enough lists to fully utilize the whole GPU.
 * With progBal we generate progressively smaller lists, which improves
 * load balancing. As we only know the current count on our own thread,
 * we will need to estimate the current total amount of i-entries.
 * As the lists get concatenated later, this estimate depends
 * both on nthread and our own thread index.
 */
static void split_sci_entry(nbnxn_pairlist_t *nbl,
                            int nsp_target_av,
                            gmx_bool progBal, float nsp_tot_est,
                            int thread, int nthread)
{
    int nsp_max;
    int cj4_start, cj4_end, j4len;
    int sci;
    int nsp, nsp_sci, nsp_cj4, nsp_cj4_e, nsp_cj4_p;

    if (progBal)
    {
        float nsp_est;

        /* Estimate the total numbers of ci's of the nblist combined
         * over all threads using the target number of ci's.
         */
        nsp_est = (nsp_tot_est*thread)/nthread + nbl->nci_tot;

        /* The first ci blocks should be larger, to avoid overhead.
         * The last ci blocks should be smaller, to improve load balancing.
         * The factor 3/2 makes the first block 3/2 times the target average
         * and ensures that the total number of blocks end up equal to
         * that of equally sized blocks of size nsp_target_av.
         */
        nsp_max = static_cast<int>(nsp_target_av*(nsp_tot_est*1.5/(nsp_est + nsp_tot_est)));
    }
    else
    {
        nsp_max = nsp_target_av;
    }

    cj4_start = nbl->sci[nbl->nsci-1].cj4_ind_start;
    cj4_end   = nbl->sci[nbl->nsci-1].cj4_ind_end;
    j4len     = cj4_end - cj4_start;

    if (j4len > 1 && j4len*c_gpuNumClusterPerCell*c_nbnxnGpuJgroupSize > nsp_max)
    {
        /* Remove the last ci entry and process the cj4's again */
        nbl->nsci -= 1;

        sci        = nbl->nsci;
        nsp        = 0;
        nsp_sci    = 0;
        nsp_cj4_e  = 0;
        nsp_cj4    = 0;
        for (int cj4 = cj4_start; cj4 < cj4_end; cj4++)
        {
            nsp_cj4_p = nsp_cj4;
            /* Count the number of cluster pairs in this cj4 group */
            nsp_cj4   = 0;
            for (int p = 0; p < c_gpuNumClusterPerCell*c_nbnxnGpuJgroupSize; p++)
            {
                nsp_cj4 += (nbl->cj4[cj4].imei[0].imask >> p) & 1;
            }

            /* If adding the current cj4 with nsp_cj4 pairs get us further
             * away from our target nsp_max, split the list before this cj4.
             */
            if (nsp > 0 && nsp_max - nsp < nsp + nsp_cj4 - nsp_max)
            {
                /* Split the list at cj4 */
                nbl->sci[sci].cj4_ind_end = cj4;
                /* Create a new sci entry */
                sci++;
                nbl->nsci++;
                if (nbl->nsci+1 > nbl->sci_nalloc)
                {
                    nb_realloc_sci(nbl, nbl->nsci+1);
                }
                nbl->sci[sci].sci           = nbl->sci[nbl->nsci-1].sci;
                nbl->sci[sci].shift         = nbl->sci[nbl->nsci-1].shift;
                nbl->sci[sci].cj4_ind_start = cj4;
                nsp_sci                     = nsp;
                nsp_cj4_e                   = nsp_cj4_p;
                nsp                         = 0;
            }
            nsp += nsp_cj4;
        }

        /* Put the remaining cj4's in the last sci entry */
        nbl->sci[sci].cj4_ind_end = cj4_end;

        /* Possibly balance out the last two sci's
         * by moving the last cj4 of the second last sci.
         */
        if (nsp_sci - nsp_cj4_e >= nsp + nsp_cj4_e)
        {
            nbl->sci[sci-1].cj4_ind_end--;
            nbl->sci[sci].cj4_ind_start--;
        }

        nbl->nsci++;
    }
}

/* Clost this super/sub list i entry */
static void close_ci_entry_supersub(nbnxn_pairlist_t *nbl,
                                    int nsp_max_av,
                                    gmx_bool progBal, float nsp_tot_est,
                                    int thread, int nthread)
{
    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    int j4len = nbl->sci[nbl->nsci].cj4_ind_end - nbl->sci[nbl->nsci].cj4_ind_start;
    if (j4len > 0)
    {
        /* We can only have complete blocks of 4 j-entries in a list,
         * so round the count up before closing.
         */
        nbl->ncj4         = (nbl->work->cj_ind + c_nbnxnGpuJgroupSize - 1)/c_nbnxnGpuJgroupSize;
        nbl->work->cj_ind = nbl->ncj4*c_nbnxnGpuJgroupSize;

        nbl->nsci++;

        if (nsp_max_av > 0)
        {
            /* Measure the size of the new entry and potentially split it */
            split_sci_entry(nbl, nsp_max_av, progBal, nsp_tot_est,
                            thread, nthread);
        }
    }
}

/* Syncs the working array before adding another grid pair to the list */
static void sync_work(nbnxn_pairlist_t *nbl)
{
    if (!nbl->bSimple)
    {
        nbl->work->cj_ind   = nbl->ncj4*c_nbnxnGpuJgroupSize;
        nbl->work->cj4_init = nbl->ncj4;
    }
}

/* Clears an nbnxn_pairlist_t data structure */
static void clear_pairlist(nbnxn_pairlist_t *nbl)
{
    nbl->nci           = 0;
    nbl->nsci          = 0;
    nbl->ncj           = 0;
    nbl->ncjInUse      = 0;
    nbl->ncj4          = 0;
    nbl->nci_tot       = 0;
    nbl->nciOuter      = -1;
    nbl->nexcl         = 1;

    nbl->work->ncj_noq = 0;
    nbl->work->ncj_hlj = 0;
}

/* Clears a group scheme pair list */
static void clear_pairlist_fep(t_nblist *nl)
{
    nl->nri = 0;
    nl->nrj = 0;
    if (nl->jindex == nullptr)
    {
        snew(nl->jindex, 1);
    }
    nl->jindex[0] = 0;
}

/* Sets a simple list i-cell bounding box, including PBC shift */
static gmx_inline void set_icell_bb_simple(const nbnxn_bb_t *bb, int ci,
                                           real shx, real shy, real shz,
                                           nbnxn_bb_t *bb_ci)
{
    bb_ci->lower[BB_X] = bb[ci].lower[BB_X] + shx;
    bb_ci->lower[BB_Y] = bb[ci].lower[BB_Y] + shy;
    bb_ci->lower[BB_Z] = bb[ci].lower[BB_Z] + shz;
    bb_ci->upper[BB_X] = bb[ci].upper[BB_X] + shx;
    bb_ci->upper[BB_Y] = bb[ci].upper[BB_Y] + shy;
    bb_ci->upper[BB_Z] = bb[ci].upper[BB_Z] + shz;
}

#if NBNXN_BBXXXX
/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
static void set_icell_bbxxxx_supersub(const float *bb, int ci,
                                      real shx, real shy, real shz,
                                      float *bb_ci)
{
    int ia = ci*(c_gpuNumClusterPerCell >> STRIDE_PBB_2LOG)*NNBSBB_XXXX;
    for (int m = 0; m < (c_gpuNumClusterPerCell >> STRIDE_PBB_2LOG)*NNBSBB_XXXX; m += NNBSBB_XXXX)
    {
        for (int i = 0; i < STRIDE_PBB; i++)
        {
            bb_ci[m+0*STRIDE_PBB+i] = bb[ia+m+0*STRIDE_PBB+i] + shx;
            bb_ci[m+1*STRIDE_PBB+i] = bb[ia+m+1*STRIDE_PBB+i] + shy;
            bb_ci[m+2*STRIDE_PBB+i] = bb[ia+m+2*STRIDE_PBB+i] + shz;
            bb_ci[m+3*STRIDE_PBB+i] = bb[ia+m+3*STRIDE_PBB+i] + shx;
            bb_ci[m+4*STRIDE_PBB+i] = bb[ia+m+4*STRIDE_PBB+i] + shy;
            bb_ci[m+5*STRIDE_PBB+i] = bb[ia+m+5*STRIDE_PBB+i] + shz;
        }
    }
}
#endif

/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
gmx_unused static void set_icell_bb_supersub(const nbnxn_bb_t *bb, int ci,
                                             real shx, real shy, real shz,
                                             nbnxn_bb_t *bb_ci)
{
    for (int i = 0; i < c_gpuNumClusterPerCell; i++)
    {
        set_icell_bb_simple(bb, ci*c_gpuNumClusterPerCell+i,
                            shx, shy, shz,
                            &bb_ci[i]);
    }
}

/* Copies PBC shifted i-cell atom coordinates x,y,z to working array */
static void icell_set_x_simple(int ci,
                               real shx, real shy, real shz,
                               int stride, const real *x,
                               nbnxn_list_work_t *work)
{
    int ia = ci*NBNXN_CPU_CLUSTER_I_SIZE;

    for (int i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        work->x_ci[i*STRIDE_XYZ+XX] = x[(ia+i)*stride+XX] + shx;
        work->x_ci[i*STRIDE_XYZ+YY] = x[(ia+i)*stride+YY] + shy;
        work->x_ci[i*STRIDE_XYZ+ZZ] = x[(ia+i)*stride+ZZ] + shz;
    }
}

/* Copies PBC shifted super-cell atom coordinates x,y,z to working array */
static void icell_set_x_supersub(int ci,
                                 real shx, real shy, real shz,
                                 int stride, const real *x,
                                 nbnxn_list_work_t *work)
{
#if !GMX_SIMD4_HAVE_REAL

    real * x_ci = work->x_ci;

    int    ia = ci*c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize;
    for (int i = 0; i < c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize; i++)
    {
        x_ci[i*DIM + XX] = x[(ia+i)*stride + XX] + shx;
        x_ci[i*DIM + YY] = x[(ia+i)*stride + YY] + shy;
        x_ci[i*DIM + ZZ] = x[(ia+i)*stride + ZZ] + shz;
    }

#else /* !GMX_SIMD4_HAVE_REAL */

    real * x_ci = work->x_ci_simd;

    for (int si = 0; si < c_gpuNumClusterPerCell; si++)
    {
        for (int i = 0; i < c_nbnxnGpuClusterSize; i += GMX_SIMD4_WIDTH)
        {
            int io = si*c_nbnxnGpuClusterSize + i;
            int ia = ci*c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize + io;
            for (int j = 0; j < GMX_SIMD4_WIDTH; j++)
            {
                x_ci[io*DIM + j + XX*GMX_SIMD4_WIDTH] = x[(ia + j)*stride + XX] + shx;
                x_ci[io*DIM + j + YY*GMX_SIMD4_WIDTH] = x[(ia + j)*stride + YY] + shy;
                x_ci[io*DIM + j + ZZ*GMX_SIMD4_WIDTH] = x[(ia + j)*stride + ZZ] + shz;
            }
        }
    }

#endif /* !GMX_SIMD4_HAVE_REAL */
}

static real minimum_subgrid_size_xy(const nbnxn_grid_t *grid)
{
    if (grid->bSimple)
    {
        return std::min(grid->sx, grid->sy);
    }
    else
    {
        return std::min(grid->sx/c_gpuNumClusterPerCellX,
                        grid->sy/c_gpuNumClusterPerCellY);
    }
}

static real effective_buffer_1x1_vs_MxN(const nbnxn_grid_t *gridi,
                                        const nbnxn_grid_t *gridj)
{
    const real eff_1x1_buffer_fac_overest = 0.1;

    /* Determine an atom-pair list cut-off buffer size for atom pairs,
     * to be added to rlist (including buffer) used for MxN.
     * This is for converting an MxN list to a 1x1 list. This means we can't
     * use the normal buffer estimate, as we have an MxN list in which
     * some atom pairs beyond rlist are missing. We want to capture
     * the beneficial effect of buffering by extra pairs just outside rlist,
     * while removing the useless pairs that are further away from rlist.
     * (Also the buffer could have been set manually not using the estimate.)
     * This buffer size is an overestimate.
     * We add 10% of the smallest grid sub-cell dimensions.
     * Note that the z-size differs per cell and we don't use this,
     * so we overestimate.
     * With PME, the 10% value gives a buffer that is somewhat larger
     * than the effective buffer with a tolerance of 0.005 kJ/mol/ps.
     * Smaller tolerances or using RF lead to a smaller effective buffer,
     * so 10% gives a safe overestimate.
     */
    return eff_1x1_buffer_fac_overest*(minimum_subgrid_size_xy(gridi) +
                                       minimum_subgrid_size_xy(gridj));
}

/* Clusters at the cut-off only increase rlist by 60% of their size */
static real nbnxn_rlist_inc_outside_fac = 0.6;

/* Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int cluster_size_j, real atom_density)
{
    int  cluster_size_i;
    real vol_inc_i, vol_inc_j;

    /* We should get this from the setup, but currently it's the same for
     * all setups, including GPUs.
     */
    cluster_size_i = NBNXN_CPU_CLUSTER_I_SIZE;

    vol_inc_i = (cluster_size_i - 1)/atom_density;
    vol_inc_j = (cluster_size_j - 1)/atom_density;

    return nbnxn_rlist_inc_outside_fac*std::cbrt(vol_inc_i + vol_inc_j);
}

/* Estimates the interaction volume^2 for non-local interactions */
static real nonlocal_vol2(const struct gmx_domdec_zones_t *zones, rvec ls, real r)
{
    real cl, ca, za;
    real vold_est;
    real vol2_est_tot;

    vol2_est_tot = 0;

    /* Here we simply add up the volumes of 1, 2 or 3 1D decomposition
     * not home interaction volume^2. As these volumes are not additive,
     * this is an overestimate, but it would only be significant in the limit
     * of small cells, where we anyhow need to split the lists into
     * as small parts as possible.
     */

    for (int z = 0; z < zones->n; z++)
    {
        if (zones->shift[z][XX] + zones->shift[z][YY] + zones->shift[z][ZZ] == 1)
        {
            cl = 0;
            ca = 1;
            za = 1;
            for (int d = 0; d < DIM; d++)
            {
                if (zones->shift[z][d] == 0)
                {
                    cl += 0.5*ls[d];
                    ca *= ls[d];
                    za *= zones->size[z].x1[d] - zones->size[z].x0[d];
                }
            }

            /* 4 octants of a sphere */
            vold_est  = 0.25*M_PI*r*r*r*r;
            /* 4 quarter pie slices on the edges */
            vold_est += 4*cl*M_PI/6.0*r*r*r;
            /* One rectangular volume on a face */
            vold_est += ca*0.5*r*r;

            vol2_est_tot += vold_est*za;
        }
    }

    return vol2_est_tot;
}

/* Estimates the average size of a full j-list for super/sub setup */
static void get_nsubpair_target(const nbnxn_search_t  nbs,
                                int                   iloc,
                                real                  rlist,
                                int                   min_ci_balanced,
                                int                  *nsubpair_target,
                                float                *nsubpair_tot_est)
{
    /* The target value of 36 seems to be the optimum for Kepler.
     * Maxwell is less sensitive to the exact value.
     */
    const int           nsubpair_target_min = 36;
    const nbnxn_grid_t *grid;
    rvec                ls;
    real                r_eff_sup, vol_est, nsp_est, nsp_est_nl;

    grid = &nbs->grid[0];

    /* We don't need to balance list sizes if:
     * - We didn't request balancing.
     * - The number of grid cells >= the number of lists requested,
     *   since we will always generate at least #cells lists.
     * - We don't have any cells, since then there won't be any lists.
     */
    if (min_ci_balanced <= 0 || grid->nc >= min_ci_balanced || grid->nc == 0)
    {
        /* nsubpair_target==0 signals no balancing */
        *nsubpair_target  = 0;
        *nsubpair_tot_est = 0;

        return;
    }

    ls[XX] = (grid->c1[XX] - grid->c0[XX])/(grid->ncx*c_gpuNumClusterPerCellX);
    ls[YY] = (grid->c1[YY] - grid->c0[YY])/(grid->ncy*c_gpuNumClusterPerCellY);
    ls[ZZ] = grid->na_c/(grid->atom_density*ls[XX]*ls[YY]);

    /* The average length of the diagonal of a sub cell */
    real diagonal = std::sqrt(ls[XX]*ls[XX] + ls[YY]*ls[YY] + ls[ZZ]*ls[ZZ]);

    /* The formulas below are a heuristic estimate of the average nsj per si*/
    r_eff_sup = rlist + nbnxn_rlist_inc_outside_fac*gmx::square((grid->na_c - 1.0)/grid->na_c)*0.5*diagonal;

    if (!nbs->DomDec || nbs->zones->n == 1)
    {
        nsp_est_nl = 0;
    }
    else
    {
        nsp_est_nl =
            gmx::square(grid->atom_density/grid->na_c)*
            nonlocal_vol2(nbs->zones, ls, r_eff_sup);
    }

    if (LOCAL_I(iloc))
    {
        /* Sub-cell interacts with itself */
        vol_est  = ls[XX]*ls[YY]*ls[ZZ];
        /* 6/2 rectangular volume on the faces */
        vol_est += (ls[XX]*ls[YY] + ls[XX]*ls[ZZ] + ls[YY]*ls[ZZ])*r_eff_sup;
        /* 12/2 quarter pie slices on the edges */
        vol_est += 2*(ls[XX] + ls[YY] + ls[ZZ])*0.25*M_PI*gmx::square(r_eff_sup);
        /* 4 octants of a sphere */
        vol_est += 0.5*4.0/3.0*M_PI*gmx::power3(r_eff_sup);

        /* Estimate the number of cluster pairs as the local number of
         * clusters times the volume they interact with times the density.
         */
        nsp_est = grid->nsubc_tot*vol_est*grid->atom_density/grid->na_c;

        /* Subtract the non-local pair count */
        nsp_est -= nsp_est_nl;

        /* For small cut-offs nsp_est will be an underesimate.
         * With DD nsp_est_nl is an overestimate so nsp_est can get negative.
         * So to avoid too small or negative nsp_est we set a minimum of
         * all cells interacting with all 3^3 direct neighbors (3^3-1)/2+1=14.
         * This might be a slight overestimate for small non-periodic groups of
         * atoms as will occur for a local domain with DD, but for small
         * groups of atoms we'll anyhow be limited by nsubpair_target_min,
         * so this overestimation will not matter.
         */
        nsp_est = std::max(nsp_est, grid->nsubc_tot*static_cast<real>(14));

        if (debug)
        {
            fprintf(debug, "nsp_est local %5.1f non-local %5.1f\n",
                    nsp_est, nsp_est_nl);
        }
    }
    else
    {
        nsp_est = nsp_est_nl;
    }

    /* Thus the (average) maximum j-list size should be as follows.
     * Since there is overhead, we shouldn't make the lists too small
     * (and we can't chop up j-groups) so we use a minimum target size of 36.
     */
    *nsubpair_target  = std::max(nsubpair_target_min,
                                 static_cast<int>(nsp_est/min_ci_balanced + 0.5));
    *nsubpair_tot_est = static_cast<int>(nsp_est);

    if (debug)
    {
        fprintf(debug, "nbl nsp estimate %.1f, nsubpair_target %d\n",
                nsp_est, *nsubpair_target);
    }
}

/* Debug list print function */
static void print_nblist_ci_cj(FILE *fp, const nbnxn_pairlist_t *nbl)
{
    for (int i = 0; i < nbl->nci; i++)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj %3d\n",
                nbl->ci[i].ci, nbl->ci[i].shift,
                nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start);

        for (int j = nbl->ci[i].cj_ind_start; j < nbl->ci[i].cj_ind_end; j++)
        {
            fprintf(fp, "  cj %5d  imask %x\n",
                    nbl->cj[j].cj,
                    nbl->cj[j].excl);
        }
    }
}

/* Debug list print function */
static void print_nblist_sci_cj(FILE *fp, const nbnxn_pairlist_t *nbl)
{
    for (int i = 0; i < nbl->nsci; i++)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d\n",
                nbl->sci[i].sci, nbl->sci[i].shift,
                nbl->sci[i].cj4_ind_end - nbl->sci[i].cj4_ind_start);

        int ncp = 0;
        for (int j4 = nbl->sci[i].cj4_ind_start; j4 < nbl->sci[i].cj4_ind_end; j4++)
        {
            for (int j = 0; j < c_nbnxnGpuJgroupSize; j++)
            {
                fprintf(fp, "  sj %5d  imask %x\n",
                        nbl->cj4[j4].cj[j],
                        nbl->cj4[j4].imei[0].imask);
                for (int si = 0; si < c_gpuNumClusterPerCell; si++)
                {
                    if (nbl->cj4[j4].imei[0].imask & (1U << (j*c_gpuNumClusterPerCell + si)))
                    {
                        ncp++;
                    }
                }
            }
        }
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d ncp %3d\n",
                nbl->sci[i].sci, nbl->sci[i].shift,
                nbl->sci[i].cj4_ind_end - nbl->sci[i].cj4_ind_start,
                ncp);
    }
}

/* Combine pair lists *nbl generated on multiple threads nblc */
static void combine_nblists(int nnbl, nbnxn_pairlist_t **nbl,
                            nbnxn_pairlist_t *nblc)
{
    int nsci, ncj4, nexcl;

    if (nblc->bSimple)
    {
        gmx_incons("combine_nblists does not support simple lists");
    }

    nsci  = nblc->nsci;
    ncj4  = nblc->ncj4;
    nexcl = nblc->nexcl;
    for (int i = 0; i < nnbl; i++)
    {
        nsci  += nbl[i]->nsci;
        ncj4  += nbl[i]->ncj4;
        nexcl += nbl[i]->nexcl;
    }

    if (nsci > nblc->sci_nalloc)
    {
        nb_realloc_sci(nblc, nsci);
    }
    if (ncj4 > nblc->cj4_nalloc)
    {
        nblc->cj4_nalloc = over_alloc_small(ncj4);
        nbnxn_realloc_void((void **)&nblc->cj4,
                           nblc->ncj4*sizeof(*nblc->cj4),
                           nblc->cj4_nalloc*sizeof(*nblc->cj4),
                           nblc->alloc, nblc->free);
    }
    if (nexcl > nblc->excl_nalloc)
    {
        nblc->excl_nalloc = over_alloc_small(nexcl);
        nbnxn_realloc_void((void **)&nblc->excl,
                           nblc->nexcl*sizeof(*nblc->excl),
                           nblc->excl_nalloc*sizeof(*nblc->excl),
                           nblc->alloc, nblc->free);
    }

    /* Each thread should copy its own data to the combined arrays,
     * as otherwise data will go back and forth between different caches.
     */
#if GMX_OPENMP && !(defined __clang_analyzer__)
    // cppcheck-suppress unreadVariable
    int nthreads = gmx_omp_nthreads_get(emntPairsearch);
#endif

#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int n = 0; n < nnbl; n++)
    {
        try
        {
            int                     sci_offset;
            int                     cj4_offset;
            int                     excl_offset;
            const nbnxn_pairlist_t *nbli;

            /* Determine the offset in the combined data for our thread */
            sci_offset  = nblc->nsci;
            cj4_offset  = nblc->ncj4;
            excl_offset = nblc->nexcl;

            for (int i = 0; i < n; i++)
            {
                sci_offset  += nbl[i]->nsci;
                cj4_offset  += nbl[i]->ncj4;
                excl_offset += nbl[i]->nexcl;
            }

            nbli = nbl[n];

            for (int i = 0; i < nbli->nsci; i++)
            {
                nblc->sci[sci_offset+i]                = nbli->sci[i];
                nblc->sci[sci_offset+i].cj4_ind_start += cj4_offset;
                nblc->sci[sci_offset+i].cj4_ind_end   += cj4_offset;
            }

            for (int j4 = 0; j4 < nbli->ncj4; j4++)
            {
                nblc->cj4[cj4_offset+j4]                   = nbli->cj4[j4];
                nblc->cj4[cj4_offset+j4].imei[0].excl_ind += excl_offset;
                nblc->cj4[cj4_offset+j4].imei[1].excl_ind += excl_offset;
            }

            for (int j4 = 0; j4 < nbli->nexcl; j4++)
            {
                nblc->excl[excl_offset+j4] = nbli->excl[j4];
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    for (int n = 0; n < nnbl; n++)
    {
        nblc->nsci    += nbl[n]->nsci;
        nblc->ncj4    += nbl[n]->ncj4;
        nblc->nci_tot += nbl[n]->nci_tot;
        nblc->nexcl   += nbl[n]->nexcl;
    }
}

static void balance_fep_lists(const nbnxn_search_t  nbs,
                              nbnxn_pairlist_set_t *nbl_lists)
{
    int       nnbl;
    int       nri_tot, nrj_tot, nrj_target;
    int       th_dest;
    t_nblist *nbld;

    nnbl = nbl_lists->nnbl;

    if (nnbl == 1)
    {
        /* Nothing to balance */
        return;
    }

    /* Count the total i-lists and pairs */
    nri_tot = 0;
    nrj_tot = 0;
    for (int th = 0; th < nnbl; th++)
    {
        nri_tot += nbl_lists->nbl_fep[th]->nri;
        nrj_tot += nbl_lists->nbl_fep[th]->nrj;
    }

    nrj_target = (nrj_tot + nnbl - 1)/nnbl;

    assert(gmx_omp_nthreads_get(emntNonbonded) == nnbl);

#pragma omp parallel for schedule(static) num_threads(nnbl)
    for (int th = 0; th < nnbl; th++)
    {
        try
        {
            t_nblist *nbl;

            nbl = nbs->work[th].nbl_fep;

            /* Note that here we allocate for the total size, instead of
             * a per-thread esimate (which is hard to obtain).
             */
            if (nri_tot > nbl->maxnri)
            {
                nbl->maxnri = over_alloc_large(nri_tot);
                reallocate_nblist(nbl);
            }
            if (nri_tot > nbl->maxnri || nrj_tot > nbl->maxnrj)
            {
                nbl->maxnrj = over_alloc_small(nrj_tot);
                srenew(nbl->jjnr, nbl->maxnrj);
                srenew(nbl->excl_fep, nbl->maxnrj);
            }

            clear_pairlist_fep(nbl);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    /* Loop over the source lists and assign and copy i-entries */
    th_dest = 0;
    nbld    = nbs->work[th_dest].nbl_fep;
    for (int th = 0; th < nnbl; th++)
    {
        t_nblist *nbls;

        nbls = nbl_lists->nbl_fep[th];

        for (int i = 0; i < nbls->nri; i++)
        {
            int nrj;

            /* The number of pairs in this i-entry */
            nrj = nbls->jindex[i+1] - nbls->jindex[i];

            /* Decide if list th_dest is too large and we should procede
             * to the next destination list.
             */
            if (th_dest+1 < nnbl && nbld->nrj > 0 &&
                nbld->nrj + nrj - nrj_target > nrj_target - nbld->nrj)
            {
                th_dest++;
                nbld = nbs->work[th_dest].nbl_fep;
            }

            nbld->iinr[nbld->nri]  = nbls->iinr[i];
            nbld->gid[nbld->nri]   = nbls->gid[i];
            nbld->shift[nbld->nri] = nbls->shift[i];

            for (int j = nbls->jindex[i]; j < nbls->jindex[i+1]; j++)
            {
                nbld->jjnr[nbld->nrj]     = nbls->jjnr[j];
                nbld->excl_fep[nbld->nrj] = nbls->excl_fep[j];
                nbld->nrj++;
            }
            nbld->nri++;
            nbld->jindex[nbld->nri] = nbld->nrj;
        }
    }

    /* Swap the list pointers */
    for (int th = 0; th < nnbl; th++)
    {
        t_nblist *nbl_tmp;

        nbl_tmp                = nbl_lists->nbl_fep[th];
        nbl_lists->nbl_fep[th] = nbs->work[th].nbl_fep;
        nbs->work[th].nbl_fep  = nbl_tmp;

        if (debug)
        {
            fprintf(debug, "nbl_fep[%d] nri %4d nrj %4d\n",
                    th,
                    nbl_lists->nbl_fep[th]->nri,
                    nbl_lists->nbl_fep[th]->nrj);
        }
    }
}

/* Returns the next ci to be processes by our thread */
static gmx_bool next_ci(const nbnxn_grid_t *grid,
                        int nth, int ci_block,
                        int *ci_x, int *ci_y,
                        int *ci_b, int *ci)
{
    (*ci_b)++;
    (*ci)++;

    if (*ci_b == ci_block)
    {
        /* Jump to the next block assigned to this task */
        *ci   += (nth - 1)*ci_block;
        *ci_b  = 0;
    }

    if (*ci >= grid->nc)
    {
        return FALSE;
    }

    while (*ci >= grid->cxy_ind[*ci_x*grid->ncy + *ci_y + 1])
    {
        *ci_y += 1;
        if (*ci_y == grid->ncy)
        {
            *ci_x += 1;
            *ci_y  = 0;
        }
    }

    return TRUE;
}

/* Returns the distance^2 for which we put cell pairs in the list
 * without checking atom pair distances. This is usually < rlist^2.
 */
static float boundingbox_only_distance2(const nbnxn_grid_t *gridi,
                                        const nbnxn_grid_t *gridj,
                                        real                rlist,
                                        gmx_bool            simple)
{
    /* If the distance between two sub-cell bounding boxes is less
     * than this distance, do not check the distance between
     * all particle pairs in the sub-cell, since then it is likely
     * that the box pair has atom pairs within the cut-off.
     * We use the nblist cut-off minus 0.5 times the average x/y diagonal
     * spacing of the sub-cells. Around 40% of the checked pairs are pruned.
     * Using more than 0.5 gains at most 0.5%.
     * If forces are calculated more than twice, the performance gain
     * in the force calculation outweighs the cost of checking.
     * Note that with subcell lists, the atom-pair distance check
     * is only performed when only 1 out of 8 sub-cells in within range,
     * this is because the GPU is much faster than the cpu.
     */
    real bbx, bby;
    real rbb2;

    bbx = 0.5*(gridi->sx + gridj->sx);
    bby = 0.5*(gridi->sy + gridj->sy);
    if (!simple)
    {
        bbx /= c_gpuNumClusterPerCellX;
        bby /= c_gpuNumClusterPerCellY;
    }

    rbb2 = std::max(0.0, rlist - 0.5*std::sqrt(bbx*bbx + bby*bby));
    rbb2 = rbb2 * rbb2;

#if !GMX_DOUBLE
    return rbb2;
#else
    return (float)((1+GMX_FLOAT_EPS)*rbb2);
#endif
}

static int get_ci_block_size(const nbnxn_grid_t *gridi,
                             gmx_bool bDomDec, int nth)
{
    const int ci_block_enum      = 5;
    const int ci_block_denom     = 11;
    const int ci_block_min_atoms = 16;
    int       ci_block;

    /* Here we decide how to distribute the blocks over the threads.
     * We use prime numbers to try to avoid that the grid size becomes
     * a multiple of the number of threads, which would lead to some
     * threads getting "inner" pairs and others getting boundary pairs,
     * which in turns will lead to load imbalance between threads.
     * Set the block size as 5/11/ntask times the average number of cells
     * in a y,z slab. This should ensure a quite uniform distribution
     * of the grid parts of the different thread along all three grid
     * zone boundaries with 3D domain decomposition. At the same time
     * the blocks will not become too small.
     */
    ci_block = (gridi->nc*ci_block_enum)/(ci_block_denom*gridi->ncx*nth);

    /* Ensure the blocks are not too small: avoids cache invalidation */
    if (ci_block*gridi->na_sc < ci_block_min_atoms)
    {
        ci_block = (ci_block_min_atoms + gridi->na_sc - 1)/gridi->na_sc;
    }

    /* Without domain decomposition
     * or with less than 3 blocks per task, divide in nth blocks.
     */
    if (!bDomDec || nth*3*ci_block > gridi->nc)
    {
        ci_block = (gridi->nc + nth - 1)/nth;
    }

    if (ci_block > 1 && (nth - 1)*ci_block >= gridi->nc)
    {
        /* Some threads have no work. Although reducing the block size
         * does not decrease the block count on the first few threads,
         * with GPUs better mixing of "upper" cells that have more empty
         * clusters results in a somewhat lower max load over all threads.
         * Without GPUs the regime of so few atoms per thread is less
         * performance relevant, but with 8-wide SIMD the same reasoning
         * applies, since the pair list uses 4 i-atom "sub-clusters".
         */
        ci_block--;
    }

    return ci_block;
}

/* Returns the number of bits to right-shift a cluster index to obtain
 * the corresponding force buffer flag index.
 */
static int getBufferFlagShift(int numAtomsPerCluster)
{
    int bufferFlagShift = 0;
    while ((numAtomsPerCluster << bufferFlagShift) < NBNXN_BUFFERFLAG_SIZE)
    {
        bufferFlagShift++;
    }

    return bufferFlagShift;
}

/* Generates the part of pair-list nbl assigned to our thread */
static void nbnxn_make_pairlist_part(const nbnxn_search_t nbs,
                                     const nbnxn_grid_t *gridi,
                                     const nbnxn_grid_t *gridj,
                                     nbnxn_search_work_t *work,
                                     const nbnxn_atomdata_t *nbat,
                                     const t_blocka &exclusions,
                                     real rlist,
                                     int nb_kernel_type,
                                     int ci_block,
                                     gmx_bool bFBufferFlag,
                                     int nsubpair_max,
                                     gmx_bool progBal,
                                     float nsubpair_tot_est,
                                     int th, int nth,
                                     nbnxn_pairlist_t *nbl,
                                     t_nblist *nbl_fep)
{
    int               na_cj_2log;
    matrix            box;
    real              rlist2, rl_fep2 = 0;
    float             rbb2;
    int               ci_b, ci, ci_x, ci_y, ci_xy, cj;
    ivec              shp;
    int               shift;
    real              shx, shy, shz;
    int               cell0_i;
    const nbnxn_bb_t *bb_i = nullptr;
#if NBNXN_BBXXXX
    const float      *pbb_i = nullptr;
#endif
    const float      *bbcz_i, *bbcz_j;
    const int        *flags_i;
    real              bx0, bx1, by0, by1, bz0, bz1;
    real              bz1_frac;
    real              d2cx, d2z, d2z_cx, d2z_cy, d2zx, d2zxy, d2xy;
    int               cxf, cxl, cyf, cyf_x, cyl;
    int               numDistanceChecks;
    int               gridi_flag_shift = 0, gridj_flag_shift = 0;
    gmx_bitmask_t    *gridj_flag       = nullptr;
    int               ncj_old_i, ncj_old_j;

    nbs_cycle_start(&work->cc[enbsCCsearch]);

    if (gridj->bSimple != nbl->bSimple)
    {
        gmx_incons("Grid incompatible with pair-list");
    }

    sync_work(nbl);
    nbl->na_sc = gridj->na_sc;
    nbl->na_ci = gridj->na_c;
    nbl->na_cj = nbnxn_kernel_to_cluster_j_size(nb_kernel_type);
    na_cj_2log = get_2log(nbl->na_cj);

    nbl->rlist  = rlist;

    if (bFBufferFlag)
    {
        /* Determine conversion of clusters to flag blocks */
        gridi_flag_shift = getBufferFlagShift(nbl->na_ci);
        gridj_flag_shift = getBufferFlagShift(nbl->na_cj);

        gridj_flag       = work->buffer_flags.flag;
    }

    copy_mat(nbs->box, box);

    rlist2 = nbl->rlist*nbl->rlist;

    if (nbs->bFEP && !nbl->bSimple)
    {
        /* Determine an atom-pair list cut-off distance for FEP atom pairs.
         * We should not simply use rlist, since then we would not have
         * the small, effective buffering of the NxN lists.
         * The buffer is on overestimate, but the resulting cost for pairs
         * beyond rlist is neglible compared to the FEP pairs within rlist.
         */
        rl_fep2 = nbl->rlist + effective_buffer_1x1_vs_MxN(gridi, gridj);

        if (debug)
        {
            fprintf(debug, "nbl_fep atom-pair rlist %f\n", rl_fep2);
        }
        rl_fep2 = rl_fep2*rl_fep2;
    }

    rbb2 = boundingbox_only_distance2(gridi, gridj, nbl->rlist, nbl->bSimple);

    if (debug)
    {
        fprintf(debug, "nbl bounding box only distance %f\n", std::sqrt(rbb2));
    }

    /* Set the shift range */
    for (int d = 0; d < DIM; d++)
    {
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        if (d >= ePBC2npbcdim(nbs->ePBC) || nbs->dd_dim[d])
        {
            shp[d] = 0;
        }
        else
        {
            if (d == XX &&
                box[XX][XX] - fabs(box[YY][XX]) - fabs(box[ZZ][XX]) < std::sqrt(rlist2))
            {
                shp[d] = 2;
            }
            else
            {
                shp[d] = 1;
            }
        }
    }

#if NBNXN_BBXXXX
    if (gridi->bSimple)
    {
        bb_i  = gridi->bb;
    }
    else
    {
        pbb_i = gridi->pbb;
    }
#else
    /* We use the normal bounding box format for both grid types */
    bb_i  = gridi->bb;
#endif
    bbcz_i  = gridi->bbcz;
    flags_i = gridi->flags;
    cell0_i = gridi->cell0;

    bbcz_j = gridj->bbcz;

    if (debug)
    {
        fprintf(debug, "nbl nc_i %d col.av. %.1f ci_block %d\n",
                gridi->nc, gridi->nc/(double)(gridi->ncx*gridi->ncy), ci_block);
    }

    numDistanceChecks = 0;

    /* Initially ci_b and ci to 1 before where we want them to start,
     * as they will both be incremented in next_ci.
     */
    ci_b = -1;
    ci   = th*ci_block - 1;
    ci_x = 0;
    ci_y = 0;
    while (next_ci(gridi, nth, ci_block, &ci_x, &ci_y, &ci_b, &ci))
    {
        if (nbl->bSimple && flags_i[ci] == 0)
        {
            continue;
        }

        ncj_old_i = nbl->ncj;

        d2cx = 0;
        if (gridj != gridi && shp[XX] == 0)
        {
            if (nbl->bSimple)
            {
                bx1 = bb_i[ci].upper[BB_X];
            }
            else
            {
                bx1 = gridi->c0[XX] + (ci_x+1)*gridi->sx;
            }
            if (bx1 < gridj->c0[XX])
            {
                d2cx = gmx::square(gridj->c0[XX] - bx1);

                if (d2cx >= rlist2)
                {
                    continue;
                }
            }
        }

        ci_xy = ci_x*gridi->ncy + ci_y;

        /* Loop over shift vectors in three dimensions */
        for (int tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
        {
            shz = tz*box[ZZ][ZZ];

            bz0 = bbcz_i[ci*NNBSBB_D  ] + shz;
            bz1 = bbcz_i[ci*NNBSBB_D+1] + shz;

            if (tz == 0)
            {
                d2z = 0;
            }
            else if (tz < 0)
            {
                d2z = gmx::square(bz1);
            }
            else
            {
                d2z = gmx::square(bz0 - box[ZZ][ZZ]);
            }

            d2z_cx = d2z + d2cx;

            if (d2z_cx >= rlist2)
            {
                continue;
            }

            bz1_frac = bz1/(gridi->cxy_ind[ci_xy+1] - gridi->cxy_ind[ci_xy]);
            if (bz1_frac < 0)
            {
                bz1_frac = 0;
            }
            /* The check with bz1_frac close to or larger than 1 comes later */

            for (int ty = -shp[YY]; ty <= shp[YY]; ty++)
            {
                shy = ty*box[YY][YY] + tz*box[ZZ][YY];

                if (nbl->bSimple)
                {
                    by0 = bb_i[ci].lower[BB_Y] + shy;
                    by1 = bb_i[ci].upper[BB_Y] + shy;
                }
                else
                {
                    by0 = gridi->c0[YY] + (ci_y  )*gridi->sy + shy;
                    by1 = gridi->c0[YY] + (ci_y+1)*gridi->sy + shy;
                }

                get_cell_range(by0, by1,
                               gridj->ncy, gridj->c0[YY], gridj->sy, gridj->inv_sy,
                               d2z_cx, rlist2,
                               &cyf, &cyl);

                if (cyf > cyl)
                {
                    continue;
                }

                d2z_cy = d2z;
                if (by1 < gridj->c0[YY])
                {
                    d2z_cy += gmx::square(gridj->c0[YY] - by1);
                }
                else if (by0 > gridj->c1[YY])
                {
                    d2z_cy += gmx::square(by0 - gridj->c1[YY]);
                }

                for (int tx = -shp[XX]; tx <= shp[XX]; tx++)
                {
                    shift = XYZ2IS(tx, ty, tz);

                    if (c_pbcShiftBackward && gridi == gridj && shift > CENTRAL)
                    {
                        continue;
                    }

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    if (nbl->bSimple)
                    {
                        bx0 = bb_i[ci].lower[BB_X] + shx;
                        bx1 = bb_i[ci].upper[BB_X] + shx;
                    }
                    else
                    {
                        bx0 = gridi->c0[XX] + (ci_x  )*gridi->sx + shx;
                        bx1 = gridi->c0[XX] + (ci_x+1)*gridi->sx + shx;
                    }

                    get_cell_range(bx0, bx1,
                                   gridj->ncx, gridj->c0[XX], gridj->sx, gridj->inv_sx,
                                   d2z_cy, rlist2,
                                   &cxf, &cxl);

                    if (cxf > cxl)
                    {
                        continue;
                    }

                    if (nbl->bSimple)
                    {
                        new_ci_entry(nbl, cell0_i+ci, shift, flags_i[ci]);
                    }
                    else
                    {
                        new_sci_entry(nbl, cell0_i+ci, shift);
                    }

                    if ((!c_pbcShiftBackward || (shift == CENTRAL &&
                                                 gridi == gridj)) &&
                        cxf < ci_x)
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    if (nbl->bSimple)
                    {
                        set_icell_bb_simple(bb_i, ci, shx, shy, shz,
                                            nbl->work->bb_ci);
                    }
                    else
                    {
#if NBNXN_BBXXXX
                        set_icell_bbxxxx_supersub(pbb_i, ci, shx, shy, shz,
                                                  nbl->work->pbb_ci);
#else
                        set_icell_bb_supersub(bb_i, ci, shx, shy, shz,
                                              nbl->work->bb_ci);
#endif
                    }

                    nbs->icell_set_x(cell0_i+ci, shx, shy, shz,
                                     nbat->xstride, nbat->x,
                                     nbl->work);

                    for (int cx = cxf; cx <= cxl; cx++)
                    {
                        d2zx = d2z;
                        if (gridj->c0[XX] + cx*gridj->sx > bx1)
                        {
                            d2zx += gmx::square(gridj->c0[XX] + cx*gridj->sx - bx1);
                        }
                        else if (gridj->c0[XX] + (cx+1)*gridj->sx < bx0)
                        {
                            d2zx += gmx::square(gridj->c0[XX] + (cx+1)*gridj->sx - bx0);
                        }

                        if (gridi == gridj &&
                            cx == 0 &&
                            (!c_pbcShiftBackward || shift == CENTRAL) &&
                            cyf < ci_y)
                        {
                            /* Leave the pairs with i > j.
                             * Skip half of y when i and j have the same x.
                             */
                            cyf_x = ci_y;
                        }
                        else
                        {
                            cyf_x = cyf;
                        }

                        for (int cy = cyf_x; cy <= cyl; cy++)
                        {
                            const int columnStart = gridj->cxy_ind[cx*gridj->ncy + cy];
                            const int columnEnd   = gridj->cxy_ind[cx*gridj->ncy + cy + 1];

                            d2zxy = d2zx;
                            if (gridj->c0[YY] + cy*gridj->sy > by1)
                            {
                                d2zxy += gmx::square(gridj->c0[YY] + cy*gridj->sy - by1);
                            }
                            else if (gridj->c0[YY] + (cy+1)*gridj->sy < by0)
                            {
                                d2zxy += gmx::square(gridj->c0[YY] + (cy+1)*gridj->sy - by0);
                            }
                            if (columnStart < columnEnd && d2zxy < rlist2)
                            {
                                /* To improve efficiency in the common case
                                 * of a homogeneous particle distribution,
                                 * we estimate the index of the middle cell
                                 * in range (midCell). We search down and up
                                 * starting from this index.
                                 *
                                 * Note that the bbcz_j array contains bounds
                                 * for i-clusters, thus for clusters of 4 atoms.
                                 * For the common case where the j-cluster size
                                 * is 8, we could step with a stride of 2,
                                 * but we do not do this because it would
                                 * complicate this code even more.
                                 */
                                int midCell = columnStart + static_cast<int>(bz1_frac*(columnEnd - columnStart));
                                if (midCell >= columnEnd)
                                {
                                    midCell = columnEnd - 1;
                                }

                                d2xy = d2zxy - d2z;

                                /* Find the lowest cell that can possibly
                                 * be within range.
                                 * Check if we hit the bottom of the grid,
                                 * if the j-cell is below the i-cell and if so,
                                 * if it is within range.
                                 */
                                int downTestCell = midCell;
                                while (downTestCell >= columnStart &&
                                       (bbcz_j[downTestCell*NNBSBB_D + 1] >= bz0 ||
                                        d2xy + gmx::square(bbcz_j[downTestCell*NNBSBB_D + 1] - bz0) < rlist2))
                                {
                                    downTestCell--;
                                }
                                int firstCell = downTestCell + 1;

                                /* Find the highest cell that can possibly
                                 * be within range.
                                 * Check if we hit the top of the grid,
                                 * if the j-cell is above the i-cell and if so,
                                 * if it is within range.
                                 */
                                int upTestCell = midCell + 1;
                                while (upTestCell < columnEnd &&
                                       (bbcz_j[upTestCell*NNBSBB_D] <= bz1 ||
                                        d2xy + gmx::square(bbcz_j[upTestCell*NNBSBB_D] - bz1) < rlist2))
                                {
                                    upTestCell++;
                                }
                                int lastCell = upTestCell - 1;

#define NBNXN_REFCODE 0
#if NBNXN_REFCODE
                                {
                                    /* Simple reference code, for debugging,
                                     * overrides the more complex code above.
                                     */
                                    firstCell = columnEnd;
                                    lastCell  = -1;
                                    for (int k = columnStart; k < columnEnd; k++)
                                    {
                                        if (d2xy + gmx::square(bbcz_j[k*NNBSBB_D + 1] - bz0) < rlist2 &&
                                            k < firstCell)
                                        {
                                            firstCell = k;
                                        }
                                        if (d2xy + gmx::square(bbcz_j[k*NNBSBB_D] - bz1) < rlist2 &&
                                            k > lastCell)
                                        {
                                            lastCell = k;
                                        }
                                    }
                                }
#endif

                                if (gridi == gridj)
                                {
                                    /* We want each atom/cell pair only once,
                                     * only use cj >= ci.
                                     */
                                    if (!c_pbcShiftBackward || shift == CENTRAL)
                                    {
                                        firstCell = std::max(firstCell, ci);
                                    }
                                }

                                if (firstCell <= lastCell)
                                {
                                    GMX_ASSERT(firstCell >= columnStart && lastCell < columnEnd, "The range should reside within the current grid column");

                                    /* For f buffer flags with simple lists */
                                    ncj_old_j = nbl->ncj;

                                    if (nbl->bSimple)
                                    {
                                        /* We have a maximum of 2 j-clusters
                                         * per i-cluster sized cell.
                                         */
                                        check_cell_list_space_simple(nbl, 2*(lastCell - firstCell + 1));
                                    }
                                    else
                                    {
                                        check_cell_list_space_supersub(nbl, lastCell - firstCell + 1);
                                    }

                                    switch (nb_kernel_type)
                                    {
                                        case nbnxnk4x4_PlainC:
                                            makeClusterListSimple(gridj,
                                                                  nbl, ci, firstCell, lastCell,
                                                                  (gridi == gridj && shift == CENTRAL),
                                                                  nbat->x,
                                                                  rlist2, rbb2,
                                                                  &numDistanceChecks);
                                            break;
#ifdef GMX_NBNXN_SIMD_4XN
                                        case nbnxnk4xN_SIMD_4xN:
                                            makeClusterListSimd4xn(gridj,
                                                                   nbl, ci, firstCell, lastCell,
                                                                   (gridi == gridj && shift == CENTRAL),
                                                                   nbat->x,
                                                                   rlist2, rbb2,
                                                                   &numDistanceChecks);
                                            break;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
                                        case nbnxnk4xN_SIMD_2xNN:
                                            makeClusterListSimd2xnn(gridj,
                                                                    nbl, ci, firstCell, lastCell,
                                                                    (gridi == gridj && shift == CENTRAL),
                                                                    nbat->x,
                                                                    rlist2, rbb2,
                                                                    &numDistanceChecks);
                                            break;
#endif
                                        case nbnxnk8x8x8_PlainC:
                                        case nbnxnk8x8x8_GPU:
                                            for (cj = firstCell; cj <= lastCell; cj++)
                                            {
                                                make_cluster_list_supersub(gridi, gridj,
                                                                           nbl, ci, cj,
                                                                           (gridi == gridj && shift == CENTRAL && ci == cj),
                                                                           nbat->xstride, nbat->x,
                                                                           rlist2, rbb2,
                                                                           &numDistanceChecks);
                                            }
                                            break;
                                    }

                                    if (bFBufferFlag && nbl->ncj > ncj_old_j)
                                    {
                                        int cbf = nbl->cj[ncj_old_j].cj >> gridj_flag_shift;
                                        int cbl = nbl->cj[nbl->ncj-1].cj >> gridj_flag_shift;
                                        for (int cb = cbf; cb <= cbl; cb++)
                                        {
                                            bitmask_init_bit(&gridj_flag[cb], th);
                                        }
                                    }

                                    nbl->ncjInUse += nbl->ncj - ncj_old_j;
                                }
                            }
                        }
                    }

                    /* Set the exclusions for this ci list */
                    if (nbl->bSimple)
                    {
                        setExclusionsForSimpleIentry(nbs,
                                                     nbl,
                                                     shift == CENTRAL && gridi == gridj,
                                                     na_cj_2log,
                                                     nbl->ci[nbl->nci],
                                                     exclusions);

                        if (nbs->bFEP)
                        {
                            make_fep_list(nbs, nbat, nbl,
                                          shift == CENTRAL && gridi == gridj,
                                          &(nbl->ci[nbl->nci]),
                                          gridi, gridj, nbl_fep);
                        }
                    }
                    else
                    {
                        setExclusionsForGpuIentry(nbs,
                                                  nbl,
                                                  shift == CENTRAL && gridi == gridj,
                                                  nbl->sci[nbl->nsci],
                                                  exclusions);

                        if (nbs->bFEP)
                        {
                            make_fep_list_supersub(nbs, nbat, nbl,
                                                   shift == CENTRAL && gridi == gridj,
                                                   &(nbl->sci[nbl->nsci]),
                                                   shx, shy, shz,
                                                   rl_fep2,
                                                   gridi, gridj, nbl_fep);
                        }
                    }

                    /* Close this ci list */
                    if (nbl->bSimple)
                    {
                        close_ci_entry_simple(nbl);
                    }
                    else
                    {
                        close_ci_entry_supersub(nbl,
                                                nsubpair_max,
                                                progBal, nsubpair_tot_est,
                                                th, nth);
                    }
                }
            }
        }

        if (bFBufferFlag && nbl->ncj > ncj_old_i)
        {
            bitmask_init_bit(&(work->buffer_flags.flag[(gridi->cell0+ci)>>gridi_flag_shift]), th);
        }
    }

    work->ndistc = numDistanceChecks;

    nbs_cycle_stop(&work->cc[enbsCCsearch]);

    GMX_ASSERT(nbl->ncjInUse == nbl->ncj || nbs->bFEP, "Without free-energy all cj pair-list entries should be in use. Note that subsequent code does not make use of the equality, this check is only here to catch bugs");

    if (debug)
    {
        fprintf(debug, "number of distance checks %d\n", numDistanceChecks);

        if (nbl->bSimple)
        {
            print_nblist_statistics_simple(debug, nbl, nbs, rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug, nbl, nbs, rlist);
        }

        if (nbs->bFEP)
        {
            fprintf(debug, "nbl FEP list pairs: %d\n", nbl_fep->nrj);
        }
    }
}

static void reduce_buffer_flags(const nbnxn_search_t        nbs,
                                int                         nsrc,
                                const nbnxn_buffer_flags_t *dest)
{
    for (int s = 0; s < nsrc; s++)
    {
        gmx_bitmask_t * flag = nbs->work[s].buffer_flags.flag;

        for (int b = 0; b < dest->nflag; b++)
        {
            bitmask_union(&(dest->flag[b]), flag[b]);
        }
    }
}

static void print_reduction_cost(const nbnxn_buffer_flags_t *flags, int nout)
{
    int           nelem, nkeep, ncopy, nred, out;
    gmx_bitmask_t mask_0;

    nelem = 0;
    nkeep = 0;
    ncopy = 0;
    nred  = 0;
    bitmask_init_bit(&mask_0, 0);
    for (int b = 0; b < flags->nflag; b++)
    {
        if (bitmask_is_equal(flags->flag[b], mask_0))
        {
            /* Only flag 0 is set, no copy of reduction required */
            nelem++;
            nkeep++;
        }
        else if (!bitmask_is_zero(flags->flag[b]))
        {
            int c = 0;
            for (out = 0; out < nout; out++)
            {
                if (bitmask_is_set(flags->flag[b], out))
                {
                    c++;
                }
            }
            nelem += c;
            if (c == 1)
            {
                ncopy++;
            }
            else
            {
                nred += c;
            }
        }
    }

    fprintf(debug, "nbnxn reduction: #flag %d #list %d elem %4.2f, keep %4.2f copy %4.2f red %4.2f\n",
            flags->nflag, nout,
            nelem/(double)(flags->nflag),
            nkeep/(double)(flags->nflag),
            ncopy/(double)(flags->nflag),
            nred/(double)(flags->nflag));
}

/* Copies the list entries from src to dest when cjStart <= *cjGlobal < cjEnd.
 * *cjGlobal is updated with the cj count in src.
 * When setFlags==true, flag bit t is set in flag for all i and j clusters.
 */
template<bool setFlags>
static void copySelectedListRange(const nbnxn_ci_t * gmx_restrict srcCi,
                                  const nbnxn_pairlist_t * gmx_restrict src,
                                  nbnxn_pairlist_t * gmx_restrict dest,
                                  gmx_bitmask_t *flag,
                                  int iFlagShift, int jFlagShift, int t)
{
    int ncj = srcCi->cj_ind_end - srcCi->cj_ind_start;

    if (dest->nci + 1 >= dest->ci_nalloc)
    {
        nb_realloc_ci(dest, dest->nci + 1);
    }
    check_cell_list_space_simple(dest, ncj);

    dest->ci[dest->nci]              = *srcCi;
    dest->ci[dest->nci].cj_ind_start = dest->ncj;
    dest->ci[dest->nci].cj_ind_end   = dest->ncj + ncj;

    if (setFlags)
    {
        bitmask_init_bit(&flag[srcCi->ci >> iFlagShift], t);
    }

    for (int j = srcCi->cj_ind_start; j < srcCi->cj_ind_end; j++)
    {
        dest->cj[dest->ncj++] = src->cj[j];

        if (setFlags)
        {
            /* NOTE: This is relatively expensive, since this
             * operation is done for all elements in the list,
             * whereas at list generation this is done only
             * once for each flag entry.
             */
            bitmask_init_bit(&flag[src->cj[j].cj >> jFlagShift], t);
        }
    }

    dest->nci++;
}

/* This routine re-balances the pairlists such that all are nearly equally
 * sized. Only whole i-entries are moved between lists. These are moved
 * between the ends of the lists, such that the buffer reduction cost should
 * not change significantly.
 * Note that all original reduction flags are currently kept. This can lead
 * to reduction of parts of the force buffer that could be avoided. But since
 * the original lists are quite balanced, this will only give minor overhead.
 */
static void rebalanceSimpleLists(int                              numLists,
                                 nbnxn_pairlist_t * const * const srcSet,
                                 nbnxn_pairlist_t               **destSet,
                                 nbnxn_search_work_t             *searchWork)
{
    int ncjTotal = 0;
    for (int s = 0; s < numLists; s++)
    {
        ncjTotal += srcSet[s]->ncjInUse;
    }
    int ncjTarget = (ncjTotal + numLists - 1)/numLists;

#pragma omp parallel num_threads(numLists)
    {
        int t       = gmx_omp_get_thread_num();

        int cjStart = ncjTarget* t;
        int cjEnd   = ncjTarget*(t + 1);

        /* The destination pair-list for task/thread t */
        nbnxn_pairlist_t *dest = destSet[t];

        clear_pairlist(dest);
        dest->bSimple = srcSet[0]->bSimple;
        dest->na_ci   = srcSet[0]->na_ci;
        dest->na_cj   = srcSet[0]->na_cj;

        /* Note that the flags in the work struct (still) contain flags
         * for all entries that are present in srcSet->nbl[t].
         */
        gmx_bitmask_t *flag       = searchWork[t].buffer_flags.flag;

        int            iFlagShift = getBufferFlagShift(dest->na_ci);
        int            jFlagShift = getBufferFlagShift(dest->na_cj);

        int            cjGlobal   = 0;
        for (int s = 0; s < numLists && cjGlobal < cjEnd; s++)
        {
            const nbnxn_pairlist_t *src = srcSet[s];

            if (cjGlobal + src->ncjInUse > cjStart)
            {
                for (int i = 0; i < src->nci && cjGlobal < cjEnd; i++)
                {
                    const nbnxn_ci_t *srcCi = &src->ci[i];
                    int               ncj   = srcCi->cj_ind_end - srcCi->cj_ind_start;
                    if (cjGlobal >= cjStart)
                    {
                        /* If the source list is not our own, we need to set
                         * extra flags (the template bool parameter).
                         */
                        if (s != t)
                        {
                            copySelectedListRange
                            <true>
                                (srcCi, src, dest,
                                flag, iFlagShift, jFlagShift, t);
                        }
                        else
                        {
                            copySelectedListRange
                            <false>
                                (srcCi, src,
                                dest, flag, iFlagShift, jFlagShift, t);
                        }
                    }
                    cjGlobal += ncj;
                }
            }
            else
            {
                cjGlobal += src->ncjInUse;
            }
        }

        dest->ncjInUse = dest->ncj;
    }

#ifndef NDEBUG
    int ncjTotalNew = 0;
    for (int s = 0; s < numLists; s++)
    {
        ncjTotalNew += destSet[s]->ncjInUse;
    }
    GMX_RELEASE_ASSERT(ncjTotalNew == ncjTotal, "The total size of the lists before and after rebalancing should match");
#endif
}

/* Returns if the pairlists are so imbalanced that it is worth rebalancing. */
static bool checkRebalanceSimpleLists(const nbnxn_pairlist_set_t *listSet)
{
    int numLists = listSet->nnbl;
    int ncjMax   = 0;
    int ncjTotal = 0;
    for (int s = 0; s < numLists; s++)
    {
        ncjMax    = std::max(ncjMax, listSet->nbl[s]->ncjInUse);
        ncjTotal += listSet->nbl[s]->ncjInUse;
    }
    if (debug)
    {
        fprintf(debug, "Pair-list ncjMax %d ncjTotal %d\n", ncjMax, ncjTotal);
    }
    /* The rebalancing adds 3% extra time to the search. Heuristically we
     * determined that under common conditions the non-bonded kernel balance
     * improvement will outweigh this when the imbalance is more than 3%.
     * But this will, obviously, depend on search vs kernel time and nstlist.
     */
    const real rebalanceTolerance = 1.03;

    return numLists*ncjMax > ncjTotal*rebalanceTolerance;
}

/* Perform a count (linear) sort to sort the smaller lists to the end.
 * This avoids load imbalance on the GPU, as large lists will be
 * scheduled and executed first and the smaller lists later.
 * Load balancing between multi-processors only happens at the end
 * and there smaller lists lead to more effective load balancing.
 * The sorting is done on the cj4 count, not on the actual pair counts.
 * Not only does this make the sort faster, but it also results in
 * better load balancing than using a list sorted on exact load.
 * This function swaps the pointer in the pair list to avoid a copy operation.
 */
static void sort_sci(nbnxn_pairlist_t *nbl)
{
    nbnxn_list_work_t *work;
    int                m, s0, s1;
    nbnxn_sci_t       *sci_sort;

    if (nbl->ncj4 <= nbl->nsci)
    {
        /* nsci = 0 or all sci have size 1, sorting won't change the order */
        return;
    }

    work = nbl->work;

    /* We will distinguish differences up to double the average */
    m = (2*nbl->ncj4)/nbl->nsci;

    if (m + 1 > work->sort_nalloc)
    {
        work->sort_nalloc = over_alloc_large(m + 1);
        srenew(work->sort, work->sort_nalloc);
    }

    if (work->sci_sort_nalloc != nbl->sci_nalloc)
    {
        work->sci_sort_nalloc = nbl->sci_nalloc;
        nbnxn_realloc_void((void **)&work->sci_sort,
                           0,
                           work->sci_sort_nalloc*sizeof(*work->sci_sort),
                           nbl->alloc, nbl->free);
    }

    /* Count the entries of each size */
    for (int i = 0; i <= m; i++)
    {
        work->sort[i] = 0;
    }
    for (int s = 0; s < nbl->nsci; s++)
    {
        int i = std::min(m, nbl->sci[s].cj4_ind_end - nbl->sci[s].cj4_ind_start);
        work->sort[i]++;
    }
    /* Calculate the offset for each count */
    s0            = work->sort[m];
    work->sort[m] = 0;
    for (int i = m - 1; i >= 0; i--)
    {
        s1            = work->sort[i];
        work->sort[i] = work->sort[i + 1] + s0;
        s0            = s1;
    }

    /* Sort entries directly into place */
    sci_sort = work->sci_sort;
    for (int s = 0; s < nbl->nsci; s++)
    {
        int i = std::min(m, nbl->sci[s].cj4_ind_end - nbl->sci[s].cj4_ind_start);
        sci_sort[work->sort[i]++] = nbl->sci[s];
    }

    /* Swap the sci pointers so we use the new, sorted list */
    work->sci_sort = nbl->sci;
    nbl->sci       = sci_sort;
}

/* Make a local or non-local pair-list, depending on iloc */
void nbnxn_make_pairlist(const nbnxn_search_t  nbs,
                         nbnxn_atomdata_t     *nbat,
                         const t_blocka       *excl,
                         real                  rlist,
                         int                   min_ci_balanced,
                         nbnxn_pairlist_set_t *nbl_list,
                         int                   iloc,
                         int                   nb_kernel_type,
                         t_nrnb               *nrnb)
{
    nbnxn_grid_t      *gridi, *gridj;
    int                nzi, zj0, zj1;
    int                nsubpair_target;
    float              nsubpair_tot_est;
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                ci_block;
    gmx_bool           CombineNBLists;
    gmx_bool           progBal;
    int                np_tot, np_noq, np_hlj, nap;

    nnbl            = nbl_list->nnbl;
    nbl             = nbl_list->nbl;
    CombineNBLists  = nbl_list->bCombined;

    if (debug)
    {
        fprintf(debug, "ns making %d nblists\n", nnbl);
    }

    nbat->bUseBufferFlags = (nbat->nout > 1);
    /* We should re-init the flags before making the first list */
    if (nbat->bUseBufferFlags && LOCAL_I(iloc))
    {
        init_buffer_flags(&nbat->buffer_flags, nbat->natoms);
    }

    if (nbl_list->bSimple)
    {
#if GMX_SIMD
        switch (nb_kernel_type)
        {
#ifdef GMX_NBNXN_SIMD_4XN
            case nbnxnk4xN_SIMD_4xN:
                nbs->icell_set_x = icell_set_x_simd_4xn;
                break;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
            case nbnxnk4xN_SIMD_2xNN:
                nbs->icell_set_x = icell_set_x_simd_2xnn;
                break;
#endif
            default:
                nbs->icell_set_x = icell_set_x_simple;
                break;
        }
#else   // GMX_SIMD
        /* MSVC 2013 complains about switch statements without case */
        nbs->icell_set_x = icell_set_x_simple;
#endif  // GMX_SIMD
    }
    else
    {
        nbs->icell_set_x = icell_set_x_supersub;
    }

    if (LOCAL_I(iloc))
    {
        /* Only zone (grid) 0 vs 0 */
        nzi = 1;
        zj0 = 0;
        zj1 = 1;
    }
    else
    {
        nzi = nbs->zones->nizone;
    }

    if (!nbl_list->bSimple && min_ci_balanced > 0)
    {
        get_nsubpair_target(nbs, iloc, rlist, min_ci_balanced,
                            &nsubpair_target, &nsubpair_tot_est);
    }
    else
    {
        nsubpair_target  = 0;
        nsubpair_tot_est = 0;
    }

    /* Clear all pair-lists */
    for (int th = 0; th < nnbl; th++)
    {
        clear_pairlist(nbl[th]);

        if (nbs->bFEP)
        {
            clear_pairlist_fep(nbl_list->nbl_fep[th]);
        }
    }

    for (int zi = 0; zi < nzi; zi++)
    {
        gridi = &nbs->grid[zi];

        if (NONLOCAL_I(iloc))
        {
            zj0 = nbs->zones->izone[zi].j0;
            zj1 = nbs->zones->izone[zi].j1;
            if (zi == 0)
            {
                zj0++;
            }
        }
        for (int zj = zj0; zj < zj1; zj++)
        {
            gridj = &nbs->grid[zj];

            if (debug)
            {
                fprintf(debug, "ns search grid %d vs %d\n", zi, zj);
            }

            nbs_cycle_start(&nbs->cc[enbsCCsearch]);

            ci_block = get_ci_block_size(gridi, nbs->DomDec, nnbl);

            /* With GPU: generate progressively smaller lists for
             * load balancing for local only or non-local with 2 zones.
             */
            progBal = (LOCAL_I(iloc) || nbs->zones->n <= 2);

#pragma omp parallel for num_threads(nnbl) schedule(static)
            for (int th = 0; th < nnbl; th++)
            {
                try
                {
                    /* Re-init the thread-local work flag data before making
                     * the first list (not an elegant conditional).
                     */
                    if (nbat->bUseBufferFlags && ((zi == 0 && zj == 0)))
                    {
                        init_buffer_flags(&nbs->work[th].buffer_flags, nbat->natoms);
                    }

                    if (CombineNBLists && th > 0)
                    {
                        clear_pairlist(nbl[th]);
                    }

                    /* Divide the i super cell equally over the nblists */
                    nbnxn_make_pairlist_part(nbs, gridi, gridj,
                                             &nbs->work[th], nbat, *excl,
                                             rlist,
                                             nb_kernel_type,
                                             ci_block,
                                             nbat->bUseBufferFlags,
                                             nsubpair_target,
                                             progBal, nsubpair_tot_est,
                                             th, nnbl,
                                             nbl[th],
                                             nbl_list->nbl_fep[th]);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
            nbs_cycle_stop(&nbs->cc[enbsCCsearch]);

            np_tot = 0;
            np_noq = 0;
            np_hlj = 0;
            for (int th = 0; th < nnbl; th++)
            {
                inc_nrnb(nrnb, eNR_NBNXN_DIST2, nbs->work[th].ndistc);

                if (nbl_list->bSimple)
                {
                    np_tot += nbl[th]->ncj;
                    np_noq += nbl[th]->work->ncj_noq;
                    np_hlj += nbl[th]->work->ncj_hlj;
                }
                else
                {
                    /* This count ignores potential subsequent pair pruning */
                    np_tot += nbl[th]->nci_tot;
                }
            }
            nap                   = nbl[0]->na_ci*nbl[0]->na_cj;
            nbl_list->natpair_ljq = (np_tot - np_noq)*nap - np_hlj*nap/2;
            nbl_list->natpair_lj  = np_noq*nap;
            nbl_list->natpair_q   = np_hlj*nap/2;

            if (CombineNBLists && nnbl > 1)
            {
                nbs_cycle_start(&nbs->cc[enbsCCcombine]);

                combine_nblists(nnbl-1, nbl+1, nbl[0]);

                nbs_cycle_stop(&nbs->cc[enbsCCcombine]);
            }
        }
    }

    if (nbl_list->bSimple)
    {
        if (nnbl > 1 && checkRebalanceSimpleLists(nbl_list))
        {
            rebalanceSimpleLists(nbl_list->nnbl, nbl_list->nbl, nbl_list->nbl_work, nbs->work);

            /* Swap the pointer of the sets of pair lists */
            nbnxn_pairlist_t **tmp = nbl_list->nbl;
            nbl_list->nbl          = nbl_list->nbl_work;
            nbl_list->nbl_work     = tmp;
        }
    }
    else
    {
        /* Sort the entries on size, large ones first */
        if (CombineNBLists || nnbl == 1)
        {
            sort_sci(nbl[0]);
        }
        else
        {
#pragma omp parallel for num_threads(nnbl) schedule(static)
            for (int th = 0; th < nnbl; th++)
            {
                try
                {
                    sort_sci(nbl[th]);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
        }
    }

    if (nbat->bUseBufferFlags)
    {
        reduce_buffer_flags(nbs, nbl_list->nnbl, &nbat->buffer_flags);
    }

    if (nbs->bFEP)
    {
        /* Balance the free-energy lists over all the threads */
        balance_fep_lists(nbs, nbl_list);
    }

    /* This is a fresh list, so not pruned, stored using ci and nci.
     * ciOuter and nciOuter are invalid at this point.
     */
    GMX_ASSERT(nbl_list->nbl[0]->nciOuter == -1, "nciOuter should have been set to -1 to signal that it is invalid");

    /* Special performance logging stuff (env.var. GMX_NBNXN_CYCLE) */
    if (LOCAL_I(iloc))
    {
        nbs->search_count++;
    }
    if (nbs->print_cycles &&
        (!nbs->DomDec || !LOCAL_I(iloc)) &&
        nbs->search_count % 100 == 0)
    {
        nbs_cycle_print(stderr, nbs);
    }

    /* If we have more than one list, they either got rebalancing (CPU)
     * or combined (GPU), so we should dump the final result to debug.
     */
    if (debug && nbl_list->nnbl > 1)
    {
        if (nbl_list->bSimple)
        {
            for (int t = 0; t < nbl_list->nnbl; t++)
            {
                print_nblist_statistics_simple(debug, nbl_list->nbl[t], nbs, rlist);
            }
        }
        else
        {
            print_nblist_statistics_supersub(debug, nbl_list->nbl[0], nbs, rlist);
        }
    }

    if (debug)
    {
        if (gmx_debug_at)
        {
            if (nbl_list->bSimple)
            {
                for (int t = 0; t < nbl_list->nnbl; t++)
                {
                    print_nblist_ci_cj(debug, nbl_list->nbl[t]);
                }
            }
            else
            {
                print_nblist_sci_cj(debug, nbl_list->nbl[0]);
            }
        }

        if (nbat->bUseBufferFlags)
        {
            print_reduction_cost(&nbat->buffer_flags, nbl_list->nnbl);
        }
    }
}

void nbnxnPrepareListForDynamicPruning(nbnxn_pairlist_set_t *listSet)
{
    /* TODO: Restructure the lists so we have actual outer and inner
     *       list objects so we can set a single pointer instead of
     *       swapping several pointers.
     */

    for (int i = 0; i < listSet->nnbl; i++)
    {
        /* The search produced a list in ci/cj.
         * Swap the list pointers so we get the outer list is ciOuter,cjOuter
         * and we can prune that to get an inner list in ci/cj.
         */
        nbnxn_pairlist_t *list = listSet->nbl[i];
        list->nciOuter         = list->nci;

        nbnxn_ci_t *ciTmp      = list->ciOuter;
        list->ciOuter          = list->ci;
        list->ci               = ciTmp;

        nbnxn_cj_t *cjTmp      = list->cjOuter;
        list->cjOuter          = list->cj;
        list->cj               = cjTmp;

        /* Signal that this inner list is currently invalid */
        list->nci              = -1;
    }
}
