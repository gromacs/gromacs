/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "pairlist.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

#include "atomdata.h"
#include "boundingboxes.h"
#include "clusterdistancekerneltype.h"
#include "gridset.h"
#include "nbnxm_geometry.h"
#include "nbnxm_simd.h"
#include "pairlistset.h"
#include "pairlistsets.h"
#include "pairlistwork.h"
#include "pairsearch.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

using BoundingBox   = Nbnxm::BoundingBox;   // TODO: Remove when refactoring this file
using BoundingBox1D = Nbnxm::BoundingBox1D; // TODO: Remove when refactoring this file

using Grid = Nbnxm::Grid; // TODO: Remove when refactoring this file

// Convience alias for partial Nbnxn namespace usage
using InteractionLocality = gmx::InteractionLocality;

/* We shift the i-particles backward for PBC.
 * This leads to more conditionals than shifting forward.
 * We do this to get more balanced pair lists.
 */
constexpr bool c_pbcShiftBackward = true;

/* Layout for the nonbonded NxN pair lists */
enum class NbnxnLayout
{
    NoSimd4x4, // i-cluster size 4, j-cluster size 4
    Simd4xN,   // i-cluster size 4, j-cluster size SIMD width
    Simd2xNN,  // i-cluster size 4, j-cluster size half SIMD width
    Gpu8x8x8   // i-cluster size 8, j-cluster size 8 + super-clustering
};

#if GMX_SIMD
/* Returns the j-cluster size */
template<NbnxnLayout layout>
static constexpr int jClusterSize()
{
    static_assert(layout == NbnxnLayout::NoSimd4x4 || layout == NbnxnLayout::Simd4xN
                          || layout == NbnxnLayout::Simd2xNN,
                  "Currently jClusterSize only supports CPU layouts");

    return layout == NbnxnLayout::Simd4xN
                   ? GMX_SIMD_REAL_WIDTH
                   : (layout == NbnxnLayout::Simd2xNN ? GMX_SIMD_REAL_WIDTH / 2 : c_nbnxnCpuIClusterSize);
}

/*! \brief Returns the j-cluster index given the i-cluster index.
 *
 * \tparam    jClusterSize      The number of atoms in a j-cluster
 * \tparam    jSubClusterIndex  The j-sub-cluster index (0/1), used when size(j-cluster) <
 * size(i-cluster) \param[in] ci                The i-cluster index
 */
template<int jClusterSize, int jSubClusterIndex>
static inline int cjFromCi(int ci)
{
    static_assert(jClusterSize == c_nbnxnCpuIClusterSize / 2 || jClusterSize == c_nbnxnCpuIClusterSize
                          || jClusterSize == c_nbnxnCpuIClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    static_assert(jSubClusterIndex == 0 || jSubClusterIndex == 1,
                  "Only sub-cluster indices 0 and 1 are supported");

    if (jClusterSize == c_nbnxnCpuIClusterSize / 2)
    {
        if (jSubClusterIndex == 0)
        {
            return ci << 1;
        }
        else
        {
            return ((ci + 1) << 1) - 1;
        }
    }
    else if (jClusterSize == c_nbnxnCpuIClusterSize)
    {
        return ci;
    }
    else
    {
        return ci >> 1;
    }
}

/*! \brief Returns the j-cluster index given the i-cluster index.
 *
 * \tparam    layout            The pair-list layout
 * \tparam    jSubClusterIndex  The j-sub-cluster index (0/1), used when size(j-cluster) <
 * size(i-cluster) \param[in] ci                The i-cluster index
 */
template<NbnxnLayout layout, int jSubClusterIndex>
static inline int cjFromCi(int ci)
{
    constexpr int clusterSize = jClusterSize<layout>();

    return cjFromCi<clusterSize, jSubClusterIndex>(ci);
}

/* Returns the nbnxn coordinate data index given the i-cluster index */
template<NbnxnLayout layout>
static inline int xIndexFromCi(int ci)
{
    constexpr int clusterSize = jClusterSize<layout>();

    static_assert(clusterSize == c_nbnxnCpuIClusterSize / 2 || clusterSize == c_nbnxnCpuIClusterSize
                          || clusterSize == c_nbnxnCpuIClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    if (clusterSize <= c_nbnxnCpuIClusterSize)
    {
        /* Coordinates are stored packed in groups of 4 */
        return ci * STRIDE_P4;
    }
    else
    {
        /* Coordinates packed in 8, i-cluster size is half the packing width */
        return (ci >> 1) * STRIDE_P8 + (ci & 1) * (c_packX8 >> 1);
    }
}

/* Returns the nbnxn coordinate data index given the j-cluster index */
template<NbnxnLayout layout>
static inline int xIndexFromCj(int cj)
{
    constexpr int clusterSize = jClusterSize<layout>();

    static_assert(clusterSize == c_nbnxnCpuIClusterSize / 2 || clusterSize == c_nbnxnCpuIClusterSize
                          || clusterSize == c_nbnxnCpuIClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    if (clusterSize == c_nbnxnCpuIClusterSize / 2)
    {
        /* Coordinates are stored packed in groups of 4 */
        return (cj >> 1) * STRIDE_P4 + (cj & 1) * (c_packX4 >> 1);
    }
    else if (clusterSize == c_nbnxnCpuIClusterSize)
    {
        /* Coordinates are stored packed in groups of 4 */
        return cj * STRIDE_P4;
    }
    else
    {
        /* Coordinates are stored packed in groups of 8 */
        return cj * STRIDE_P8;
    }
}
#endif // GMX_SIMD


void nbnxn_init_pairlist_fep(t_nblist* nl)
{
    nl->type      = GMX_NBLIST_INTERACTION_FREE_ENERGY;
    nl->igeometry = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
    /* The interaction functions are set in the free energy kernel fuction */
    nl->ivdw     = -1;
    nl->ivdwmod  = -1;
    nl->ielec    = -1;
    nl->ielecmod = -1;

    nl->maxnri   = 0;
    nl->maxnrj   = 0;
    nl->nri      = 0;
    nl->nrj      = 0;
    nl->iinr     = nullptr;
    nl->gid      = nullptr;
    nl->shift    = nullptr;
    nl->jindex   = nullptr;
    nl->jjnr     = nullptr;
    nl->excl_fep = nullptr;
}

static void init_buffer_flags(nbnxn_buffer_flags_t* flags, int natoms)
{
    flags->nflag = (natoms + NBNXN_BUFFERFLAG_SIZE - 1) / NBNXN_BUFFERFLAG_SIZE;
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

/* Returns the pair-list cutoff between a bounding box and a grid cell given an atom-to-atom pair-list cutoff
 *
 * Given a cutoff distance between atoms, this functions returns the cutoff
 * distance2 between a bounding box of a group of atoms and a grid cell.
 * Since atoms can be geometrically outside of the cell they have been
 * assigned to (when atom groups instead of individual atoms are assigned
 * to cells), this distance returned can be larger than the input.
 */
static real listRangeForBoundingBoxToGridCell(real rlist, const Grid::Dimensions& gridDims)
{
    return rlist + gridDims.maxAtomGroupRadius;
}
/* Returns the pair-list cutoff between a grid cells given an atom-to-atom pair-list cutoff
 *
 * Given a cutoff distance between atoms, this functions returns the cutoff
 * distance2 between two grid cells.
 * Since atoms can be geometrically outside of the cell they have been
 * assigned to (when atom groups instead of individual atoms are assigned
 * to cells), this distance returned can be larger than the input.
 */
static real listRangeForGridCellToGridCell(real                    rlist,
                                           const Grid::Dimensions& iGridDims,
                                           const Grid::Dimensions& jGridDims)
{
    return rlist + iGridDims.maxAtomGroupRadius + jGridDims.maxAtomGroupRadius;
}

/* Determines the cell range along one dimension that
 * the bounding box b0 - b1 sees.
 */
template<int dim>
static void
get_cell_range(real b0, real b1, const Grid::Dimensions& jGridDims, real d2, real rlist, int* cf, int* cl)
{
    real listRangeBBToCell2 = gmx::square(listRangeForBoundingBoxToGridCell(rlist, jGridDims));
    real distanceInCells    = (b0 - jGridDims.lowerCorner[dim]) * jGridDims.invCellSize[dim];
    *cf                     = std::max(static_cast<int>(distanceInCells), 0);

    while (*cf > 0
           && d2 + gmx::square((b0 - jGridDims.lowerCorner[dim]) - (*cf - 1 + 1) * jGridDims.cellSize[dim])
                      < listRangeBBToCell2)
    {
        (*cf)--;
    }

    *cl = std::min(static_cast<int>((b1 - jGridDims.lowerCorner[dim]) * jGridDims.invCellSize[dim]),
                   jGridDims.numCells[dim] - 1);
    while (*cl < jGridDims.numCells[dim] - 1
           && d2 + gmx::square((*cl + 1) * jGridDims.cellSize[dim] - (b1 - jGridDims.lowerCorner[dim]))
                      < listRangeBBToCell2)
    {
        (*cl)++;
    }
}

/* Reference code calculating the distance^2 between two bounding boxes */
/*
   static float box_dist2(float bx0, float bx1, float by0,
                       float by1, float bz0, float bz1,
                       const BoundingBox *bb)
   {
    float d2;
    float dl, dh, dm, dm0;

    d2 = 0;

    dl  = bx0 - bb->upper.x;
    dh  = bb->lower.x - bx1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = by0 - bb->upper.y;
    dh  = bb->lower.y - by1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    dl  = bz0 - bb->upper.z;
    dh  = bb->lower.z - bz1;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0f);
    d2 += dm0*dm0;

    return d2;
   }
 */

#if !NBNXN_SEARCH_BB_SIMD4

/*! \brief Plain C code calculating the distance^2 between two bounding boxes in xyz0 format
 *
 * \param[in] bb_i  First bounding box
 * \param[in] bb_j  Second bounding box
 */
static float clusterBoundingBoxDistance2(const BoundingBox& bb_i, const BoundingBox& bb_j)
{
    float dl  = bb_i.lower.x - bb_j.upper.x;
    float dh  = bb_j.lower.x - bb_i.upper.x;
    float dm  = std::max(dl, dh);
    float dm0 = std::max(dm, 0.0F);
    float d2  = dm0 * dm0;

    dl  = bb_i.lower.y - bb_j.upper.y;
    dh  = bb_j.lower.y - bb_i.upper.y;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0F);
    d2 += dm0 * dm0;

    dl  = bb_i.lower.z - bb_j.upper.z;
    dh  = bb_j.lower.z - bb_i.upper.z;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0F);
    d2 += dm0 * dm0;

    return d2;
}

#else /* NBNXN_SEARCH_BB_SIMD4 */

/*! \brief 4-wide SIMD code calculating the distance^2 between two bounding boxes in xyz0 format
 *
 * \param[in] bb_i  First bounding box, should be aligned for 4-wide SIMD
 * \param[in] bb_j  Second bounding box, should be aligned for 4-wide SIMD
 */
static float clusterBoundingBoxDistance2(const BoundingBox& bb_i, const BoundingBox& bb_j)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    const Simd4Float bb_i_S0 = load4(bb_i.lower.ptr());
    const Simd4Float bb_i_S1 = load4(bb_i.upper.ptr());
    const Simd4Float bb_j_S0 = load4(bb_j.lower.ptr());
    const Simd4Float bb_j_S1 = load4(bb_j.upper.ptr());

    const Simd4Float dl_S = bb_i_S0 - bb_j_S1;
    const Simd4Float dh_S = bb_j_S0 - bb_i_S1;

    const Simd4Float dm_S  = max(dl_S, dh_S);
    const Simd4Float dm0_S = max(dm_S, simd4SetZeroF());

    return dotProduct(dm0_S, dm0_S);
}

/* Calculate bb bounding distances of bb_i[si,...,si+3] and store them in d2 */
template<int boundingBoxStart>
static inline void gmx_simdcall clusterBoundingBoxDistance2_xxxx_simd4_inner(const float*     bb_i,
                                                                             float*           d2,
                                                                             const Simd4Float xj_l,
                                                                             const Simd4Float yj_l,
                                                                             const Simd4Float zj_l,
                                                                             const Simd4Float xj_h,
                                                                             const Simd4Float yj_h,
                                                                             const Simd4Float zj_h)
{
    constexpr int stride = c_packedBoundingBoxesDimSize;

    const int shi = boundingBoxStart * Nbnxm::c_numBoundingBoxBounds1D * DIM;

    const Simd4Float zero = setZero();

    const Simd4Float xi_l = load4(bb_i + shi + 0 * stride);
    const Simd4Float yi_l = load4(bb_i + shi + 1 * stride);
    const Simd4Float zi_l = load4(bb_i + shi + 2 * stride);
    const Simd4Float xi_h = load4(bb_i + shi + 3 * stride);
    const Simd4Float yi_h = load4(bb_i + shi + 4 * stride);
    const Simd4Float zi_h = load4(bb_i + shi + 5 * stride);

    const Simd4Float dx_0 = xi_l - xj_h;
    const Simd4Float dy_0 = yi_l - yj_h;
    const Simd4Float dz_0 = zi_l - zj_h;

    const Simd4Float dx_1 = xj_l - xi_h;
    const Simd4Float dy_1 = yj_l - yi_h;
    const Simd4Float dz_1 = zj_l - zi_h;

    const Simd4Float mx = max(dx_0, dx_1);
    const Simd4Float my = max(dy_0, dy_1);
    const Simd4Float mz = max(dz_0, dz_1);

    const Simd4Float m0x = max(mx, zero);
    const Simd4Float m0y = max(my, zero);
    const Simd4Float m0z = max(mz, zero);

    const Simd4Float d2x = m0x * m0x;
    const Simd4Float d2y = m0y * m0y;
    const Simd4Float d2z = m0z * m0z;

    const Simd4Float d2s = d2x + d2y;
    const Simd4Float d2t = d2s + d2z;

    store4(d2 + boundingBoxStart, d2t);
}

/* 4-wide SIMD code for nsi bb distances for bb format xxxxyyyyzzzz */
static void clusterBoundingBoxDistance2_xxxx_simd4(const float* bb_j, const int nsi, const float* bb_i, float* d2)
{
    constexpr int stride = c_packedBoundingBoxesDimSize;

    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    const Simd4Float xj_l = Simd4Float(bb_j[0 * stride]);
    const Simd4Float yj_l = Simd4Float(bb_j[1 * stride]);
    const Simd4Float zj_l = Simd4Float(bb_j[2 * stride]);
    const Simd4Float xj_h = Simd4Float(bb_j[3 * stride]);
    const Simd4Float yj_h = Simd4Float(bb_j[4 * stride]);
    const Simd4Float zj_h = Simd4Float(bb_j[5 * stride]);

    /* Here we "loop" over si (0,stride) from 0 to nsi with step stride.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    clusterBoundingBoxDistance2_xxxx_simd4_inner<0>(bb_i, d2, xj_l, yj_l, zj_l, xj_h, yj_h, zj_h);
    if (stride < nsi)
    {
        clusterBoundingBoxDistance2_xxxx_simd4_inner<stride>(bb_i, d2, xj_l, yj_l, zj_l, xj_h, yj_h, zj_h);
    }
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */


/* Returns if any atom pair from two clusters is within distance sqrt(rlist2) */
static inline gmx_bool
clusterpair_in_range(const NbnxnPairlistGpuWork& work, int si, int csj, int stride, const real* x_j, real rlist2)
{
#if !GMX_SIMD4_HAVE_REAL

    /* Plain C version.
     * All coordinates are stored as xyzxyz...
     */

    const real* x_i = work.iSuperClusterData.x.data();

    for (int i = 0; i < c_nbnxnGpuClusterSize; i++)
    {
        int i0 = (si * c_nbnxnGpuClusterSize + i) * DIM;
        for (int j = 0; j < c_nbnxnGpuClusterSize; j++)
        {
            int j0 = (csj * c_nbnxnGpuClusterSize + j) * stride;

            real d2 = gmx::square(x_i[i0] - x_j[j0]) + gmx::square(x_i[i0 + 1] - x_j[j0 + 1])
                      + gmx::square(x_i[i0 + 2] - x_j[j0 + 2]);

            if (d2 < rlist2)
            {
                return TRUE;
            }
        }
    }

    return FALSE;

#else /* !GMX_SIMD4_HAVE_REAL */

    /* 4-wide SIMD version.
     * The coordinates x_i are stored as xxxxyyyy..., x_j is stored xyzxyz...
     * Using 8-wide AVX(2) is not faster on Intel Sandy Bridge and Haswell.
     */
    static_assert(c_nbnxnGpuClusterSize == 8 || c_nbnxnGpuClusterSize == 4,
                  "A cluster is hard-coded to 4/8 atoms.");

    Simd4Real rc2_S = Simd4Real(rlist2);

    const real* x_i = work.iSuperClusterData.xSimd.data();

    int       dim_stride = c_nbnxnGpuClusterSize * DIM;
    Simd4Real ix_S0      = load4(x_i + si * dim_stride + 0 * GMX_SIMD4_WIDTH);
    Simd4Real iy_S0      = load4(x_i + si * dim_stride + 1 * GMX_SIMD4_WIDTH);
    Simd4Real iz_S0      = load4(x_i + si * dim_stride + 2 * GMX_SIMD4_WIDTH);

    Simd4Real ix_S1, iy_S1, iz_S1;
    if (c_nbnxnGpuClusterSize == 8)
    {
        ix_S1 = load4(x_i + si * dim_stride + 3 * GMX_SIMD4_WIDTH);
        iy_S1 = load4(x_i + si * dim_stride + 4 * GMX_SIMD4_WIDTH);
        iz_S1 = load4(x_i + si * dim_stride + 5 * GMX_SIMD4_WIDTH);
    }
    /* We loop from the outer to the inner particles to maximize
     * the chance that we find a pair in range quickly and return.
     */
    int j0 = csj * c_nbnxnGpuClusterSize;
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

        jx0_S = Simd4Real(x_j[j0 * stride + 0]);
        jy0_S = Simd4Real(x_j[j0 * stride + 1]);
        jz0_S = Simd4Real(x_j[j0 * stride + 2]);

        jx1_S = Simd4Real(x_j[j1 * stride + 0]);
        jy1_S = Simd4Real(x_j[j1 * stride + 1]);
        jz1_S = Simd4Real(x_j[j1 * stride + 2]);

        /* Calculate distance */
        dx_S0 = ix_S0 - jx0_S;
        dy_S0 = iy_S0 - jy0_S;
        dz_S0 = iz_S0 - jz0_S;
        dx_S2 = ix_S0 - jx1_S;
        dy_S2 = iy_S0 - jy1_S;
        dz_S2 = iz_S0 - jz1_S;
        if (c_nbnxnGpuClusterSize == 8)
        {
            dx_S1 = ix_S1 - jx0_S;
            dy_S1 = iy_S1 - jy0_S;
            dz_S1 = iz_S1 - jz0_S;
            dx_S3 = ix_S1 - jx1_S;
            dy_S3 = iy_S1 - jy1_S;
            dz_S3 = iz_S1 - jz1_S;
        }

        /* rsq = dx*dx+dy*dy+dz*dz */
        rsq_S0 = norm2(dx_S0, dy_S0, dz_S0);
        rsq_S2 = norm2(dx_S2, dy_S2, dz_S2);
        if (c_nbnxnGpuClusterSize == 8)
        {
            rsq_S1 = norm2(dx_S1, dy_S1, dz_S1);
            rsq_S3 = norm2(dx_S3, dy_S3, dz_S3);
        }

        wco_S0 = (rsq_S0 < rc2_S);
        wco_S2 = (rsq_S2 < rc2_S);
        if (c_nbnxnGpuClusterSize == 8)
        {
            wco_S1 = (rsq_S1 < rc2_S);
            wco_S3 = (rsq_S3 < rc2_S);
        }
        if (c_nbnxnGpuClusterSize == 8)
        {
            wco_any_S01 = wco_S0 || wco_S1;
            wco_any_S23 = wco_S2 || wco_S3;
            wco_any_S   = wco_any_S01 || wco_any_S23;
        }
        else
        {
            wco_any_S = wco_S0 || wco_S2;
        }

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
static inline int nblCj(gmx::ArrayRef<const nbnxn_cj_t> cjList, int cjIndex)
{
    return cjList[cjIndex].cj;
}

/* Returns the j-cluster index for index cjIndex in a cj4 list */
static inline int nblCj(gmx::ArrayRef<const nbnxn_cj4_t> cj4List, int cjIndex)
{
    return cj4List[cjIndex / c_nbnxnGpuJgroupSize].cj[cjIndex & (c_nbnxnGpuJgroupSize - 1)];
}

/* Returns the i-interaction mask of the j sub-cell for index cj_ind */
static unsigned int nbl_imask0(const NbnxnPairlistGpu* nbl, int cj_ind)
{
    return nbl->cj4[cj_ind / c_nbnxnGpuJgroupSize].imei[0].imask;
}

NbnxnPairlistCpu::NbnxnPairlistCpu() :
    na_ci(c_nbnxnCpuIClusterSize),
    na_cj(0),
    rlist(0),
    ncjInUse(0),
    nci_tot(0),
    work(std::make_unique<NbnxnPairlistCpuWork>())
{
}

NbnxnPairlistGpu::NbnxnPairlistGpu(gmx::PinningPolicy pinningPolicy) :
    na_ci(c_nbnxnGpuClusterSize),
    na_cj(c_nbnxnGpuClusterSize),
    na_sc(c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize),
    rlist(0),
    sci({}, { pinningPolicy }),
    cj4({}, { pinningPolicy }),
    excl({}, { pinningPolicy }),
    nci_tot(0),
    work(std::make_unique<NbnxnPairlistGpuWork>())
{
    static_assert(c_nbnxnGpuNumClusterPerSupercluster == c_gpuNumClusterPerCell,
                  "The search code assumes that the a super-cluster matches a search grid cell");

    static_assert(sizeof(cj4[0].imei[0].imask) * 8 >= c_nbnxnGpuJgroupSize * c_gpuNumClusterPerCell,
                  "The i super-cluster cluster interaction mask does not contain a sufficient "
                  "number of bits");

    static_assert(sizeof(excl[0]) * 8 >= c_nbnxnGpuJgroupSize * c_gpuNumClusterPerCell,
                  "The GPU exclusion mask does not contain a sufficient number of bits");

    // We always want a first entry without any exclusions
    excl.resize(1);
}

// TODO: Move to pairlistset.cpp
PairlistSet::PairlistSet(const InteractionLocality locality, const PairlistParams& pairlistParams) :
    locality_(locality),
    params_(pairlistParams)
{
    isCpuType_ = (params_.pairlistType == PairlistType::Simple4x2
                  || params_.pairlistType == PairlistType::Simple4x4
                  || params_.pairlistType == PairlistType::Simple4x8);
    // Currently GPU lists are always combined
    combineLists_ = !isCpuType_;

    const int numLists = gmx_omp_nthreads_get(emntNonbonded);

    if (!combineLists_ && numLists > NBNXN_BUFFERFLAG_MAX_THREADS)
    {
        gmx_fatal(FARGS,
                  "%d OpenMP threads were requested. Since the non-bonded force buffer reduction "
                  "is prohibitively slow with more than %d threads, we do not allow this. Use %d "
                  "or less OpenMP threads.",
                  numLists, NBNXN_BUFFERFLAG_MAX_THREADS, NBNXN_BUFFERFLAG_MAX_THREADS);
    }

    if (isCpuType_)
    {
        cpuLists_.resize(numLists);
        if (numLists > 1)
        {
            cpuListsWork_.resize(numLists);
        }
    }
    else
    {
        /* Only list 0 is used on the GPU, use normal allocation for i>0 */
        gpuLists_.emplace_back(gmx::PinningPolicy::PinnedIfSupported);
        /* Lists 0 to numLists are use for constructing lists in parallel
         * on the CPU using numLists threads (and then merged into list 0).
         */
        for (int i = 1; i < numLists; i++)
        {
            gpuLists_.emplace_back(gmx::PinningPolicy::CannotBePinned);
        }
    }
    if (params_.haveFep)
    {
        fepLists_.resize(numLists);

        /* Execute in order to avoid memory interleaving between threads */
#pragma omp parallel for num_threads(numLists) schedule(static)
        for (int i = 0; i < numLists; i++)
        {
            try
            {
                /* We used to allocate all normal lists locally on each thread
                 * as well. The question is if allocating the object on the
                 * master thread (but all contained list memory thread local)
                 * impacts performance.
                 */
                fepLists_[i] = std::make_unique<t_nblist>();
                nbnxn_init_pairlist_fep(fepLists_[i].get());
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
    }
}

/* Print statistics of a pair list, used for debug output */
static void print_nblist_statistics(FILE*                   fp,
                                    const NbnxnPairlistCpu& nbl,
                                    const Nbnxm::GridSet&   gridSet,
                                    const real              rl)
{
    const Grid&             grid = gridSet.grids()[0];
    const Grid::Dimensions& dims = grid.dimensions();

    fprintf(fp, "nbl nci %zu ncj %d\n", nbl.ci.size(), nbl.ncjInUse);
    const int numAtomsJCluster = grid.geometry().numAtomsJCluster;
    const double numAtomsPerCell = nbl.ncjInUse / static_cast<double>(grid.numCells()) * numAtomsJCluster;
    fprintf(fp, "nbl na_cj %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n", nbl.na_cj, rl,
            nbl.ncjInUse, nbl.ncjInUse / static_cast<double>(grid.numCells()), numAtomsPerCell,
            numAtomsPerCell
                    / (0.5 * 4.0 / 3.0 * M_PI * rl * rl * rl * grid.numCells() * numAtomsJCluster
                       / (dims.gridSize[XX] * dims.gridSize[YY] * dims.gridSize[ZZ])));

    fprintf(fp, "nbl average j cell list length %.1f\n",
            0.25 * nbl.ncjInUse / std::max(static_cast<double>(nbl.ci.size()), 1.0));

    int cs[SHIFTS] = { 0 };
    int npexcl     = 0;
    for (const nbnxn_ci_t& ciEntry : nbl.ci)
    {
        cs[ciEntry.shift & NBNXN_CI_SHIFT] += ciEntry.cj_ind_end - ciEntry.cj_ind_start;

        int j = ciEntry.cj_ind_start;
        while (j < ciEntry.cj_ind_end && nbl.cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            npexcl++;
            j++;
        }
    }
    fprintf(fp, "nbl cell pairs, total: %zu excl: %d %.1f%%\n", nbl.cj.size(), npexcl,
            100 * npexcl / std::max(static_cast<double>(nbl.cj.size()), 1.0));
    for (int s = 0; s < SHIFTS; s++)
    {
        if (cs[s] > 0)
        {
            fprintf(fp, "nbl shift %2d ncj %3d\n", s, cs[s]);
        }
    }
}

/* Print statistics of a pair lists, used for debug output */
static void print_nblist_statistics(FILE*                   fp,
                                    const NbnxnPairlistGpu& nbl,
                                    const Nbnxm::GridSet&   gridSet,
                                    const real              rl)
{
    const Grid&             grid = gridSet.grids()[0];
    const Grid::Dimensions& dims = grid.dimensions();

    fprintf(fp, "nbl nsci %zu ncj4 %zu nsi %d excl4 %zu\n", nbl.sci.size(), nbl.cj4.size(),
            nbl.nci_tot, nbl.excl.size());
    const int numAtomsCluster = grid.geometry().numAtomsICluster;
    const double numAtomsPerCell = nbl.nci_tot / static_cast<double>(grid.numClusters()) * numAtomsCluster;
    fprintf(fp, "nbl na_c %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n", nbl.na_ci, rl,
            nbl.nci_tot, nbl.nci_tot / static_cast<double>(grid.numClusters()), numAtomsPerCell,
            numAtomsPerCell
                    / (0.5 * 4.0 / 3.0 * M_PI * rl * rl * rl * grid.numClusters() * numAtomsCluster
                       / (dims.gridSize[XX] * dims.gridSize[YY] * dims.gridSize[ZZ])));

    double sum_nsp                       = 0;
    double sum_nsp2                      = 0;
    int    nsp_max                       = 0;
    int    c[c_gpuNumClusterPerCell + 1] = { 0 };
    for (const nbnxn_sci_t& sci : nbl.sci)
    {
        int nsp = 0;
        for (int j4 = sci.cj4_ind_start; j4 < sci.cj4_ind_end; j4++)
        {
            for (int j = 0; j < c_nbnxnGpuJgroupSize; j++)
            {
                int b = 0;
                for (int si = 0; si < c_gpuNumClusterPerCell; si++)
                {
                    if (nbl.cj4[j4].imei[0].imask & (1U << (j * c_gpuNumClusterPerCell + si)))
                    {
                        b++;
                    }
                }
                nsp += b;
                c[b]++;
            }
        }
        sum_nsp += nsp;
        sum_nsp2 += nsp * nsp;
        nsp_max = std::max(nsp_max, nsp);
    }
    if (!nbl.sci.empty())
    {
        sum_nsp /= nbl.sci.size();
        sum_nsp2 /= nbl.sci.size();
    }
    fprintf(fp, "nbl #cluster-pairs: av %.1f stddev %.1f max %d\n", sum_nsp,
            std::sqrt(sum_nsp2 - sum_nsp * sum_nsp), nsp_max);

    if (!nbl.cj4.empty())
    {
        for (int b = 0; b <= c_gpuNumClusterPerCell; b++)
        {
            fprintf(fp, "nbl j-list #i-subcell %d %7d %4.1f\n", b, c[b],
                    100.0 * c[b] / size_t{ nbl.cj4.size() * c_nbnxnGpuJgroupSize });
        }
    }
}

/* Returns a reference to the exclusion mask for j-cluster-group \p cj4 and warp \p warp
 * Generates a new exclusion entry when the j-cluster-group uses
 * the default all-interaction mask at call time, so the returned mask
 * can be modified when needed.
 */
static nbnxn_excl_t& get_exclusion_mask(NbnxnPairlistGpu* nbl, int cj4, int warp)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* No exclusions set, make a new list entry */
        const size_t oldSize = nbl->excl.size();
        GMX_ASSERT(oldSize >= 1, "We should always have entry [0]");
        /* Add entry with default values: no exclusions */
        nbl->excl.resize(oldSize + 1);
        nbl->cj4[cj4].imei[warp].excl_ind = oldSize;
    }

    return nbl->excl[nbl->cj4[cj4].imei[warp].excl_ind];
}

/* Sets self exclusions and excludes half of the double pairs in the self cluster-pair \p nbl->cj4[cj4Index].cj[jOffsetInGroup]
 *
 * \param[in,out] nbl             The cluster pair list
 * \param[in]     cj4Index        The j-cluster group index into \p nbl->cj4
 * \param[in]     jOffsetInGroup  The j-entry offset in \p nbl->cj4[cj4Index]
 * \param[in]     iClusterInCell  The i-cluster index in the cell
 */
static void setSelfAndNewtonExclusionsGpu(NbnxnPairlistGpu* nbl,
                                          const int         cj4Index,
                                          const int         jOffsetInGroup,
                                          const int         iClusterInCell)
{
    constexpr int numJatomsPerPart = c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit;

    /* The exclusions are stored separately for each part of the split */
    for (int part = 0; part < c_nbnxnGpuClusterpairSplit; part++)
    {
        const int jOffset = part * numJatomsPerPart;
        /* Make a new exclusion mask entry for each part, if we don't already have one yet */
        nbnxn_excl_t& excl = get_exclusion_mask(nbl, cj4Index, part);

        /* Set all bits with j-index <= i-index */
        for (int jIndexInPart = 0; jIndexInPart < numJatomsPerPart; jIndexInPart++)
        {
            for (int i = jOffset + jIndexInPart; i < c_nbnxnGpuClusterSize; i++)
            {
                excl.pair[jIndexInPart * c_nbnxnGpuClusterSize + i] &=
                        ~(1U << (jOffsetInGroup * c_gpuNumClusterPerCell + iClusterInCell));
            }
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
    return (rdiag && ci * 2 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_0
                                  : (rdiag && ci * 2 + 1 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_1
                                                               : NBNXN_INTERACTION_MASK_ALL));
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=4 */
gmx_unused static unsigned int get_imask_simd_j4(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=8 */
gmx_unused static unsigned int get_imask_simd_j8(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj * 2 ? NBNXN_INTERACTION_MASK_DIAG_J8_0
                                  : (rdiag && ci == cj * 2 + 1 ? NBNXN_INTERACTION_MASK_DIAG_J8_1
                                                               : NBNXN_INTERACTION_MASK_ALL));
}

#if GMX_SIMD
#    if GMX_SIMD_REAL_WIDTH == 2
#        define get_imask_simd_4xn get_imask_simd_j2
#    endif
#    if GMX_SIMD_REAL_WIDTH == 4
#        define get_imask_simd_4xn get_imask_simd_j4
#    endif
#    if GMX_SIMD_REAL_WIDTH == 8
#        define get_imask_simd_4xn get_imask_simd_j8
#        define get_imask_simd_2xnn get_imask_simd_j4
#    endif
#    if GMX_SIMD_REAL_WIDTH == 16
#        define get_imask_simd_2xnn get_imask_simd_j8
#    endif
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
static void makeClusterListSimple(const Grid&       jGrid,
                                  NbnxnPairlistCpu* nbl,
                                  int               icluster,
                                  int               jclusterFirst,
                                  int               jclusterLast,
                                  bool              excludeSubDiagonal,
                                  const real* gmx_restrict x_j,
                                  real                     rlist2,
                                  float                    rbb2,
                                  int* gmx_restrict numDistanceChecks)
{
    const BoundingBox* gmx_restrict bb_ci = nbl->work->iClusterData.bb.data();
    const real* gmx_restrict x_ci         = nbl->work->iClusterData.x.data();

    gmx_bool InRange;

    InRange = FALSE;
    while (!InRange && jclusterFirst <= jclusterLast)
    {
        real d2 = clusterBoundingBoxDistance2(bb_ci[0], jGrid.jBoundingBoxes()[jclusterFirst]);
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
            int cjf_gl = jGrid.cellOffset() + jclusterFirst;
            for (int i = 0; i < c_nbnxnCpuIClusterSize && !InRange; i++)
            {
                for (int j = 0; j < c_nbnxnCpuIClusterSize; j++)
                {
                    InRange =
                            InRange
                            || (gmx::square(x_ci[i * STRIDE_XYZ + XX]
                                            - x_j[(cjf_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + XX])
                                        + gmx::square(x_ci[i * STRIDE_XYZ + YY]
                                                      - x_j[(cjf_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + YY])
                                        + gmx::square(x_ci[i * STRIDE_XYZ + ZZ]
                                                      - x_j[(cjf_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + ZZ])
                                < rlist2);
                }
            }
            *numDistanceChecks += c_nbnxnCpuIClusterSize * c_nbnxnCpuIClusterSize;
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
        real d2 = clusterBoundingBoxDistance2(bb_ci[0], jGrid.jBoundingBoxes()[jclusterLast]);
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
            int cjl_gl = jGrid.cellOffset() + jclusterLast;
            for (int i = 0; i < c_nbnxnCpuIClusterSize && !InRange; i++)
            {
                for (int j = 0; j < c_nbnxnCpuIClusterSize; j++)
                {
                    InRange =
                            InRange
                            || (gmx::square(x_ci[i * STRIDE_XYZ + XX]
                                            - x_j[(cjl_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + XX])
                                        + gmx::square(x_ci[i * STRIDE_XYZ + YY]
                                                      - x_j[(cjl_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + YY])
                                        + gmx::square(x_ci[i * STRIDE_XYZ + ZZ]
                                                      - x_j[(cjl_gl * c_nbnxnCpuIClusterSize + j) * STRIDE_XYZ + ZZ])
                                < rlist2);
                }
            }
            *numDistanceChecks += c_nbnxnCpuIClusterSize * c_nbnxnCpuIClusterSize;
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
            nbnxn_cj_t cjEntry;
            cjEntry.cj   = jGrid.cellOffset() + jcluster;
            cjEntry.excl = get_imask(excludeSubDiagonal, icluster, jcluster);
            nbl->cj.push_back(cjEntry);
        }
        /* Increase the closing index in the i list */
        nbl->ci.back().cj_ind_end = nbl->cj.size();
    }
}

#ifdef GMX_NBNXN_SIMD_4XN
#    include "pairlist_simd_4xm.h"
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
#    include "pairlist_simd_2xmm.h"
#endif

/* Plain C or SIMD4 code for making a pair list of super-cell sci vs scj.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_supersub(const Grid&       iGrid,
                                       const Grid&       jGrid,
                                       NbnxnPairlistGpu* nbl,
                                       const int         sci,
                                       const int         scj,
                                       const bool        excludeSubDiagonal,
                                       const int         stride,
                                       const real*       x,
                                       const real        rlist2,
                                       const float       rbb2,
                                       int*              numDistanceChecks)
{
    NbnxnPairlistGpuWork& work = *nbl->work;

#if NBNXN_BBXXXX
    const float* pbb_ci = work.iSuperClusterData.bbPacked.data();
#else
    const BoundingBox* bb_ci = work.iSuperClusterData.bb.data();
#endif

    assert(c_nbnxnGpuClusterSize == iGrid.geometry().numAtomsICluster);
    assert(c_nbnxnGpuClusterSize == jGrid.geometry().numAtomsICluster);

    /* We generate the pairlist mainly based on bounding-box distances
     * and do atom pair distance based pruning on the GPU.
     * Only if a j-group contains a single cluster-pair, we try to prune
     * that pair based on atom distances on the CPU to avoid empty j-groups.
     */
#define PRUNE_LIST_CPU_ONE 1
#define PRUNE_LIST_CPU_ALL 0

#if PRUNE_LIST_CPU_ONE
    int ci_last = -1;
#endif

    float* d2l = work.distanceBuffer.data();

    for (int subc = 0; subc < jGrid.numClustersPerCell()[scj]; subc++)
    {
        const int cj4_ind   = work.cj_ind / c_nbnxnGpuJgroupSize;
        const int cj_offset = work.cj_ind - cj4_ind * c_nbnxnGpuJgroupSize;
        const int cj        = scj * c_gpuNumClusterPerCell + subc;

        const int cj_gl = jGrid.cellOffset() * c_gpuNumClusterPerCell + cj;

        int ci1;
        if (excludeSubDiagonal && sci == scj)
        {
            ci1 = subc + 1;
        }
        else
        {
            ci1 = iGrid.numClustersPerCell()[sci];
        }

#if NBNXN_BBXXXX
        /* Determine all ci1 bb distances in one call with SIMD4 */
        const int offset = packedBoundingBoxesIndex(cj) + (cj & (c_packedBoundingBoxesDimSize - 1));
        clusterBoundingBoxDistance2_xxxx_simd4(jGrid.packedBoundingBoxes().data() + offset, ci1,
                                               pbb_ci, d2l);
        *numDistanceChecks += c_nbnxnGpuClusterSize * 2;
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
            d2l[ci] = clusterBoundingBoxDistance2(bb_ci[ci], jGrid.jBoundingBoxes()[cj]);
            *numDistanceChecks += 2;
#endif
            float d2 = d2l[ci];

#if PRUNE_LIST_CPU_ALL
            /* Check if the distance is within the distance where
             * we use only the bounding box distance rbb,
             * or within the cut-off and there is at least one atom pair
             * within the cut-off. This check is very costly.
             */
            *numDistanceChecks += c_nbnxnGpuClusterSize * c_nbnxnGpuClusterSize;
            if (d2 < rbb2 || (d2 < rlist2 && clusterpair_in_range(work, ci, cj_gl, stride, x, rlist2)))
#else
            /* Check if the distance between the two bounding boxes
             * in within the pair-list cut-off.
             */
            if (d2 < rlist2)
#endif
            {
                /* Flag this i-subcell to be taken into account */
                imask |= (1U << (cj_offset * c_gpuNumClusterPerCell + ci));

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
        if (npair == 1 && d2l[ci_last] >= rbb2
            && !clusterpair_in_range(work, ci_last, cj_gl, stride, x, rlist2))
        {
            imask &= ~(1U << (cj_offset * c_gpuNumClusterPerCell + ci_last));
            npair--;
        }
#endif

        if (npair > 0)
        {
            /* We have at least one cluster pair: add a j-entry */
            if (static_cast<size_t>(cj4_ind) == nbl->cj4.size())
            {
                nbl->cj4.resize(nbl->cj4.size() + 1);
            }
            nbnxn_cj4_t* cj4 = &nbl->cj4[cj4_ind];

            cj4->cj[cj_offset] = cj_gl;

            /* Set the exclusions for the ci==sj entry.
             * Here we don't bother to check if this entry is actually flagged,
             * as it will nearly always be in the list.
             */
            if (excludeSubDiagonal && sci == scj)
            {
                setSelfAndNewtonExclusionsGpu(nbl, cj4_ind, cj_offset, subc);
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
            nbl->sci.back().cj4_ind_end =
                    (nbl->work->cj_ind + c_nbnxnGpuJgroupSize - 1) / c_nbnxnGpuJgroupSize;
        }
    }
}

/* Returns how many contiguous j-clusters we have starting in the i-list */
template<typename CjListType>
static int numContiguousJClusters(const int                       cjIndexStart,
                                  const int                       cjIndexEnd,
                                  gmx::ArrayRef<const CjListType> cjList)
{
    const int firstJCluster = nblCj(cjList, cjIndexStart);

    int numContiguous = 0;

    while (cjIndexStart + numContiguous < cjIndexEnd
           && nblCj(cjList, cjIndexStart + numContiguous) == firstJCluster + numContiguous)
    {
        numContiguous++;
    }

    return numContiguous;
}

/*! \internal
 * \brief Helper struct for efficient searching for excluded atoms in a j-list
 */
struct JListRanges
{
    /*! \brief Constructs a j-list range from \p cjList with the given index range */
    template<typename CjListType>
    JListRanges(int cjIndexStart, int cjIndexEnd, gmx::ArrayRef<const CjListType> cjList);

    int cjIndexStart; //!< The start index in the j-list
    int cjIndexEnd;   //!< The end index in the j-list
    int cjFirst;      //!< The j-cluster with index cjIndexStart
    int cjLast;       //!< The j-cluster with index cjIndexEnd-1
    int numDirect; //!< Up to cjIndexStart+numDirect the j-clusters are cjFirst + the index offset
};

#ifndef DOXYGEN
template<typename CjListType>
JListRanges::JListRanges(int cjIndexStart, int cjIndexEnd, gmx::ArrayRef<const CjListType> cjList) :
    cjIndexStart(cjIndexStart),
    cjIndexEnd(cjIndexEnd)
{
    GMX_ASSERT(cjIndexEnd > cjIndexStart, "JListRanges should only be called with non-empty lists");

    cjFirst = nblCj(cjList, cjIndexStart);
    cjLast  = nblCj(cjList, cjIndexEnd - 1);

    /* Determine how many contiguous j-cells we have starting
     * from the first i-cell. This number can be used to directly
     * calculate j-cell indices for excluded atoms.
     */
    numDirect = numContiguousJClusters(cjIndexStart, cjIndexEnd, cjList);
}
#endif // !DOXYGEN

/* Return the index of \p jCluster in the given range or -1 when not present
 *
 * Note: This code is executed very often and therefore performance is
 *       important. It should be inlined and fully optimized.
 */
template<typename CjListType>
static inline int findJClusterInJList(int                             jCluster,
                                      const JListRanges&              ranges,
                                      gmx::ArrayRef<const CjListType> cjList)
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
        index          = -1;
        int rangeStart = ranges.cjIndexStart + ranges.numDirect;
        int rangeEnd   = ranges.cjIndexEnd;
        int rangeMiddle;
        while (index == -1 && rangeStart < rangeEnd)
        {
            rangeMiddle = (rangeStart + rangeEnd) >> 1;

            const int clusterMiddle = nblCj(cjList, rangeMiddle);

            if (jCluster == clusterMiddle)
            {
                index = rangeMiddle;
            }
            else if (jCluster < clusterMiddle)
            {
                rangeEnd = rangeMiddle;
            }
            else
            {
                rangeStart = rangeMiddle + 1;
            }
        }
    }

    return index;
}

// TODO: Get rid of the two functions below by renaming sci to ci (or something better)

/* Return the i-entry in the list we are currently operating on */
static nbnxn_ci_t* getOpenIEntry(NbnxnPairlistCpu* nbl)
{
    return &nbl->ci.back();
}

/* Return the i-entry in the list we are currently operating on */
static nbnxn_sci_t* getOpenIEntry(NbnxnPairlistGpu* nbl)
{
    return &nbl->sci.back();
}

/* Set all atom-pair exclusions for a simple type list i-entry
 *
 * Set all atom-pair exclusions from the topology stored in exclusions
 * as masks in the pair-list for simple list entry iEntry.
 */
static void setExclusionsForIEntry(const Nbnxm::GridSet& gridSet,
                                   NbnxnPairlistCpu*     nbl,
                                   gmx_bool              diagRemoved,
                                   int                   na_cj_2log,
                                   const nbnxn_ci_t&     iEntry,
                                   const t_blocka&       exclusions)
{
    if (iEntry.cj_ind_end == iEntry.cj_ind_start)
    {
        /* Empty list: no exclusions */
        return;
    }

    const JListRanges ranges(iEntry.cj_ind_start, iEntry.cj_ind_end, gmx::makeConstArrayRef(nbl->cj));

    const int iCluster = iEntry.ci;

    gmx::ArrayRef<const int> cell        = gridSet.cells();
    gmx::ArrayRef<const int> atomIndices = gridSet.atomIndices();

    /* Loop over the atoms in the i-cluster */
    for (int i = 0; i < nbl->na_ci; i++)
    {
        const int iIndex = iCluster * nbl->na_ci + i;
        const int iAtom  = atomIndices[iIndex];
        if (iAtom >= 0)
        {
            /* Loop over the topology-based exclusions for this i-atom */
            for (int exclIndex = exclusions.index[iAtom]; exclIndex < exclusions.index[iAtom + 1];
                 exclIndex++)
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
                            findJClusterInJList(jCluster, ranges, gmx::makeConstArrayRef(nbl->cj));

                    if (index >= 0)
                    {
                        /* We found an exclusion, clear the corresponding
                         * interaction bit.
                         */
                        const int innerJ = jIndex - (jCluster << na_cj_2log);

                        nbl->cj[index].excl &= ~(1U << ((i << na_cj_2log) + innerJ));
                    }
                }
            }
        }
    }
}

/* Add a new i-entry to the FEP list and copy the i-properties */
static inline void fep_list_new_nri_copy(t_nblist* nlist)
{
    /* Add a new i-entry */
    nlist->nri++;

    assert(nlist->nri < nlist->maxnri);

    /* Duplicate the last i-entry, except for jindex, which continues */
    nlist->iinr[nlist->nri]   = nlist->iinr[nlist->nri - 1];
    nlist->shift[nlist->nri]  = nlist->shift[nlist->nri - 1];
    nlist->gid[nlist->nri]    = nlist->gid[nlist->nri - 1];
    nlist->jindex[nlist->nri] = nlist->nrj;
}

/* Rellocate FEP list for size nl->maxnri, TODO: replace by C++ */
static void reallocate_nblist(t_nblist* nl)
{
    if (gmx_debug_at)
    {
        fprintf(debug,
                "reallocating neigborlist (ielec=%d, ivdw=%d, igeometry=%d, type=%d), maxnri=%d\n",
                nl->ielec, nl->ivdw, nl->igeometry, nl->type, nl->maxnri);
    }
    srenew(nl->iinr, nl->maxnri);
    srenew(nl->gid, nl->maxnri);
    srenew(nl->shift, nl->maxnri);
    srenew(nl->jindex, nl->maxnri + 1);
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
static void make_fep_list(gmx::ArrayRef<const int> atomIndices,
                          const nbnxn_atomdata_t*  nbat,
                          NbnxnPairlistCpu*        nbl,
                          gmx_bool                 bDiagRemoved,
                          nbnxn_ci_t*              nbl_ci,
                          real gmx_unused shx,
                          real gmx_unused shy,
                          real gmx_unused shz,
                          real gmx_unused rlist_fep2,
                          const Grid&     iGrid,
                          const Grid&     jGrid,
                          t_nblist*       nlist)
{
    int      ci, cj_ind_start, cj_ind_end, cja, cjr;
    int      nri_max;
    int      gid_i = 0, gid_j, gid;
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
    nri_max = nbl->na_ci * nbl->na_cj * (cj_ind_end - cj_ind_start);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    const int numAtomsJCluster = jGrid.geometry().numAtomsJCluster;

    const nbnxn_atomdata_t::Params& nbatParams = nbat->params();

    const int ngid = nbatParams.nenergrp;

    /* TODO: Consider adding a check in grompp and changing this to an assert */
    const int numBitsInEnergyGroupIdsForAtomsInJCluster = sizeof(gid_cj) * 8;
    if (ngid * numAtomsJCluster > numBitsInEnergyGroupIdsForAtomsInJCluster)
    {
        gmx_fatal(FARGS,
                  "The Verlet scheme with %dx%d kernels and free-energy only supports up to %zu "
                  "energy groups",
                  iGrid.geometry().numAtomsICluster, numAtomsJCluster,
                  (sizeof(gid_cj) * 8) / numAtomsJCluster);
    }

    egp_shift = nbatParams.neg_2log;
    egp_mask  = (1 << egp_shift) - 1;

    /* Loop over the atoms in the i sub-cell */
    bFEP_i_all = TRUE;
    for (int i = 0; i < nbl->na_ci; i++)
    {
        ind_i = ci * nbl->na_ci + i;
        ai    = atomIndices[ind_i];
        if (ai >= 0)
        {
            nri                    = nlist->nri;
            nlist->jindex[nri + 1] = nlist->jindex[nri];
            nlist->iinr[nri]       = ai;
            /* The actual energy group pair index is set later */
            nlist->gid[nri]   = 0;
            nlist->shift[nri] = nbl_ci->shift & NBNXN_CI_SHIFT;

            bFEP_i = iGrid.atomIsPerturbed(ci - iGrid.cellOffset(), i);

            bFEP_i_all = bFEP_i_all && bFEP_i;

            if (nlist->nrj + (cj_ind_end - cj_ind_start) * nbl->na_cj > nlist->maxnrj)
            {
                nlist->maxnrj = over_alloc_small(nlist->nrj + (cj_ind_end - cj_ind_start) * nbl->na_cj);
                srenew(nlist->jjnr, nlist->maxnrj);
                srenew(nlist->excl_fep, nlist->maxnrj);
            }

            if (ngid > 1)
            {
                gid_i = (nbatParams.energrp[ci] >> (egp_shift * i)) & egp_mask;
            }

            for (int cj_ind = cj_ind_start; cj_ind < cj_ind_end; cj_ind++)
            {
                unsigned int fep_cj;

                cja = nbl->cj[cj_ind].cj;

                if (numAtomsJCluster == jGrid.geometry().numAtomsICluster)
                {
                    cjr    = cja - jGrid.cellOffset();
                    fep_cj = jGrid.fepBits(cjr);
                    if (ngid > 1)
                    {
                        gid_cj = nbatParams.energrp[cja];
                    }
                }
                else if (2 * numAtomsJCluster == jGrid.geometry().numAtomsICluster)
                {
                    cjr = cja - jGrid.cellOffset() * 2;
                    /* Extract half of the ci fep/energrp mask */
                    fep_cj = (jGrid.fepBits(cjr >> 1) >> ((cjr & 1) * numAtomsJCluster))
                             & ((1 << numAtomsJCluster) - 1);
                    if (ngid > 1)
                    {
                        gid_cj = nbatParams.energrp[cja >> 1] >> ((cja & 1) * numAtomsJCluster * egp_shift)
                                 & ((1 << (numAtomsJCluster * egp_shift)) - 1);
                    }
                }
                else
                {
                    cjr = cja - (jGrid.cellOffset() >> 1);
                    /* Combine two ci fep masks/energrp */
                    fep_cj = jGrid.fepBits(cjr * 2)
                             + (jGrid.fepBits(cjr * 2 + 1) << jGrid.geometry().numAtomsICluster);
                    if (ngid > 1)
                    {
                        gid_cj = nbatParams.energrp[cja * 2]
                                 + (nbatParams.energrp[cja * 2 + 1]
                                    << (jGrid.geometry().numAtomsICluster * egp_shift));
                    }
                }

                if (bFEP_i || fep_cj != 0)
                {
                    for (int j = 0; j < nbl->na_cj; j++)
                    {
                        /* Is this interaction perturbed and not excluded? */
                        ind_j = cja * nbl->na_cj + j;
                        aj    = atomIndices[ind_j];
                        if (aj >= 0 && (bFEP_i || (fep_cj & (1 << j))) && (!bDiagRemoved || ind_j >= ind_i))
                        {
                            if (ngid > 1)
                            {
                                gid_j = (gid_cj >> (j * egp_shift)) & egp_mask;
                                gid   = GID(gid_i, gid_j, ngid);

                                if (nlist->nrj > nlist->jindex[nri] && nlist->gid[nri] != gid)
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
                            nlist->jjnr[nlist->nrj] = aj;
                            nlist->excl_fep[nlist->nrj] =
                                    (nbl->cj[cj_ind].excl >> (i * nbl->na_cj + j)) & 1;
                            nlist->nrj++;

                            /* Exclude it from the normal list.
                             * Note that the charge has been set to zero,
                             * but we need to avoid 0/0, as perturbed atoms
                             * can be on top of each other.
                             */
                            nbl->cj[cj_ind].excl &= ~(1U << (i * nbl->na_cj + j));
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
        nbl->ncjInUse -= cj_ind_end - cj_ind_start;
    }
}

/* Return the index of atom a within a cluster */
static inline int cj_mod_cj4(int cj)
{
    return cj & (c_nbnxnGpuJgroupSize - 1);
}

/* Convert a j-cluster to a cj4 group */
static inline int cj_to_cj4(int cj)
{
    return cj / c_nbnxnGpuJgroupSize;
}

/* Return the index of an j-atom within a warp */
static inline int a_mod_wj(int a)
{
    return a & (c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit - 1);
}

/* As make_fep_list above, but for super/sub lists. */
static void make_fep_list(gmx::ArrayRef<const int> atomIndices,
                          const nbnxn_atomdata_t*  nbat,
                          NbnxnPairlistGpu*        nbl,
                          gmx_bool                 bDiagRemoved,
                          const nbnxn_sci_t*       nbl_sci,
                          real                     shx,
                          real                     shy,
                          real                     shz,
                          real                     rlist_fep2,
                          const Grid&              iGrid,
                          const Grid&              jGrid,
                          t_nblist*                nlist)
{
    int                nri_max;
    int                c_abs;
    int                ind_i, ind_j, ai, aj;
    int                nri;
    gmx_bool           bFEP_i;
    real               xi, yi, zi;
    const nbnxn_cj4_t* cj4;

    const int numJClusterGroups = nbl_sci->numJClusterGroups();
    if (numJClusterGroups == 0)
    {
        /* Empty list */
        return;
    }

    const int sci = nbl_sci->sci;

    const int cj4_ind_start = nbl_sci->cj4_ind_start;
    const int cj4_ind_end   = nbl_sci->cj4_ind_end;

    /* Here we process one super-cell, max #atoms na_sc, versus a list
     * cj4 entries, each with max c_nbnxnGpuJgroupSize cj's, each
     * of size na_cj atoms.
     * On the GPU we don't support energy groups (yet).
     * So for each of the na_sc i-atoms, we need max one FEP list
     * for each max_nrj_fep j-atoms.
     */
    nri_max = nbl->na_sc * nbl->na_cj * (1 + (numJClusterGroups * c_nbnxnGpuJgroupSize) / max_nrj_fep);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    /* Loop over the atoms in the i super-cluster */
    for (int c = 0; c < c_gpuNumClusterPerCell; c++)
    {
        c_abs = sci * c_gpuNumClusterPerCell + c;

        for (int i = 0; i < nbl->na_ci; i++)
        {
            ind_i = c_abs * nbl->na_ci + i;
            ai    = atomIndices[ind_i];
            if (ai >= 0)
            {
                nri                    = nlist->nri;
                nlist->jindex[nri + 1] = nlist->jindex[nri];
                nlist->iinr[nri]       = ai;
                /* With GPUs, energy groups are not supported */
                nlist->gid[nri]   = 0;
                nlist->shift[nri] = nbl_sci->shift & NBNXN_CI_SHIFT;

                bFEP_i = iGrid.atomIsPerturbed(c_abs - iGrid.cellOffset() * c_gpuNumClusterPerCell, i);

                xi = nbat->x()[ind_i * nbat->xstride + XX] + shx;
                yi = nbat->x()[ind_i * nbat->xstride + YY] + shy;
                zi = nbat->x()[ind_i * nbat->xstride + ZZ] + shz;

                const int nrjMax = nlist->nrj + numJClusterGroups * c_nbnxnGpuJgroupSize * nbl->na_cj;
                if (nrjMax > nlist->maxnrj)
                {
                    nlist->maxnrj = over_alloc_small(nrjMax);
                    srenew(nlist->jjnr, nlist->maxnrj);
                    srenew(nlist->excl_fep, nlist->maxnrj);
                }

                for (int cj4_ind = cj4_ind_start; cj4_ind < cj4_ind_end; cj4_ind++)
                {
                    cj4 = &nbl->cj4[cj4_ind];

                    for (int gcj = 0; gcj < c_nbnxnGpuJgroupSize; gcj++)
                    {
                        if ((cj4->imei[0].imask & (1U << (gcj * c_gpuNumClusterPerCell + c))) == 0)
                        {
                            /* Skip this ci for this cj */
                            continue;
                        }

                        const int cjr = cj4->cj[gcj] - jGrid.cellOffset() * c_gpuNumClusterPerCell;

                        if (bFEP_i || jGrid.clusterIsPerturbed(cjr))
                        {
                            for (int j = 0; j < nbl->na_cj; j++)
                            {
                                /* Is this interaction perturbed and not excluded? */
                                ind_j = (jGrid.cellOffset() * c_gpuNumClusterPerCell + cjr) * nbl->na_cj + j;
                                aj = atomIndices[ind_j];
                                if (aj >= 0 && (bFEP_i || jGrid.atomIsPerturbed(cjr, j))
                                    && (!bDiagRemoved || ind_j >= ind_i))
                                {
                                    int          excl_pair;
                                    unsigned int excl_bit;
                                    real         dx, dy, dz;

                                    const int jHalf =
                                            j / (c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit);
                                    nbnxn_excl_t& excl = get_exclusion_mask(nbl, cj4_ind, jHalf);

                                    excl_pair = a_mod_wj(j) * nbl->na_ci + i;
                                    excl_bit  = (1U << (gcj * c_gpuNumClusterPerCell + c));

                                    dx = nbat->x()[ind_j * nbat->xstride + XX] - xi;
                                    dy = nbat->x()[ind_j * nbat->xstride + YY] - yi;
                                    dz = nbat->x()[ind_j * nbat->xstride + ZZ] - zi;

                                    /* The unpruned GPU list has more than 2/3
                                     * of the atom pairs beyond rlist. Using
                                     * this list will cause a lot of overhead
                                     * in the CPU FEP kernels, especially
                                     * relative to the fast GPU kernels.
                                     * So we prune the FEP list here.
                                     */
                                    if (dx * dx + dy * dy + dz * dz < rlist_fep2)
                                    {
                                        if (nlist->nrj - nlist->jindex[nri] >= max_nrj_fep)
                                        {
                                            fep_list_new_nri_copy(nlist);
                                            nri = nlist->nri;
                                        }

                                        /* Add it to the FEP list */
                                        nlist->jjnr[nlist->nrj] = aj;
                                        nlist->excl_fep[nlist->nrj] =
                                                (excl.pair[excl_pair] & excl_bit) ? 1 : 0;
                                        nlist->nrj++;
                                    }

                                    /* Exclude it from the normal list.
                                     * Note that the charge and LJ parameters have
                                     * been set to zero, but we need to avoid 0/0,
                                     * as perturbed atoms can be on top of each other.
                                     */
                                    excl.pair[excl_pair] &= ~excl_bit;
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
static void setExclusionsForIEntry(const Nbnxm::GridSet& gridSet,
                                   NbnxnPairlistGpu*     nbl,
                                   gmx_bool              diagRemoved,
                                   int gmx_unused     na_cj_2log,
                                   const nbnxn_sci_t& iEntry,
                                   const t_blocka&    exclusions)
{
    if (iEntry.numJClusterGroups() == 0)
    {
        /* Empty list */
        return;
    }

    /* Set the search ranges using start and end j-cluster indices.
     * Note that here we can not use cj4_ind_end, since the last cj4
     * can be only partially filled, so we use cj_ind.
     */
    const JListRanges ranges(iEntry.cj4_ind_start * c_nbnxnGpuJgroupSize, nbl->work->cj_ind,
                             gmx::makeConstArrayRef(nbl->cj4));

    GMX_ASSERT(nbl->na_ci == c_nbnxnGpuClusterSize, "na_ci should match the GPU cluster size");
    constexpr int c_clusterSize      = c_nbnxnGpuClusterSize;
    constexpr int c_superClusterSize = c_nbnxnGpuNumClusterPerSupercluster * c_nbnxnGpuClusterSize;

    const int iSuperCluster = iEntry.sci;

    gmx::ArrayRef<const int> atomIndices = gridSet.atomIndices();
    gmx::ArrayRef<const int> cell        = gridSet.cells();

    /* Loop over the atoms in the i super-cluster */
    for (int i = 0; i < c_superClusterSize; i++)
    {
        const int iIndex = iSuperCluster * c_superClusterSize + i;
        const int iAtom  = atomIndices[iIndex];
        if (iAtom >= 0)
        {
            const int iCluster = i / c_clusterSize;

            /* Loop over the topology-based exclusions for this i-atom */
            for (int exclIndex = exclusions.index[iAtom]; exclIndex < exclusions.index[iAtom + 1];
                 exclIndex++)
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
                if (diagRemoved && jIndex <= iSuperCluster * nbl->na_sc + i)
                {
                    continue;
                }

                const int jCluster = jIndex / c_clusterSize;

                /* Check whether the cluster is in our list? */
                if (jCluster >= ranges.cjFirst && jCluster <= ranges.cjLast)
                {
                    const int index =
                            findJClusterInJList(jCluster, ranges, gmx::makeConstArrayRef(nbl->cj4));

                    if (index >= 0)
                    {
                        /* We found an exclusion, clear the corresponding
                         * interaction bit.
                         */
                        const unsigned int pairMask =
                                (1U << (cj_mod_cj4(index) * c_gpuNumClusterPerCell + iCluster));
                        /* Check if the i-cluster interacts with the j-cluster */
                        if (nbl_imask0(nbl, index) & pairMask)
                        {
                            const int innerI = (i & (c_clusterSize - 1));
                            const int innerJ = (jIndex & (c_clusterSize - 1));

                            /* Determine which j-half (CUDA warp) we are in */
                            const int jHalf = innerJ / (c_clusterSize / c_nbnxnGpuClusterpairSplit);

                            nbnxn_excl_t& interactionMask =
                                    get_exclusion_mask(nbl, cj_to_cj4(index), jHalf);

                            interactionMask.pair[a_mod_wj(innerJ) * c_clusterSize + innerI] &= ~pairMask;
                        }
                    }
                }
            }
        }
    }
}

/* Make a new ci entry at the back of nbl->ci */
static void addNewIEntry(NbnxnPairlistCpu* nbl, int ci, int shift, int flags)
{
    nbnxn_ci_t ciEntry;
    ciEntry.ci    = ci;
    ciEntry.shift = shift;
    /* Store the interaction flags along with the shift */
    ciEntry.shift |= flags;
    ciEntry.cj_ind_start = nbl->cj.size();
    ciEntry.cj_ind_end   = nbl->cj.size();
    nbl->ci.push_back(ciEntry);
}

/* Make a new sci entry at index nbl->nsci */
static void addNewIEntry(NbnxnPairlistGpu* nbl, int sci, int shift, int gmx_unused flags)
{
    nbnxn_sci_t sciEntry;
    sciEntry.sci           = sci;
    sciEntry.shift         = shift;
    sciEntry.cj4_ind_start = nbl->cj4.size();
    sciEntry.cj4_ind_end   = nbl->cj4.size();

    nbl->sci.push_back(sciEntry);
}

/* Sort the simple j-list cj on exclusions.
 * Entries with exclusions will all be sorted to the beginning of the list.
 */
static void sort_cj_excl(nbnxn_cj_t* cj, int ncj, NbnxnPairlistCpuWork* work)
{
    work->cj.resize(ncj);

    /* Make a list of the j-cells involving exclusions */
    int jnew = 0;
    for (int j = 0; j < ncj; j++)
    {
        if (cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) || (jnew == 1 && cj[0].excl != NBNXN_INTERACTION_MASK_ALL)))
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
static void closeIEntry(NbnxnPairlistCpu* nbl,
                        int gmx_unused sp_max_av,
                        gmx_bool gmx_unused progBal,
                        float gmx_unused nsp_tot_est,
                        int gmx_unused thread,
                        int gmx_unused nthread)
{
    nbnxn_ci_t& ciEntry = nbl->ci.back();

    /* All content of the new ci entry have already been filled correctly,
     * we only need to sort and increase counts or remove the entry when empty.
     */
    const int jlen = ciEntry.cj_ind_end - ciEntry.cj_ind_start;
    if (jlen > 0)
    {
        sort_cj_excl(nbl->cj.data() + ciEntry.cj_ind_start, jlen, nbl->work.get());

        /* The counts below are used for non-bonded pair/flop counts
         * and should therefore match the available kernel setups.
         */
        if (!(ciEntry.shift & NBNXN_CI_DO_COUL(0)))
        {
            nbl->work->ncj_noq += jlen;
        }
        else if ((ciEntry.shift & NBNXN_CI_HALF_LJ(0)) || !(ciEntry.shift & NBNXN_CI_DO_LJ(0)))
        {
            nbl->work->ncj_hlj += jlen;
        }
    }
    else
    {
        /* Entry is empty: remove it  */
        nbl->ci.pop_back();
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
static void split_sci_entry(NbnxnPairlistGpu* nbl,
                            int               nsp_target_av,
                            gmx_bool          progBal,
                            float             nsp_tot_est,
                            int               thread,
                            int               nthread)
{
    int nsp_max;

    if (progBal)
    {
        float nsp_est;

        /* Estimate the total numbers of ci's of the nblist combined
         * over all threads using the target number of ci's.
         */
        nsp_est = (nsp_tot_est * thread) / nthread + nbl->nci_tot;

        /* The first ci blocks should be larger, to avoid overhead.
         * The last ci blocks should be smaller, to improve load balancing.
         * The factor 3/2 makes the first block 3/2 times the target average
         * and ensures that the total number of blocks end up equal to
         * that of equally sized blocks of size nsp_target_av.
         */
        nsp_max = static_cast<int>(nsp_target_av * (nsp_tot_est * 1.5 / (nsp_est + nsp_tot_est)));
    }
    else
    {
        nsp_max = nsp_target_av;
    }

    const int cj4_start = nbl->sci.back().cj4_ind_start;
    const int cj4_end   = nbl->sci.back().cj4_ind_end;
    const int j4len     = cj4_end - cj4_start;

    if (j4len > 1 && j4len * c_gpuNumClusterPerCell * c_nbnxnGpuJgroupSize > nsp_max)
    {
        /* Modify the last ci entry and process the cj4's again */

        int nsp       = 0;
        int nsp_sci   = 0;
        int nsp_cj4_e = 0;
        int nsp_cj4   = 0;
        for (int cj4 = cj4_start; cj4 < cj4_end; cj4++)
        {
            int nsp_cj4_p = nsp_cj4;
            /* Count the number of cluster pairs in this cj4 group */
            nsp_cj4 = 0;
            for (int p = 0; p < c_gpuNumClusterPerCell * c_nbnxnGpuJgroupSize; p++)
            {
                nsp_cj4 += (nbl->cj4[cj4].imei[0].imask >> p) & 1;
            }

            /* If adding the current cj4 with nsp_cj4 pairs get us further
             * away from our target nsp_max, split the list before this cj4.
             */
            if (nsp > 0 && nsp_max - nsp < nsp + nsp_cj4 - nsp_max)
            {
                /* Split the list at cj4 */
                nbl->sci.back().cj4_ind_end = cj4;
                /* Create a new sci entry */
                nbnxn_sci_t sciNew;
                sciNew.sci           = nbl->sci.back().sci;
                sciNew.shift         = nbl->sci.back().shift;
                sciNew.cj4_ind_start = cj4;
                nbl->sci.push_back(sciNew);

                nsp_sci   = nsp;
                nsp_cj4_e = nsp_cj4_p;
                nsp       = 0;
            }
            nsp += nsp_cj4;
        }

        /* Put the remaining cj4's in the last sci entry */
        nbl->sci.back().cj4_ind_end = cj4_end;

        /* Possibly balance out the last two sci's
         * by moving the last cj4 of the second last sci.
         */
        if (nsp_sci - nsp_cj4_e >= nsp + nsp_cj4_e)
        {
            GMX_ASSERT(nbl->sci.size() >= 2, "We expect at least two elements");
            nbl->sci[nbl->sci.size() - 2].cj4_ind_end--;
            nbl->sci[nbl->sci.size() - 1].cj4_ind_start--;
        }
    }
}

/* Clost this super/sub list i entry */
static void closeIEntry(NbnxnPairlistGpu* nbl, int nsp_max_av, gmx_bool progBal, float nsp_tot_est, int thread, int nthread)
{
    nbnxn_sci_t& sciEntry = *getOpenIEntry(nbl);

    /* All content of the new ci entry have already been filled correctly,
     * we only need to, potentially, split or remove the entry when empty.
     */
    int j4len = sciEntry.numJClusterGroups();
    if (j4len > 0)
    {
        /* We can only have complete blocks of 4 j-entries in a list,
         * so round the count up before closing.
         */
        int ncj4          = (nbl->work->cj_ind + c_nbnxnGpuJgroupSize - 1) / c_nbnxnGpuJgroupSize;
        nbl->work->cj_ind = ncj4 * c_nbnxnGpuJgroupSize;

        if (nsp_max_av > 0)
        {
            /* Measure the size of the new entry and potentially split it */
            split_sci_entry(nbl, nsp_max_av, progBal, nsp_tot_est, thread, nthread);
        }
    }
    else
    {
        /* Entry is empty: remove it  */
        nbl->sci.pop_back();
    }
}

/* Syncs the working array before adding another grid pair to the GPU list */
static void sync_work(NbnxnPairlistCpu gmx_unused* nbl) {}

/* Syncs the working array before adding another grid pair to the GPU list */
static void sync_work(NbnxnPairlistGpu* nbl)
{
    nbl->work->cj_ind = nbl->cj4.size() * c_nbnxnGpuJgroupSize;
}

/* Clears an NbnxnPairlistCpu data structure */
static void clear_pairlist(NbnxnPairlistCpu* nbl)
{
    nbl->ci.clear();
    nbl->cj.clear();
    nbl->ncjInUse = 0;
    nbl->nci_tot  = 0;
    nbl->ciOuter.clear();
    nbl->cjOuter.clear();

    nbl->work->ncj_noq = 0;
    nbl->work->ncj_hlj = 0;
}

/* Clears an NbnxnPairlistGpu data structure */
static void clear_pairlist(NbnxnPairlistGpu* nbl)
{
    nbl->sci.clear();
    nbl->cj4.clear();
    nbl->excl.resize(1);
    nbl->nci_tot = 0;
}

/* Clears an atom-atom-style pair list */
static void clear_pairlist_fep(t_nblist* nl)
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
static inline void
set_icell_bb_simple(gmx::ArrayRef<const BoundingBox> bb, int ci, real shx, real shy, real shz, BoundingBox* bb_ci)
{
    bb_ci->lower.x = bb[ci].lower.x + shx;
    bb_ci->lower.y = bb[ci].lower.y + shy;
    bb_ci->lower.z = bb[ci].lower.z + shz;
    bb_ci->upper.x = bb[ci].upper.x + shx;
    bb_ci->upper.y = bb[ci].upper.y + shy;
    bb_ci->upper.z = bb[ci].upper.z + shz;
}

/* Sets a simple list i-cell bounding box, including PBC shift */
static inline void set_icell_bb(const Grid& iGrid, int ci, real shx, real shy, real shz, NbnxnPairlistCpuWork* work)
{
    set_icell_bb_simple(iGrid.iBoundingBoxes(), ci, shx, shy, shz, &work->iClusterData.bb[0]);
}

#if NBNXN_BBXXXX
/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
static void set_icell_bbxxxx_supersub(gmx::ArrayRef<const float> bb, int ci, real shx, real shy, real shz, float* bb_ci)
{
    constexpr int cellBBStride = packedBoundingBoxesIndex(c_gpuNumClusterPerCell);
    constexpr int pbbStride    = c_packedBoundingBoxesDimSize;
    const int     ia           = ci * cellBBStride;
    for (int m = 0; m < cellBBStride; m += c_packedBoundingBoxesSize)
    {
        for (int i = 0; i < pbbStride; i++)
        {
            bb_ci[m + 0 * pbbStride + i] = bb[ia + m + 0 * pbbStride + i] + shx;
            bb_ci[m + 1 * pbbStride + i] = bb[ia + m + 1 * pbbStride + i] + shy;
            bb_ci[m + 2 * pbbStride + i] = bb[ia + m + 2 * pbbStride + i] + shz;
            bb_ci[m + 3 * pbbStride + i] = bb[ia + m + 3 * pbbStride + i] + shx;
            bb_ci[m + 4 * pbbStride + i] = bb[ia + m + 4 * pbbStride + i] + shy;
            bb_ci[m + 5 * pbbStride + i] = bb[ia + m + 5 * pbbStride + i] + shz;
        }
    }
}
#endif

/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
gmx_unused static void set_icell_bb_supersub(gmx::ArrayRef<const BoundingBox> bb,
                                             int                              ci,
                                             real                             shx,
                                             real                             shy,
                                             real                             shz,
                                             BoundingBox*                     bb_ci)
{
    for (int i = 0; i < c_gpuNumClusterPerCell; i++)
    {
        set_icell_bb_simple(bb, ci * c_gpuNumClusterPerCell + i, shx, shy, shz, &bb_ci[i]);
    }
}

/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
gmx_unused static void set_icell_bb(const Grid& iGrid, int ci, real shx, real shy, real shz, NbnxnPairlistGpuWork* work)
{
#if NBNXN_BBXXXX
    set_icell_bbxxxx_supersub(iGrid.packedBoundingBoxes(), ci, shx, shy, shz,
                              work->iSuperClusterData.bbPacked.data());
#else
    set_icell_bb_supersub(iGrid.iBoundingBoxes(), ci, shx, shy, shz, work->iSuperClusterData.bb.data());
#endif
}

/* Copies PBC shifted i-cell atom coordinates x,y,z to working array */
static void icell_set_x_simple(int                                 ci,
                               real                                shx,
                               real                                shy,
                               real                                shz,
                               int                                 stride,
                               const real*                         x,
                               NbnxnPairlistCpuWork::IClusterData* iClusterData)
{
    const int ia = ci * c_nbnxnCpuIClusterSize;

    for (int i = 0; i < c_nbnxnCpuIClusterSize; i++)
    {
        iClusterData->x[i * STRIDE_XYZ + XX] = x[(ia + i) * stride + XX] + shx;
        iClusterData->x[i * STRIDE_XYZ + YY] = x[(ia + i) * stride + YY] + shy;
        iClusterData->x[i * STRIDE_XYZ + ZZ] = x[(ia + i) * stride + ZZ] + shz;
    }
}

static void icell_set_x(int                             ci,
                        real                            shx,
                        real                            shy,
                        real                            shz,
                        int                             stride,
                        const real*                     x,
                        const ClusterDistanceKernelType kernelType,
                        NbnxnPairlistCpuWork*           work)
{
    switch (kernelType)
    {
#if GMX_SIMD
#    ifdef GMX_NBNXN_SIMD_4XN
        case ClusterDistanceKernelType::CpuSimd_4xM:
            icell_set_x_simd_4xn(ci, shx, shy, shz, stride, x, work);
            break;
#    endif
#    ifdef GMX_NBNXN_SIMD_2XNN
        case ClusterDistanceKernelType::CpuSimd_2xMM:
            icell_set_x_simd_2xnn(ci, shx, shy, shz, stride, x, work);
            break;
#    endif
#endif
        case ClusterDistanceKernelType::CpuPlainC:
            icell_set_x_simple(ci, shx, shy, shz, stride, x, &work->iClusterData);
            break;
        default: GMX_ASSERT(false, "Unhandled case"); break;
    }
}

/* Copies PBC shifted super-cell atom coordinates x,y,z to working array */
static void icell_set_x(int                       ci,
                        real                      shx,
                        real                      shy,
                        real                      shz,
                        int                       stride,
                        const real*               x,
                        ClusterDistanceKernelType gmx_unused kernelType,
                        NbnxnPairlistGpuWork*                work)
{
#if !GMX_SIMD4_HAVE_REAL

    real* x_ci = work->iSuperClusterData.x.data();

    int ia = ci * c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize;
    for (int i = 0; i < c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize; i++)
    {
        x_ci[i * DIM + XX] = x[(ia + i) * stride + XX] + shx;
        x_ci[i * DIM + YY] = x[(ia + i) * stride + YY] + shy;
        x_ci[i * DIM + ZZ] = x[(ia + i) * stride + ZZ] + shz;
    }

#else /* !GMX_SIMD4_HAVE_REAL */

    real* x_ci = work->iSuperClusterData.xSimd.data();

    for (int si = 0; si < c_gpuNumClusterPerCell; si++)
    {
        for (int i = 0; i < c_nbnxnGpuClusterSize; i += GMX_SIMD4_WIDTH)
        {
            int io = si * c_nbnxnGpuClusterSize + i;
            int ia = ci * c_gpuNumClusterPerCell * c_nbnxnGpuClusterSize + io;
            for (int j = 0; j < GMX_SIMD4_WIDTH; j++)
            {
                x_ci[io * DIM + j + XX * GMX_SIMD4_WIDTH] = x[(ia + j) * stride + XX] + shx;
                x_ci[io * DIM + j + YY * GMX_SIMD4_WIDTH] = x[(ia + j) * stride + YY] + shy;
                x_ci[io * DIM + j + ZZ * GMX_SIMD4_WIDTH] = x[(ia + j) * stride + ZZ] + shz;
            }
        }
    }

#endif /* !GMX_SIMD4_HAVE_REAL */
}

static real minimum_subgrid_size_xy(const Grid& grid)
{
    const Grid::Dimensions& dims = grid.dimensions();

    if (grid.geometry().isSimple)
    {
        return std::min(dims.cellSize[XX], dims.cellSize[YY]);
    }
    else
    {
        return std::min(dims.cellSize[XX] / c_gpuNumClusterPerCellX,
                        dims.cellSize[YY] / c_gpuNumClusterPerCellY);
    }
}

static real effective_buffer_1x1_vs_MxN(const Grid& iGrid, const Grid& jGrid)
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
    return eff_1x1_buffer_fac_overest * (minimum_subgrid_size_xy(iGrid) + minimum_subgrid_size_xy(jGrid));
}

/* Estimates the interaction volume^2 for non-local interactions */
static real nonlocal_vol2(const struct gmx_domdec_zones_t* zones, const rvec ls, real r)
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
                    cl += 0.5 * ls[d];
                    ca *= ls[d];
                    za *= zones->size[z].x1[d] - zones->size[z].x0[d];
                }
            }

            /* 4 octants of a sphere */
            vold_est = 0.25 * M_PI * r * r * r * r;
            /* 4 quarter pie slices on the edges */
            vold_est += 4 * cl * M_PI / 6.0 * r * r * r;
            /* One rectangular volume on a face */
            vold_est += ca * 0.5 * r * r;

            vol2_est_tot += vold_est * za;
        }
    }

    return vol2_est_tot;
}

/* Estimates the average size of a full j-list for super/sub setup */
static void get_nsubpair_target(const Nbnxm::GridSet&     gridSet,
                                const InteractionLocality iloc,
                                const real                rlist,
                                const int                 min_ci_balanced,
                                int*                      nsubpair_target,
                                float*                    nsubpair_tot_est)
{
    /* The target value of 36 seems to be the optimum for Kepler.
     * Maxwell is less sensitive to the exact value.
     */
    const int nsubpair_target_min = 36;
    real      r_eff_sup, vol_est, nsp_est, nsp_est_nl;

    const Grid& grid = gridSet.grids()[0];

    /* We don't need to balance list sizes if:
     * - We didn't request balancing.
     * - The number of grid cells >= the number of lists requested,
     *   since we will always generate at least #cells lists.
     * - We don't have any cells, since then there won't be any lists.
     */
    if (min_ci_balanced <= 0 || grid.numCells() >= min_ci_balanced || grid.numCells() == 0)
    {
        /* nsubpair_target==0 signals no balancing */
        *nsubpair_target  = 0;
        *nsubpair_tot_est = 0;

        return;
    }

    gmx::RVec               ls;
    const int               numAtomsCluster = grid.geometry().numAtomsICluster;
    const Grid::Dimensions& dims            = grid.dimensions();

    ls[XX] = dims.cellSize[XX] / c_gpuNumClusterPerCellX;
    ls[YY] = dims.cellSize[YY] / c_gpuNumClusterPerCellY;
    ls[ZZ] = numAtomsCluster / (dims.atomDensity * ls[XX] * ls[YY]);

    /* The formulas below are a heuristic estimate of the average nsj per si*/
    r_eff_sup = rlist + nbnxn_get_rlist_effective_inc(numAtomsCluster, ls);

    if (!gridSet.domainSetup().haveMultipleDomains || gridSet.domainSetup().zones->n == 1)
    {
        nsp_est_nl = 0;
    }
    else
    {
        nsp_est_nl = gmx::square(dims.atomDensity / numAtomsCluster)
                     * nonlocal_vol2(gridSet.domainSetup().zones, ls, r_eff_sup);
    }

    if (iloc == InteractionLocality::Local)
    {
        /* Sub-cell interacts with itself */
        vol_est = ls[XX] * ls[YY] * ls[ZZ];
        /* 6/2 rectangular volume on the faces */
        vol_est += (ls[XX] * ls[YY] + ls[XX] * ls[ZZ] + ls[YY] * ls[ZZ]) * r_eff_sup;
        /* 12/2 quarter pie slices on the edges */
        vol_est += 2 * (ls[XX] + ls[YY] + ls[ZZ]) * 0.25 * M_PI * gmx::square(r_eff_sup);
        /* 4 octants of a sphere */
        vol_est += 0.5 * 4.0 / 3.0 * M_PI * gmx::power3(r_eff_sup);

        /* Estimate the number of cluster pairs as the local number of
         * clusters times the volume they interact with times the density.
         */
        nsp_est = grid.numClusters() * vol_est * dims.atomDensity / numAtomsCluster;

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
        nsp_est = std::max(nsp_est, grid.numClusters() * 14._real);

        if (debug)
        {
            fprintf(debug, "nsp_est local %5.1f non-local %5.1f\n", nsp_est, nsp_est_nl);
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
    *nsubpair_target  = std::max(nsubpair_target_min, roundToInt(nsp_est / min_ci_balanced));
    *nsubpair_tot_est = nsp_est;

    if (debug)
    {
        fprintf(debug, "nbl nsp estimate %.1f, nsubpair_target %d\n", nsp_est, *nsubpair_target);
    }
}

/* Debug list print function */
static void print_nblist_ci_cj(FILE* fp, const NbnxnPairlistCpu& nbl)
{
    for (const nbnxn_ci_t& ciEntry : nbl.ci)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj %3d\n", ciEntry.ci, ciEntry.shift,
                ciEntry.cj_ind_end - ciEntry.cj_ind_start);

        for (int j = ciEntry.cj_ind_start; j < ciEntry.cj_ind_end; j++)
        {
            fprintf(fp, "  cj %5d  imask %x\n", nbl.cj[j].cj, nbl.cj[j].excl);
        }
    }
}

/* Debug list print function */
static void print_nblist_sci_cj(FILE* fp, const NbnxnPairlistGpu& nbl)
{
    for (const nbnxn_sci_t& sci : nbl.sci)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d\n", sci.sci, sci.shift, sci.numJClusterGroups());

        int ncp = 0;
        for (int j4 = sci.cj4_ind_start; j4 < sci.cj4_ind_end; j4++)
        {
            for (int j = 0; j < c_nbnxnGpuJgroupSize; j++)
            {
                fprintf(fp, "  sj %5d  imask %x\n", nbl.cj4[j4].cj[j], nbl.cj4[j4].imei[0].imask);
                for (int si = 0; si < c_gpuNumClusterPerCell; si++)
                {
                    if (nbl.cj4[j4].imei[0].imask & (1U << (j * c_gpuNumClusterPerCell + si)))
                    {
                        ncp++;
                    }
                }
            }
        }
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d ncp %3d\n", sci.sci, sci.shift,
                sci.numJClusterGroups(), ncp);
    }
}

/* Combine pair lists *nbl generated on multiple threads nblc */
static void combine_nblists(gmx::ArrayRef<const NbnxnPairlistGpu> nbls, NbnxnPairlistGpu* nblc)
{
    int nsci  = nblc->sci.size();
    int ncj4  = nblc->cj4.size();
    int nexcl = nblc->excl.size();
    for (auto& nbl : nbls)
    {
        nsci += nbl.sci.size();
        ncj4 += nbl.cj4.size();
        nexcl += nbl.excl.size();
    }

    /* Resize with the final, combined size, so we can fill in parallel */
    /* NOTE: For better performance we should use default initialization */
    nblc->sci.resize(nsci);
    nblc->cj4.resize(ncj4);
    nblc->excl.resize(nexcl);

    /* Each thread should copy its own data to the combined arrays,
     * as otherwise data will go back and forth between different caches.
     */
    const int gmx_unused nthreads = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (gmx::index n = 0; n < nbls.ssize(); n++)
    {
        try
        {
            /* Determine the offset in the combined data for our thread.
             * Note that the original sizes in nblc are lost.
             */
            int sci_offset  = nsci;
            int cj4_offset  = ncj4;
            int excl_offset = nexcl;

            for (gmx::index i = n; i < nbls.ssize(); i++)
            {
                sci_offset -= nbls[i].sci.size();
                cj4_offset -= nbls[i].cj4.size();
                excl_offset -= nbls[i].excl.size();
            }

            const NbnxnPairlistGpu& nbli = nbls[n];

            for (size_t i = 0; i < nbli.sci.size(); i++)
            {
                nblc->sci[sci_offset + i] = nbli.sci[i];
                nblc->sci[sci_offset + i].cj4_ind_start += cj4_offset;
                nblc->sci[sci_offset + i].cj4_ind_end += cj4_offset;
            }

            for (size_t j4 = 0; j4 < nbli.cj4.size(); j4++)
            {
                nblc->cj4[cj4_offset + j4] = nbli.cj4[j4];
                nblc->cj4[cj4_offset + j4].imei[0].excl_ind += excl_offset;
                nblc->cj4[cj4_offset + j4].imei[1].excl_ind += excl_offset;
            }

            for (size_t j4 = 0; j4 < nbli.excl.size(); j4++)
            {
                nblc->excl[excl_offset + j4] = nbli.excl[j4];
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    for (auto& nbl : nbls)
    {
        nblc->nci_tot += nbl.nci_tot;
    }
}

static void balance_fep_lists(gmx::ArrayRef<std::unique_ptr<t_nblist>> fepLists,
                              gmx::ArrayRef<PairsearchWork>            work)
{
    const int numLists = fepLists.ssize();

    if (numLists == 1)
    {
        /* Nothing to balance */
        return;
    }

    /* Count the total i-lists and pairs */
    int nri_tot = 0;
    int nrj_tot = 0;
    for (const auto& list : fepLists)
    {
        nri_tot += list->nri;
        nrj_tot += list->nrj;
    }

    const int nrj_target = (nrj_tot + numLists - 1) / numLists;

    GMX_ASSERT(gmx_omp_nthreads_get(emntNonbonded) == numLists,
               "We should have as many work objects as FEP lists");

#pragma omp parallel for schedule(static) num_threads(numLists)
    for (int th = 0; th < numLists; th++)
    {
        try
        {
            t_nblist* nbl = work[th].nbl_fep.get();

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
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    /* Loop over the source lists and assign and copy i-entries */
    int       th_dest = 0;
    t_nblist* nbld    = work[th_dest].nbl_fep.get();
    for (int th = 0; th < numLists; th++)
    {
        const t_nblist* nbls = fepLists[th].get();

        for (int i = 0; i < nbls->nri; i++)
        {
            int nrj;

            /* The number of pairs in this i-entry */
            nrj = nbls->jindex[i + 1] - nbls->jindex[i];

            /* Decide if list th_dest is too large and we should procede
             * to the next destination list.
             */
            if (th_dest + 1 < numLists && nbld->nrj > 0
                && nbld->nrj + nrj - nrj_target > nrj_target - nbld->nrj)
            {
                th_dest++;
                nbld = work[th_dest].nbl_fep.get();
            }

            nbld->iinr[nbld->nri]  = nbls->iinr[i];
            nbld->gid[nbld->nri]   = nbls->gid[i];
            nbld->shift[nbld->nri] = nbls->shift[i];

            for (int j = nbls->jindex[i]; j < nbls->jindex[i + 1]; j++)
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
    for (int th = 0; th < numLists; th++)
    {
        fepLists[th].swap(work[th].nbl_fep);

        if (debug)
        {
            fprintf(debug, "nbl_fep[%d] nri %4d nrj %4d\n", th, fepLists[th]->nri, fepLists[th]->nrj);
        }
    }
}

/* Returns the next ci to be processes by our thread */
static gmx_bool next_ci(const Grid& grid, int nth, int ci_block, int* ci_x, int* ci_y, int* ci_b, int* ci)
{
    (*ci_b)++;
    (*ci)++;

    if (*ci_b == ci_block)
    {
        /* Jump to the next block assigned to this task */
        *ci += (nth - 1) * ci_block;
        *ci_b = 0;
    }

    if (*ci >= grid.numCells())
    {
        return FALSE;
    }

    while (*ci >= grid.firstCellInColumn(*ci_x * grid.dimensions().numCells[YY] + *ci_y + 1))
    {
        *ci_y += 1;
        if (*ci_y == grid.dimensions().numCells[YY])
        {
            *ci_x += 1;
            *ci_y = 0;
        }
    }

    return TRUE;
}

/* Returns the distance^2 for which we put cell pairs in the list
 * without checking atom pair distances. This is usually < rlist^2.
 */
static float boundingbox_only_distance2(const Grid::Dimensions& iGridDims,
                                        const Grid::Dimensions& jGridDims,
                                        real                    rlist,
                                        gmx_bool                simple)
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

    bbx = 0.5 * (iGridDims.cellSize[XX] + jGridDims.cellSize[XX]);
    bby = 0.5 * (iGridDims.cellSize[YY] + jGridDims.cellSize[YY]);
    if (!simple)
    {
        bbx /= c_gpuNumClusterPerCellX;
        bby /= c_gpuNumClusterPerCellY;
    }

    rbb2 = std::max(0.0, rlist - 0.5 * std::sqrt(bbx * bbx + bby * bby));
    rbb2 = rbb2 * rbb2;

#if !GMX_DOUBLE
    return rbb2;
#else
    return (float)((1 + GMX_FLOAT_EPS) * rbb2);
#endif
}

static int get_ci_block_size(const Grid& iGrid, const bool haveMultipleDomains, const int numLists)
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
    GMX_ASSERT(iGrid.dimensions().numCells[XX] > 0, "Grid can't be empty");
    GMX_ASSERT(numLists > 0, "We need at least one list");
    ci_block = (iGrid.numCells() * ci_block_enum)
               / (ci_block_denom * iGrid.dimensions().numCells[XX] * numLists);

    const int numAtomsPerCell = iGrid.geometry().numAtomsPerCell;

    /* Ensure the blocks are not too small: avoids cache invalidation */
    if (ci_block * numAtomsPerCell < ci_block_min_atoms)
    {
        ci_block = (ci_block_min_atoms + numAtomsPerCell - 1) / numAtomsPerCell;
    }

    /* Without domain decomposition
     * or with less than 3 blocks per task, divide in nth blocks.
     */
    if (!haveMultipleDomains || numLists * 3 * ci_block > iGrid.numCells())
    {
        ci_block = (iGrid.numCells() + numLists - 1) / numLists;
    }

    if (ci_block > 1 && (numLists - 1) * ci_block >= iGrid.numCells())
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

static bool pairlistIsSimple(const NbnxnPairlistCpu gmx_unused& pairlist)
{
    return true;
}

static bool pairlistIsSimple(const NbnxnPairlistGpu gmx_unused& pairlist)
{
    return false;
}

static void makeClusterListWrapper(NbnxnPairlistCpu* nbl,
                                   const Grid gmx_unused&          iGrid,
                                   const int                       ci,
                                   const Grid&                     jGrid,
                                   const int                       firstCell,
                                   const int                       lastCell,
                                   const bool                      excludeSubDiagonal,
                                   const nbnxn_atomdata_t*         nbat,
                                   const real                      rlist2,
                                   const real                      rbb2,
                                   const ClusterDistanceKernelType kernelType,
                                   int*                            numDistanceChecks)
{
    switch (kernelType)
    {
        case ClusterDistanceKernelType::CpuPlainC:
            makeClusterListSimple(jGrid, nbl, ci, firstCell, lastCell, excludeSubDiagonal,
                                  nbat->x().data(), rlist2, rbb2, numDistanceChecks);
            break;
#ifdef GMX_NBNXN_SIMD_4XN
        case ClusterDistanceKernelType::CpuSimd_4xM:
            makeClusterListSimd4xn(jGrid, nbl, ci, firstCell, lastCell, excludeSubDiagonal,
                                   nbat->x().data(), rlist2, rbb2, numDistanceChecks);
            break;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        case ClusterDistanceKernelType::CpuSimd_2xMM:
            makeClusterListSimd2xnn(jGrid, nbl, ci, firstCell, lastCell, excludeSubDiagonal,
                                    nbat->x().data(), rlist2, rbb2, numDistanceChecks);
            break;
#endif
        default: GMX_ASSERT(false, "Unhandled kernel type");
    }
}

static void makeClusterListWrapper(NbnxnPairlistGpu* nbl,
                                   const Grid& gmx_unused    iGrid,
                                   const int                 ci,
                                   const Grid&               jGrid,
                                   const int                 firstCell,
                                   const int                 lastCell,
                                   const bool                excludeSubDiagonal,
                                   const nbnxn_atomdata_t*   nbat,
                                   const real                rlist2,
                                   const real                rbb2,
                                   ClusterDistanceKernelType gmx_unused kernelType,
                                   int*                                 numDistanceChecks)
{
    for (int cj = firstCell; cj <= lastCell; cj++)
    {
        make_cluster_list_supersub(iGrid, jGrid, nbl, ci, cj, excludeSubDiagonal, nbat->xstride,
                                   nbat->x().data(), rlist2, rbb2, numDistanceChecks);
    }
}

static int getNumSimpleJClustersInList(const NbnxnPairlistCpu& nbl)
{
    return nbl.cj.size();
}

static int getNumSimpleJClustersInList(const gmx_unused NbnxnPairlistGpu& nbl)
{
    return 0;
}

static void incrementNumSimpleJClustersInList(NbnxnPairlistCpu* nbl, int ncj_old_j)
{
    nbl->ncjInUse += nbl->cj.size();
    nbl->ncjInUse -= ncj_old_j;
}

static void incrementNumSimpleJClustersInList(NbnxnPairlistGpu gmx_unused* nbl, int gmx_unused ncj_old_j)
{
}

static void checkListSizeConsistency(const NbnxnPairlistCpu& nbl, const bool haveFreeEnergy)
{
    GMX_RELEASE_ASSERT(static_cast<size_t>(nbl.ncjInUse) == nbl.cj.size() || haveFreeEnergy,
                       "Without free-energy all cj pair-list entries should be in use. "
                       "Note that subsequent code does not make use of the equality, "
                       "this check is only here to catch bugs");
}

static void checkListSizeConsistency(const NbnxnPairlistGpu gmx_unused& nbl, bool gmx_unused haveFreeEnergy)
{
    /* We currently can not check consistency here */
}

/* Set the buffer flags for newly added entries in the list */
static void setBufferFlags(const NbnxnPairlistCpu& nbl,
                           const int               ncj_old_j,
                           const int               gridj_flag_shift,
                           gmx_bitmask_t*          gridj_flag,
                           const int               th)
{
    if (gmx::ssize(nbl.cj) > ncj_old_j)
    {
        int cbFirst = nbl.cj[ncj_old_j].cj >> gridj_flag_shift;
        int cbLast  = nbl.cj.back().cj >> gridj_flag_shift;
        for (int cb = cbFirst; cb <= cbLast; cb++)
        {
            bitmask_init_bit(&gridj_flag[cb], th);
        }
    }
}

static void setBufferFlags(const NbnxnPairlistGpu gmx_unused& nbl,
                           int gmx_unused ncj_old_j,
                           int gmx_unused gridj_flag_shift,
                           gmx_bitmask_t gmx_unused* gridj_flag,
                           int gmx_unused th)
{
    GMX_ASSERT(false, "This function should never be called");
}

/* Generates the part of pair-list nbl assigned to our thread */
template<typename T>
static void nbnxn_make_pairlist_part(const Nbnxm::GridSet&   gridSet,
                                     const Grid&             iGrid,
                                     const Grid&             jGrid,
                                     PairsearchWork*         work,
                                     const nbnxn_atomdata_t* nbat,
                                     const t_blocka&         exclusions,
                                     real                    rlist,
                                     const PairlistType      pairlistType,
                                     int                     ci_block,
                                     gmx_bool                bFBufferFlag,
                                     int                     nsubpair_max,
                                     gmx_bool                progBal,
                                     float                   nsubpair_tot_est,
                                     int                     th,
                                     int                     nth,
                                     T*                      nbl,
                                     t_nblist*               nbl_fep)
{
    int            na_cj_2log;
    matrix         box;
    real           rl_fep2 = 0;
    float          rbb2;
    int            ci_b, ci, ci_x, ci_y, ci_xy;
    ivec           shp;
    real           bx0, bx1, by0, by1, bz0, bz1;
    real           bz1_frac;
    real           d2cx, d2z, d2z_cx, d2z_cy, d2zx, d2zxy, d2xy;
    int            cxf, cxl, cyf, cyf_x, cyl;
    int            numDistanceChecks;
    int            gridi_flag_shift = 0, gridj_flag_shift = 0;
    gmx_bitmask_t* gridj_flag = nullptr;
    int            ncj_old_i, ncj_old_j;

    if (jGrid.geometry().isSimple != pairlistIsSimple(*nbl)
        || iGrid.geometry().isSimple != pairlistIsSimple(*nbl))
    {
        gmx_incons("Grid incompatible with pair-list");
    }

    sync_work(nbl);
    GMX_ASSERT(nbl->na_ci == jGrid.geometry().numAtomsICluster,
               "The cluster sizes in the list and grid should match");
    nbl->na_cj = JClusterSizePerListType[pairlistType];
    na_cj_2log = get_2log(nbl->na_cj);

    nbl->rlist = rlist;

    if (bFBufferFlag)
    {
        /* Determine conversion of clusters to flag blocks */
        gridi_flag_shift = getBufferFlagShift(nbl->na_ci);
        gridj_flag_shift = getBufferFlagShift(nbl->na_cj);

        gridj_flag = work->buffer_flags.flag;
    }

    gridSet.getBox(box);

    const bool haveFep = gridSet.haveFep();

    const real rlist2 = nbl->rlist * nbl->rlist;

    // Select the cluster pair distance kernel type
    const ClusterDistanceKernelType kernelType = getClusterDistanceKernelType(pairlistType, *nbat);

    if (haveFep && !pairlistIsSimple(*nbl))
    {
        /* Determine an atom-pair list cut-off distance for FEP atom pairs.
         * We should not simply use rlist, since then we would not have
         * the small, effective buffering of the NxN lists.
         * The buffer is on overestimate, but the resulting cost for pairs
         * beyond rlist is neglible compared to the FEP pairs within rlist.
         */
        rl_fep2 = nbl->rlist + effective_buffer_1x1_vs_MxN(iGrid, jGrid);

        if (debug)
        {
            fprintf(debug, "nbl_fep atom-pair rlist %f\n", rl_fep2);
        }
        rl_fep2 = rl_fep2 * rl_fep2;
    }

    const Grid::Dimensions& iGridDims = iGrid.dimensions();
    const Grid::Dimensions& jGridDims = jGrid.dimensions();

    rbb2 = boundingbox_only_distance2(iGridDims, jGridDims, nbl->rlist, pairlistIsSimple(*nbl));

    if (debug)
    {
        fprintf(debug, "nbl bounding box only distance %f\n", std::sqrt(rbb2));
    }

    const bool isIntraGridList = (&iGrid == &jGrid);

    /* Set the shift range */
    for (int d = 0; d < DIM; d++)
    {
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        if (d >= ePBC2npbcdim(gridSet.domainSetup().ePBC)
            || gridSet.domainSetup().haveMultipleDomainsPerDim[d])
        {
            shp[d] = 0;
        }
        else
        {
            const real listRangeCellToCell =
                    listRangeForGridCellToGridCell(rlist, iGrid.dimensions(), jGrid.dimensions());
            if (d == XX && box[XX][XX] - fabs(box[YY][XX]) - fabs(box[ZZ][XX]) < listRangeCellToCell)
            {
                shp[d] = 2;
            }
            else
            {
                shp[d] = 1;
            }
        }
    }
    const bool                       bSimple = pairlistIsSimple(*nbl);
    gmx::ArrayRef<const BoundingBox> bb_i;
#if NBNXN_BBXXXX
    gmx::ArrayRef<const float> pbb_i;
    if (bSimple)
    {
        bb_i = iGrid.iBoundingBoxes();
    }
    else
    {
        pbb_i = iGrid.packedBoundingBoxes();
    }
#else
    /* We use the normal bounding box format for both grid types */
    bb_i = iGrid.iBoundingBoxes();
#endif
    gmx::ArrayRef<const BoundingBox1D> bbcz_i  = iGrid.zBoundingBoxes();
    gmx::ArrayRef<const int>           flags_i = iGrid.clusterFlags();
    gmx::ArrayRef<const BoundingBox1D> bbcz_j  = jGrid.zBoundingBoxes();
    int                                cell0_i = iGrid.cellOffset();

    if (debug)
    {
        fprintf(debug, "nbl nc_i %d col.av. %.1f ci_block %d\n", iGrid.numCells(),
                iGrid.numCells() / static_cast<double>(iGrid.numColumns()), ci_block);
    }

    numDistanceChecks = 0;

    const real listRangeBBToJCell2 =
            gmx::square(listRangeForBoundingBoxToGridCell(rlist, jGrid.dimensions()));

    /* Initially ci_b and ci to 1 before where we want them to start,
     * as they will both be incremented in next_ci.
     */
    ci_b = -1;
    ci   = th * ci_block - 1;
    ci_x = 0;
    ci_y = 0;
    while (next_ci(iGrid, nth, ci_block, &ci_x, &ci_y, &ci_b, &ci))
    {
        if (bSimple && flags_i[ci] == 0)
        {
            continue;
        }
        ncj_old_i = getNumSimpleJClustersInList(*nbl);

        d2cx = 0;
        if (!isIntraGridList && shp[XX] == 0)
        {
            if (bSimple)
            {
                bx1 = bb_i[ci].upper.x;
            }
            else
            {
                bx1 = iGridDims.lowerCorner[XX] + (real(ci_x) + 1) * iGridDims.cellSize[XX];
            }
            if (bx1 < jGridDims.lowerCorner[XX])
            {
                d2cx = gmx::square(jGridDims.lowerCorner[XX] - bx1);

                if (d2cx >= listRangeBBToJCell2)
                {
                    continue;
                }
            }
        }

        ci_xy = ci_x * iGridDims.numCells[YY] + ci_y;

        /* Loop over shift vectors in three dimensions */
        for (int tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
        {
            const real shz = real(tz) * box[ZZ][ZZ];

            bz0 = bbcz_i[ci].lower + shz;
            bz1 = bbcz_i[ci].upper + shz;

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

            bz1_frac = bz1 / real(iGrid.numCellsInColumn(ci_xy));
            if (bz1_frac < 0)
            {
                bz1_frac = 0;
            }
            /* The check with bz1_frac close to or larger than 1 comes later */

            for (int ty = -shp[YY]; ty <= shp[YY]; ty++)
            {
                const real shy = real(ty) * box[YY][YY] + real(tz) * box[ZZ][YY];

                if (bSimple)
                {
                    by0 = bb_i[ci].lower.y + shy;
                    by1 = bb_i[ci].upper.y + shy;
                }
                else
                {
                    by0 = iGridDims.lowerCorner[YY] + (real(ci_y)) * iGridDims.cellSize[YY] + shy;
                    by1 = iGridDims.lowerCorner[YY] + (real(ci_y) + 1) * iGridDims.cellSize[YY] + shy;
                }

                get_cell_range<YY>(by0, by1, jGridDims, d2z_cx, rlist, &cyf, &cyl);

                if (cyf > cyl)
                {
                    continue;
                }

                d2z_cy = d2z;
                if (by1 < jGridDims.lowerCorner[YY])
                {
                    d2z_cy += gmx::square(jGridDims.lowerCorner[YY] - by1);
                }
                else if (by0 > jGridDims.upperCorner[YY])
                {
                    d2z_cy += gmx::square(by0 - jGridDims.upperCorner[YY]);
                }

                for (int tx = -shp[XX]; tx <= shp[XX]; tx++)
                {
                    const int shift = XYZ2IS(tx, ty, tz);

                    const bool excludeSubDiagonal = (isIntraGridList && shift == CENTRAL);

                    if (c_pbcShiftBackward && isIntraGridList && shift > CENTRAL)
                    {
                        continue;
                    }

                    const real shx =
                            real(tx) * box[XX][XX] + real(ty) * box[YY][XX] + real(tz) * box[ZZ][XX];

                    if (bSimple)
                    {
                        bx0 = bb_i[ci].lower.x + shx;
                        bx1 = bb_i[ci].upper.x + shx;
                    }
                    else
                    {
                        bx0 = iGridDims.lowerCorner[XX] + (real(ci_x)) * iGridDims.cellSize[XX] + shx;
                        bx1 = iGridDims.lowerCorner[XX] + (real(ci_x) + 1) * iGridDims.cellSize[XX] + shx;
                    }

                    get_cell_range<XX>(bx0, bx1, jGridDims, d2z_cy, rlist, &cxf, &cxl);

                    if (cxf > cxl)
                    {
                        continue;
                    }

                    addNewIEntry(nbl, cell0_i + ci, shift, flags_i[ci]);

                    if ((!c_pbcShiftBackward || excludeSubDiagonal) && cxf < ci_x)
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    set_icell_bb(iGrid, ci, shx, shy, shz, nbl->work.get());

                    icell_set_x(cell0_i + ci, shx, shy, shz, nbat->xstride, nbat->x().data(),
                                kernelType, nbl->work.get());

                    for (int cx = cxf; cx <= cxl; cx++)
                    {
                        const real cx_real = cx;
                        d2zx               = d2z;
                        if (jGridDims.lowerCorner[XX] + cx_real * jGridDims.cellSize[XX] > bx1)
                        {
                            d2zx += gmx::square(jGridDims.lowerCorner[XX]
                                                + cx_real * jGridDims.cellSize[XX] - bx1);
                        }
                        else if (jGridDims.lowerCorner[XX] + (cx_real + 1) * jGridDims.cellSize[XX] < bx0)
                        {
                            d2zx += gmx::square(jGridDims.lowerCorner[XX]
                                                + (cx_real + 1) * jGridDims.cellSize[XX] - bx0);
                        }

                        if (isIntraGridList && cx == 0 && (!c_pbcShiftBackward || shift == CENTRAL)
                            && cyf < ci_y)
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
                            const int columnStart =
                                    jGrid.firstCellInColumn(cx * jGridDims.numCells[YY] + cy);
                            const int columnEnd =
                                    jGrid.firstCellInColumn(cx * jGridDims.numCells[YY] + cy + 1);

                            const real cy_real = cy;
                            d2zxy              = d2zx;
                            if (jGridDims.lowerCorner[YY] + cy_real * jGridDims.cellSize[YY] > by1)
                            {
                                d2zxy += gmx::square(jGridDims.lowerCorner[YY]
                                                     + cy_real * jGridDims.cellSize[YY] - by1);
                            }
                            else if (jGridDims.lowerCorner[YY] + (cy_real + 1) * jGridDims.cellSize[YY] < by0)
                            {
                                d2zxy += gmx::square(jGridDims.lowerCorner[YY]
                                                     + (cy_real + 1) * jGridDims.cellSize[YY] - by0);
                            }
                            if (columnStart < columnEnd && d2zxy < listRangeBBToJCell2)
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
                                int midCell =
                                        columnStart
                                        + static_cast<int>(bz1_frac
                                                           * static_cast<real>(columnEnd - columnStart));
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
                                while (downTestCell >= columnStart
                                       && (bbcz_j[downTestCell].upper >= bz0
                                           || d2xy + gmx::square(bbcz_j[downTestCell].upper - bz0) < rlist2))
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
                                while (upTestCell < columnEnd
                                       && (bbcz_j[upTestCell].lower <= bz1
                                           || d2xy + gmx::square(bbcz_j[upTestCell].lower - bz1) < rlist2))
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
                                        if (d2xy + gmx::square(bbcz_j[k * NNBSBB_D + 1] - bz0) < rlist2
                                            && k < firstCell)
                                        {
                                            firstCell = k;
                                        }
                                        if (d2xy + gmx::square(bbcz_j[k * NNBSBB_D] - bz1) < rlist2
                                            && k > lastCell)
                                        {
                                            lastCell = k;
                                        }
                                    }
                                }
#endif

                                if (isIntraGridList)
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
                                    GMX_ASSERT(firstCell >= columnStart && lastCell < columnEnd,
                                               "The range should reside within the current grid "
                                               "column");

                                    /* For f buffer flags with simple lists */
                                    ncj_old_j = getNumSimpleJClustersInList(*nbl);

                                    makeClusterListWrapper(nbl, iGrid, ci, jGrid, firstCell, lastCell,
                                                           excludeSubDiagonal, nbat, rlist2, rbb2,
                                                           kernelType, &numDistanceChecks);

                                    if (bFBufferFlag)
                                    {
                                        setBufferFlags(*nbl, ncj_old_j, gridj_flag_shift, gridj_flag, th);
                                    }

                                    incrementNumSimpleJClustersInList(nbl, ncj_old_j);
                                }
                            }
                        }
                    }

                    /* Set the exclusions for this ci list */
                    setExclusionsForIEntry(gridSet, nbl, excludeSubDiagonal, na_cj_2log,
                                           *getOpenIEntry(nbl), exclusions);

                    if (haveFep)
                    {
                        make_fep_list(gridSet.atomIndices(), nbat, nbl, excludeSubDiagonal,
                                      getOpenIEntry(nbl), shx, shy, shz, rl_fep2, iGrid, jGrid, nbl_fep);
                    }

                    /* Close this ci list */
                    closeIEntry(nbl, nsubpair_max, progBal, nsubpair_tot_est, th, nth);
                }
            }
        }

        if (bFBufferFlag && getNumSimpleJClustersInList(*nbl) > ncj_old_i)
        {
            bitmask_init_bit(&(work->buffer_flags.flag[(iGrid.cellOffset() + ci) >> gridi_flag_shift]), th);
        }
    }

    work->ndistc = numDistanceChecks;

    checkListSizeConsistency(*nbl, haveFep);

    if (debug)
    {
        fprintf(debug, "number of distance checks %d\n", numDistanceChecks);

        print_nblist_statistics(debug, *nbl, gridSet, rlist);

        if (haveFep)
        {
            fprintf(debug, "nbl FEP list pairs: %d\n", nbl_fep->nrj);
        }
    }
}

static void reduce_buffer_flags(gmx::ArrayRef<PairsearchWork> searchWork,
                                int                           nsrc,
                                const nbnxn_buffer_flags_t*   dest)
{
    for (int s = 0; s < nsrc; s++)
    {
        gmx_bitmask_t* flag = searchWork[s].buffer_flags.flag;

        for (int b = 0; b < dest->nflag; b++)
        {
            bitmask_union(&(dest->flag[b]), flag[b]);
        }
    }
}

static void print_reduction_cost(const nbnxn_buffer_flags_t* flags, int nout)
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

    fprintf(debug,
            "nbnxn reduction: #flag %d #list %d elem %4.2f, keep %4.2f copy %4.2f red %4.2f\n",
            flags->nflag, nout, nelem / static_cast<double>(flags->nflag),
            nkeep / static_cast<double>(flags->nflag), ncopy / static_cast<double>(flags->nflag),
            nred / static_cast<double>(flags->nflag));
}

/* Copies the list entries from src to dest when cjStart <= *cjGlobal < cjEnd.
 * *cjGlobal is updated with the cj count in src.
 * When setFlags==true, flag bit t is set in flag for all i and j clusters.
 */
template<bool setFlags>
static void copySelectedListRange(const nbnxn_ci_t* gmx_restrict srcCi,
                                  const NbnxnPairlistCpu* gmx_restrict src,
                                  NbnxnPairlistCpu* gmx_restrict dest,
                                  gmx_bitmask_t*                 flag,
                                  int                            iFlagShift,
                                  int                            jFlagShift,
                                  int                            t)
{
    const int ncj = srcCi->cj_ind_end - srcCi->cj_ind_start;

    dest->ci.push_back(*srcCi);
    dest->ci.back().cj_ind_start = dest->cj.size();
    dest->ci.back().cj_ind_end   = dest->ci.back().cj_ind_start + ncj;

    if (setFlags)
    {
        bitmask_init_bit(&flag[srcCi->ci >> iFlagShift], t);
    }

    for (int j = srcCi->cj_ind_start; j < srcCi->cj_ind_end; j++)
    {
        dest->cj.push_back(src->cj[j]);

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
}

#if defined(__GNUC__) && !defined(__clang__) && !defined(__ICC) && __GNUC__ == 7
/* Avoid gcc 7 avx512 loop vectorization bug (actually only needed with -mavx512f) */
#    pragma GCC push_options
#    pragma GCC optimize("no-tree-vectorize")
#endif

/* Returns the number of cluster pairs that are in use summed over all lists */
static int countClusterpairs(gmx::ArrayRef<const NbnxnPairlistCpu> pairlists)
{
    /* gcc 7 with -mavx512f can miss the contributions of 16 consecutive
     * elements to the sum calculated in this loop. Above we have disabled
     * loop vectorization to avoid this bug.
     */
    int ncjTotal = 0;
    for (const auto& pairlist : pairlists)
    {
        ncjTotal += pairlist.ncjInUse;
    }
    return ncjTotal;
}

#if defined(__GNUC__) && !defined(__clang__) && !defined(__ICC) && __GNUC__ == 7
#    pragma GCC pop_options
#endif

/* This routine re-balances the pairlists such that all are nearly equally
 * sized. Only whole i-entries are moved between lists. These are moved
 * between the ends of the lists, such that the buffer reduction cost should
 * not change significantly.
 * Note that all original reduction flags are currently kept. This can lead
 * to reduction of parts of the force buffer that could be avoided. But since
 * the original lists are quite balanced, this will only give minor overhead.
 */
static void rebalanceSimpleLists(gmx::ArrayRef<const NbnxnPairlistCpu> srcSet,
                                 gmx::ArrayRef<NbnxnPairlistCpu>       destSet,
                                 gmx::ArrayRef<PairsearchWork>         searchWork)
{
    const int ncjTotal  = countClusterpairs(srcSet);
    const int numLists  = srcSet.ssize();
    const int ncjTarget = (ncjTotal + numLists - 1) / numLists;

#pragma omp parallel num_threads(numLists)
    {
        int t = gmx_omp_get_thread_num();

        int cjStart = ncjTarget * t;
        int cjEnd   = ncjTarget * (t + 1);

        /* The destination pair-list for task/thread t */
        NbnxnPairlistCpu& dest = destSet[t];

        clear_pairlist(&dest);
        dest.na_cj = srcSet[0].na_cj;

        /* Note that the flags in the work struct (still) contain flags
         * for all entries that are present in srcSet->nbl[t].
         */
        gmx_bitmask_t* flag = searchWork[t].buffer_flags.flag;

        int iFlagShift = getBufferFlagShift(dest.na_ci);
        int jFlagShift = getBufferFlagShift(dest.na_cj);

        int cjGlobal = 0;
        for (int s = 0; s < numLists && cjGlobal < cjEnd; s++)
        {
            const NbnxnPairlistCpu* src = &srcSet[s];

            if (cjGlobal + src->ncjInUse > cjStart)
            {
                for (gmx::index i = 0; i < gmx::ssize(src->ci) && cjGlobal < cjEnd; i++)
                {
                    const nbnxn_ci_t* srcCi = &src->ci[i];
                    int               ncj   = srcCi->cj_ind_end - srcCi->cj_ind_start;
                    if (cjGlobal >= cjStart)
                    {
                        /* If the source list is not our own, we need to set
                         * extra flags (the template bool parameter).
                         */
                        if (s != t)
                        {
                            copySelectedListRange<true>(srcCi, src, &dest, flag, iFlagShift, jFlagShift, t);
                        }
                        else
                        {
                            copySelectedListRange<false>(srcCi, src, &dest, flag, iFlagShift,
                                                         jFlagShift, t);
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

        dest.ncjInUse = dest.cj.size();
    }

#ifndef NDEBUG
    const int ncjTotalNew = countClusterpairs(destSet);
    GMX_RELEASE_ASSERT(ncjTotalNew == ncjTotal,
                       "The total size of the lists before and after rebalancing should match");
#endif
}

/* Returns if the pairlists are so imbalanced that it is worth rebalancing. */
static bool checkRebalanceSimpleLists(gmx::ArrayRef<const NbnxnPairlistCpu> lists)
{
    int numLists = lists.ssize();
    int ncjMax   = 0;
    int ncjTotal = 0;
    for (int s = 0; s < numLists; s++)
    {
        ncjMax = std::max(ncjMax, lists[s].ncjInUse);
        ncjTotal += lists[s].ncjInUse;
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

    return real(numLists * ncjMax) > real(ncjTotal) * rebalanceTolerance;
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
static void sort_sci(NbnxnPairlistGpu* nbl)
{
    if (nbl->cj4.size() <= nbl->sci.size())
    {
        /* nsci = 0 or all sci have size 1, sorting won't change the order */
        return;
    }

    NbnxnPairlistGpuWork& work = *nbl->work;

    /* We will distinguish differences up to double the average */
    const int m = static_cast<int>((2 * ssize(nbl->cj4)) / ssize(nbl->sci));

    /* Resize work.sci_sort so we can sort into it */
    work.sci_sort.resize(nbl->sci.size());

    std::vector<int>& sort = work.sortBuffer;
    /* Set up m + 1 entries in sort, initialized at 0 */
    sort.clear();
    sort.resize(m + 1, 0);
    /* Count the entries of each size */
    for (const nbnxn_sci_t& sci : nbl->sci)
    {
        int i = std::min(m, sci.numJClusterGroups());
        sort[i]++;
    }
    /* Calculate the offset for each count */
    int s0  = sort[m];
    sort[m] = 0;
    for (gmx::index i = m - 1; i >= 0; i--)
    {
        int s1  = sort[i];
        sort[i] = sort[i + 1] + s0;
        s0      = s1;
    }

    /* Sort entries directly into place */
    gmx::ArrayRef<nbnxn_sci_t> sci_sort = work.sci_sort;
    for (const nbnxn_sci_t& sci : nbl->sci)
    {
        int i               = std::min(m, sci.numJClusterGroups());
        sci_sort[sort[i]++] = sci;
    }

    /* Swap the sci pointers so we use the new, sorted list */
    std::swap(nbl->sci, work.sci_sort);
}

/* Returns the i-zone range for pairlist construction for the give locality */
static Range<int> getIZoneRange(const Nbnxm::GridSet::DomainSetup& domainSetup,
                                const InteractionLocality          locality)
{
    if (domainSetup.doTestParticleInsertion)
    {
        /* With TPI we do grid 1, the inserted molecule, versus grid 0, the rest */
        return { 1, 2 };
    }
    else if (locality == InteractionLocality::Local)
    {
        /* Local: only zone (grid) 0 vs 0 */
        return { 0, 1 };
    }
    else
    {
        /* Non-local: we need all i-zones */
        return { 0, int(domainSetup.zones->iZones.size()) };
    }
}

/* Returns the j-zone range for pairlist construction for the give locality and i-zone */
static Range<int> getJZoneRange(const gmx_domdec_zones_t& ddZones,
                                const InteractionLocality locality,
                                const int                 iZone)
{
    if (locality == InteractionLocality::Local)
    {
        /* Local: zone 0 vs 0 or with TPI 1 vs 0 */
        return { 0, 1 };
    }
    else if (iZone == 0)
    {
        /* Non-local: we need to avoid the local (zone 0 vs 0) interactions */
        return { 1, *ddZones.iZones[iZone].jZoneRange.end() };
    }
    else
    {
        /* Non-local with non-local i-zone: use all j-zones */
        return ddZones.iZones[iZone].jZoneRange;
    }
}

//! Prepares CPU lists produced by the search for dynamic pruning
static void prepareListsForDynamicPruning(gmx::ArrayRef<NbnxnPairlistCpu> lists);

void PairlistSet::constructPairlists(const Nbnxm::GridSet&         gridSet,
                                     gmx::ArrayRef<PairsearchWork> searchWork,
                                     nbnxn_atomdata_t*             nbat,
                                     const t_blocka*               excl,
                                     const int                     minimumIlistCountForGpuBalancing,
                                     t_nrnb*                       nrnb,
                                     SearchCycleCounting*          searchCycleCounting)
{
    const real rlist = params_.rlistOuter;

    int      nsubpair_target;
    float    nsubpair_tot_est;
    int      ci_block;
    gmx_bool progBal;
    int      np_tot, np_noq, np_hlj, nap;

    const int numLists = (isCpuType_ ? cpuLists_.size() : gpuLists_.size());

    if (debug)
    {
        fprintf(debug, "ns making %d nblists\n", numLists);
    }

    nbat->bUseBufferFlags = (nbat->out.size() > 1);
    /* We should re-init the flags before making the first list */
    if (nbat->bUseBufferFlags && locality_ == InteractionLocality::Local)
    {
        init_buffer_flags(&nbat->buffer_flags, nbat->numAtoms());
    }

    if (!isCpuType_ && minimumIlistCountForGpuBalancing > 0)
    {
        get_nsubpair_target(gridSet, locality_, rlist, minimumIlistCountForGpuBalancing,
                            &nsubpair_target, &nsubpair_tot_est);
    }
    else
    {
        nsubpair_target  = 0;
        nsubpair_tot_est = 0;
    }

    /* Clear all pair-lists */
    for (int th = 0; th < numLists; th++)
    {
        if (isCpuType_)
        {
            clear_pairlist(&cpuLists_[th]);
        }
        else
        {
            clear_pairlist(&gpuLists_[th]);
        }

        if (params_.haveFep)
        {
            clear_pairlist_fep(fepLists_[th].get());
        }
    }

    const gmx_domdec_zones_t& ddZones = *gridSet.domainSetup().zones;

    const auto iZoneRange = getIZoneRange(gridSet.domainSetup(), locality_);

    for (const int iZone : iZoneRange)
    {
        const Grid& iGrid = gridSet.grids()[iZone];

        const auto jZoneRange = getJZoneRange(ddZones, locality_, iZone);

        for (int jZone : jZoneRange)
        {
            const Grid& jGrid = gridSet.grids()[jZone];

            if (debug)
            {
                fprintf(debug, "ns search grid %d vs %d\n", iZone, jZone);
            }

            searchCycleCounting->start(enbsCCsearch);

            ci_block = get_ci_block_size(iGrid, gridSet.domainSetup().haveMultipleDomains, numLists);

            /* With GPU: generate progressively smaller lists for
             * load balancing for local only or non-local with 2 zones.
             */
            progBal = (locality_ == InteractionLocality::Local || ddZones.n <= 2);

#pragma omp parallel for num_threads(numLists) schedule(static)
            for (int th = 0; th < numLists; th++)
            {
                try
                {
                    /* Re-init the thread-local work flag data before making
                     * the first list (not an elegant conditional).
                     */
                    if (nbat->bUseBufferFlags && (iZone == 0 && jZone == 0))
                    {
                        init_buffer_flags(&searchWork[th].buffer_flags, nbat->numAtoms());
                    }

                    if (combineLists_ && th > 0)
                    {
                        GMX_ASSERT(!isCpuType_, "Can only combine GPU lists");

                        clear_pairlist(&gpuLists_[th]);
                    }

                    PairsearchWork& work = searchWork[th];

                    work.cycleCounter.start();

                    t_nblist* fepListPtr = (fepLists_.empty() ? nullptr : fepLists_[th].get());

                    /* Divide the i cells equally over the pairlists */
                    if (isCpuType_)
                    {
                        nbnxn_make_pairlist_part(gridSet, iGrid, jGrid, &work, nbat, *excl, rlist,
                                                 params_.pairlistType, ci_block, nbat->bUseBufferFlags,
                                                 nsubpair_target, progBal, nsubpair_tot_est, th,
                                                 numLists, &cpuLists_[th], fepListPtr);
                    }
                    else
                    {
                        nbnxn_make_pairlist_part(gridSet, iGrid, jGrid, &work, nbat, *excl, rlist,
                                                 params_.pairlistType, ci_block, nbat->bUseBufferFlags,
                                                 nsubpair_target, progBal, nsubpair_tot_est, th,
                                                 numLists, &gpuLists_[th], fepListPtr);
                    }

                    work.cycleCounter.stop();
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }
            searchCycleCounting->stop(enbsCCsearch);

            np_tot = 0;
            np_noq = 0;
            np_hlj = 0;
            for (int th = 0; th < numLists; th++)
            {
                inc_nrnb(nrnb, eNR_NBNXN_DIST2, searchWork[th].ndistc);

                if (isCpuType_)
                {
                    const NbnxnPairlistCpu& nbl = cpuLists_[th];
                    np_tot += nbl.cj.size();
                    np_noq += nbl.work->ncj_noq;
                    np_hlj += nbl.work->ncj_hlj;
                }
                else
                {
                    const NbnxnPairlistGpu& nbl = gpuLists_[th];
                    /* This count ignores potential subsequent pair pruning */
                    np_tot += nbl.nci_tot;
                }
            }
            if (isCpuType_)
            {
                nap = cpuLists_[0].na_ci * cpuLists_[0].na_cj;
            }
            else
            {
                nap = gmx::square(gpuLists_[0].na_ci);
            }
            natpair_ljq_ = (np_tot - np_noq) * nap - np_hlj * nap / 2;
            natpair_lj_  = np_noq * nap;
            natpair_q_   = np_hlj * nap / 2;

            if (combineLists_ && numLists > 1)
            {
                GMX_ASSERT(!isCpuType_, "Can only combine GPU lists");

                searchCycleCounting->start(enbsCCcombine);

                combine_nblists(gmx::constArrayRefFromArray(&gpuLists_[1], numLists - 1), &gpuLists_[0]);

                searchCycleCounting->stop(enbsCCcombine);
            }
        }
    }

    if (isCpuType_)
    {
        if (numLists > 1 && checkRebalanceSimpleLists(cpuLists_))
        {
            rebalanceSimpleLists(cpuLists_, cpuListsWork_, searchWork);

            /* Swap the sets of pair lists */
            cpuLists_.swap(cpuListsWork_);
        }
    }
    else
    {
        /* Sort the entries on size, large ones first */
        if (combineLists_ || gpuLists_.size() == 1)
        {
            sort_sci(&gpuLists_[0]);
        }
        else
        {
#pragma omp parallel for num_threads(numLists) schedule(static)
            for (int th = 0; th < numLists; th++)
            {
                try
                {
                    sort_sci(&gpuLists_[th]);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }
        }
    }

    if (nbat->bUseBufferFlags)
    {
        reduce_buffer_flags(searchWork, numLists, &nbat->buffer_flags);
    }

    if (gridSet.haveFep())
    {
        /* Balance the free-energy lists over all the threads */
        balance_fep_lists(fepLists_, searchWork);
    }

    if (isCpuType_)
    {
        /* This is a fresh list, so not pruned, stored using ci.
         * ciOuter is invalid at this point.
         */
        GMX_ASSERT(cpuLists_[0].ciOuter.empty(), "ciOuter is invalid so it should be empty");
    }

    /* If we have more than one list, they either got rebalancing (CPU)
     * or combined (GPU), so we should dump the final result to debug.
     */
    if (debug)
    {
        if (isCpuType_ && cpuLists_.size() > 1)
        {
            for (auto& cpuList : cpuLists_)
            {
                print_nblist_statistics(debug, cpuList, gridSet, rlist);
            }
        }
        else if (!isCpuType_ && gpuLists_.size() > 1)
        {
            print_nblist_statistics(debug, gpuLists_[0], gridSet, rlist);
        }
    }

    if (debug)
    {
        if (gmx_debug_at)
        {
            if (isCpuType_)
            {
                for (auto& cpuList : cpuLists_)
                {
                    print_nblist_ci_cj(debug, cpuList);
                }
            }
            else
            {
                print_nblist_sci_cj(debug, gpuLists_[0]);
            }
        }

        if (nbat->bUseBufferFlags)
        {
            print_reduction_cost(&nbat->buffer_flags, numLists);
        }
    }

    if (params_.useDynamicPruning && isCpuType_)
    {
        prepareListsForDynamicPruning(cpuLists_);
    }
}

void PairlistSets::construct(const InteractionLocality iLocality,
                             PairSearch*               pairSearch,
                             nbnxn_atomdata_t*         nbat,
                             const t_blocka*           excl,
                             const int64_t             step,
                             t_nrnb*                   nrnb)
{
    pairlistSet(iLocality).constructPairlists(pairSearch->gridSet(), pairSearch->work(), nbat, excl,
                                              minimumIlistCountForGpuBalancing_, nrnb,
                                              &pairSearch->cycleCounting_);

    if (iLocality == InteractionLocality::Local)
    {
        outerListCreationStep_ = step;
    }
    else
    {
        GMX_RELEASE_ASSERT(outerListCreationStep_ == step,
                           "Outer list should be created at the same step as the inner list");
    }

    /* Special performance logging stuff (env.var. GMX_NBNXN_CYCLE) */
    if (iLocality == InteractionLocality::Local)
    {
        pairSearch->cycleCounting_.searchCount_++;
    }
    if (pairSearch->cycleCounting_.recordCycles_
        && (!pairSearch->gridSet().domainSetup().haveMultipleDomains || iLocality == InteractionLocality::NonLocal)
        && pairSearch->cycleCounting_.searchCount_ % 100 == 0)
    {
        pairSearch->cycleCounting_.printCycles(stderr, pairSearch->work());
    }
}

void nonbonded_verlet_t::constructPairlist(const InteractionLocality iLocality,
                                           const t_blocka*           excl,
                                           int64_t                   step,
                                           t_nrnb*                   nrnb)
{
    pairlistSets_->construct(iLocality, pairSearch_.get(), nbat.get(), excl, step, nrnb);

    if (useGpu())
    {
        /* Launch the transfer of the pairlist to the GPU.
         *
         * NOTE: The launch overhead is currently not timed separately
         */
        Nbnxm::gpu_init_pairlist(gpu_nbv, pairlistSets().pairlistSet(iLocality).gpuList(), iLocality);
    }
}

static void prepareListsForDynamicPruning(gmx::ArrayRef<NbnxnPairlistCpu> lists)
{
    /* TODO: Restructure the lists so we have actual outer and inner
     *       list objects so we can set a single pointer instead of
     *       swapping several pointers.
     */

    for (auto& list : lists)
    {
        /* The search produced a list in ci/cj.
         * Swap the list pointers so we get the outer list is ciOuter,cjOuter
         * and we can prune that to get an inner list in ci/cj.
         */
        GMX_RELEASE_ASSERT(list.ciOuter.empty() && list.cjOuter.empty(),
                           "The outer lists should be empty before preparation");

        std::swap(list.ci, list.ciOuter);
        std::swap(list.cj, list.cjOuter);
    }
}
