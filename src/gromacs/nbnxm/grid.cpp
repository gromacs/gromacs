/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * \brief
 * Implements the Grid class.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "grid.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/allocator.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/vec.h"

#include "boundingbox.h"
#include "boundingbox_simd.h"
#include "gridsetdata.h"
#include "nbnxm.h"
#include "nbnxm_geometry.h"
#include "pairlist.h"
#include "pairlistparams.h"

namespace gmx
{

Grid::Geometry::Geometry(const PairlistType pairlistType) :
    isSimple_(pairlistType != PairlistType::Hierarchical8x8x8),
    numAtomsICluster_(IClusterSizePerListType[pairlistType]),
    numAtomsJCluster_(JClusterSizePerListType[pairlistType]),
    numAtomsPerCell_((isSimple_ ? 1 : sc_gpuNumClusterPerCell(pairlistType)) * numAtomsICluster_),
    numAtomsICluster2Log_(get_2log(numAtomsICluster_)),
    pairlistType_(pairlistType)
{
}

namespace
{
/*! \brief Named constant
 *
 * We want the allocator for these two HostVectors to be propagated on
 * copy construction of Grid during the construction of the vector of
 * Grids in GridSet. */
constexpr bool allocatorShouldPropagateDuringCopyConstruction = true;
} // namespace

Grid::Grid(const PairlistType pairlistType, const int ddZone, const bool& haveFep, PinningPolicy pinningPolicy) :
    geometry_(pairlistType),
    ddZone_(ddZone),
    cxy_na_(gmx::HostAllocationPolicy(pinningPolicy, allocatorShouldPropagateDuringCopyConstruction)),
    cxy_ind_(gmx::HostAllocationPolicy(pinningPolicy, allocatorShouldPropagateDuringCopyConstruction)),
    haveFep_(haveFep)
{
}

RVec Grid::averageCellSize() const
{
    const PairlistType    layoutType      = geometry().pairlistType_;
    const int             numAtomsCluster = geometry().numAtomsICluster_;
    const GridDimensions& dims            = dimensions();

    RVec ls;
    ls[XX] = dims.cellSize[XX] / sc_gpuNumClusterPerCellX(layoutType);
    ls[YY] = dims.cellSize[YY] / sc_gpuNumClusterPerCellY(layoutType);
    ls[ZZ] = numAtomsCluster / (dims.atomDensity * ls[XX] * ls[YY]);

    return ls;
}

/*! \brief Returns the atom density (> 0) of a rectangular grid */
static real gridAtomDensity(int numAtoms, const RVec& gridBoundingBoxSize)
{
    if (numAtoms == 0)
    {
        /* To avoid zero density we use a minimum of 1 atom */
        numAtoms = 1;
    }

    return static_cast<real>(numAtoms)
           / (gridBoundingBoxSize[XX] * gridBoundingBoxSize[YY] * gridBoundingBoxSize[ZZ]);
}

//! \brief Get approximate dimensions of each cell. Returns the length along X and Y.
static std::array<real, DIM - 1> getTargetCellLength(const Grid::Geometry& geometry, const real atomDensity)
{
    if (geometry.isSimple_)
    {
        /* To minimize the zero interactions, we should make
         * the largest of the i/j cell cubic.
         */
        int numAtomsInCell = std::max(geometry.numAtomsICluster_, geometry.numAtomsJCluster_);

        /* Approximately cubic cells */
        real tlen = std::cbrt(numAtomsInCell / atomDensity);
        return { tlen, tlen };
    }
    else
    {
        /* Approximately cubic sub cells */
        real tlen = std::cbrt(geometry.numAtomsICluster_ / atomDensity);
        return { tlen * sc_gpuNumClusterPerCellX(geometry.pairlistType_),
                 tlen * sc_gpuNumClusterPerCellY(geometry.pairlistType_) };
    }
}

static int getMaxNumCells(const Grid::Geometry& geometry, const int numAtoms, const int numColumns)
{
    if (geometry.numAtomsJCluster_ <= geometry.numAtomsICluster_)
    {
        return numAtoms / geometry.numAtomsPerCell_ + numColumns;
    }
    else
    {
        return numAtoms / geometry.numAtomsPerCell_
               + numColumns * geometry.numAtomsJCluster_ / geometry.numAtomsICluster_;
    }
}

//! Returns whether \p ddZone is the home zone
static bool isHomeZone(const int ddZone)
{
    return ddZone == 0;
}

void Grid::resizeBoundingBoxesAndFlags(const int maxNumCells)
{
    if (!geometry_.isSimple_)
    {
        numClusters_.resize(maxNumCells);
    }
    bbcz_.resize(maxNumCells);

    /* This resize also zeros the contents, this avoid possible
     * floating exceptions in SIMD with the unused bb elements.
     */
    if (geometry_.isSimple_)
    {
        bb_.resize(maxNumCells);
    }
    else
    {
        if (sc_boundingBoxCornersAsQuadruplets(geometry_.pairlistType_))
        {
            pbb_.resize(packedBoundingBoxesIndex(
                    maxNumCells * sc_gpuNumClusterPerCell(geometry_.pairlistType_)));
        }
        else
        {
            bb_.resize(maxNumCells * sc_gpuNumClusterPerCell(geometry_.pairlistType_));
        }
    }

    if (geometry_.numAtomsJCluster_ == geometry_.numAtomsICluster_)
    {
        bbj_ = bb_;
    }
    else
    {
        GMX_ASSERT(geometry_.isSimple_, "Only CPU lists should have different i/j cluster sizes");

        bbjStorage_.resize(maxNumCells * geometry_.numAtomsICluster_ / geometry_.numAtomsJCluster_);
        bbj_ = bbjStorage_;
    }

    flags_.resize(maxNumCells);
    if (haveFep_)
    {
        fep_.resize(maxNumCells * geometry_.numAtomsPerCell_ / geometry_.numAtomsICluster_);
    }
}

void Grid::setDimensions(const int   ddZone,
                         const int   numAtomsTotal,
                         const int   numAtomsWithoutFillers,
                         const RVec& lowerCorner,
                         const RVec& upperCorner,
                         real*       atomDensity,
                         const real  maxAtomGroupRadius)
{
    dimensions_.lowerCorner = lowerCorner;
    dimensions_.upperCorner = upperCorner;

    /* We allow passing lowerCorner=upperCorner, in which case we need to
     * create a finite sized bounding box to avoid division by zero.
     * We use a minimum size such that the volume fits in float with some
     * margin for computing and using the atom number density.
     */
    constexpr real c_minimumGridSize = 1e-10;
    for (int d = 0; d < DIM; d++)
    {
        GMX_ASSERT(upperCorner[d] >= lowerCorner[d],
                   "Upper corner should be larger than the lower corner");
        if (dimensions_.upperCorner[d] - dimensions_.lowerCorner[d] < c_minimumGridSize)
        {
            /* Ensure we apply a correction to the bounding box */
            real correction =
                    std::max(std::abs(lowerCorner[d]) * GMX_REAL_EPS, 0.5_real * c_minimumGridSize);
            dimensions_.lowerCorner[d] -= correction;
            dimensions_.upperCorner[d] += correction;
        }
    }

    const RVec gridBoundingBoxSize = dimensions_.upperCorner - dimensions_.lowerCorner;

    dimensions_.gridSize = gridBoundingBoxSize;

    /* For the home zone we compute the density when not set (=-1) or when =0 */
    GMX_ASSERT(atomDensity, "atomDensity cannot be nullptr");
    if (isHomeZone(ddZone) && *atomDensity <= 0)
    {
        *atomDensity = gridAtomDensity(numAtomsWithoutFillers, gridBoundingBoxSize);
    }

    dimensions_.atomDensity        = *atomDensity;
    dimensions_.maxAtomGroupRadius = maxAtomGroupRadius;

    if (numAtomsWithoutFillers > 0)
    {
        GMX_ASSERT(*atomDensity > 0, "With one or more atoms, the density should be positive");

        /* target cell length */
        const std::array<real, DIM - 1> tlen = getTargetCellLength(geometry_, *atomDensity);

        /* We round ncx and ncy down, because we get less cell pairs
         * in the pairlist when the fixed cell dimensions (x,y) are
         * larger than the variable one (z) than the other way around.
         */
        dimensions_.numCells[XX] = std::max(1, static_cast<int>(gridBoundingBoxSize[XX] / tlen[XX]));
        dimensions_.numCells[YY] = std::max(1, static_cast<int>(gridBoundingBoxSize[YY] / tlen[YY]));
    }
    else
    {
        dimensions_.numCells[XX] = 1;
        dimensions_.numCells[YY] = 1;
    }

    for (int d = 0; d < DIM - 1; d++)
    {
        dimensions_.cellSize[d]    = gridBoundingBoxSize[d] / dimensions_.numCells[d];
        dimensions_.invCellSize[d] = 1 / dimensions_.cellSize[d];
    }

    if (ddZone > 0)
    {
        /* This is a non-home zone, add an extra row of cells
         * for particles communicated for bonded interactions.
         * These can be beyond the cut-off. It doesn't matter where
         * they end up on the grid, but for performance it's better
         * if they don't end up in cells that can be within cut-off range.
         */
        dimensions_.numCells[XX]++;
        dimensions_.numCells[YY]++;
    }

    /* We need one additional cell entry for particles moved by DD */
    cxy_na_.resize(numColumns() + 1);
    cxy_ind_.resize(numColumns() + 2);

    /* Worst case scenario of 1 atom in each last cell */
    const int maxNumCells = getMaxNumCells(geometry_, numAtomsTotal, numColumns());

    resizeBoundingBoxesAndFlags(maxNumCells);
}

/* We need to sort particles in grid columns on z-coordinate.
 * As particle are very often distributed homogeneously, we use a sorting
 * algorithm similar to pigeonhole sort. We multiply the z-coordinate
 * by a factor, cast to an int and try to store in that hole. If the hole
 * is full, we move this or another particle. A second pass is needed to make
 * contiguous elements. SORT_GRID_OVERSIZE is the ratio of holes to particles.
 * 4 is the optimal value for homogeneous particle distribution and allows
 * for an O(#particles) sort up till distributions were all particles are
 * concentrated in 1/4 of the space. No NlogN fallback is implemented,
 * as it can be expensive to detect inhomogeneous particle distributions.
 */
/*! \brief Ratio of grid cells to atoms */
static constexpr int c_sortGridRatio = 4;
/*! \brief Maximum ratio of holes used, in the worst case all particles
 * end up in the last hole and we need num. atoms extra holes at the end.
 */
static constexpr int c_sortGridMaxSizeFactor = c_sortGridRatio + 1;

/*! \brief Sorts particle index a on coordinates x along dim.
 *
 * Backwards tells if we want decreasing iso increasing coordinates.
 * h0 is the minimum of the coordinate range.
 * invh is the 1/length of the sorting range.
 * n_per_h (>=n) is the expected average number of particles per 1/invh
 * sort is the sorting work array.
 * sort should have a size of at least n_per_h*c_sortGridRatio + n,
 * or easier, allocate at least n*c_sortGridMaxSizeFactor elements.
 */
static void sort_atoms(int                  dim,
                       gmx_bool             Backwards,
                       int gmx_unused       dd_zone,
                       bool gmx_unused      relevantAtomsAreWithinGridBounds,
                       int*                 a,
                       int                  n,
                       ArrayRef<const RVec> x,
                       real                 h0,
                       real                 invh,
                       int                  n_per_h,
                       ArrayRef<int>        sort)
{
    if (n <= 1)
    {
        /* Nothing to do */
        return;
    }

    GMX_ASSERT(n <= n_per_h, "We require n <= n_per_h");

    /* Transform the inverse range height into the inverse hole height */
    invh *= n_per_h * c_sortGridRatio;

    /* Set nsort to the maximum possible number of holes used.
     * In worst case all n elements end up in the last bin.
     */
    int nsort = n_per_h * c_sortGridRatio + n;

    /* Determine the index range used, so we can limit it for the second pass */
    int zi_min = INT_MAX;
    int zi_max = -1;

    /* Sort the particles using a simple index sort */
    for (int i = 0; i < n; i++)
    {
        /* The cast takes care of float-point rounding effects below zero.
         * This code assumes particles are less than 1/SORT_GRID_OVERSIZE
         * times the box height out of the box.
         */
        int zi = static_cast<int>((x[a[i]][dim] - h0) * invh);

#ifndef NDEBUG
        /* As we can have rounding effect, we use > iso >= here */
        if (relevantAtomsAreWithinGridBounds && (zi < 0 || (dd_zone == 0 && zi > n_per_h * c_sortGridRatio)))
        {
            gmx_fatal(FARGS,
                      "(int)((x[%d][%c]=%f - %f)*%f) = %d, not in 0 - %d*%d\n",
                      a[i],
                      'x' + dim,
                      x[a[i]][dim],
                      h0,
                      invh,
                      zi,
                      n_per_h,
                      c_sortGridRatio);
        }
#endif
        if (zi < 0)
        {
            zi = 0;
        }

        /* In a non-local domain, particles communicated for bonded interactions
         * can be far beyond the grid size, which is set by the non-bonded
         * cut-off distance. We sort such particles into the last cell.
         */
        if (zi > n_per_h * c_sortGridRatio)
        {
            zi = n_per_h * c_sortGridRatio;
        }

        /* Ideally this particle should go in sort cell zi,
         * but that might already be in use,
         * in that case find the first empty cell higher up
         */
        if (sort[zi] < 0)
        {
            sort[zi] = a[i];
            zi_min   = std::min(zi_min, zi);
            zi_max   = std::max(zi_max, zi);
        }
        else
        {
            /* We have multiple atoms in the same sorting slot.
             * Sort on real z for minimal bounding box size.
             * There is an extra check for identical z to ensure
             * well-defined output order, independent of input order
             * to ensure binary reproducibility after restarts.
             */
            while (sort[zi] >= 0
                   && (x[a[i]][dim] > x[sort[zi]][dim]
                       || (x[a[i]][dim] == x[sort[zi]][dim] && a[i] > sort[zi])))
            {
                zi++;
            }

            if (sort[zi] >= 0)
            {
                /* Shift all elements by one slot until we find an empty slot */
                int cp  = sort[zi];
                int zim = zi + 1;
                while (sort[zim] >= 0)
                {
                    int tmp   = sort[zim];
                    sort[zim] = cp;
                    cp        = tmp;
                    zim++;
                }
                sort[zim] = cp;
                zi_max    = std::max(zi_max, zim);
            }
            sort[zi] = a[i];
            zi_max   = std::max(zi_max, zi);
        }
    }

    int c = 0;
    if (!Backwards)
    {
        for (int zi = 0; zi < nsort; zi++)
        {
            if (sort[zi] >= 0)
            {
                a[c++]   = sort[zi];
                sort[zi] = -1;
            }
        }
    }
    else
    {
        for (int zi = zi_max; zi >= zi_min; zi--)
        {
            if (sort[zi] >= 0)
            {
                a[c++]   = sort[zi];
                sort[zi] = -1;
            }
        }
    }
    if (c < n)
    {
        gmx_incons("Lost particles while sorting");
    }
}

#if GMX_DOUBLE
//! Returns double up to one least significant float bit smaller than x
static double R2F_D(const float x)
{
    return static_cast<float>(x >= 0 ? (1 - GMX_FLOAT_EPS) * x : (1 + GMX_FLOAT_EPS) * x);
}
//! Returns double up to one least significant float bit larger than x
static double R2F_U(const float x)
{
    return static_cast<float>(x >= 0 ? (1 + GMX_FLOAT_EPS) * x : (1 - GMX_FLOAT_EPS) * x);
}
#else
//! Returns x
static float R2F_D(const float x)
{
    return x;
}
//! Returns x
static float R2F_U(const float x)
{
    return x;
}
#endif

//! Computes the bounding box for na coordinates in order x,y,z, bb order xyz0
static void calc_bounding_box(const int numAtoms, const int stride, const real* x, BoundingBox* bb)
{
    int  i  = 0;
    real xl = x[i + XX];
    real xh = x[i + XX];
    real yl = x[i + YY];
    real yh = x[i + YY];
    real zl = x[i + ZZ];
    real zh = x[i + ZZ];
    i += stride;
    for (int j = 1; j < numAtoms; j++)
    {
        xl = std::min(xl, x[i + XX]);
        xh = std::max(xh, x[i + XX]);
        yl = std::min(yl, x[i + YY]);
        yh = std::max(yh, x[i + YY]);
        zl = std::min(zl, x[i + ZZ]);
        zh = std::max(zh, x[i + ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb->lower.x = R2F_D(xl);
    bb->lower.y = R2F_D(yl);
    bb->lower.z = R2F_D(zl);
    bb->upper.x = R2F_U(xh);
    bb->upper.y = R2F_U(yh);
    bb->upper.z = R2F_U(zh);
}

/*! \brief Computes the bounding box for packed coordinates
 * \tparam     packSize  The pack size for the coordinates, also the number of atoms per cell
 * \param[in]  numAtoms  The actual number of atoms in this cell
 * \param[in]  x         Packed coodinates
 * \param[out] bb        Pointer to the bounding box
 */
template<int packSize>
static void calcBoundingBoxXPacked(const int numAtoms, const real* x, BoundingBox* bb)
{
    real xl = x[XX * packSize];
    real xh = x[XX * packSize];
    real yl = x[YY * packSize];
    real yh = x[YY * packSize];
    real zl = x[ZZ * packSize];
    real zh = x[ZZ * packSize];
    for (int j = 1; j < numAtoms; j++)
    {
        xl = std::min(xl, x[j + XX * packSize]);
        xh = std::max(xh, x[j + XX * packSize]);
        yl = std::min(yl, x[j + YY * packSize]);
        yh = std::max(yh, x[j + YY * packSize]);
        zl = std::min(zl, x[j + ZZ * packSize]);
        zh = std::max(zh, x[j + ZZ * packSize]);
    }
    /* Note: possible double to float conversion here */
    bb->lower.x = R2F_D(xl);
    bb->lower.y = R2F_D(yl);
    bb->lower.z = R2F_D(zl);
    bb->upper.x = R2F_U(xh);
    bb->upper.y = R2F_U(yh);
    bb->upper.z = R2F_U(zh);
}

/*! \brief Computes the whole plus half bounding boxes for packed coordinates
 *
 * \tparam     packSize  The pack size for the coordinates, also the number of atoms per cell
 * \param[in]  numAtoms  The actual number of atoms in this cell
 * \param[in]  x         Packed coodinates
 * \param[out] bb        Pointer to the bounding box for the whole cell
 * \param[out] bbj       Pointer to the bounding boxes for the two halves of the cell
 */
template<int packSize>
gmx_unused static void calcBoundingBoxHalves(const int numAtoms, const real* x, BoundingBox* bb, BoundingBox* bbj)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    constexpr int halfPackSize = packSize / 2;

    calcBoundingBoxXPacked<packSize>(std::min(numAtoms, halfPackSize), x, bbj);

    if (numAtoms > halfPackSize)
    {
        calcBoundingBoxXPacked<packSize>(
                std::min(numAtoms - halfPackSize, halfPackSize), x + halfPackSize, bbj + 1);
    }
    else
    {
        /* Set the "empty" bounding box to the same as the first one,
         * so we don't need to treat special cases in the rest of the code.
         */
#if NBNXN_SEARCH_BB_SIMD4
        store4(bbj[1].lower.ptr(), load4(bbj[0].lower.ptr()));
        store4(bbj[1].upper.ptr(), load4(bbj[0].upper.ptr()));
#else
        bbj[1] = bbj[0];
#endif
    }

#if NBNXN_SEARCH_BB_SIMD4
    store4(bb->lower.ptr(), min(load4(bbj[0].lower.ptr()), load4(bbj[1].lower.ptr())));
    store4(bb->upper.ptr(), max(load4(bbj[0].upper.ptr()), load4(bbj[1].upper.ptr())));
#else
    {
        bb->lower = BoundingBox::Corner::min(bbj[0].lower, bbj[1].lower);
        bb->upper = BoundingBox::Corner::max(bbj[0].upper, bbj[1].upper);
    }
#endif
}

/*! \brief Computes the bounding box for na coordinates in order xyz, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx(int na, int stride, const real* x, float* bb)
{
    int  i  = 0;
    real xl = x[i + XX];
    real xh = x[i + XX];
    real yl = x[i + YY];
    real yh = x[i + YY];
    real zl = x[i + ZZ];
    real zh = x[i + ZZ];
    i += stride;
    for (int j = 1; j < na; j++)
    {
        xl = std::min(xl, x[i + XX]);
        xh = std::max(xh, x[i + XX]);
        yl = std::min(yl, x[i + YY]);
        yh = std::max(yh, x[i + YY]);
        zl = std::min(zl, x[i + ZZ]);
        zh = std::max(zh, x[i + ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb[0 * c_packedBoundingBoxesDimSize] = R2F_D(xl);
    bb[1 * c_packedBoundingBoxesDimSize] = R2F_D(yl);
    bb[2 * c_packedBoundingBoxesDimSize] = R2F_D(zl);
    bb[3 * c_packedBoundingBoxesDimSize] = R2F_U(xh);
    bb[4 * c_packedBoundingBoxesDimSize] = R2F_U(yh);
    bb[5 * c_packedBoundingBoxesDimSize] = R2F_U(zh);
}

#if NBNXN_SEARCH_SIMD4_FLOAT_X_BB

/*! \brief Computes the bounding box for na coordinates in order xyz?, bb order xyz0 */
static void calc_bounding_box_simd4(int na, const float* x, BoundingBox* bb)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    static_assert(sizeof(BoundingBox::Corner) == GMX_SIMD4_WIDTH * sizeof(float),
                  "The Corner struct should hold exactly 4 floats");

    Simd4Float bb_0_S = load4(x);
    Simd4Float bb_1_S = bb_0_S;

    for (int i = 1; i < na; i++)
    {
        Simd4Float x_S = load4(x + i * GMX_SIMD4_WIDTH);
        bb_0_S         = min(bb_0_S, x_S);
        bb_1_S         = max(bb_1_S, x_S);
    }

    store4(bb->lower.ptr(), bb_0_S);
    store4(bb->upper.ptr(), bb_1_S);
}

/*! \brief Computes the bounding box for na coordinates in order xyz?, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx_simd4(int na, const float* x, real* bb)
{
    alignas(GMX_SIMD_ALIGNMENT) BoundingBox bbWorkAligned;

    calc_bounding_box_simd4(na, x, &bbWorkAligned);

    bb[0 * c_packedBoundingBoxesDimSize] = bbWorkAligned.lower.x;
    bb[1 * c_packedBoundingBoxesDimSize] = bbWorkAligned.lower.y;
    bb[2 * c_packedBoundingBoxesDimSize] = bbWorkAligned.lower.z;
    bb[3 * c_packedBoundingBoxesDimSize] = bbWorkAligned.upper.x;
    bb[4 * c_packedBoundingBoxesDimSize] = bbWorkAligned.upper.y;
    bb[5 * c_packedBoundingBoxesDimSize] = bbWorkAligned.upper.z;
}

#endif /* NBNXN_SEARCH_SIMD4_FLOAT_X_BB */


/*! \brief Combines pairs of consecutive bounding boxes */
static void combine_bounding_box_pairs(const Grid& grid, ArrayRef<const BoundingBox> bb, ArrayRef<BoundingBox> bbj)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    for (int i = 0; i < grid.numColumns(); i++)
    {
        /* Starting bb in a column is expected to be 2-aligned */
        const int sc2 = grid.firstCellInColumn(i) >> 1;
        /* For odd numbers skip the last bb here */
        const int nc2 = (grid.numAtomsInColumn(i) + 3) >> (2 + 1);
        for (int c2 = sc2; c2 < sc2 + nc2; c2++)
        {
#if NBNXN_SEARCH_BB_SIMD4
            Simd4Float min_S, max_S;

            min_S = min(load4(bb[c2 * 2 + 0].lower.ptr()), load4(bb[c2 * 2 + 1].lower.ptr()));
            max_S = max(load4(bb[c2 * 2 + 0].upper.ptr()), load4(bb[c2 * 2 + 1].upper.ptr()));
            store4(bbj[c2].lower.ptr(), min_S);
            store4(bbj[c2].upper.ptr(), max_S);
#else
            bbj[c2].lower = BoundingBox::Corner::min(bb[c2 * 2 + 0].lower, bb[c2 * 2 + 1].lower);
            bbj[c2].upper = BoundingBox::Corner::max(bb[c2 * 2 + 0].upper, bb[c2 * 2 + 1].upper);
#endif
        }
        if (((grid.numAtomsInColumn(i) + 3) >> 2) & 1)
        {
            /* The bb count in this column is odd: duplicate the last bb */
            int c2        = sc2 + nc2;
            bbj[c2].lower = bb[c2 * 2].lower;
            bbj[c2].upper = bb[c2 * 2].upper;
        }
    }
}


/*! \brief Prints the average bb size, used for debug output */
static void print_bbsizes_simple(FILE* fp, const Grid& grid)
{
    dvec       ba             = { 0 };
    const auto iBoundingBoxes = grid.iBoundingBoxes();
    for (int c = 0; c < grid.numCells(); c++)
    {
        const BoundingBox& bb = iBoundingBoxes[c];
        ba[XX] += bb.upper.x - bb.lower.x;
        ba[YY] += bb.upper.y - bb.lower.y;
        ba[ZZ] += bb.upper.z - bb.lower.z;
    }
    if (grid.numCells() > 0)
    {
        dsvmul(1.0 / grid.numCells(), ba, ba);
    }

    const GridDimensions& dims = grid.dimensions();
    real                  avgCellSizeZ =
            (dims.atomDensity > 0 ? grid.geometry().numAtomsICluster_
                                            / (dims.atomDensity * dims.cellSize[XX] * dims.cellSize[YY])
                                  : 0.0);

    fprintf(fp,
            "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            dims.cellSize[XX],
            dims.cellSize[YY],
            avgCellSizeZ,
            ba[XX],
            ba[YY],
            ba[ZZ],
            ba[XX] * dims.invCellSize[XX],
            ba[YY] * dims.invCellSize[YY],
            dims.atomDensity > 0 ? ba[ZZ] / avgCellSizeZ : 0.0);
}

/*! \brief Prints the average bb size, used for debug output */
static void print_bbsizes_supersub(FILE* fp, const Grid& grid)
{
    dvec ba;

    clear_dvec(ba);
    int        ns         = 0;
    const auto layoutType = grid.geometry().pairlistType_;
    for (int c = 0; c < grid.numCells(); c++)
    {
        if (sc_boundingBoxCornersAsQuadruplets(layoutType))
        {
            for (int s = 0; s < grid.numClustersPerCell()[c]; s += c_packedBoundingBoxesDimSize)
            {
                int cs_w = (c * sc_gpuNumClusterPerCell(layoutType) + s) / c_packedBoundingBoxesDimSize;
                auto boundingBoxes = grid.packedBoundingBoxes().subArray(
                        cs_w * c_packedBoundingBoxesSize, c_packedBoundingBoxesSize);
                for (int i = 0; i < c_packedBoundingBoxesDimSize; i++)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        ba[d] += boundingBoxes[(DIM + d) * c_packedBoundingBoxesDimSize + i]
                                 - boundingBoxes[(0 + d) * c_packedBoundingBoxesDimSize + i];
                    }
                }
            }
        }
        else
        {
            for (int s = 0; s < grid.numClustersPerCell()[c]; s++)
            {
                const int index = c * sc_gpuNumClusterPerCell(layoutType) + s;
                const ArrayRef<const BoundingBox> iboundingBoxes = grid.iBoundingBoxes();
                const auto&                       bb             = iboundingBoxes[index];
                ba[XX] += bb.upper.x - bb.lower.x;
                ba[YY] += bb.upper.y - bb.lower.y;
                ba[ZZ] += bb.upper.z - bb.lower.z;
            }
        }
        ns += grid.numClustersPerCell()[c];
    }
    dsvmul(1.0 / ns, ba, ba);

    const GridDimensions& dims = grid.dimensions();
    const real            avgClusterSizeZ =
            (dims.atomDensity > 0 ? grid.geometry().numAtomsPerCell_
                                            / (dims.atomDensity * dims.cellSize[XX] * dims.cellSize[YY]
                                               * sc_gpuNumClusterPerCellZ(layoutType))
                                  : 0.0);

    fprintf(fp,
            "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            dims.cellSize[XX] / sc_gpuNumClusterPerCellX(layoutType),
            dims.cellSize[YY] / sc_gpuNumClusterPerCellY(layoutType),
            avgClusterSizeZ,
            ba[XX],
            ba[YY],
            ba[ZZ],
            ba[XX] * sc_gpuNumClusterPerCellX(layoutType) * dims.invCellSize[XX],
            ba[YY] * sc_gpuNumClusterPerCellY(layoutType) * dims.invCellSize[YY],
            dims.atomDensity > 0 ? ba[ZZ] / avgClusterSizeZ : 0.0);
}

/*!\brief Set non-bonded interaction flags for the current cluster.
 *
 * Sorts atoms on LJ coefficients: !=0 first, ==0 at the end.
 */
static void sort_cluster_on_flag(int                     numAtomsInCluster,
                                 int                     atomStart,
                                 int                     atomEnd,
                                 ArrayRef<const int32_t> atomInfo,
                                 ArrayRef<int>           order,
                                 int*                    flags)
{
    constexpr int c_maxNumAtomsInCluster = 8;
    int           sort1[c_maxNumAtomsInCluster];
    int           sort2[c_maxNumAtomsInCluster];

    GMX_ASSERT(numAtomsInCluster <= c_maxNumAtomsInCluster,
               "Need to increase c_maxNumAtomsInCluster to support larger clusters");

    *flags = 0;

    int subc = 0;
    for (int s = atomStart; s < atomEnd; s += numAtomsInCluster)
    {
        /* Make lists for this (sub-)cell on atoms with and without LJ */
        int  n1       = 0;
        int  n2       = 0;
        bool haveQ    = false;
        int  a_lj_max = -1;
        for (int a = s; a < std::min(s + numAtomsInCluster, atomEnd); a++)
        {
            haveQ = haveQ || bool(atomInfo[order[a]] & sc_atomInfo_HasCharge);

            if (atomInfo[order[a]] & (sc_atomInfo_HasVdw))
            {
                sort1[n1++] = order[a];
                a_lj_max    = a;
            }
            else
            {
                sort2[n2++] = order[a];
            }
        }

        /* If we don't have atoms with LJ, there's nothing to sort */
        if (n1 > 0)
        {
            *flags |= NBNXN_CI_DO_LJ(subc);

            if (2 * n1 <= numAtomsInCluster)
            {
                /* Only sort when strictly necessary. Ordering particles
                 * Ordering particles can lead to less accurate summation
                 * due to rounding, both for LJ and Coulomb interactions.
                 */
                if (2 * (a_lj_max - s) >= numAtomsInCluster)
                {
                    for (int i = 0; i < n1; i++)
                    {
                        order[atomStart + i] = sort1[i];
                    }
                    for (int j = 0; j < n2; j++)
                    {
                        order[atomStart + n1 + j] = sort2[j];
                    }
                }

                *flags |= NBNXN_CI_HALF_LJ(subc);
            }
        }
        if (haveQ)
        {
            *flags |= NBNXN_CI_DO_COUL(subc);
        }
        subc++;
    }
}

/*! \brief Fill a pair search cell with atoms.
 *
 * Potentially sorts atoms and sets the interaction flags.
 */
void Grid::fillCell(GridSetData*            gridSetData,
                    nbnxn_atomdata_t*       nbat,
                    int                     atomStart,
                    int                     atomEnd,
                    ArrayRef<const int32_t> atomInfo,
                    ArrayRef<const RVec>    x)
{
    const int numAtoms = atomEnd - atomStart;

    const ArrayRef<int>& cells       = gridSetData->cells;
    const ArrayRef<int>& atomIndices = gridSetData->atomIndices;

    if (geometry_.isSimple_)
    {
        /* Note that non-local grids are already sorted.
         * Then sort_cluster_on_flag will only set the flags and the sorting
         * will not affect the atom order.
         */
        sort_cluster_on_flag(geometry_.numAtomsICluster_,
                             atomStart,
                             atomEnd,
                             atomInfo,
                             atomIndices,
                             flags_.data() + atomToCluster(atomStart) - cellOffset_);
    }

    if (haveFep_)
    {
        const auto layoutType = geometry_.pairlistType_;
        /* Set the fep flag for perturbed atoms in this (sub-)cell */

        /* The grid-local cluster/(sub-)cell index */
        int cell = atomToCluster(atomStart)
                   - cellOffset_ * (geometry_.isSimple_ ? 1 : sc_gpuNumClusterPerCell(layoutType));
        fep_[cell] = 0;
        for (int at = atomStart; at < atomEnd; at++)
        {
            if (atomIndices[at] >= 0 && (atomInfo[atomIndices[at]] & sc_atomInfo_FreeEnergyPerturbation))
            {
                fep_[cell] |= (1 << (at - atomStart));
            }
        }
    }

    /* Now we have sorted the atoms, set the cell indices */
    for (int at = atomStart; at < atomEnd; at++)
    {
        cells[atomIndices[at]] = at;
    }

    copy_rvec_to_nbat_real(atomIndices.data() + atomStart,
                           numAtoms,
                           geometry_.numAtomsICluster_,
                           as_rvec_array(x.data()),
                           nbat->XFormat,
                           nbat->x().data(),
                           atomStart);

    if (nbat->XFormat == nbatX4)
    {
        /* Store the bounding boxes as xyz.xyz. */
        size_t       offset = atomToCluster(atomStart - cellOffset_ * geometry_.numAtomsICluster_);
        BoundingBox* bb_ptr = bb_.data() + offset;

#if GMX_SIMD
        if (2 * geometry_.numAtomsJCluster_ == geometry_.numAtomsICluster_)
        {
            calcBoundingBoxHalves<c_packX4>(numAtoms,
                                            nbat->x().data() + atom_to_x_index<c_packX4>(atomStart),
                                            bb_ptr,
                                            bbj_.data() + offset * 2);
        }
        else
#endif
        {
            calcBoundingBoxXPacked<c_packX4>(
                    numAtoms, nbat->x().data() + atom_to_x_index<c_packX4>(atomStart), bb_ptr);
        }
    }
    else if (nbat->XFormat == nbatX8)
    {
        /* Store the bounding boxes as xyz.xyz. */
        size_t       offset = atomToCluster(atomStart - cellOffset_ * geometry_.numAtomsICluster_);
        BoundingBox* bb_ptr = bb_.data() + offset;

        calcBoundingBoxXPacked<c_packX8>(
                numAtoms, nbat->x().data() + atom_to_x_index<c_packX8>(atomStart), bb_ptr);
    }
    else if (!geometry_.isSimple_ && sc_boundingBoxCornersAsQuadruplets(geometry_.pairlistType_))
    {
        /* Store the bounding boxes in a format convenient
         * for SIMD4 calculations: xxxxyyyyzzzz...
         */
        const int clusterIndex = ((atomStart - cellOffset_ * geometry_.numAtomsPerCell_)
                                  >> geometry_.numAtomsICluster2Log_);
        float*    pbb_ptr      = pbb_.data() + packedBoundingBoxesIndex(clusterIndex)
                         + (clusterIndex & (c_packedBoundingBoxesDimSize - 1));

#if NBNXN_SEARCH_SIMD4_FLOAT_X_BB
        if (nbat->XFormat == nbatXYZQ)
        {
            calc_bounding_box_xxxx_simd4(numAtoms, nbat->x().data() + atomStart * nbat->xstride, pbb_ptr);
        }
        else
#endif
        {
            calc_bounding_box_xxxx(
                    numAtoms, nbat->xstride, nbat->x().data() + atomStart * nbat->xstride, pbb_ptr);
        }
        if (gmx_debug_at)
        {
            fprintf(debug,
                    "cell %4d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    atomToCluster(atomStart),
                    pbb_ptr[0 * c_packedBoundingBoxesDimSize],
                    pbb_ptr[3 * c_packedBoundingBoxesDimSize],
                    pbb_ptr[1 * c_packedBoundingBoxesDimSize],
                    pbb_ptr[4 * c_packedBoundingBoxesDimSize],
                    pbb_ptr[2 * c_packedBoundingBoxesDimSize],
                    pbb_ptr[5 * c_packedBoundingBoxesDimSize]);
        }
    }
    else
    {
        /* Store the bounding boxes as xyz.xyz. */
        BoundingBox* bb_ptr =
                bb_.data() + atomToCluster(atomStart - cellOffset_ * geometry_.numAtomsPerCell_);

        calc_bounding_box(numAtoms, nbat->xstride, nbat->x().data() + atomStart * nbat->xstride, bb_ptr);

        if (gmx_debug_at)
        {
            int bbo = atomToCluster(atomStart - cellOffset_ * geometry_.numAtomsPerCell_);
            fprintf(debug,
                    "cell %4d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    atomToCluster(atomStart),
                    bb_[bbo].lower.x,
                    bb_[bbo].lower.y,
                    bb_[bbo].lower.z,
                    bb_[bbo].upper.x,
                    bb_[bbo].upper.y,
                    bb_[bbo].upper.z);
        }
    }
}

void Grid::sortColumnsCpuGeometry(GridSetData*            gridSetData,
                                  int                     dd_zone,
                                  ArrayRef<const int32_t> atomInfo,
                                  ArrayRef<const RVec>    x,
                                  nbnxn_atomdata_t*       nbat,
                                  const Range<int>        columnRange,
                                  ArrayRef<int>           sort_work)
{
    if (debug)
    {
        fprintf(debug,
                "cell_offset %d sorting columns %d - %d\n",
                cellOffset_,
                *columnRange.begin(),
                *columnRange.end());
    }

    const bool relevantAtomsAreWithinGridBounds = (dimensions_.maxAtomGroupRadius == 0);

    const int numAtomsPerCell = geometry_.numAtomsPerCell_;

    /* Sort the atoms within each x,y column in 3 dimensions */
    for (int cxy : columnRange)
    {
        const int numAtoms   = numAtomsInColumn(cxy);
        const int numCellsZ  = cxy_ind_[cxy + 1] - cxy_ind_[cxy];
        const int atomOffset = firstAtomInColumn(cxy);

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,
                   FALSE,
                   dd_zone,
                   relevantAtomsAreWithinGridBounds,
                   gridSetData->atomIndices.data() + atomOffset,
                   numAtoms,
                   x,
                   dimensions_.lowerCorner[ZZ],
                   1.0 / dimensions_.gridSize[ZZ],
                   numCellsZ * numAtomsPerCell,
                   sort_work);

        /* Fill the ncz cells in this column */
        const int firstCell  = firstCellInColumn(cxy);
        int       cellFilled = firstCell;
        for (int cellZ = 0; cellZ < numCellsZ; cellZ++)
        {
            const int cell = firstCell + cellZ;

            const int atomOffsetCell       = atomOffset + cellZ * numAtomsPerCell;
            const int numAtomsLeftInColumn = std::max(numAtoms - (atomOffsetCell - atomOffset), 0);
            const int numAtomsCell         = std::min(numAtomsPerCell, numAtomsLeftInColumn);

            fillCell(gridSetData, nbat, atomOffsetCell, atomOffsetCell + numAtomsCell, atomInfo, x);

            /* This copy to bbcz is not really necessary.
             * But it allows to use the same grid search code
             * for the simple and supersub cell setups.
             */
            if (numAtomsCell > 0)
            {
                cellFilled = cell;
            }
            bbcz_[cell].lower = bb_[cellFilled].lower.z;
            bbcz_[cell].upper = bb_[cellFilled].upper.z;
        }

        /* Set the unused atom indices to -1 */
        GMX_ASSERT(gmx::ssize(gridSetData->atomIndices) >= atomOffset + numCellsZ * numAtomsPerCell,
                   "We need sufficient space in atomIndices");
        for (int ind = numAtoms; ind < numCellsZ * numAtomsPerCell; ind++)
        {
            gridSetData->atomIndices[atomOffset + ind] = -1;
        }
    }
}

/* Spatially sort the atoms within one grid column */
void Grid::sortColumnsGpuGeometry(GridSetData*            gridSetData,
                                  int                     dd_zone,
                                  ArrayRef<const int32_t> atomInfo,
                                  ArrayRef<const RVec>    x,
                                  nbnxn_atomdata_t*       nbat,
                                  const Range<int>        columnRange,
                                  ArrayRef<int>           sort_work)
{
    if (debug)
    {
        fprintf(debug,
                "cell_offset %d sorting columns %d - %d\n",
                cellOffset_,
                *columnRange.begin(),
                *columnRange.end());
    }

    const auto layoutType                       = geometry_.pairlistType_;
    const bool relevantAtomsAreWithinGridBounds = (dimensions_.maxAtomGroupRadius == 0);

    const int numAtomsPerCell = geometry_.numAtomsPerCell_;

    const int subdiv_x = geometry_.numAtomsICluster_;
    const int subdiv_y = sc_gpuNumClusterPerCellX(layoutType) * subdiv_x;
    const int subdiv_z = sc_gpuNumClusterPerCellY(layoutType) * subdiv_y;

    /* Extract the atom index array that will be filled here */
    const ArrayRef<int>& atomIndices = gridSetData->atomIndices;

    /* Sort the atoms within each x,y column in 3 dimensions.
     * Loop over all columns on the x/y grid.
     */
    for (int cxy : columnRange)
    {
        const int gridX = cxy / dimensions_.numCells[YY];
        const int gridY = cxy - gridX * dimensions_.numCells[YY];

        const int numAtomsInColumn = cxy_na_[cxy];
        const int numCellsInColumn = cxy_ind_[cxy + 1] - cxy_ind_[cxy];
        const int atomOffset       = firstAtomInColumn(cxy);

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,
                   FALSE,
                   dd_zone,
                   relevantAtomsAreWithinGridBounds,
                   atomIndices.data() + atomOffset,
                   numAtomsInColumn,
                   x,
                   dimensions_.lowerCorner[ZZ],
                   1.0 / dimensions_.gridSize[ZZ],
                   numCellsInColumn * numAtomsPerCell,
                   sort_work);

        /* This loop goes over the cells and clusters along z at once */
        for (int sub_z = 0; sub_z < numCellsInColumn * sc_gpuNumClusterPerCellZ(layoutType); sub_z++)
        {
            const int atomOffsetZ = atomOffset + sub_z * subdiv_z;
            const int numAtomsZ = std::min(subdiv_z, numAtomsInColumn - (atomOffsetZ - atomOffset));
            int       cz        = -1;
            /* We have already sorted on z */

            if (sub_z % sc_gpuNumClusterPerCellZ(layoutType) == 0)
            {
                cz             = sub_z / sc_gpuNumClusterPerCellZ(layoutType);
                const int cell = cxy_ind_[cxy] + cz;

                /* The number of atoms in this cell/super-cluster */
                const int numAtoms =
                        std::min(numAtomsPerCell, numAtomsInColumn - (atomOffsetZ - atomOffset));

                numClusters_[cell] = std::min(sc_gpuNumClusterPerCell(layoutType),
                                              divideRoundUp(numAtoms, geometry_.numAtomsICluster_));

                /* Store the z-boundaries of the bounding box of the cell */
                bbcz_[cell].lower = x[atomIndices[atomOffsetZ]][ZZ];
                bbcz_[cell].upper = x[atomIndices[atomOffsetZ + numAtoms - 1]][ZZ];
            }

            if (sc_gpuNumClusterPerCellY(layoutType) > 1)
            {
                /* Sort the atoms along y */
                sort_atoms(YY,
                           (sub_z & 1) != 0,
                           dd_zone,
                           relevantAtomsAreWithinGridBounds,
                           atomIndices.data() + atomOffsetZ,
                           numAtomsZ,
                           x,
                           dimensions_.lowerCorner[YY] + gridY * dimensions_.cellSize[YY],
                           dimensions_.invCellSize[YY],
                           subdiv_z,
                           sort_work);
            }

            for (int sub_y = 0; sub_y < sc_gpuNumClusterPerCellY(layoutType); sub_y++)
            {
                const int atomOffsetY = atomOffsetZ + sub_y * subdiv_y;
                const int numAtomsY = std::min(subdiv_y, numAtomsInColumn - (atomOffsetY - atomOffset));

                if (sc_gpuNumClusterPerCellX(layoutType) > 1)
                {
                    /* Sort the atoms along x */
                    sort_atoms(XX,
                               ((cz * sc_gpuNumClusterPerCellY(layoutType) + sub_y) & 1) != 0,
                               dd_zone,
                               relevantAtomsAreWithinGridBounds,
                               atomIndices.data() + atomOffsetY,
                               numAtomsY,
                               x,
                               dimensions_.lowerCorner[XX] + gridX * dimensions_.cellSize[XX],
                               dimensions_.invCellSize[XX],
                               subdiv_y,
                               sort_work);
                }

                for (int sub_x = 0; sub_x < sc_gpuNumClusterPerCellX(layoutType); sub_x++)
                {
                    const int atomOffsetX = atomOffsetY + sub_x * subdiv_x;
                    const int numAtomsX =
                            std::min(subdiv_x, numAtomsInColumn - (atomOffsetX - atomOffset));

                    fillCell(gridSetData, nbat, atomOffsetX, atomOffsetX + numAtomsX, atomInfo, x);
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for (int ind = numAtomsInColumn; ind < numCellsInColumn * numAtomsPerCell; ind++)
        {
            atomIndices[atomOffset + ind] = -1;
        }
    }
}

/*! \brief Sets the cell index in the cell array for atom \p atomIndex and increments the atom count for the grid column */
static void setCellAndAtomCount(ArrayRef<int> cell, int cellIndex, ArrayRef<int> cxy_na, int atomIndex)
{
    cell[atomIndex] = cellIndex;
    cxy_na[cellIndex] += 1;
}

void Grid::calcColumnIndices(const GridDimensions&  gridDims,
                             const UpdateGroupsCog* updateGroupsCog,
                             const Range<int>       atomRange,
                             ArrayRef<const RVec>   x,
                             const int              dd_zone,
                             const int*             move,
                             const int              thread,
                             const int              nthread,
                             ArrayRef<int>          cell,
                             ArrayRef<int>          cxy_na)
{
    const int numColumns = gridDims.numCells[XX] * gridDims.numCells[YY];

    /* We add one extra cell for particles which moved during DD */
    for (int i = 0; i < numColumns + 1; i++)
    {
        cxy_na[i] = 0;
    }

    // Use Index to avoid overflow of int with large atom and thread counts
    const Index rangeSize     = atomRange.size();
    const int   taskAtomStart = *atomRange.begin() + ((thread + 0) * rangeSize) / nthread;
    const int   taskAtomEnd   = *atomRange.begin() + ((thread + 1) * rangeSize) / nthread;

    if (dd_zone == 0)
    {
        /* Home zone */
        for (int i = taskAtomStart; i < taskAtomEnd; i++)
        {
            if (move == nullptr || move[i] >= 0)
            {

                const RVec& coord = (updateGroupsCog ? updateGroupsCog->cogForAtom(i) : x[i]);

                /* We need to be careful with rounding,
                 * particles might be a few bits outside the local zone.
                 * The int cast takes care of the lower bound,
                 * we will explicitly take care of the upper bound.
                 */
                int cx = static_cast<int>((coord[XX] - gridDims.lowerCorner[XX])
                                          * gridDims.invCellSize[XX]);
                int cy = static_cast<int>((coord[YY] - gridDims.lowerCorner[YY])
                                          * gridDims.invCellSize[YY]);

#ifndef NDEBUG
                if (cx < 0 || cx > gridDims.numCells[XX] || cy < 0 || cy > gridDims.numCells[YY])
                {
                    gmx_fatal(FARGS,
                              "grid cell cx %d cy %d out of range (max %d %d)\n"
                              "atom %f %f %f, grid->c0 %f %f",
                              cx,
                              cy,
                              gridDims.numCells[XX],
                              gridDims.numCells[YY],
                              x[i][XX],
                              x[i][YY],
                              x[i][ZZ],
                              gridDims.lowerCorner[XX],
                              gridDims.lowerCorner[YY]);
                }
#endif
                /* Take care of potential rounding issues */
                cx = std::min(cx, gridDims.numCells[XX] - 1);
                cy = std::min(cy, gridDims.numCells[YY] - 1);

                /* For the moment cell will contain only the, grid local,
                 * x and y indices, not z.
                 */
                setCellAndAtomCount(cell, cx * gridDims.numCells[YY] + cy, cxy_na, i);
            }
            else
            {
                /* Put this moved particle after the end of the grid,
                 * so we can process it later without using conditionals.
                 */
                setCellAndAtomCount(cell, numColumns, cxy_na, i);
            }
        }
    }
    else
    {
        /* Non-home zone */
        for (int i = taskAtomStart; i < taskAtomEnd; i++)
        {
            int cx = static_cast<int>((x[i][XX] - gridDims.lowerCorner[XX]) * gridDims.invCellSize[XX]);
            int cy = static_cast<int>((x[i][YY] - gridDims.lowerCorner[YY]) * gridDims.invCellSize[YY]);

            /* For non-home zones there could be particles outside
             * the non-bonded cut-off range, which have been communicated
             * for bonded interactions only. For the result it doesn't
             * matter where these end up on the grid. For performance
             * we put them in an extra row at the border.
             */
            cx = std::max(cx, 0);
            cx = std::min(cx, gridDims.numCells[XX] - 1);
            cy = std::max(cy, 0);
            cy = std::min(cy, gridDims.numCells[YY] - 1);

            /* For the moment cell will contain only the, grid local,
             * x and y indices, not z.
             */
            setCellAndAtomCount(cell, cx * gridDims.numCells[YY] + cy, cxy_na, i);
        }
    }
}

/*! \brief Resizes grid and atom data which depend on the number of cells */
static void resizeForNumberOfCells(const int         numNbnxnAtoms,
                                   const int         numAtomsMoved,
                                   const int         ddZone,
                                   GridSetData*      gridSetData,
                                   nbnxn_atomdata_t* nbat)
{
    /* Note: gridSetData->cellIndices was already resized before */

    /* To avoid conditionals we store the moved particles at the end of a,
     * make sure we have enough space.
     */
    gridSetData->atomIndices.resize(numNbnxnAtoms + numAtomsMoved);

    /* Make space in nbat for storing the atom coordinates */
    nbat->resizeCoordinateBuffer(numNbnxnAtoms, ddZone);
}

void Grid::setCellIndices(int                  ddZone,
                          int                  cellOffset,
                          GridSetData*         gridSetData,
                          ArrayRef<GridWork>   gridWork,
                          const Range<int>     atomRange,
                          ArrayRef<const int>  atomInfo,
                          ArrayRef<const RVec> x,
                          nbnxn_atomdata_t*    nbat)
{
    cellOffset_ = cellOffset;

    const int nthread = gmx_omp_nthreads_get(ModuleMultiThread::Pairsearch);

    const int numAtomsPerCell = geometry_.numAtomsPerCell_;

    /* Make the cell index as a function of x and y */
    int ncz_max = 0;
    int ncz     = 0;
    cxy_ind_[0] = 0;
    for (int i = 0; i < numColumns() + 1; i++)
    {
        /* We set ncz_max at the beginning of the loop iso at the end
         * to skip i=grid->ncx*grid->numCells[YY] which are moved particles
         * that do not need to be ordered on the grid.
         */
        if (ncz > ncz_max)
        {
            ncz_max = ncz;
        }
        int cxy_na_i = gridWork[0].numAtomsPerColumn[i];
        for (int thread = 1; thread < nthread; thread++)
        {
            cxy_na_i += gridWork[thread].numAtomsPerColumn[i];
        }
        ncz = divideRoundUp(cxy_na_i, numAtomsPerCell);
        if (geometry_.numAtomsJCluster_ == 2 * numAtomsPerCell)
        {
            /* Make the number of cell a multiple of 2 */
            ncz = (ncz + 1) & ~1;
        }
        cxy_ind_[i + 1] = cxy_ind_[i] + ncz;
        /* Clear cxy_na_, so we can reuse the array below */
        cxy_na_[i] = 0;
    }
    numCellsTotal_     = cxy_ind_[numColumns()] - cxy_ind_[0];
    numCellsColumnMax_ = ncz_max;

    /* Resize grid and atom data which depend on the number of cells.
     * Note that only the home zone (temporarily) contains moved atoms. The halo zones are
     * set up after moving atoms between home zones and removing them from the local grid.
     */
    int numAtomsMoved = 0;
    if (isHomeZone(ddZone))
    {
        // Note that atomRange can include filler particles, but an overestimate is fine
        numAtomsMoved = (cxy_ind_[numColumns() + 1] - cxy_ind_[numColumns()]) * numAtomsPerCell;
    }
    resizeForNumberOfCells(atomIndexEnd(), numAtomsMoved, ddZone, gridSetData, nbat);

    if (debug)
    {
        fprintf(debug,
                "ns na_sc %d na_c %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                numAtomsPerCell,
                geometry_.numAtomsICluster_,
                numCellsTotal_,
                dimensions_.numCells[XX],
                dimensions_.numCells[YY],
                numCellsTotal_ / (static_cast<double>(numColumns())),
                ncz_max);
        if (gmx_debug_at)
        {
            int i = 0;
            for (int cy = 0; cy < dimensions_.numCells[YY]; cy++)
            {
                for (int cx = 0; cx < dimensions_.numCells[XX]; cx++)
                {
                    fprintf(debug, " %2d", cxy_ind_[i + 1] - cxy_ind_[i]);
                    i++;
                }
                fprintf(debug, "\n");
            }
        }
    }

    /* Make sure the work array for sorting is large enough */
    const int worstCaseSortBufferSize = ncz_max * numAtomsPerCell * c_sortGridMaxSizeFactor;
    if (worstCaseSortBufferSize > Index(gridWork[0].sortBuffer.size()))
    {
        for (GridWork& work : gridWork)
        {
            /* Elements not in use should be -1 */
            work.sortBuffer.resize(worstCaseSortBufferSize, -1);
        }
    }

    /* Now we know the dimensions we can fill the grid.
     * This is the first, unsorted fill. We sort the columns after this.
     */
    ArrayRef<int> cells       = gridSetData->cells;
    ArrayRef<int> atomIndices = gridSetData->atomIndices;
    for (int i : atomRange)
    {
        /* At this point nbs->cell contains the local grid x,y indices */
        const int cxy                                        = cells[i];
        atomIndices[firstAtomInColumn(cxy) + cxy_na_[cxy]++] = i;
    }

    if (isHomeZone(ddZone))
    {
        /* Set the cell indices for the moved particles */
        int n0 = numCellsTotal_ * numAtomsPerCell;
        int n1 = numCellsTotal_ * numAtomsPerCell + cxy_na_[numColumns()];
        for (int i = n0; i < n1; i++)
        {
            cells[atomIndices[i]] = i;
        }
    }

    /* Sort the super-cell columns along z into the sub-cells. */
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            Range<int> columnRange(((thread + 0) * numColumns()) / nthread,
                                   ((thread + 1) * numColumns()) / nthread);
            if (geometry_.isSimple_)
            {
                sortColumnsCpuGeometry(
                        gridSetData, ddZone, atomInfo, x, nbat, columnRange, gridWork[thread].sortBuffer);
            }
            else
            {
                sortColumnsGpuGeometry(
                        gridSetData, ddZone, atomInfo, x, nbat, columnRange, gridWork[thread].sortBuffer);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    if (geometry_.isSimple_ && geometry_.numAtomsJCluster_ == 2 * numAtomsPerCell)
    {
        combine_bounding_box_pairs(*this, bb_, bbj_);
    }

    if (!geometry_.isSimple_)
    {
        numClustersTotal_ = 0;
        for (int i = 0; i < numCellsTotal_; i++)
        {
            numClustersTotal_ += numClusters_[i];
        }
    }

    if (debug)
    {
        if (geometry_.isSimple_)
        {
            print_bbsizes_simple(debug, *this);
        }
        else
        {
            fprintf(debug,
                    "ns non-zero sub-cells: %d average atoms %.2f\n",
                    numClustersTotal_,
                    atomRange.size() / static_cast<double>(numClustersTotal_));

            print_bbsizes_supersub(debug, *this);
        }
    }
}

real generateAndFill2DGrid(Grid*                  grid,
                           ArrayRef<GridWork>     gridWork,
                           HostVector<int>*       cells,
                           const rvec             lowerCorner,
                           const rvec             upperCorner,
                           const UpdateGroupsCog* updateGroupsCog,
                           const Range<int>       atomRange,
                           const int              numAtomsWithoutFillers,
                           real*                  atomDensity,
                           const real             maxAtomGroupRadius,
                           ArrayRef<const RVec>   x,
                           const int              ddZone,
                           const int*             move,
                           const bool             computeGridDensityRatio)
{
    GMX_ASSERT(numAtomsWithoutFillers <= atomRange.end(),
               "The real atoms are a subset of atomRange");

    grid->setDimensions(
            ddZone, atomRange.size(), numAtomsWithoutFillers, lowerCorner, upperCorner, atomDensity, maxAtomGroupRadius);

    for (GridWork& work : gridWork)
    {
        work.numAtomsPerColumn.resize(grid->numColumns() + 1);
    }

    /* Make space for the new cell indices */
    cells->resize(*atomRange.end());

    const int nthread = gmx_omp_nthreads_get(ModuleMultiThread::Pairsearch);
    GMX_ASSERT(nthread > 0, "We expect the OpenMP thread count to be set");

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            Grid::calcColumnIndices(grid->dimensions(),
                                    updateGroupsCog,
                                    atomRange,
                                    x,
                                    ddZone,
                                    move,
                                    thread,
                                    nthread,
                                    *cells,
                                    gridWork[thread].numAtomsPerColumn);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    real gridDensityRatio = 0;

    if (computeGridDensityRatio)
    {
        GMX_RELEASE_ASSERT(
                numAtomsWithoutFillers > 0,
                "Number of atoms should be > 0 when density ratio computing is requested");

        // Compute the effective density ratio on the current grid
        int64_t sumAtomsInColumnSquared = 0;
        for (int i = 0; i < grid->numColumns(); i++)
        {
            int64_t numAtomsInColumn = 0;
            for (int thread = 0; thread < nthread; thread++)
            {
                numAtomsInColumn += gridWork[thread].numAtomsPerColumn[i];
            }
            sumAtomsInColumnSquared += square(numAtomsInColumn);
        }
        // The effective density divided by the uniform density
        gridDensityRatio = sumAtomsInColumnSquared * grid->numColumns()
                           / gmx::square(real(numAtomsWithoutFillers));
        if (debug)
        {
            fprintf(debug, "ns grid effective density ratio %f\n", gridDensityRatio);
        }
    }

    return gridDensityRatio;
}

void Grid::setNonLocalGrid(const int                           ddZone,
                           const GridDimensions&               dimensions,
                           ArrayRef<const std::pair<int, int>> columns,
                           const int                           cellOffset,
                           ArrayRef<const int32_t>             atomInfo,
                           ArrayRef<const RVec>                x,
                           GridSetData*                        gridSetData,
                           nbnxn_atomdata_t*                   nbat)
{
    ddZone_ = ddZone;

    dimensions_ = dimensions;

    // Clear the grid
    cxy_na_.clear();
    cxy_ind_.resize(1);

    cellOffset_ = cellOffset;

    const int numAtomsPerCell = geometry_.numAtomsPerCell_;

    // Set the cluster counts for the columns in the range of communicated columns
    int lastColumnIndex = -1;
    int cellIndex       = 0;
    numCellsColumnMax_  = 0;
    for (const auto& columnInfo : columns)
    {
        const int columnIndex = columnInfo.first;

        GMX_ASSERT(columnIndex > lastColumnIndex, "The input columns should be ordered");

        // When we skipped columns this fills those with appropriate values
        cxy_na_.resize(columnIndex, 0);
        cxy_ind_.resize(columnIndex + 1, cellIndex);
        numClusters_.resize(columnIndex, 0);

        const int numCellsInColumn = columnInfo.second;
        const int numAtomsInColumn = numCellsInColumn * numAtomsPerCell;
        GMX_ASSERT(numAtomsInColumn % geometry().numAtomsICluster_ == 0
                           && numAtomsInColumn % geometry().numAtomsJCluster_ == 0,
                   "The number of cell in a column should be a multiple of the cluster sizes");

        cellIndex += numCellsInColumn;
        cxy_na_.push_back(numAtomsInColumn);
        cxy_ind_.push_back(cellIndex);

        numCellsColumnMax_ = std::max(numCellsColumnMax_, numAtomsInColumn);

        lastColumnIndex = columnIndex;
    }

    GMX_ASSERT(gmx::ssize(cxy_ind_) == lastColumnIndex + 2,
               "cxy_ind should have lastColumnIndex + 1 entries");

    // Small optimization: re-compute the cell count along x from the current number of columns
    dimensions_.numCells[XX] = gmx::divideRoundUp(lastColumnIndex + 1, dimensions_.numCells[YY]);

    // Set the data for the remaining, empty, columns
    cxy_na_.resize(numColumns(), 0);
    cxy_ind_.resize(numColumns() + 1, cellIndex);
    numClusters_.resize(numColumns(), 0);

    numCellsTotal_ = cxy_ind_.back() - cxy_ind_[0];

    // Note that cells and atomIndices are currently (scaled) identity mappings
    const int numAtomsOverAllGrids = (cellOffset_ + numCellsTotal_) * numAtomsPerCell;
    gridSetData->cells.resize(numAtomsOverAllGrids);
    gridSetData->atomIndices.resize(numAtomsOverAllGrids);
    resizeBoundingBoxesAndFlags(numCellsTotal_);

    // Set the bounding boxes and the interaction flags for all cells in the grid
    numClusters_.resize(numCellsTotal_);
    // Prepare a non-member OpenMP reduction variable
    int numClustersTotal = 0;

    const int gmx_unused numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Pairsearch);
#pragma omp parallel for num_threads(numThreads) reduction(+ : numClustersTotal) schedule(static)
    for (int cell = 0; cell < numCellsTotal_; cell++)
    {
        const int atomOffsetCell = (cellOffset_ + cell) * numAtomsPerCell;

        // This should be removed when we remove the (identity) atomIndices
        for (int a = atomOffsetCell; a < atomOffsetCell + numAtomsPerCell; a++)
        {
            gridSetData->atomIndices[a] = a;
        }

        // Do not pass filler particles to fillCell()
        int atomEnd = atomOffsetCell + numAtomsPerCell;
        while (atomEnd > atomOffsetCell && atomInfo[atomEnd - 1] == gmx::sc_atomInfo_IsFillerParticle)
        {
            atomEnd--;
        }

        int numClustersInCell;

        if (geometry_.isSimple_)
        {
            fillCell(gridSetData, nbat, atomOffsetCell, atomEnd, atomInfo, x);

            int nonEmptyCell;
            if (atomEnd > atomOffsetCell)
            {
                numClustersInCell = 1;
                nonEmptyCell      = cell;
            }
            else
            {
                numClustersInCell = 0;
                // We have a cell with only filler particles, copy the previous bounding box
                nonEmptyCell = cell - 1;
            }
            bbcz_[cell].lower = bb_[nonEmptyCell].lower.z;
            bbcz_[cell].upper = bb_[nonEmptyCell].upper.z;
        }
        else
        {
            const int numAtomsPerCluster = geometry_.numAtomsICluster_;

            numClustersInCell = divideRoundUp(atomEnd - atomOffsetCell, numAtomsPerCluster);

            bbcz_[cell].lower = x[atomOffsetCell][ZZ];
            bbcz_[cell].upper = x[atomOffsetCell][ZZ];
            for (int c = 0; c < numClustersInCell; c++)
            {
                const int atomClusterStart = atomOffsetCell + c * numAtomsPerCluster;
                const int atomClusterEnd = std::min(atomClusterStart + numAtomsPerCluster, atomEnd);
                fillCell(gridSetData, nbat, atomClusterStart, atomClusterEnd, atomInfo, x);

#if GMX_DOUBLE
                GMX_RELEASE_ASSERT(false,
                                   "GPU bounding boxes not (yet) handled with double precision");
#else
                for (int a = atomClusterStart; a < atomClusterEnd; a++)
                {
                    bbcz_[cell].lower = std::min(bbcz_[cell].lower, x[a][ZZ]);
                    bbcz_[cell].upper = std::max(bbcz_[cell].upper, x[a][ZZ]);
                }
#endif
            }

            GMX_ASSERT(bbcz_[cell].upper >= bbcz_[cell].lower,
                       "Upper bound should be >= lower bound");
        }

        numClusters_[cell] = numClustersInCell;
        numClustersTotal += numClustersInCell;
    }

    // Store the reduced value in the member variable. The separate reduction
    // variable will not be needed when GROMACS requires OpenMP 5.1.
    numClustersTotal_ = numClustersTotal;

    if (geometry_.isSimple_ && geometry_.numAtomsJCluster_ == 2 * numAtomsPerCell)
    {
        combine_bounding_box_pairs(*this, bb_, bbj_);
    }
}

} // namespace gmx
