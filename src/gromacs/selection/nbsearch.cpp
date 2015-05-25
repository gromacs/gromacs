/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements neighborhood searching for analysis (from nbsearch.h).
 *
 * High-level overview of the algorithm is at \ref page_analysisnbsearch.
 *
 * \todo
 * The grid implementation could still be optimized in several different ways:
 *   - Pruning grid cells from the search list if they are completely outside
 *     the sphere that is being considered.
 *   - A better heuristic could be added for falling back to simple loops for a
 *     small number of reference particles.
 *   - A better heuristic for selecting the grid size.
 *   - A multi-level grid implementation could be used to be able to use small
 *     grids for short cutoffs with very inhomogeneous particle distributions
 *     without a memory cost.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "nbsearch.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <vector>

#include "thread_mpi/mutex.h"

#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/position.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \brief
 * Computes the bounding box for a set of positions.
 *
 * \param[in]  posCount Number of positions in \p x.
 * \param[in]  x        Positions to compute the bounding box for.
 * \param[out] origin   Origin of the bounding box.
 * \param[out] size     Size of the bounding box.
 */
void computeBoundingBox(int posCount, const rvec x[], rvec origin, rvec size)
{
    rvec maxBound;
    copy_rvec(x[0], origin);
    copy_rvec(x[0], maxBound);
    for (int i = 1; i < posCount; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            if (origin[d] > x[i][d])
            {
                origin[d] = x[i][d];
            }
            if (maxBound[d] < x[i][d])
            {
                maxBound[d] = x[i][d];
            }
        }
    }
    rvec_sub(maxBound, origin, size);
}

}   // namespace

namespace internal
{

/********************************************************************
 * Implementation class declarations
 */

class AnalysisNeighborhoodSearchImpl
{
    public:
        typedef AnalysisNeighborhoodPairSearch::ImplPointer
            PairSearchImplPointer;
        typedef std::vector<PairSearchImplPointer> PairSearchList;
        typedef std::vector<std::vector<int> > CellList;

        explicit AnalysisNeighborhoodSearchImpl(real cutoff);
        ~AnalysisNeighborhoodSearchImpl();

        /*! \brief
         * Initializes the search with a given box and reference positions.
         *
         * \param[in] mode            Search mode to use.
         * \param[in] bXY             Whether to use 2D searching.
         * \param[in] excls           Exclusions.
         * \param[in] pbc             PBC information.
         * \param[in] positions       Set of reference positions.
         */
        void init(AnalysisNeighborhood::SearchMode     mode,
                  bool                                 bXY,
                  const t_blocka                      *excls,
                  const t_pbc                         *pbc,
                  const AnalysisNeighborhoodPositions &positions);
        PairSearchImplPointer getPairSearch();

        real cutoffSquared() const { return cutoff2_; }
        bool usesGridSearch() const { return bGrid_; }

    private:
        /*! \brief
         * Checks the efficiency and possibility of doing grid-based searching.
         *
         * \param[in] bForce  If `true`, grid search will be forced if possible.
         * \returns   `false` if grid search is not suitable.
         */
        bool checkGridSearchEfficiency(bool bForce);
        /*! \brief
         * Determines a suitable grid size and sets up the cells.
         *
         * \param[in] box          Box vectors (should not have zero vectors).
         * \param[in] bSingleCell  If `true`, the corresponding dimension will
         *     be forced to use a single cell.
         * \param[in] posCount     Number of positions that will be put on the
         *     grid.
         * \returns   `false` if grid search is not suitable.
         */
        bool initGridCells(const matrix box, bool bSingleCell[DIM],
                           int posCount);
        /*! \brief
         * Sets ua a search grid for a given box.
         *
         * \param[in] pbc      Information about the box.
         * \param[in] posCount Number of positions in \p x.
         * \param[in] x        Reference positions that will be put on the grid.
         * \param[in] bForce   If `true`, grid searching will be used if at all
         *     possible, even if a simple search might give better performance.
         * \returns   `false` if grid search is not suitable.
         */
        bool initGrid(const t_pbc &pbc, int posCount, const rvec x[], bool bForce);
        /*! \brief
         * Maps a point into a grid cell.
         *
         * \param[in]  x    Point to map.
         * \param[out] cell Fractional cell coordinates of \p x on the grid.
         * \param[out] xout Coordinates to use.
         *
         * \p xout will be within the rectangular unit cell in dimensions where
         * the grid is periodic.  For other dimensions, both \p xout and
         * \p cell can be outside the grid/unit cell.
         */
        void mapPointToGridCell(const rvec x, rvec cell, rvec xout) const;
        /*! \brief
         * Calculates linear index of a grid cell.
         *
         * \param[in]  cell Cell indices (must be within the grid).
         * \returns    Linear index of \p cell.
         */
        int getGridCellIndex(const ivec cell) const;
        /*! \brief
         * Adds an index into a grid cell.
         *
         * \param[in]  cell Fractional cell coordinates into which \p i should
         *     be added.
         * \param[in]  i    Index to add.
         *
         * \p cell should satisfy the conditions that \p mapPointToGridCell()
         * produces.
         */
        void addToGridCell(const rvec cell, int i);
        /*! \brief
         * Initializes a cell pair loop for a dimension.
         *
         * \param[in]     centerCell Fractional cell coordiates of the particle
         *     for which pairs are being searched.
         * \param[in,out] cell       Current/initial cell to loop over.
         * \param[in,out] upperBound Last cell to loop over.
         * \param[in]     dim        Dimension to initialize in this call.
         *
         * Initializes `cell[dim]` and `upperBound[dim]` for looping over
         * neighbors of a particle at position given by \p centerCell.
         * If 'dim != ZZ`, `cell[d]` (`d > dim`) set the plane/row of cells
         * for which the loop is initialized.  The loop should then go from
         * `cell[dim]` until `upperBound[dim]`, inclusive.
         * `cell[d]` with `d < dim` or `upperBound[d]` with `d != dim` are not
         * modified by this function.
         *
         * `cell` and `upperBound` may be outside the grid for periodic
         * dimensions and need to be shifted separately: to simplify the
         * looping, the range is always (roughly) symmetric around the value in
         * `centerCell`.
         */
        void initCellRange(const rvec centerCell, ivec cell,
                           ivec upperBound, int dim) const;
        /*! \brief
         * Advances cell pair loop to the next cell.
         *
         * \param[in]     centerCell Fractional cell coordiates of the particle
         *     for which pairs are being searched.
         * \param[in,out] cell       Current (in)/next (out) cell in the loop.
         * \param[in,out] upperBound Last cell in the loop for each dimension.
         */
        bool nextCell(const rvec centerCell, ivec cell, ivec upperBound) const;
        /*! \brief
         * Calculates the index and shift of a grid cell during looping.
         *
         * \param[in]  cell       Unshifted cell index.
         * \param[out] shift      Shift to apply to get the periodic distance
         *     for distances between the cells.
         * \returns    Grid cell index corresponding to `cell`.
         */
        int shiftCell(const ivec cell, rvec shift) const;

        //! Whether to try grid searching.
        bool                    bTryGrid_;
        //! The cutoff.
        real                    cutoff_;
        //! The cutoff squared.
        real                    cutoff2_;
        //! Whether to do searching in XY plane only.
        bool                    bXY_;

        //! Number of reference points for the current frame.
        int                     nref_;
        //! Reference point positions.
        const rvec             *xref_;
        //! Reference position exclusion IDs.
        const int              *refExclusionIds_;
        //! Reference position indices (NULL if no indices).
        const int              *refIndices_;
        //! Exclusions.
        const t_blocka         *excls_;
        //! PBC data.
        t_pbc                   pbc_;

        //! Whether grid searching is actually used for the current positions.
        bool                    bGrid_;
        //! false if the box is rectangular.
        bool                    bTric_;
        //! Whether the grid is periodic in a dimension.
        bool                    bGridPBC_[DIM];
        //! Array for storing in-unit-cell reference positions.
        std::vector<RVec>       xrefAlloc_;
        //! Origin of the grid (zero for periodic dimensions).
        rvec                    gridOrigin_;
        //! Size of a single grid cell.
        rvec                    cellSize_;
        //! Inverse of \p cellSize_. Zero for dimensions where grid is not used.
        rvec                    invCellSize_;
        /*! \brief
         * Shift in cell coordinates (for triclinic boxes) in X when crossing
         * the Z periodic boundary.
         */
        real                    cellShiftZX_;
        /*! \brief
         * Shift in cell coordinates (for triclinic boxes) in Y when crossing
         * the Z periodic boundary.
         */
        real                    cellShiftZY_;
        /*! \brief
         * Shift in cell coordinates (for triclinic boxes) in X when crossing
         * the Y periodic boundary.
         */
        real                    cellShiftYX_;
        //! Number of cells along each dimension.
        ivec                    ncelldim_;
        //! Data structure to hold the grid cell contents.
        CellList                cells_;

        tMPI::mutex             createPairSearchMutex_;
        PairSearchList          pairSearchList_;

        friend class AnalysisNeighborhoodPairSearchImpl;

        GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisNeighborhoodSearchImpl);
};

class AnalysisNeighborhoodPairSearchImpl
{
    public:
        explicit AnalysisNeighborhoodPairSearchImpl(const AnalysisNeighborhoodSearchImpl &search)
            : search_(search)
        {
            testPosCount_     = 0;
            testPositions_    = NULL;
            testExclusionIds_ = NULL;
            testIndices_      = NULL;
            nexcl_            = 0;
            excl_             = NULL;
            clear_rvec(xtest_);
            clear_rvec(testcell_);
            clear_ivec(currCell_);
            clear_ivec(cellBound_);
            reset(-1);
        }

        //! Initializes a search to find reference positions neighboring \p x.
        void startSearch(const AnalysisNeighborhoodPositions &positions);
        //! Searches for the next neighbor.
        template <class Action>
        bool searchNext(Action action);
        //! Initializes a pair representing the pair found by searchNext().
        void initFoundPair(AnalysisNeighborhoodPair *pair) const;
        //! Advances to the next test position, skipping any remaining pairs.
        void nextTestPosition();

    private:
        //! Clears the loop indices.
        void reset(int testIndex);
        //! Checks whether a reference positiong should be excluded.
        bool isExcluded(int j);

        //! Parent search object.
        const AnalysisNeighborhoodSearchImpl   &search_;
        //! Number of test positions.
        int                                     testPosCount_;
        //! Reference to the test positions.
        const rvec                             *testPositions_;
        //! Reference to the test exclusion indices.
        const int                              *testExclusionIds_;
        //! Reference to the test position indices.
        const int                              *testIndices_;
        //! Number of excluded reference positions for current test particle.
        int                                     nexcl_;
        //! Exclusions for current test particle.
        const int                              *excl_;
        //! Index of the currently active test position in \p testPositions_.
        int                                     testIndex_;
        //! Stores test position during a pair loop.
        rvec                                    xtest_;
        //! Stores the previous returned position during a pair loop.
        int                                     previ_;
        //! Stores the pair distance corresponding to previ_;
        real                                    prevr2_;
        //! Stores the shortest distance vector corresponding to previ_;
        rvec                                    prevdx_;
        //! Stores the current exclusion index during loops.
        int                                     exclind_;
        //! Stores the fractional test particle cell location during loops.
        rvec                                    testcell_;
        //! Stores the current cell during pair loops.
        ivec                                    currCell_;
        //! Stores the current loop upper bounds for each dimension during pair loops.
        ivec                                    cellBound_;
        //! Stores the index within the current cell during pair loops.
        int                                     prevcai_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisNeighborhoodPairSearchImpl);
};

/********************************************************************
 * AnalysisNeighborhoodSearchImpl
 */

AnalysisNeighborhoodSearchImpl::AnalysisNeighborhoodSearchImpl(real cutoff)
{
    bTryGrid_       = true;
    cutoff_         = cutoff;
    if (cutoff_ <= 0)
    {
        cutoff_     = cutoff2_ = GMX_REAL_MAX;
        bTryGrid_   = false;
    }
    else
    {
        cutoff2_        = sqr(cutoff_);
    }
    bXY_             = false;
    nref_            = 0;
    xref_            = NULL;
    refExclusionIds_ = NULL;
    refIndices_      = NULL;
    std::memset(&pbc_, 0, sizeof(pbc_));

    bGrid_          = false;
    bTric_          = false;
    bGridPBC_[XX]   = true;
    bGridPBC_[YY]   = true;
    bGridPBC_[ZZ]   = true;

    clear_rvec(gridOrigin_);
    clear_rvec(cellSize_);
    clear_rvec(invCellSize_);
    clear_ivec(ncelldim_);
}

AnalysisNeighborhoodSearchImpl::~AnalysisNeighborhoodSearchImpl()
{
    PairSearchList::const_iterator i;
    for (i = pairSearchList_.begin(); i != pairSearchList_.end(); ++i)
    {
        GMX_RELEASE_ASSERT(i->unique(),
                           "Dangling AnalysisNeighborhoodPairSearch reference");
    }
}

AnalysisNeighborhoodSearchImpl::PairSearchImplPointer
AnalysisNeighborhoodSearchImpl::getPairSearch()
{
    tMPI::lock_guard<tMPI::mutex> lock(createPairSearchMutex_);
    // TODO: Consider whether this needs to/can be faster, e.g., by keeping a
    // separate pool of unused search objects.
    PairSearchList::const_iterator i;
    for (i = pairSearchList_.begin(); i != pairSearchList_.end(); ++i)
    {
        if (i->unique())
        {
            return *i;
        }
    }
    PairSearchImplPointer pairSearch(new AnalysisNeighborhoodPairSearchImpl(*this));
    pairSearchList_.push_back(pairSearch);
    return pairSearch;
}

bool AnalysisNeighborhoodSearchImpl::checkGridSearchEfficiency(bool bForce)
{
    // Find the extent of the sphere in cells.
    ivec  range;
    for (int dd = 0; dd < DIM; ++dd)
    {
        range[dd] = static_cast<int>(ceil(cutoff_ * invCellSize_[dd]));
    }

    // Calculate the fraction of cell pairs that need to be searched,
    // and check that the cutoff is not too large for periodic dimensions.
    real coveredCells = 1.0;
    for (int dd = 0; dd < DIM; ++dd)
    {
        const int cellCount    = ncelldim_[dd];
        const int coveredCount = 2 * range[dd] + 1;
        if (bGridPBC_[dd])
        {
            if (coveredCount > cellCount)
            {
                // Cutoff is too close to half the box size for grid searching
                // (it is not possible to find a single shift for every pair of
                // grid cells).
                return false;
            }
            coveredCells *= coveredCount;
        }
        else
        {
            if (range[dd] >= cellCount - 1)
            {
                range[dd]     = cellCount - 1;
                coveredCells *= cellCount;
            }
            else if (coveredCount > cellCount)
            {
                // The sum of range+1, range+2, ..., range+N/2, ... range+1.
                coveredCells *= range[dd] +
                    static_cast<real>((cellCount + 1)/2 * (cellCount/2 + 1)) / cellCount;
            }
            else
            {
                // The sum of range+1, ..., 2*range+1, ..., 2*range+1, ... range+1.
                coveredCells *= coveredCount -
                    static_cast<real>(range[dd] * (range[dd] + 1)) / cellCount;
            }
        }
    }
    // Magic constant that would need tuning for optimal performance:
    // Don't do grid searching if nearly all cell pairs would anyways need to
    // be looped through.
    const int totalCellCount = ncelldim_[XX] * ncelldim_[YY] * ncelldim_[ZZ];
    if (!bForce && coveredCells >= 0.5 * totalCellCount)
    {
        return false;
    }
    return true;
}

bool AnalysisNeighborhoodSearchImpl::initGridCells(
        const matrix box, bool bSingleCell[DIM], int posCount)
{
    // Determine the size of cubes where there are on average 10 positions.
    // The loop takes care of cases where some of the box edges are shorter
    // than the the desired cube size; in such cases, a single grid cell is
    // used in these dimensions, and the cube size is determined only from the
    // larger box vectors.  Such boxes should be rare, but the bounding box
    // approach can result in very flat boxes with certain types of selections
    // (e.g., for interfacial systems or for small number of atoms).
    real targetsize   = 0.0;
    int  prevDimCount = 4;
    while (true)
    {
        real volume   = 1.0;
        int  dimCount = 3;
        for (int dd = 0; dd < DIM; ++dd)
        {
            const real boxSize = box[dd][dd];
            if (boxSize < targetsize)
            {
                bSingleCell[dd] = true;
                if (bGridPBC_[dd])
                {
                    return false;
                }
            }
            if (bSingleCell[dd])
            {
                --dimCount;
            }
            else
            {
                volume *= boxSize;
            }
        }
        if (dimCount == 0 || dimCount == prevDimCount)
        {
            break;
        }
        targetsize   = pow(volume * 10 / posCount, static_cast<real>(1./dimCount));
        prevDimCount = dimCount;
    }

    int totalCellCount = 1;
    for (int dd = 0; dd < DIM; ++dd)
    {
        int cellCount;
        if (bSingleCell[dd])
        {
            cellCount = 1;
        }
        else
        {
            cellCount = std::max(1, static_cast<int>(box[dd][dd] / targetsize));
            // TODO: If the cell count is one or two, it would be better to
            // just fall back to bSingleCell[dd] = true, and leave the rest to
            // the efficiency check later.
            if (bGridPBC_[dd] && cellCount < 3)
            {
                return false;
            }
        }
        totalCellCount *= cellCount;
        ncelldim_[dd]   = cellCount;
    }
    if (totalCellCount <= 3)
    {
        return false;
    }
    // Never decrease the size of the cell vector to avoid reallocating
    // memory for the nested vectors.  The actual size of the vector is not
    // used outside this function.
    if (cells_.size() < static_cast<size_t>(totalCellCount))
    {
        cells_.resize(totalCellCount);
    }
    for (int ci = 0; ci < totalCellCount; ++ci)
    {
        cells_[ci].clear();
    }
    return true;
}

bool AnalysisNeighborhoodSearchImpl::initGrid(
        const t_pbc &pbc, int posCount, const rvec x[], bool bForce)
{
    if (posCount == 0)
    {
        return false;
    }

    switch (pbc.ePBC)
    {
        case epbcNONE:
            bGridPBC_[XX] = false;
            bGridPBC_[YY] = false;
            bGridPBC_[ZZ] = false;
            break;
        case epbcXY:
            bGridPBC_[XX] = true;
            bGridPBC_[YY] = true;
            bGridPBC_[ZZ] = false;
            break;
        case epbcXYZ:
            bGridPBC_[XX] = true;
            bGridPBC_[YY] = true;
            bGridPBC_[ZZ] = true;
            break;
        default:
            // Grid searching not supported for now with screw.
            return false;
    }

    bool   bSingleCell[DIM] = {false, false, bXY_};
    matrix box;
    copy_mat(pbc.box, box);
    // TODO: In principle, we could use the bounding box for periodic
    // dimensions as well if the bounding box is sufficiently far from the box
    // edges.
    rvec   origin, boundingBoxSize;
    computeBoundingBox(posCount, x, origin, boundingBoxSize);
    clear_rvec(gridOrigin_);
    for (int dd = 0; dd < DIM; ++dd)
    {
        if (!bGridPBC_[dd] && !bSingleCell[dd])
        {
            gridOrigin_[dd] = origin[dd];
            clear_rvec(box[dd]);
            box[dd][dd] = boundingBoxSize[dd];
        }
        // TODO: In case the zero vector comes from the bounding box, this does
        // not lead to a very efficient grid search, but that should be rare.
        if (box[dd][dd] <= 0.0)
        {
            GMX_ASSERT(!bGridPBC_[dd], "Periodic box vector is zero");
            bSingleCell[dd] = true;
            clear_rvec(box[dd]);
            box[dd][dd] = 1.0;
        }
    }

    if (!initGridCells(box, bSingleCell, posCount))
    {
        return false;
    }

    bTric_ = TRICLINIC(pbc.box);
    for (int dd = 0; dd < DIM; ++dd)
    {
        cellSize_[dd] = box[dd][dd] / ncelldim_[dd];
        if (bSingleCell[dd])
        {
            invCellSize_[dd] = 0.0;
        }
        else
        {
            invCellSize_[dd] = 1.0 / cellSize_[dd];
        }
    }
    if (bTric_)
    {
        cellShiftZY_ = box[ZZ][YY] * invCellSize_[YY];
        cellShiftZX_ = box[ZZ][XX] * invCellSize_[XX];
        cellShiftYX_ = box[YY][XX] * invCellSize_[XX];
    }
    return checkGridSearchEfficiency(bForce);
}

void AnalysisNeighborhoodSearchImpl::mapPointToGridCell(const rvec x,
                                                        rvec       cell,
                                                        rvec       xout) const
{
    rvec xtmp;
    rvec_sub(x, gridOrigin_, xtmp);
    // The reverse order is necessary for triclinic cells: shifting in Z may
    // modify also X and Y, and shifting in Y may modify X, so the mapping to
    // a rectangular grid needs to be done in this order.
    for (int dd = DIM - 1; dd >= 0; --dd)
    {
        real cellIndex = xtmp[dd] * invCellSize_[dd];
        if (bGridPBC_[dd])
        {
            const real cellCount = ncelldim_[dd];
            while (cellIndex < 0)
            {
                cellIndex += cellCount;
                rvec_inc(xtmp, pbc_.box[dd]);
            }
            while (cellIndex >= cellCount)
            {
                cellIndex -= cellCount;
                rvec_dec(xtmp, pbc_.box[dd]);
            }
        }
        cell[dd] = cellIndex;
    }
    copy_rvec(xtmp, xout);
}

int AnalysisNeighborhoodSearchImpl::getGridCellIndex(const ivec cell) const
{
    GMX_ASSERT(cell[XX] >= 0 && cell[XX] < ncelldim_[XX],
               "Grid cell X index out of range");
    GMX_ASSERT(cell[YY] >= 0 && cell[YY] < ncelldim_[YY],
               "Grid cell Y index out of range");
    GMX_ASSERT(cell[ZZ] >= 0 && cell[ZZ] < ncelldim_[ZZ],
               "Grid cell Z index out of range");
    return cell[XX] + cell[YY] * ncelldim_[XX]
           + cell[ZZ] * ncelldim_[XX] * ncelldim_[YY];
}

void AnalysisNeighborhoodSearchImpl::addToGridCell(const rvec cell, int i)
{
    ivec icell;
    for (int dd = 0; dd < DIM; ++dd)
    {
        int cellIndex = static_cast<int>(floor(cell[dd]));
        if (!bGridPBC_[dd])
        {
            const int cellCount = ncelldim_[dd];
            if (cellIndex < 0)
            {
                cellIndex = 0;
            }
            else if (cellIndex >= cellCount)
            {
                cellIndex = cellCount - 1;
            }
        }
        icell[dd] = cellIndex;
    }
    const int ci = getGridCellIndex(icell);
    cells_[ci].push_back(i);
}

void AnalysisNeighborhoodSearchImpl::initCellRange(
        const rvec centerCell, ivec currCell, ivec upperBound, int dim) const
{
    // TODO: Prune off cells that are completely outside the cutoff.
    const real range       = cutoff_ * invCellSize_[dim];
    real       startOffset = centerCell[dim] - range;
    real       endOffset   = centerCell[dim] + range;
    if (bTric_)
    {
        switch (dim)
        {
            case ZZ:
                break;
            case YY:
                if (currCell[ZZ] < 0)
                {
                    startOffset += cellShiftZY_;
                    endOffset   += cellShiftZY_;
                }
                else if (currCell[ZZ] >= ncelldim_[ZZ])
                {
                    startOffset -= cellShiftZY_;
                    endOffset   -= cellShiftZY_;
                }
                break;
            case XX:
                if (currCell[ZZ] < 0)
                {
                    startOffset += cellShiftZX_;
                    endOffset   += cellShiftZX_;
                }
                else if (currCell[ZZ] >= ncelldim_[ZZ])
                {
                    startOffset -= cellShiftZX_;
                    endOffset   -= cellShiftZX_;
                }
                if (currCell[YY] < 0)
                {
                    startOffset += cellShiftYX_;
                    endOffset   += cellShiftYX_;
                }
                else if (currCell[YY] >= ncelldim_[YY])
                {
                    startOffset -= cellShiftYX_;
                    endOffset   -= cellShiftYX_;
                }
                break;
        }
    }
    // For non-periodic dimensions, clamp to the actual grid edges.
    if (!bGridPBC_[dim])
    {
        // If endOffset < 0 or startOffset > N, these may cause the whole
        // test position/grid plane/grid row to be skipped.
        if (startOffset < 0)
        {
            startOffset = 0;
        }
        const int cellCount = ncelldim_[dim];
        if (endOffset > cellCount - 1)
        {
            endOffset = cellCount - 1;
        }
    }
    currCell[dim]   = static_cast<int>(floor(startOffset));
    upperBound[dim] = static_cast<int>(floor(endOffset));
}

bool AnalysisNeighborhoodSearchImpl::nextCell(
        const rvec centerCell, ivec cell, ivec upperBound) const
{
    int dim = 0;
    while (dim < DIM)
    {
next:
        ++cell[dim];
        if (cell[dim] > upperBound[dim])
        {
            ++dim;
            continue;
        }
        for (int d = dim - 1; d >= 0; --d)
        {
            initCellRange(centerCell, cell, upperBound, d);
            if (cell[d] > upperBound[d])
            {
                dim = d + 1;
                goto next;
            }
        }
        return true;
    }
    return false;
}

int AnalysisNeighborhoodSearchImpl::shiftCell(const ivec cell, rvec shift) const
{
    ivec shiftedCell;
    copy_ivec(cell, shiftedCell);

    clear_rvec(shift);
    for (int d = 0; d < DIM; ++d)
    {
        const int cellCount = ncelldim_[d];
        if (bGridPBC_[d])
        {
            // A single shift may not be sufficient if the cell must be shifted
            // in more than one dimension, although for each individual
            // dimension it would be.
            while (shiftedCell[d] < 0)
            {
                shiftedCell[d] += cellCount;
                rvec_inc(shift, pbc_.box[d]);
            }
            while (shiftedCell[d] >= cellCount)
            {
                shiftedCell[d] -= cellCount;
                rvec_dec(shift, pbc_.box[d]);
            }
        }
    }

    return getGridCellIndex(shiftedCell);
}

void AnalysisNeighborhoodSearchImpl::init(
        AnalysisNeighborhood::SearchMode     mode,
        bool                                 bXY,
        const t_blocka                      *excls,
        const t_pbc                         *pbc,
        const AnalysisNeighborhoodPositions &positions)
{
    GMX_RELEASE_ASSERT(positions.index_ == -1,
                       "Individual indexed positions not supported as reference");
    bXY_ = bXY;
    if (bXY_ && pbc != NULL && pbc->ePBC != epbcNONE)
    {
        if (pbc->ePBC != epbcXY && pbc->ePBC != epbcXYZ)
        {
            std::string message =
                formatString("Computations in the XY plane are not supported with PBC type '%s'",
                             EPBC(pbc->ePBC));
            GMX_THROW(NotImplementedError(message));
        }
        if (pbc->ePBC == epbcXYZ &&
            (std::fabs(pbc->box[ZZ][XX]) > GMX_REAL_EPS*pbc->box[ZZ][ZZ] ||
             std::fabs(pbc->box[ZZ][YY]) > GMX_REAL_EPS*pbc->box[ZZ][ZZ]))
        {
            GMX_THROW(NotImplementedError("Computations in the XY plane are not supported when the last box vector is not parallel to the Z axis"));
        }
        // Use a single grid cell in Z direction.
        matrix box;
        copy_mat(pbc->box, box);
        clear_rvec(box[ZZ]);
        set_pbc(&pbc_, epbcXY, box);
    }
    else if (pbc != NULL)
    {
        pbc_ = *pbc;
    }
    else
    {
        pbc_.ePBC = epbcNONE;
        clear_mat(pbc_.box);
    }
    nref_ = positions.count_;
    if (mode == AnalysisNeighborhood::eSearchMode_Simple)
    {
        bGrid_ = false;
    }
    else if (bTryGrid_)
    {
        bGrid_ = initGrid(pbc_, positions.count_, positions.x_,
                          mode == AnalysisNeighborhood::eSearchMode_Grid);
    }
    refIndices_ = positions.indices_;
    if (bGrid_)
    {
        xrefAlloc_.resize(nref_);
        xref_ = as_rvec_array(&xrefAlloc_[0]);

        for (int i = 0; i < nref_; ++i)
        {
            const int ii = (refIndices_ != NULL) ? refIndices_[i] : i;
            rvec      refcell;
            mapPointToGridCell(positions.x_[ii], refcell, xrefAlloc_[i]);
            addToGridCell(refcell, i);
        }
    }
    else if (refIndices_ != NULL)
    {
        xrefAlloc_.resize(nref_);
        xref_ = as_rvec_array(&xrefAlloc_[0]);
        for (int i = 0; i < nref_; ++i)
        {
            copy_rvec(positions.x_[refIndices_[i]], xrefAlloc_[i]);
        }
    }
    else
    {
        xref_ = positions.x_;
    }
    excls_           = excls;
    refExclusionIds_ = NULL;
    if (excls != NULL)
    {
        // TODO: Check that the IDs are ascending, or remove the limitation.
        refExclusionIds_ = positions.exclusionIds_;
        GMX_RELEASE_ASSERT(refExclusionIds_ != NULL,
                           "Exclusion IDs must be set for reference positions "
                           "when exclusions are enabled");
    }
}

/********************************************************************
 * AnalysisNeighborhoodPairSearchImpl
 */

void AnalysisNeighborhoodPairSearchImpl::reset(int testIndex)
{
    testIndex_ = testIndex;
    if (testIndex_ >= 0 && testIndex_ < testPosCount_)
    {
        const int index =
            (testIndices_ != NULL ? testIndices_[testIndex] : testIndex);
        if (search_.bGrid_)
        {
            search_.mapPointToGridCell(testPositions_[index], testcell_, xtest_);
            search_.initCellRange(testcell_, currCell_, cellBound_, ZZ);
            search_.initCellRange(testcell_, currCell_, cellBound_, YY);
            search_.initCellRange(testcell_, currCell_, cellBound_, XX);
        }
        else
        {
            copy_rvec(testPositions_[index], xtest_);
        }
        if (search_.excls_ != NULL)
        {
            const int exclIndex  = testExclusionIds_[index];
            if (exclIndex < search_.excls_->nr)
            {
                const int startIndex = search_.excls_->index[exclIndex];
                nexcl_ = search_.excls_->index[exclIndex + 1] - startIndex;
                excl_  = &search_.excls_->a[startIndex];
            }
            else
            {
                nexcl_ = 0;
                excl_  = NULL;
            }
        }
    }
    previ_     = -1;
    prevr2_    = 0.0;
    clear_rvec(prevdx_);
    exclind_   = 0;
    prevcai_   = -1;
}

void AnalysisNeighborhoodPairSearchImpl::nextTestPosition()
{
    if (testIndex_ < testPosCount_)
    {
        ++testIndex_;
        reset(testIndex_);
    }
}

bool AnalysisNeighborhoodPairSearchImpl::isExcluded(int j)
{
    if (exclind_ < nexcl_)
    {
        const int index =
            (search_.refIndices_ != NULL ? search_.refIndices_[j] : j);
        const int refId = search_.refExclusionIds_[index];
        while (exclind_ < nexcl_ && excl_[exclind_] < refId)
        {
            ++exclind_;
        }
        if (exclind_ < nexcl_ && refId == excl_[exclind_])
        {
            ++exclind_;
            return true;
        }
    }
    return false;
}

void AnalysisNeighborhoodPairSearchImpl::startSearch(
        const AnalysisNeighborhoodPositions &positions)
{
    testPosCount_     = positions.count_;
    testPositions_    = positions.x_;
    testExclusionIds_ = positions.exclusionIds_;
    testIndices_      = positions.indices_;
    GMX_RELEASE_ASSERT(search_.excls_ == NULL || testExclusionIds_ != NULL,
                       "Exclusion IDs must be set when exclusions are enabled");
    if (positions.index_ < 0)
    {
        reset(0);
    }
    else
    {
        // Somewhat of a hack: setup the array such that only the last position
        // will be used.
        testPosCount_ = positions.index_ + 1;
        reset(positions.index_);
    }
}

template <class Action>
bool AnalysisNeighborhoodPairSearchImpl::searchNext(Action action)
{
    while (testIndex_ < testPosCount_)
    {
        if (search_.bGrid_)
        {
            int cai = prevcai_ + 1;

            do
            {
                rvec      shift;
                const int ci       = search_.shiftCell(currCell_, shift);
                const int cellSize = static_cast<int>(search_.cells_[ci].size());
                for (; cai < cellSize; ++cai)
                {
                    const int i = search_.cells_[ci][cai];
                    if (isExcluded(i))
                    {
                        continue;
                    }
                    rvec       dx;
                    rvec_sub(search_.xref_[i], xtest_, dx);
                    rvec_sub(dx, shift, dx);
                    const real r2
                        = search_.bXY_
                            ? dx[XX]*dx[XX] + dx[YY]*dx[YY]
                            : norm2(dx);
                    if (r2 <= search_.cutoff2_)
                    {
                        if (action(i, r2, dx))
                        {
                            prevcai_ = cai;
                            previ_   = i;
                            prevr2_  = r2;
                            copy_rvec(dx, prevdx_);
                            return true;
                        }
                    }
                }
                exclind_ = 0;
                cai      = 0;
            }
            while (search_.nextCell(testcell_, currCell_, cellBound_));
        }
        else
        {
            for (int i = previ_ + 1; i < search_.nref_; ++i)
            {
                if (isExcluded(i))
                {
                    continue;
                }
                rvec dx;
                if (search_.pbc_.ePBC != epbcNONE)
                {
                    pbc_dx(&search_.pbc_, search_.xref_[i], xtest_, dx);
                }
                else
                {
                    rvec_sub(search_.xref_[i], xtest_, dx);
                }
                const real r2
                    = search_.bXY_
                        ? dx[XX]*dx[XX] + dx[YY]*dx[YY]
                        : norm2(dx);
                if (r2 <= search_.cutoff2_)
                {
                    if (action(i, r2, dx))
                    {
                        previ_  = i;
                        prevr2_ = r2;
                        copy_rvec(dx, prevdx_);
                        return true;
                    }
                }
            }
        }
        nextTestPosition();
    }
    return false;
}

void AnalysisNeighborhoodPairSearchImpl::initFoundPair(
        AnalysisNeighborhoodPair *pair) const
{
    if (previ_ < 0)
    {
        *pair = AnalysisNeighborhoodPair();
    }
    else
    {
        *pair = AnalysisNeighborhoodPair(previ_, testIndex_, prevr2_, prevdx_);
    }
}

}   // namespace internal

namespace
{

/*! \brief
 * Search action to find the next neighbor.
 *
 * Used as the action for AnalysisNeighborhoodPairSearchImpl::searchNext() to
 * find the next neighbor.
 *
 * Simply breaks the loop on the first found neighbor.
 */
bool withinAction(int /*i*/, real /*r2*/, const rvec /* dx */)
{
    return true;
}

/*! \brief
 * Search action find the minimum distance.
 *
 * Used as the action for AnalysisNeighborhoodPairSearchImpl::searchNext() to
 * find the nearest neighbor.
 *
 * With this action, AnalysisNeighborhoodPairSearchImpl::searchNext() always
 * returns false, and the output is put into the variables passed by pointer
 * into the constructor.  If no neighbors are found, the output variables are
 * not modified, i.e., the caller must initialize them.
 */
class MindistAction
{
    public:
        /*! \brief
         * Initializes the action with given output locations.
         *
         * \param[out] closestPoint Index of the closest reference location.
         * \param[out] minDist2     Minimum distance squared.
         * \param[out] dx           Shortest distance vector.
         *
         * The constructor call does not modify the pointed values, but only
         * stores the pointers for later use.
         * See the class description for additional semantics.
         */
        MindistAction(int *closestPoint, real *minDist2, rvec *dx)
            : closestPoint_(*closestPoint), minDist2_(*minDist2), dx_(*dx)
        {
        }

        //! Processes a neighbor to find the nearest point.
        bool operator()(int i, real r2, const rvec dx)
        {
            if (r2 < minDist2_)
            {
                closestPoint_ = i;
                minDist2_     = r2;
                copy_rvec(dx, dx_);
            }
            return false;
        }

    private:
        int     &closestPoint_;
        real    &minDist2_;
        rvec    &dx_;

        GMX_DISALLOW_ASSIGN(MindistAction);
};

}   // namespace

/********************************************************************
 * AnalysisNeighborhood::Impl
 */

class AnalysisNeighborhood::Impl
{
    public:
        typedef AnalysisNeighborhoodSearch::ImplPointer SearchImplPointer;
        typedef std::vector<SearchImplPointer> SearchList;

        Impl()
            : cutoff_(0), excls_(NULL), mode_(eSearchMode_Automatic), bXY_(false)
        {
        }
        ~Impl()
        {
            SearchList::const_iterator i;
            for (i = searchList_.begin(); i != searchList_.end(); ++i)
            {
                GMX_RELEASE_ASSERT(i->unique(),
                                   "Dangling AnalysisNeighborhoodSearch reference");
            }
        }

        SearchImplPointer getSearch();

        tMPI::mutex             createSearchMutex_;
        SearchList              searchList_;
        real                    cutoff_;
        const t_blocka         *excls_;
        SearchMode              mode_;
        bool                    bXY_;
};

AnalysisNeighborhood::Impl::SearchImplPointer
AnalysisNeighborhood::Impl::getSearch()
{
    tMPI::lock_guard<tMPI::mutex> lock(createSearchMutex_);
    // TODO: Consider whether this needs to/can be faster, e.g., by keeping a
    // separate pool of unused search objects.
    SearchList::const_iterator i;
    for (i = searchList_.begin(); i != searchList_.end(); ++i)
    {
        if (i->unique())
        {
            return *i;
        }
    }
    SearchImplPointer search(new internal::AnalysisNeighborhoodSearchImpl(cutoff_));
    searchList_.push_back(search);
    return search;
}

/********************************************************************
 * AnalysisNeighborhood
 */

AnalysisNeighborhood::AnalysisNeighborhood()
    : impl_(new Impl)
{
}

AnalysisNeighborhood::~AnalysisNeighborhood()
{
}

void AnalysisNeighborhood::setCutoff(real cutoff)
{
    GMX_RELEASE_ASSERT(impl_->searchList_.empty(),
                       "Changing the cutoff after initSearch() not currently supported");
    impl_->cutoff_ = cutoff;
}

void AnalysisNeighborhood::setXYMode(bool bXY)
{
    impl_->bXY_ = bXY;
}

void AnalysisNeighborhood::setTopologyExclusions(const t_blocka *excls)
{
    GMX_RELEASE_ASSERT(impl_->searchList_.empty(),
                       "Changing the exclusions after initSearch() not currently supported");
    impl_->excls_ = excls;
}

void AnalysisNeighborhood::setMode(SearchMode mode)
{
    impl_->mode_ = mode;
}

AnalysisNeighborhood::SearchMode AnalysisNeighborhood::mode() const
{
    return impl_->mode_;
}

AnalysisNeighborhoodSearch
AnalysisNeighborhood::initSearch(const t_pbc                         *pbc,
                                 const AnalysisNeighborhoodPositions &positions)
{
    Impl::SearchImplPointer search(impl_->getSearch());
    search->init(mode(), impl_->bXY_, impl_->excls_,
                 pbc, positions);
    return AnalysisNeighborhoodSearch(search);
}

/********************************************************************
 * AnalysisNeighborhoodSearch
 */

AnalysisNeighborhoodSearch::AnalysisNeighborhoodSearch()
{
}

AnalysisNeighborhoodSearch::AnalysisNeighborhoodSearch(const ImplPointer &impl)
    : impl_(impl)
{
}

void AnalysisNeighborhoodSearch::reset()
{
    impl_.reset();
}

AnalysisNeighborhood::SearchMode AnalysisNeighborhoodSearch::mode() const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    return (impl_->usesGridSearch()
            ? AnalysisNeighborhood::eSearchMode_Grid
            : AnalysisNeighborhood::eSearchMode_Simple);
}

bool AnalysisNeighborhoodSearch::isWithin(
        const AnalysisNeighborhoodPositions &positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    return pairSearch.searchNext(&withinAction);
}

real AnalysisNeighborhoodSearch::minimumDistance(
        const AnalysisNeighborhoodPositions &positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    real          minDist2     = impl_->cutoffSquared();
    int           closestPoint = -1;
    rvec          dx           = {0.0, 0.0, 0.0};
    MindistAction action(&closestPoint, &minDist2, &dx);
    (void)pairSearch.searchNext(action);
    return sqrt(minDist2);
}

AnalysisNeighborhoodPair
AnalysisNeighborhoodSearch::nearestPoint(
        const AnalysisNeighborhoodPositions &positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    real          minDist2     = impl_->cutoffSquared();
    int           closestPoint = -1;
    rvec          dx           = {0.0, 0.0, 0.0};
    MindistAction action(&closestPoint, &minDist2, &dx);
    (void)pairSearch.searchNext(action);
    return AnalysisNeighborhoodPair(closestPoint, 0, minDist2, dx);
}

AnalysisNeighborhoodPairSearch
AnalysisNeighborhoodSearch::startPairSearch(
        const AnalysisNeighborhoodPositions &positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    Impl::PairSearchImplPointer pairSearch(impl_->getPairSearch());
    pairSearch->startSearch(positions);
    return AnalysisNeighborhoodPairSearch(pairSearch);
}

/********************************************************************
 * AnalysisNeighborhoodPairSearch
 */

AnalysisNeighborhoodPairSearch::AnalysisNeighborhoodPairSearch(
        const ImplPointer &impl)
    : impl_(impl)
{
}

bool AnalysisNeighborhoodPairSearch::findNextPair(AnalysisNeighborhoodPair *pair)
{
    bool bFound = impl_->searchNext(&withinAction);
    impl_->initFoundPair(pair);
    return bFound;
}

void AnalysisNeighborhoodPairSearch::skipRemainingPairsForTestPosition()
{
    impl_->nextTestPosition();
}

} // namespace gmx
