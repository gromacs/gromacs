/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
 * \brief
 * Implements neighborhood searching for analysis (from nbsearch.h).
 *
 * High-level overview of the algorithm is at \ref page_analysisnbsearch.
 *
 * \todo
 * The grid implementation could still be optimized in several different ways:
 *   - A better heuristic for selecting the grid size or falling back to a
 *     simple all-pairs search.
 *   - A multi-level grid implementation could be used to be able to use small
 *     grids for short cutoffs with very inhomogeneous particle distributions
 *     without a memory cost.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/nbsearch.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/position.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"
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

} // namespace

namespace internal
{

/********************************************************************
 * Implementation class declarations
 */

class AnalysisNeighborhoodSearchImpl
{
public:
    typedef AnalysisNeighborhoodPairSearch::ImplPointer PairSearchImplPointer;
    typedef std::vector<PairSearchImplPointer>          PairSearchList;
    typedef std::vector<std::vector<int>>               CellList;

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
    void                  init(AnalysisNeighborhood::SearchMode     mode,
                               bool                                 bXY,
                               const ListOfLists<int>*              excls,
                               const t_pbc*                         pbc,
                               const AnalysisNeighborhoodPositions& positions);
    PairSearchImplPointer getPairSearch();

    real cutoffSquared() const { return cutoff2_; }
    bool usesGridSearch() const { return bGrid_; }

private:
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
    bool initGridCells(const matrix box, bool bSingleCell[DIM], int posCount);
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
    bool initGrid(const t_pbc& pbc, int posCount, const rvec x[], bool bForce);
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
     * Calculates linear index of a grid cell from fractional coordinates.
     *
     * \param[in]  cell Cell indices (must be within the grid).
     * \returns    Linear index of \p cell.
     */
    int getGridCellIndex(const rvec cell) const;
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
    void initCellRange(const rvec centerCell, ivec cell, ivec upperBound, int dim) const;
    /*! \brief
     * Computes the extent of the cutoff sphere on a particular cell edge.
     *
     * \param[in]     centerCell Fractional cell coordiates of the particle
     *     for which pairs are being searched.
     * \param[in]     cell       Current cell (for dimensions `>dim`).
     * \param[in]     dim        Dimension to compute in this call.
     * \returns       Fractional extent of the cutoff sphere when looping
     *    over cells in dimension `dim`, for `cell[d]` (`d > dim`).
     *
     * Input parameters are as for initCellRange(), except that if `cell`
     * is over a periodic boundary from `centerCell`, triclinic shifts
     * should have been applied to `centerCell` X/Y components.
     */
    real computeCutoffExtent(RVec centerCell, const ivec cell, int dim) const;
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
    bool bTryGrid_;
    //! The cutoff.
    real cutoff_;
    //! The cutoff squared.
    real cutoff2_;
    //! Whether to do searching in XY plane only.
    bool bXY_;

    //! Number of reference points for the current frame.
    int nref_;
    //! Reference point positions.
    const rvec* xref_;
    //! Reference position exclusion IDs.
    const int* refExclusionIds_;
    //! Reference position indices (NULL if no indices).
    const int* refIndices_;
    //! Exclusions.
    const ListOfLists<int>* excls_;
    //! PBC data.
    t_pbc pbc_;

    //! Whether grid searching is actually used for the current positions.
    bool bGrid_;
    //! false if the box is rectangular.
    bool bTric_;
    //! Whether the grid is periodic in a dimension.
    bool bGridPBC_[DIM];
    //! Array for storing in-unit-cell reference positions.
    std::vector<RVec> xrefAlloc_;
    //! Origin of the grid (zero for periodic dimensions).
    rvec gridOrigin_;
    //! Size of a single grid cell.
    rvec cellSize_;
    //! Inverse of \p cellSize_. Zero for dimensions where grid is not used.
    rvec invCellSize_;
    /*! \brief
     * Shift in cell coordinates (for triclinic boxes) in X when crossing
     * the Z periodic boundary.
     */
    real cellShiftZX_;
    /*! \brief
     * Shift in cell coordinates (for triclinic boxes) in Y when crossing
     * the Z periodic boundary.
     */
    real cellShiftZY_;
    /*! \brief
     * Shift in cell coordinates (for triclinic boxes) in X when crossing
     * the Y periodic boundary.
     */
    real cellShiftYX_;
    //! Number of cells along each dimension.
    ivec ncelldim_;
    //! Data structure to hold the grid cell contents.
    CellList cells_;

    std::mutex     createPairSearchMutex_;
    PairSearchList pairSearchList_;

    friend class AnalysisNeighborhoodPairSearchImpl;

    GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisNeighborhoodSearchImpl);
};

class AnalysisNeighborhoodPairSearchImpl
{
public:
    explicit AnalysisNeighborhoodPairSearchImpl(const AnalysisNeighborhoodSearchImpl& search) :
        search_(search)
    {
        selfSearchMode_   = false;
        testPosCount_     = 0;
        testPositions_    = nullptr;
        testExclusionIds_ = nullptr;
        testIndices_      = nullptr;
        clear_rvec(xtest_);
        clear_rvec(testcell_);
        clear_ivec(currCell_);
        clear_ivec(cellBound_);
        reset(-1);
    }

    //! Initializes a search to find reference positions neighboring \p x.
    void startSearch(const AnalysisNeighborhoodPositions& positions);
    //! Initializes a search to find reference position pairs.
    void startSelfSearch();
    //! Searches for the next neighbor.
    template<class Action>
    bool searchNext(Action action);
    //! Initializes a pair representing the pair found by searchNext().
    void initFoundPair(AnalysisNeighborhoodPair* pair) const;
    //! Advances to the next test position, skipping any remaining pairs.
    void nextTestPosition();

private:
    //! Clears the loop indices.
    void reset(int testIndex);
    //! Checks whether a reference positiong should be excluded.
    bool isExcluded(int j);

    //! Parent search object.
    const AnalysisNeighborhoodSearchImpl& search_;
    //! Whether we are searching for ref-ref pairs.
    bool selfSearchMode_;
    //! Number of test positions.
    int testPosCount_;
    //! Reference to the test positions.
    const rvec* testPositions_;
    //! Reference to the test exclusion indices.
    const int* testExclusionIds_;
    //! Reference to the test position indices.
    const int* testIndices_;
    //! Exclusions for current test particle.
    ArrayRef<const int> excl_;
    //! Index of the currently active test position in \p testPositions_.
    int testIndex_;
    //! Stores test position during a pair loop.
    rvec xtest_;
    //! Stores the previous returned position during a pair loop.
    int previ_;
    //! Stores the pair distance corresponding to previ_.
    real prevr2_;
    //! Stores the shortest distance vector corresponding to previ_.
    rvec prevdx_;
    //! Stores the current exclusion index during loops.
    int exclind_;
    //! Stores the fractional test particle cell location during loops.
    rvec testcell_;
    //! Stores the cell index corresponding to testcell_.
    int testCellIndex_;
    //! Stores the current cell during pair loops.
    ivec currCell_;
    //! Stores the current loop upper bounds for each dimension during pair loops.
    ivec cellBound_;
    //! Stores the index within the current cell during pair loops.
    int prevcai_;

    GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisNeighborhoodPairSearchImpl);
};

/********************************************************************
 * AnalysisNeighborhoodSearchImpl
 */

AnalysisNeighborhoodSearchImpl::AnalysisNeighborhoodSearchImpl(real cutoff)
{
    bTryGrid_ = true;
    cutoff_   = cutoff;
    if (cutoff_ <= 0)
    {
        cutoff_ = cutoff2_ = GMX_REAL_MAX;
        bTryGrid_          = false;
    }
    else
    {
        cutoff2_ = gmx::square(cutoff_);
    }
    bXY_             = false;
    nref_            = 0;
    xref_            = nullptr;
    refExclusionIds_ = nullptr;
    refIndices_      = nullptr;
    std::memset(&pbc_, 0, sizeof(pbc_));

    bGrid_        = false;
    bTric_        = false;
    bGridPBC_[XX] = true;
    bGridPBC_[YY] = true;
    bGridPBC_[ZZ] = true;

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
        GMX_RELEASE_ASSERT(i->use_count() == 1, "Dangling AnalysisNeighborhoodPairSearch reference");
    }
}

AnalysisNeighborhoodSearchImpl::PairSearchImplPointer AnalysisNeighborhoodSearchImpl::getPairSearch()
{
    std::lock_guard<std::mutex> lock(createPairSearchMutex_);
    // TODO: Consider whether this needs to/can be faster, e.g., by keeping a
    // separate pool of unused search objects.
    PairSearchList::const_iterator i;
    for (i = pairSearchList_.begin(); i != pairSearchList_.end(); ++i)
    {
        if (i->use_count() == 1)
        {
            return *i;
        }
    }
    PairSearchImplPointer pairSearch(new AnalysisNeighborhoodPairSearchImpl(*this));
    pairSearchList_.push_back(pairSearch);
    return pairSearch;
}

bool AnalysisNeighborhoodSearchImpl::initGridCells(const matrix box, bool bSingleCell[DIM], int posCount)
{
    // Determine the size of cubes where there are on average 10 positions.
    // The loop takes care of cases where some of the box edges are shorter
    // than the desired cube size; in such cases, a single grid cell is
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
                    // TODO: Consider if a fallback would be possible/better.
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
        targetsize   = std::pow(volume * 10 / posCount, static_cast<real>(1. / dimCount));
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
            // TODO: If the cell count is one or two, it could be better to
            // just fall back to bSingleCell[dd] = true.
            if (bGridPBC_[dd] && cellCount < 3)
            {
                return false;
            }
        }
        totalCellCount *= cellCount;
        ncelldim_[dd] = cellCount;
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

bool AnalysisNeighborhoodSearchImpl::initGrid(const t_pbc& pbc, int posCount, const rvec x[], bool bForce)
{
    if (posCount == 0)
    {
        return false;
    }

    // TODO: Use this again (can be useful when tuning initGridCells()),
    // or remove throughout.
    GMX_UNUSED_VALUE(bForce);

    switch (pbc.pbcType)
    {
        case PbcType::No:
            bGridPBC_[XX] = false;
            bGridPBC_[YY] = false;
            bGridPBC_[ZZ] = false;
            break;
        case PbcType::XY:
            bGridPBC_[XX] = true;
            bGridPBC_[YY] = true;
            bGridPBC_[ZZ] = false;
            break;
        case PbcType::Xyz:
            bGridPBC_[XX] = true;
            bGridPBC_[YY] = true;
            bGridPBC_[ZZ] = true;
            break;
        default:
            // Grid searching not supported for now with screw.
            return false;
    }

    bool   bSingleCell[DIM] = { false, false, bXY_ };
    matrix box;
    copy_mat(pbc.box, box);
    // TODO: In principle, we could use the bounding box for periodic
    // dimensions as well if the bounding box is sufficiently far from the box
    // edges.
    rvec origin, boundingBoxSize;
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
            // TODO: It could be better to avoid this when determining the cell
            // size, but this can still remain here as a fallback to avoid
            // incorrect results.
            if (std::ceil(2 * cutoff_ * invCellSize_[dd]) >= ncelldim_[dd])
            {
                // Cutoff is too close to half the box size for grid searching
                // (it is not possible to find a single shift for every pair of
                // grid cells).
                return false;
            }
        }
    }
    if (bTric_)
    {
        cellShiftZY_ = box[ZZ][YY] * invCellSize_[YY];
        cellShiftZX_ = box[ZZ][XX] * invCellSize_[XX];
        cellShiftYX_ = box[YY][XX] * invCellSize_[XX];
    }
    return true;
}

void AnalysisNeighborhoodSearchImpl::mapPointToGridCell(const rvec x, rvec cell, rvec xout) const
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
    GMX_ASSERT(cell[XX] >= 0 && cell[XX] < ncelldim_[XX], "Grid cell X index out of range");
    GMX_ASSERT(cell[YY] >= 0 && cell[YY] < ncelldim_[YY], "Grid cell Y index out of range");
    GMX_ASSERT(cell[ZZ] >= 0 && cell[ZZ] < ncelldim_[ZZ], "Grid cell Z index out of range");
    return cell[XX] + cell[YY] * ncelldim_[XX] + cell[ZZ] * ncelldim_[XX] * ncelldim_[YY];
}

int AnalysisNeighborhoodSearchImpl::getGridCellIndex(const rvec cell) const
{
    ivec icell;
    for (int dd = 0; dd < DIM; ++dd)
    {
        int cellIndex = static_cast<int>(std::floor(cell[dd]));
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
    return getGridCellIndex(icell);
}

void AnalysisNeighborhoodSearchImpl::addToGridCell(const rvec cell, int i)
{
    const int ci = getGridCellIndex(cell);
    cells_[ci].push_back(i);
}

void AnalysisNeighborhoodSearchImpl::initCellRange(const rvec centerCell, ivec currCell, ivec upperBound, int dim) const
{
    RVec shiftedCenter(centerCell);
    // Shift the center to the cell coordinates of currCell, so that
    // computeCutoffExtent() can assume simple rectangular grid.
    if (bTric_)
    {
        if (dim == XX)
        {
            if (currCell[ZZ] < 0)
            {
                shiftedCenter[XX] += cellShiftZX_;
            }
            else if (currCell[ZZ] >= ncelldim_[ZZ])
            {
                shiftedCenter[XX] -= cellShiftZX_;
            }
            if (currCell[YY] < 0)
            {
                shiftedCenter[XX] += cellShiftYX_;
            }
            else if (currCell[YY] >= ncelldim_[YY])
            {
                shiftedCenter[XX] -= cellShiftYX_;
            }
        }
        if (dim == XX || dim == YY)
        {
            if (currCell[ZZ] < 0)
            {
                shiftedCenter[YY] += cellShiftZY_;
            }
            else if (currCell[ZZ] >= ncelldim_[ZZ])
            {
                shiftedCenter[YY] -= cellShiftZY_;
            }
        }
    }
    const real range       = computeCutoffExtent(shiftedCenter, currCell, dim) * invCellSize_[dim];
    real       startOffset = shiftedCenter[dim] - range;
    real       endOffset   = shiftedCenter[dim] + range;
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
    currCell[dim]   = static_cast<int>(std::floor(startOffset));
    upperBound[dim] = static_cast<int>(std::floor(endOffset));
}

real AnalysisNeighborhoodSearchImpl::computeCutoffExtent(const RVec centerCell, const ivec cell, int dim) const
{
    if (dim == ZZ)
    {
        return cutoff_;
    }

    real dist2 = 0;
    for (int d = dim + 1; d < DIM; ++d)
    {
        real dimDist = cell[d] - centerCell[d];
        if (dimDist < -1)
        {
            dimDist += 1;
        }
        else if (dimDist <= 0)
        {
            continue;
        }
        dist2 += dimDist * dimDist * cellSize_[d] * cellSize_[d];
    }
    if (dist2 >= cutoff2_)
    {
        return 0;
    }
    return std::sqrt(cutoff2_ - dist2);
}

bool AnalysisNeighborhoodSearchImpl::nextCell(const rvec centerCell, ivec cell, ivec upperBound) const
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

void AnalysisNeighborhoodSearchImpl::init(AnalysisNeighborhood::SearchMode     mode,
                                          bool                                 bXY,
                                          const ListOfLists<int>*              excls,
                                          const t_pbc*                         pbc,
                                          const AnalysisNeighborhoodPositions& positions)
{
    GMX_RELEASE_ASSERT(positions.index_ == -1,
                       "Individual indexed positions not supported as reference");
    bXY_ = bXY;
    if (bXY_ && pbc != nullptr && pbc->pbcType != PbcType::No)
    {
        if (pbc->pbcType != PbcType::XY && pbc->pbcType != PbcType::Xyz)
        {
            std::string message = formatString(
                    "Computations in the XY plane are not supported with PBC type '%s'",
                    c_pbcTypeNames[pbc->pbcType].c_str());
            GMX_THROW(NotImplementedError(message));
        }
        if (pbc->pbcType == PbcType::Xyz
            && (std::fabs(pbc->box[ZZ][XX]) > GMX_REAL_EPS * pbc->box[ZZ][ZZ]
                || std::fabs(pbc->box[ZZ][YY]) > GMX_REAL_EPS * pbc->box[ZZ][ZZ]))
        {
            GMX_THROW(
                    NotImplementedError("Computations in the XY plane are not supported when the "
                                        "last box vector is not parallel to the Z axis"));
        }
        // Use a single grid cell in Z direction.
        matrix box;
        copy_mat(pbc->box, box);
        clear_rvec(box[ZZ]);
        set_pbc(&pbc_, PbcType::XY, box);
    }
    else if (pbc != nullptr)
    {
        pbc_ = *pbc;
    }
    else
    {
        pbc_.pbcType = PbcType::No;
        clear_mat(pbc_.box);
    }
    nref_ = positions.count_;
    if (mode == AnalysisNeighborhood::eSearchMode_Simple)
    {
        bGrid_ = false;
    }
    else if (bTryGrid_)
    {
        bGrid_ = initGrid(
                pbc_, positions.count_, positions.x_, mode == AnalysisNeighborhood::eSearchMode_Grid);
    }
    refIndices_ = positions.indices_;
    if (bGrid_)
    {
        xrefAlloc_.resize(nref_);
        xref_ = as_rvec_array(xrefAlloc_.data());

        for (int i = 0; i < nref_; ++i)
        {
            const int ii = (refIndices_ != nullptr) ? refIndices_[i] : i;
            rvec      refcell;
            mapPointToGridCell(positions.x_[ii], refcell, xrefAlloc_[i]);
            addToGridCell(refcell, i);
        }
    }
    else if (refIndices_ != nullptr)
    {
        xrefAlloc_.resize(nref_);
        xref_ = as_rvec_array(xrefAlloc_.data());
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
    refExclusionIds_ = nullptr;
    if (excls != nullptr)
    {
        // TODO: Check that the IDs are ascending, or remove the limitation.
        refExclusionIds_ = positions.exclusionIds_;
        GMX_RELEASE_ASSERT(refExclusionIds_ != nullptr,
                           "Exclusion IDs must be set for reference positions "
                           "when exclusions are enabled");
    }
}

/********************************************************************
 * AnalysisNeighborhoodPairSearchImpl
 */

void AnalysisNeighborhoodPairSearchImpl::reset(int testIndex)
{
    testIndex_     = testIndex;
    testCellIndex_ = -1;
    previ_         = -1;
    prevr2_        = 0.0;
    clear_rvec(prevdx_);
    exclind_ = 0;
    prevcai_ = -1;
    if (testIndex_ >= 0 && testIndex_ < testPosCount_)
    {
        const int index = (testIndices_ != nullptr ? testIndices_[testIndex] : testIndex);
        if (search_.bGrid_)
        {
            search_.mapPointToGridCell(testPositions_[index], testcell_, xtest_);
            search_.initCellRange(testcell_, currCell_, cellBound_, ZZ);
            search_.initCellRange(testcell_, currCell_, cellBound_, YY);
            search_.initCellRange(testcell_, currCell_, cellBound_, XX);
            if (selfSearchMode_)
            {
                testCellIndex_ = search_.getGridCellIndex(testcell_);
            }
        }
        else
        {
            copy_rvec(testPositions_[index], xtest_);
            if (selfSearchMode_)
            {
                previ_ = testIndex_;
            }
        }
        if (search_.excls_ != nullptr)
        {
            const int exclIndex = testExclusionIds_[index];
            if (exclIndex < search_.excls_->ssize())
            {
                excl_ = (*search_.excls_)[exclIndex];
            }
            else
            {
                excl_ = ArrayRef<const int>();
            }
        }
    }
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
    const int nexcl = excl_.ssize();
    if (exclind_ < nexcl)
    {
        const int index = (search_.refIndices_ != nullptr ? search_.refIndices_[j] : j);
        const int refId = search_.refExclusionIds_[index];
        while (exclind_ < nexcl && excl_[exclind_] < refId)
        {
            ++exclind_;
        }
        if (exclind_ < nexcl && refId == excl_[exclind_])
        {
            ++exclind_;
            return true;
        }
    }
    return false;
}

void AnalysisNeighborhoodPairSearchImpl::startSearch(const AnalysisNeighborhoodPositions& positions)
{
    selfSearchMode_   = false;
    testPosCount_     = positions.count_;
    testPositions_    = positions.x_;
    testExclusionIds_ = positions.exclusionIds_;
    testIndices_      = positions.indices_;
    GMX_RELEASE_ASSERT(search_.excls_ == nullptr || testExclusionIds_ != nullptr,
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

void AnalysisNeighborhoodPairSearchImpl::startSelfSearch()
{
    selfSearchMode_   = true;
    testPosCount_     = search_.nref_;
    testPositions_    = search_.xref_;
    testExclusionIds_ = search_.refExclusionIds_;
    testIndices_      = search_.refIndices_;
    GMX_RELEASE_ASSERT(search_.excls_ == nullptr || testIndices_ == nullptr,
                       "Exclusion IDs not implemented with indexed ref positions");
    reset(0);
}

template<class Action>
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
                const int ci = search_.shiftCell(currCell_, shift);
                if (selfSearchMode_ && ci > testCellIndex_)
                {
                    continue;
                }
                const int cellSize = gmx::ssize(search_.cells_[ci]);
                for (; cai < cellSize; ++cai)
                {
                    const int i = search_.cells_[ci][cai];
                    if (selfSearchMode_ && ci == testCellIndex_ && i >= testIndex_)
                    {
                        continue;
                    }
                    if (isExcluded(i))
                    {
                        continue;
                    }
                    rvec dx;
                    rvec_sub(search_.xref_[i], xtest_, dx);
                    rvec_sub(dx, shift, dx);
                    const real r2 = search_.bXY_ ? dx[XX] * dx[XX] + dx[YY] * dx[YY] : norm2(dx);
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
            } while (search_.nextCell(testcell_, currCell_, cellBound_));
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
                if (search_.pbc_.pbcType != PbcType::No)
                {
                    pbc_dx(&search_.pbc_, search_.xref_[i], xtest_, dx);
                }
                else
                {
                    rvec_sub(search_.xref_[i], xtest_, dx);
                }
                const real r2 = search_.bXY_ ? dx[XX] * dx[XX] + dx[YY] * dx[YY] : norm2(dx);
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

void AnalysisNeighborhoodPairSearchImpl::initFoundPair(AnalysisNeighborhoodPair* pair) const
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

} // namespace internal

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
    MindistAction(int* closestPoint, real* minDist2, rvec* dx) // NOLINT(readability-non-const-parameter)
        :
        closestPoint_(*closestPoint), minDist2_(*minDist2), dx_(*dx)
    {
    }
    //! Copies the action.
    MindistAction(const MindistAction&) = default;

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
    int&  closestPoint_;
    real& minDist2_;
    rvec& dx_;

    GMX_DISALLOW_ASSIGN(MindistAction);
};

} // namespace

/********************************************************************
 * AnalysisNeighborhood::Impl
 */

class AnalysisNeighborhood::Impl
{
public:
    typedef AnalysisNeighborhoodSearch::ImplPointer SearchImplPointer;
    typedef std::vector<SearchImplPointer>          SearchList;

    Impl() : cutoff_(0), excls_(nullptr), mode_(eSearchMode_Automatic), bXY_(false) {}
    ~Impl()
    {
        SearchList::const_iterator i;
        for (i = searchList_.begin(); i != searchList_.end(); ++i)
        {
            GMX_RELEASE_ASSERT(i->use_count() == 1, "Dangling AnalysisNeighborhoodSearch reference");
        }
    }

    SearchImplPointer getSearch();

    std::mutex              createSearchMutex_;
    SearchList              searchList_;
    real                    cutoff_;
    const ListOfLists<int>* excls_;
    SearchMode              mode_;
    bool                    bXY_;
};

AnalysisNeighborhood::Impl::SearchImplPointer AnalysisNeighborhood::Impl::getSearch()
{
    std::lock_guard<std::mutex> lock(createSearchMutex_);
    // TODO: Consider whether this needs to/can be faster, e.g., by keeping a
    // separate pool of unused search objects.
    SearchList::const_iterator i;
    for (i = searchList_.begin(); i != searchList_.end(); ++i)
    {
        if (i->use_count() == 1)
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

AnalysisNeighborhood::AnalysisNeighborhood() : impl_(new Impl) {}

AnalysisNeighborhood::~AnalysisNeighborhood() {}

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

void AnalysisNeighborhood::setTopologyExclusions(const ListOfLists<int>* excls)
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

AnalysisNeighborhoodSearch AnalysisNeighborhood::initSearch(const t_pbc* pbc,
                                                            const AnalysisNeighborhoodPositions& positions)
{
    Impl::SearchImplPointer search(impl_->getSearch());
    search->init(mode(), impl_->bXY_, impl_->excls_, pbc, positions);
    return AnalysisNeighborhoodSearch(search);
}

/********************************************************************
 * AnalysisNeighborhoodSearch
 */

AnalysisNeighborhoodSearch::AnalysisNeighborhoodSearch() {}

AnalysisNeighborhoodSearch::AnalysisNeighborhoodSearch(const ImplPointer& impl) : impl_(impl) {}

void AnalysisNeighborhoodSearch::reset()
{
    impl_.reset();
}

AnalysisNeighborhood::SearchMode AnalysisNeighborhoodSearch::mode() const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    return (impl_->usesGridSearch() ? AnalysisNeighborhood::eSearchMode_Grid
                                    : AnalysisNeighborhood::eSearchMode_Simple);
}

bool AnalysisNeighborhoodSearch::isWithin(const AnalysisNeighborhoodPositions& positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    return pairSearch.searchNext(&withinAction);
}

real AnalysisNeighborhoodSearch::minimumDistance(const AnalysisNeighborhoodPositions& positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    real          minDist2     = impl_->cutoffSquared();
    int           closestPoint = -1;
    rvec          dx           = { 0.0, 0.0, 0.0 };
    MindistAction action(&closestPoint, &minDist2, &dx);
    (void)pairSearch.searchNext(action);
    return std::sqrt(minDist2);
}

AnalysisNeighborhoodPair AnalysisNeighborhoodSearch::nearestPoint(const AnalysisNeighborhoodPositions& positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    internal::AnalysisNeighborhoodPairSearchImpl pairSearch(*impl_);
    pairSearch.startSearch(positions);
    real          minDist2     = impl_->cutoffSquared();
    int           closestPoint = -1;
    rvec          dx           = { 0.0, 0.0, 0.0 };
    MindistAction action(&closestPoint, &minDist2, &dx);
    (void)pairSearch.searchNext(action);
    return AnalysisNeighborhoodPair(closestPoint, 0, minDist2, dx);
}

AnalysisNeighborhoodPairSearch AnalysisNeighborhoodSearch::startSelfPairSearch() const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    Impl::PairSearchImplPointer pairSearch(impl_->getPairSearch());
    pairSearch->startSelfSearch();
    return AnalysisNeighborhoodPairSearch(pairSearch);
}

AnalysisNeighborhoodPairSearch
AnalysisNeighborhoodSearch::startPairSearch(const AnalysisNeighborhoodPositions& positions) const
{
    GMX_RELEASE_ASSERT(impl_, "Accessing an invalid search object");
    Impl::PairSearchImplPointer pairSearch(impl_->getPairSearch());
    pairSearch->startSearch(positions);
    return AnalysisNeighborhoodPairSearch(pairSearch);
}

/********************************************************************
 * AnalysisNeighborhoodPairSearch
 */

AnalysisNeighborhoodPairSearch::AnalysisNeighborhoodPairSearch(const ImplPointer& impl) :
    impl_(impl)
{
}

bool AnalysisNeighborhoodPairSearch::findNextPair(AnalysisNeighborhoodPair* pair)
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
