/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * \todo
 * The grid implementation could still be optimized in several different ways:
 *   - Triclinic grid cells are not the most efficient shape, but make PBC
 *     handling easier.
 *   - Precalculating the required PBC shift for a pair of cells outside the
 *     inner loop. After this is done, it should be quite straightforward to
 *     move to rectangular cells.
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
#include "gromacs/selection/nbsearch.h"

#include <math.h>

#include <algorithm>
#include <vector>

#include "thread_mpi/mutex.h"

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"

#include "gromacs/selection/position.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

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
         * \param[in]     mode      Search mode to use.
         * \param[in]     pbc       PBC information.
         * \param[in]     positions Set of reference positions.
         */
        void init(AnalysisNeighborhood::SearchMode     mode,
                  const t_pbc                         *pbc,
                  const AnalysisNeighborhoodPositions &positions);
        PairSearchImplPointer getPairSearch();

        real cutoffSquared() const { return cutoff2_; }
        bool usesGridSearch() const { return bGrid_; }

    private:
        //! Calculates offsets to neighboring grid cells that should be considered.
        void initGridCellNeighborList();
        /*! \brief
         * Determines a suitable grid size and sets up the cells.
         *
         * \param[in]     pbc  Information about the box.
         * \returns  false if grid search is not suitable.
         */
        bool initGridCells(const t_pbc *pbc);
        /*! \brief
         * Sets ua a search grid for a given box.
         *
         * \param[in]     pbc  Information about the box.
         * \returns  false if grid search is not suitable.
         */
        bool initGrid(const t_pbc *pbc);
        /*! \brief
         * Maps a point into a grid cell.
         *
         * \param[in]  x    Point to map.
         * \param[out] cell Indices of the grid cell in which \p x lies.
         * \param[out] xout Coordinates to use
         *     (will be within the triclinic unit cell).
         */
        void mapPointToGridCell(const rvec x, ivec cell, rvec xout) const;
        /*! \brief
         * Calculates linear index of a grid cell.
         *
         * \param[in]  cell Cell indices.
         * \returns    Linear index of \p cell.
         */
        int getGridCellIndex(const ivec cell) const;
        /*! \brief
         * Adds an index into a grid cell.
         *
         * \param[in]     cell Cell into which \p i should be added.
         * \param[in]     i    Index to add.
         */
        void addToGridCell(const ivec cell, int i);

        //! Whether to try grid searching.
        bool                    bTryGrid_;
        //! The cutoff.
        real                    cutoff_;
        //! The cutoff squared.
        real                    cutoff2_;

        //! Number of reference points for the current frame.
        int                     nref_;
        //! Reference point positions.
        const rvec             *xref_;
        //! Reference position ids (NULL if not available).
        const int              *refid_;
        //! PBC data.
        t_pbc                  *pbc_;

        //! Number of excluded reference positions for current test particle.
        int                     nexcl_;
        //! Exclusions for current test particle.
        int                    *excl_;

        //! Whether grid searching is actually used for the current positions.
        bool                    bGrid_;
        //! Array allocated for storing in-unit-cell reference positions.
        rvec                   *xref_alloc_;
        //! Allocation count for xref_alloc.
        int                     xref_nalloc_;
        //! false if the box is rectangular.
        bool                    bTric_;
        //! Box vectors of a single grid cell.
        matrix                  cellbox_;
        //! The reciprocal cell vectors as columns; the inverse of \p cellbox.
        matrix                  recipcell_;
        //! Number of cells along each dimension.
        ivec                    ncelldim_;
        //! Data structure to hold the grid cell contents.
        CellList                cells_;
        //! Number of neighboring cells to consider.
        int                     ngridnb_;
        //! Offsets of the neighboring cells to consider.
        ivec                   *gnboffs_;
        //! Allocation count for \p gnboffs.
        int                     gnboffs_nalloc_;

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
            clear_rvec(xtest_);
            clear_ivec(testcell_);
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
        //! Reference to the test positions.
        ConstArrayRef<rvec>                     testPositions_;
        //! Index of the currently active test position in \p testPositions_.
        int                                     testIndex_;
        //! Stores test position during a pair loop.
        rvec                                    xtest_;
        //! Stores the previous returned position during a pair loop.
        int                                     previ_;
        //! Stores the current exclusion index during loops.
        int                                     exclind_;
        //! Stores the test particle cell index during loops.
        ivec                                    testcell_;
        //! Stores the current cell neighbor index during pair loops.
        int                                     prevnbi_;
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
        cutoff_     = GMX_REAL_MAX;
        bTryGrid_   = false;
    }
    cutoff2_        = sqr(cutoff_);

    nref_           = 0;
    xref_           = NULL;
    refid_          = NULL;
    pbc_            = NULL;

    nexcl_          = 0;
    excl_           = NULL;

    bGrid_          = false;

    xref_alloc_     = NULL;
    xref_nalloc_    = 0;
    bTric_          = false;
    clear_mat(cellbox_);
    clear_mat(recipcell_);
    clear_ivec(ncelldim_);

    ngridnb_        = 0;
    gnboffs_        = NULL;
    gnboffs_nalloc_ = 0;
}

AnalysisNeighborhoodSearchImpl::~AnalysisNeighborhoodSearchImpl()
{
    PairSearchList::const_iterator i;
    for (i = pairSearchList_.begin(); i != pairSearchList_.end(); ++i)
    {
        GMX_RELEASE_ASSERT(i->unique(),
                           "Dangling AnalysisNeighborhoodPairSearch reference");
    }
    sfree(xref_alloc_);
    sfree(gnboffs_);
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

void AnalysisNeighborhoodSearchImpl::initGridCellNeighborList()
{
    int   maxx, maxy, maxz;
    real  rvnorm;

    /* Find the extent of the sphere in triclinic coordinates */
    maxz   = static_cast<int>(cutoff_ * recipcell_[ZZ][ZZ]) + 1;
    rvnorm = sqrt(sqr(recipcell_[YY][YY]) + sqr(recipcell_[ZZ][YY]));
    maxy   = static_cast<int>(cutoff_ * rvnorm) + 1;
    rvnorm = sqrt(sqr(recipcell_[XX][XX]) + sqr(recipcell_[YY][XX])
                  + sqr(recipcell_[ZZ][XX]));
    maxx   = static_cast<int>(cutoff_ * rvnorm) + 1;

    /* Calculate the number of cells and reallocate if necessary */
    ngridnb_ = (2 * maxx + 1) * (2 * maxy + 1) * (2 * maxz + 1);
    if (gnboffs_nalloc_ < ngridnb_)
    {
        gnboffs_nalloc_ = ngridnb_;
        srenew(gnboffs_, gnboffs_nalloc_);
    }

    /* Store the whole cube */
    /* TODO: Prune off corners that are not needed */
    int i = 0;
    for (int x = -maxx; x <= maxx; ++x)
    {
        for (int y = -maxy; y <= maxy; ++y)
        {
            for (int z = -maxz; z <= maxz; ++z)
            {
                gnboffs_[i][XX] = x;
                gnboffs_[i][YY] = y;
                gnboffs_[i][ZZ] = z;
                ++i;
            }
        }
    }
}

bool AnalysisNeighborhoodSearchImpl::initGridCells(const t_pbc *pbc)
{
    const real targetsize =
        pow(pbc->box[XX][XX] * pbc->box[YY][YY] * pbc->box[ZZ][ZZ]
            * 10 / nref_, static_cast<real>(1./3.));

    int cellCount = 1;
    for (int dd = 0; dd < DIM; ++dd)
    {
        ncelldim_[dd] = static_cast<int>(pbc->box[dd][dd] / targetsize);
        cellCount    *= ncelldim_[dd];
        if (ncelldim_[dd] < 3)
        {
            return false;
        }
    }
    // Never decrease the size of the cell vector to avoid reallocating
    // memory for the nested vectors.  The actual size of the vector is not
    // used outside this function.
    if (cells_.size() < static_cast<size_t>(cellCount))
    {
        cells_.resize(cellCount);
    }
    for (int ci = 0; ci < cellCount; ++ci)
    {
        cells_[ci].clear();
    }
    return true;
}

bool AnalysisNeighborhoodSearchImpl::initGrid(const t_pbc *pbc)
{
    /* TODO: This check could be improved. */
    if (0.5*pbc->max_cutoff2 < cutoff2_)
    {
        return false;
    }

    if (!initGridCells(pbc))
    {
        return false;
    }

    bTric_ = TRICLINIC(pbc->box);
    if (bTric_)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            svmul(1.0 / ncelldim_[dd], pbc->box[dd], cellbox_[dd]);
        }
        m_inv_ur0(cellbox_, recipcell_);
    }
    else
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            cellbox_[dd][dd]   = pbc->box[dd][dd] / ncelldim_[dd];
            recipcell_[dd][dd] = 1.0 / cellbox_[dd][dd];
        }
    }
    initGridCellNeighborList();
    return true;
}

void AnalysisNeighborhoodSearchImpl::mapPointToGridCell(const rvec x,
                                                        ivec       cell,
                                                        rvec       xout) const
{
    rvec xtmp;
    copy_rvec(x, xtmp);
    if (bTric_)
    {
        rvec tx;
        tmvmul_ur0(recipcell_, xtmp, tx);
        for (int dd = 0; dd < DIM; ++dd)
        {
            const int cellCount = ncelldim_[dd];
            int       cellIndex = static_cast<int>(floor(tx[dd]));
            while (cellIndex < 0)
            {
                cellIndex += cellCount;
                rvec_add(xtmp, pbc_->box[dd], xtmp);
            }
            while (cellIndex >= cellCount)
            {
                cellIndex -= cellCount;
                rvec_sub(xtmp, pbc_->box[dd], xtmp);
            }
            cell[dd] = cellIndex;
        }
    }
    else
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            const int cellCount = ncelldim_[dd];
            int       cellIndex = static_cast<int>(floor(xtmp[dd] * recipcell_[dd][dd]));
            while (cellIndex < 0)
            {
                cellIndex += cellCount;
                xtmp[dd]  += pbc_->box[dd][dd];
            }
            while (cellIndex >= cellCount)
            {
                cellIndex -= cellCount;
                xtmp[dd]  -= pbc_->box[dd][dd];
            }
            cell[dd] = cellIndex;
        }
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

void AnalysisNeighborhoodSearchImpl::addToGridCell(const ivec cell, int i)
{
    const int ci = getGridCellIndex(cell);
    cells_[ci].push_back(i);
}

void AnalysisNeighborhoodSearchImpl::init(
        AnalysisNeighborhood::SearchMode     mode,
        const t_pbc                         *pbc,
        const AnalysisNeighborhoodPositions &positions)
{
    GMX_RELEASE_ASSERT(positions.index_ == -1,
                       "Individual indexed positions not supported as reference");
    pbc_  = const_cast<t_pbc *>(pbc);
    nref_ = positions.count_;
    // TODO: Consider whether it would be possible to support grid searching in
    // more cases.
    if (mode == AnalysisNeighborhood::eSearchMode_Simple
        || pbc_ == NULL || pbc_->ePBC != epbcXYZ)
    {
        bGrid_ = false;
    }
    else if (bTryGrid_)
    {
        // TODO: Actually implement forcing eSearchMode_Grid
        bGrid_ = initGrid(pbc_);
    }
    if (bGrid_)
    {
        if (xref_nalloc_ < nref_)
        {
            srenew(xref_alloc_, nref_);
            xref_nalloc_ = nref_;
        }
        xref_ = xref_alloc_;

        for (int i = 0; i < nref_; ++i)
        {
            ivec refcell;
            mapPointToGridCell(positions.x_[i], refcell, xref_alloc_[i]);
            addToGridCell(refcell, i);
        }
    }
    else
    {
        xref_ = positions.x_;
    }
    // TODO: Once exclusions are supported, this may need to be initialized.
    refid_ = NULL;
}

#if 0
/*! \brief
 * Sets the exclusions for the next neighborhood search.
 *
 * \param[in,out] d     Neighborhood search data structure.
 * \param[in]     nexcl Number of reference positions to exclude from next
 *      search.
 * \param[in]     excl  Indices of reference positions to exclude.
 *
 * The set exclusions remain in effect until the next call of this function.
 */
void
gmx_ana_nbsearch_set_excl(gmx_ana_nbsearch_t *d, int nexcl, int excl[])
{

    d->nexcl = nexcl;
    d->excl  = excl;
}
#endif

/********************************************************************
 * AnalysisNeighborhoodPairSearchImpl
 */

void AnalysisNeighborhoodPairSearchImpl::reset(int testIndex)
{
    testIndex_ = testIndex;
    if (testIndex_ >= 0 && testIndex_ < static_cast<int>(testPositions_.size()))
    {
        if (search_.bGrid_)
        {
            search_.mapPointToGridCell(testPositions_[testIndex], testcell_, xtest_);
        }
        else
        {
            copy_rvec(testPositions_[testIndex_], xtest_);
        }
    }
    previ_     = -1;
    exclind_   = 0;
    prevnbi_   = 0;
    prevcai_   = -1;
}

void AnalysisNeighborhoodPairSearchImpl::nextTestPosition()
{
    if (testIndex_ < static_cast<int>(testPositions_.size()))
    {
        ++testIndex_;
        reset(testIndex_);
    }
}

bool AnalysisNeighborhoodPairSearchImpl::isExcluded(int j)
{
    if (exclind_ < search_.nexcl_)
    {
        if (search_.refid_)
        {
            while (exclind_ < search_.nexcl_
                   && search_.excl_[exclind_] < search_.refid_[j])
            {
                ++exclind_;
            }
            if (exclind_ < search_.nexcl_
                && search_.refid_[j] == search_.excl_[exclind_])
            {
                ++exclind_;
                return true;
            }
        }
        else
        {
            while (search_.bGrid_ && exclind_ < search_.nexcl_
                   && search_.excl_[exclind_] < j)
            {
                ++exclind_;
            }
            if (search_.excl_[exclind_] == j)
            {
                ++exclind_;
                return true;
            }
        }
    }
    return false;
}

void AnalysisNeighborhoodPairSearchImpl::startSearch(
        const AnalysisNeighborhoodPositions &positions)
{
    if (positions.index_ < 0)
    {
        testPositions_ = constArrayRefFromArray<rvec>(positions.x_, positions.count_);
        reset(0);
    }
    else
    {
        // Somewhat of a hack: setup the array such that only the last position
        // will be used.
        testPositions_ = constArrayRefFromArray<rvec>(positions.x_, positions.index_ + 1);
        reset(positions.index_);
    }
}

template <class Action>
bool AnalysisNeighborhoodPairSearchImpl::searchNext(Action action)
{
    while (testIndex_ < static_cast<int>(testPositions_.size()))
    {
        if (search_.bGrid_)
        {
            int nbi = prevnbi_;
            int cai = prevcai_ + 1;

            for (; nbi < search_.ngridnb_; ++nbi)
            {
                ivec cell;

                ivec_add(testcell_, search_.gnboffs_[nbi], cell);
                cell[XX] = (cell[XX] + search_.ncelldim_[XX]) % search_.ncelldim_[XX];
                cell[YY] = (cell[YY] + search_.ncelldim_[YY]) % search_.ncelldim_[YY];
                cell[ZZ] = (cell[ZZ] + search_.ncelldim_[ZZ]) % search_.ncelldim_[ZZ];

                const int ci       = search_.getGridCellIndex(cell);
                const int cellSize = static_cast<int>(search_.cells_[ci].size());
                /* TODO: Calculate the required PBC shift outside the inner loop */
                for (; cai < cellSize; ++cai)
                {
                    const int i = search_.cells_[ci][cai];
                    if (isExcluded(i))
                    {
                        continue;
                    }
                    rvec       dx;
                    pbc_dx_aiuc(search_.pbc_, xtest_, search_.xref_[i], dx);
                    const real r2 = norm2(dx);
                    if (r2 <= search_.cutoff2_)
                    {
                        if (action(i, r2))
                        {
                            prevnbi_ = nbi;
                            prevcai_ = cai;
                            previ_   = i;
                            return true;
                        }
                    }
                }
                exclind_ = 0;
                cai      = 0;
            }
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
                if (search_.pbc_)
                {
                    pbc_dx(search_.pbc_, xtest_, search_.xref_[i], dx);
                }
                else
                {
                    rvec_sub(xtest_, search_.xref_[i], dx);
                }
                const real r2 = norm2(dx);
                if (r2 <= search_.cutoff2_)
                {
                    if (action(i, r2))
                    {
                        previ_ = i;
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
        *pair = AnalysisNeighborhoodPair(previ_, testIndex_);
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
bool withinAction(int /*i*/, real /*r2*/)
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
         *
         * The constructor call does not modify the pointed values, but only
         * stores the pointers for later use.
         * See the class description for additional semantics.
         */
        MindistAction(int *closestPoint, real *minDist2)
            : closestPoint_(*closestPoint), minDist2_(*minDist2)
        {
        }

        //! Processes a neighbor to find the nearest point.
        bool operator()(int i, real r2)
        {
            if (r2 < minDist2_)
            {
                closestPoint_ = i;
                minDist2_     = r2;
            }
            return false;
        }

    private:
        int     &closestPoint_;
        real    &minDist2_;

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

        Impl() : cutoff_(0), mode_(eSearchMode_Automatic)
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
        SearchMode              mode_;
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
    search->init(mode(), pbc, positions);
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
    MindistAction action(&closestPoint, &minDist2);
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
    MindistAction action(&closestPoint, &minDist2);
    (void)pairSearch.searchNext(action);
    return AnalysisNeighborhoodPair(closestPoint, 0);
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
