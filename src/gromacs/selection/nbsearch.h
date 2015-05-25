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
/*! \file
 * \brief API for neighborhood searching for analysis.
 *
 * The main part of the API is the class gmx::AnalysisNeighborhood.
 * See \ref page_analysisnbsearch for an overview.
 *
 * The classes within this file can be used independently of the other parts
 * of the selection module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_NBSEARCH_H
#define GMX_SELECTION_NBSEARCH_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

struct t_blocka;
struct t_pbc;

namespace gmx
{

namespace internal
{
class AnalysisNeighborhoodSearchImpl;
class AnalysisNeighborhoodPairSearchImpl;
};

class AnalysisNeighborhoodSearch;
class AnalysisNeighborhoodPairSearch;

/*! \brief
 * Input positions for neighborhood searching.
 *
 * This class supports uniformly specifying sets of positions for various
 * methods in the analysis neighborhood searching classes
 * (AnalysisNeighborhood and AnalysisNeighborhoodSearch).
 *
 * Note that copies are not made: only a reference to the positions passed to
 * the constructors are kept.  The caller is responsible to ensure that those
 * positions remain in scope as long as the neighborhood search object requires
 * access to them.
 *
 * Also note that in addition to constructors here, Selection and
 * SelectionPosition provide conversions operators to this type.  It is done
 * this way to not introduce a cyclic dependency between the selection code and
 * the neighborhood search code, which in turn allows splitting this search
 * code into a separate lower-level module if desired at some point.
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class AnalysisNeighborhoodPositions
{
    public:
        /*! \brief
         * Initializes positions from a single position vector.
         *
         * For positions initialized this way, AnalysisNeighborhoodPair always
         * returns zero in the corresponding index.
         *
         * This constructor is not explicit to allow directly passing an rvec
         * to methods that accept positions.
         */
        AnalysisNeighborhoodPositions(const rvec &x)
            : count_(1), index_(-1), x_(&x), exclusionIds_(NULL), indices_(NULL)
        {
        }
        /*! \brief
         * Initializes positions from an array of position vectors.
         */
        AnalysisNeighborhoodPositions(const rvec x[], int count)
            : count_(count), index_(-1), x_(x), exclusionIds_(NULL), indices_(NULL)
        {
        }
        /*! \brief
         * Initializes positions from a vector of position vectors.
         */
        AnalysisNeighborhoodPositions(const std::vector<RVec> &x)
            : count_(x.size()), index_(-1), x_(as_rvec_array(&x[0])),
              exclusionIds_(NULL), indices_(NULL)
        {
        }

        /*! \brief
         * Sets indices to use for mapping exclusions to these positions.
         *
         * The exclusion IDs can always be set, but they are ignored unless
         * actual exclusions have been set with
         * AnalysisNeighborhood::setTopologyExclusions().
         */
        AnalysisNeighborhoodPositions &
        exclusionIds(ConstArrayRef<int> ids)
        {
            GMX_ASSERT(static_cast<int>(ids.size()) == count_,
                       "Exclusion id array should match the number of positions");
            exclusionIds_ = ids.data();
            return *this;
        }
        /*! \brief
         * Sets indices that select a subset of all positions from the array.
         *
         * If called, selected positions from the array of positions passed to
         * the constructor is used instead of the whole array.
         * All returned indices from AnalysisNeighborhoodPair objects are
         * indices to the \p indices array passed here.
         */
        AnalysisNeighborhoodPositions &
        indexed(ConstArrayRef<int> indices)
        {
            count_   = indices.size();
            indices_ = indices.data();
            return *this;
        }

        /*! \brief
         * Selects a single position to use from an array.
         *
         * If called, a single position from the array of positions passed to
         * the constructor is used instead of the whole array.
         * In contrast to the AnalysisNeighborhoodPositions(const rvec &)
         * constructor, AnalysisNeighborhoodPair objects return \p index
         * instead of zero.
         *
         * If used together with indexed(), \p index references the index array
         * passed to indexed() instead of the position array.
         */
        AnalysisNeighborhoodPositions &selectSingleFromArray(int index)
        {
            GMX_ASSERT(index >= 0 && index < count_, "Invalid position index");
            index_ = index;
            return *this;
        }

    private:
        int                     count_;
        int                     index_;
        const rvec             *x_;
        const int              *exclusionIds_;
        const int              *indices_;

        //! To access the positions for initialization.
        friend class internal::AnalysisNeighborhoodSearchImpl;
        //! To access the positions for initialization.
        friend class internal::AnalysisNeighborhoodPairSearchImpl;
};

/*! \brief
 * Neighborhood searching for analysis tools.
 *
 * See \ref page_analysisnbsearch for an overview.
 *
 * To use the search, create an object of this type, call setCutoff() to
 * initialize it, and then repeatedly call initSearch() to start a search with
 * different sets of reference positions.  For each set of reference positions,
 * use methods in the returned AnalysisNeighborhoodSearch to find the reference
 * positions that are within the given cutoff from a provided position.
 *
 * initSearch() is thread-safe and can be called from multiple threads.  Each
 * call returns a different instance of the search object that can be used
 * independently of the others.  The returned AnalysisNeighborhoodSearch
 * objects are also thread-safe, and can be used concurrently from multiple
 * threads.  It is also possible to create multiple concurrent searches within
 * a single thread.
 *
 * \todo
 * Generalize the exclusion machinery to make it easier to use for other cases
 * than atom-atom exclusions from the topology.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class AnalysisNeighborhood
{
    public:
        //! Searching algorithm to use.
        enum SearchMode
        {
            //! Select algorithm based on heuristic efficiency considerations.
            eSearchMode_Automatic,
            //! Use a simple loop over all pairs.
            eSearchMode_Simple,
            //! Use grid-based searching whenever possible.
            eSearchMode_Grid
        };

        //! Creates an uninitialized neighborhood search.
        AnalysisNeighborhood();
        ~AnalysisNeighborhood();

        /*! \brief
         * Sets cutoff distance for the neighborhood searching.
         *
         * \param[in]  cutoff Cutoff distance for the search
         *   (<=0 stands for no cutoff).
         *
         * Currently, can only be called before the first call to initSearch().
         * If this method is not called, no cutoff is used in the searches.
         *
         * Does not throw.
         */
        void setCutoff(real cutoff);
        /*! \brief
         * Sets the search to only happen in the XY plane.
         *
         * Z component of the coordinates is not used in the searching,
         * and returned distances are computed in the XY plane.
         * Only boxes with the third box vector parallel to the Z axis are
         * currently implemented.
         *
         * Does not throw.
         */
        void setXYMode(bool bXY);
        /*! \brief
         * Sets atom exclusions from a topology.
         *
         * The \p excls structure specifies the exclusions from test positions
         * to reference positions, i.e., a block starting at `excls->index[i]`
         * specifies the exclusions for test position `i`, and the indices in
         * `excls->a` are indices of the reference positions.  If `excls->nr`
         * is smaller than a test position id, then such test positions do not
         * have any exclusions.
         * It is assumed that the indices within a block of indices in
         * `excls->a` is ascending.
         *
         * Does not throw.
         *
         * \see AnalysisNeighborhoodPositions::exclusionIds()
         */
        void setTopologyExclusions(const t_blocka *excls);
        /*! \brief
         * Sets the algorithm to use for searching.
         *
         * \param[in] mode  Search mode to use.
         *
         * Note that if \p mode is \ref eSearchMode_Grid, it is still only a
         * suggestion: grid-based searching may not be possible with the
         * provided input, in which case a simple search is still used.
         * This is mainly useful for testing purposes to force a mode.
         *
         * Does not throw.
         */
        void setMode(SearchMode mode);
        //! Returns the currently active search mode.
        SearchMode mode() const;

        /*! \brief
         * Initializes neighborhood search for a set of positions.
         *
         * \param[in] pbc        PBC information for the frame.
         * \param[in] positions  Set of reference positions to use.
         * \returns   Search object that can be used to find positions from
         *      \p x within the given cutoff.
         * \throws    std::bad_alloc if out of memory.
         *
         * Currently, the input positions cannot use
         * AnalysisNeighborhoodPositions::selectSingleFromArray().
         */
        AnalysisNeighborhoodSearch
        initSearch(const t_pbc                         *pbc,
                   const AnalysisNeighborhoodPositions &positions);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \brief
 * Value type to represent a pair of positions found in neighborhood searching.
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class AnalysisNeighborhoodPair
{
    public:
        //! Initializes an invalid pair.
        AnalysisNeighborhoodPair() : refIndex_(-1), testIndex_(0), distance2_(0.0)
        {
            clear_rvec(dx_);
        }
        //! Initializes a pair object with the given data.
        AnalysisNeighborhoodPair(int refIndex, int testIndex, real distance2,
                                 const rvec dx)
            : refIndex_(refIndex), testIndex_(testIndex), distance2_(distance2)
        {
            copy_rvec(dx, dx_);
        }

        /*! \brief
         * Whether this pair is valid.
         *
         * If isValid() returns false, other methods should not be called.
         */
        bool isValid() const { return refIndex_ >= 0; }

        /*! \brief
         * Returns the index of the reference position in the pair.
         *
         * This index is always the index into the position array provided to
         * AnalysisNeighborhood::initSearch().
         */
        int refIndex() const
        {
            GMX_ASSERT(isValid(), "Accessing invalid object");
            return refIndex_;
        }
        /*! \brief
         * Returns the index of the test position in the pair.
         *
         * The contents of this index depends on the context (method call) that
         * produces the pair.
         * If there was no array in the call, this index is zero.
         */
        int testIndex() const
        {
            GMX_ASSERT(isValid(), "Accessing invalid object");
            return testIndex_;
        }
        /*! \brief
         * Returns the squared distance between the pair of positions.
         */
        real distance2() const
        {
            GMX_ASSERT(isValid(), "Accessing invalid object");
            return distance2_;
        }
        /*! \brief
         * Returns the shortest vector between the pair of positions.
         *
         * The vector is from the test position to the reference position.
         */
        const rvec &dx() const
        {
            GMX_ASSERT(isValid(), "Accessing invalid object");
            return dx_;
        }

    private:
        int                     refIndex_;
        int                     testIndex_;
        real                    distance2_;
        rvec                    dx_;
};

/*! \brief
 * Initialized neighborhood search with a fixed set of reference positions.
 *
 * An instance of this class is obtained through
 * AnalysisNeighborhood::initSearch(), and can be used to do multiple searches
 * against the provided set of reference positions.
 * It is possible to create concurrent pair searches (including from different
 * threads), as well as call other methods in this class while a pair search is
 * in progress.
 *
 * This class works like a pointer: copies of it point to the same search.
 * In general, avoid creating copies, and only use the copy/assignment support
 * for moving the variable around.  With C++11, this class would best be
 * movable.
 *
 * Methods in this class do not throw unless otherwise indicated.
 *
 * \todo
 * Make it such that reset() is not necessary to call in code that repeatedly
 * assigns the result of AnalysisNeighborhood::initSearch() to the same
 * variable (see sm_distance.cpp).
 *
 * \todo
 * Consider merging nearestPoint() and minimumDistance() by adding the distance
 * to AnalysisNeighborhoodPair.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class AnalysisNeighborhoodSearch
{
    public:
        /*! \brief
         * Internal short-hand type for a pointer to the implementation class.
         *
         * shared_ptr is used here to automatically keep a reference count to
         * track whether an implementation class is still used outside the
         * AnalysisNeighborhood object.  Ownership currently always stays with
         * AnalysisNeighborhood; it always keeps one instance of the pointer.
         */
        typedef boost::shared_ptr<internal::AnalysisNeighborhoodSearchImpl>
            ImplPointer;

        /*! \brief
         * Initializes an invalid search.
         *
         * Such an object cannot be used for searching.  It needs to be
         * assigned a value from AnalysisNeighborhood::initSearch() before it
         * can be used.  Provided to allow declaring a variable to hold the
         * search before calling AnalysisNeighborhood::initSearch().
         */
        AnalysisNeighborhoodSearch();
        /*! \brief
         * Internally initialize the search.
         *
         * Used to implement AnalysisNeighborhood::initSearch().
         * Cannot be called from user code.
         */
        explicit AnalysisNeighborhoodSearch(const ImplPointer &impl);

        /*! \brief
         * Clears this search.
         *
         * Equivalent to \c "*this = AnalysisNeighborhoodSearch();".
         * Currently, this is necessary to avoid unnecessary memory allocation
         * if the previous search variable is still in scope when you want to
         * call AnalysisNeighborhood::initSearch() again.
         */
        void reset();

        /*! \brief
         * Returns the searching algorithm that this search is using.
         *
         * The return value is never AnalysisNeighborhood::eSearchMode_Automatic.
         */
        AnalysisNeighborhood::SearchMode mode() const;

        /*! \brief
         * Check whether a point is within a neighborhood.
         *
         * \param[in] positions  Set of test positions to use.
         * \returns   true if any of the test positions is within the cutoff of
         *     any reference position.
         */
        bool isWithin(const AnalysisNeighborhoodPositions &positions) const;
        /*! \brief
         * Calculates the minimum distance from the reference points.
         *
         * \param[in] positions  Set of test positions to use.
         * \returns   The distance to the nearest reference position, or the
         *     cutoff value if there are no reference positions within the
         *     cutoff.
         */
        real minimumDistance(const AnalysisNeighborhoodPositions &positions) const;
        /*! \brief
         * Finds the closest reference point.
         *
         * \param[in] positions  Set of test positions to use.
         * \returns   The reference index identifies the reference position
         *     that is closest to the test positions.
         *     The test index identifies the test position that is closest to
         *     the provided test position.  The returned pair is invalid if
         *     no reference position is within the cutoff.
         */
        AnalysisNeighborhoodPair
        nearestPoint(const AnalysisNeighborhoodPositions &positions) const;

        /*! \brief
         * Start a search to find reference positions within a cutoff.
         *
         * \param[in] positions  Set of test positions to use.
         * \returns   Initialized search object to loop through all reference
         *     positions within the configured cutoff.
         * \throws    std::bad_alloc if out of memory.
         */
        AnalysisNeighborhoodPairSearch
        startPairSearch(const AnalysisNeighborhoodPositions &positions) const;

    private:
        typedef internal::AnalysisNeighborhoodSearchImpl Impl;

        ImplPointer             impl_;
};

/*! \brief
 * Initialized neighborhood pair search with a fixed set of positions.
 *
 * This class is used to loop through pairs of neighbors within the cutoff
 * provided to AnalysisNeighborhood.  The following code demonstrates its use:
 * \code
   gmx::AnalysisNeighborhood       nb;
   nb.setCutoff(cutoff);
   gmx::AnalysisNeighborhoodPositions refPos(xref, nref);
   gmx::AnalysisNeighborhoodSearch search = nb.initSearch(pbc, refPos);
   gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(selection);
   gmx::AnalysisNeighborhoodPair pair;
   while (pairSearch.findNextPair(&pair))
   {
       // <do something for each found pair the information in pair>
   }
 * \endcode
 *
 * It is not possible to use a single search object from multiple threads
 * concurrently.
 *
 * This class works like a pointer: copies of it point to the same search.
 * In general, avoid creating copies, and only use the copy/assignment support
 * for moving the variable around.  With C++11, this class would best be
 * movable.
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class AnalysisNeighborhoodPairSearch
{
    public:
        /*! \brief
         * Internal short-hand type for a pointer to the implementation class.
         *
         * See AnalysisNeighborhoodSearch::ImplPointer for rationale of using
         * shared_ptr and ownership semantics.
         */
        typedef boost::shared_ptr<internal::AnalysisNeighborhoodPairSearchImpl>
            ImplPointer;

        /*! \brief
         * Internally initialize the search.
         *
         * Used to implement AnalysisNeighborhoodSearch::startPairSearch().
         * Cannot be called from user code.
         */
        explicit AnalysisNeighborhoodPairSearch(const ImplPointer &impl);

        /*! \brief
         * Finds the next pair within the cutoff.
         *
         * \param[out] pair  Information about the found pair.
         * \returns    false if there were no more pairs.
         *
         * If the method returns false, \p pair will be invalid.
         *
         * \see AnalysisNeighborhoodPair
         * \see AnalysisNeighborhoodSearch::startPairSearch()
         */
        bool findNextPair(AnalysisNeighborhoodPair *pair);
        /*! \brief
         * Skip remaining pairs for a test position in the search.
         *
         * When called after findNextPair(), makes subsequent calls to
         * findNextPair() skip any pairs that have the same test position as
         * that previously returned.
         * This is useful if the caller wants to search whether any reference
         * position within the cutoff satisfies some condition.  This method
         * can be used to skip remaining pairs after the first such position
         * has been found if the remaining pairs would not have an effect on
         * the outcome.
         */
        void skipRemainingPairsForTestPosition();

    private:
        ImplPointer             impl_;
};

} // namespace gmx

#endif
