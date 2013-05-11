/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * See the class documentation for usage.
 *
 * The classes within this file can be used independently of the other parts
 * of the library.
 * The library also uses the classes internally.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_NBSEARCH_H
#define GMX_SELECTION_NBSEARCH_H

#include <boost/shared_ptr.hpp>

#include "../legacyheaders/typedefs.h"
#include "../utility/common.h"
#include "../utility/gmxassert.h"

#include "indexutil.h"

struct gmx_ana_pos_t;

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
 * Neighborhood searching for analysis tools.
 *
 * This class implements neighborhood searching routines for analysis tools.
 * The emphasis is in flexibility and ease of use; one main driver is to have
 * a common implementation of grid-based searching to avoid replicating this in
 * multiple tools (and to make more tools take advantage of the significant
 * performance improvement this allows).
 *
 * To use the search, create an object of this type, call setCutoff() to
 * initialize it, and then repeatedly call initSearch() to start a search with
 * different sets of reference positions.  For each set of reference positions,
 * use methods in the returned AnalysisNeighborhoodSearch to find the reference
 * positions that are within the given cutoff from a provided position.
 *
 * initSearch() is thread-safe and can be called from multiple threads.  Each
 * call returns a different instance of the search object that can be used
 * independently of the others.  However, the returned search objects can only
 * be used within a single thread.  It is also possible to create multiple
 * concurrent searches within a single thread.
 *
 * \todo
 * Support for exclusions.
 * The 4.5/4.6 C API had very low-level support for exclusions, which was not
 * very convenient to use, and hadn't been tested much.  The internal code that
 * it used to do the exclusion during the search itself is still there, but it
 * needs more thought on what would be a convenient way to initialize it.
 * Can be implemented once there is need for it in some calling code.
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
         * Set cutoff distance for the neighborhood searching.
         *
         * \param[in]  cutoff Cutoff distance for the search
         *   (<=0 stands for no cutoff).
         *
         * Currently, cannot be called after initSearch() has been called.
         * If this method is not called, no cutoff is used in the searches.
         *
         * Does not throw.
         */
        void setCutoff(real cutoff);
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
         * \param[in] pbc PBC information for the frame.
         * \param[in] n   Number of reference positions for the frame.
         * \param[in] x   \p n reference positions for the frame.
         * \returns   Search object that can be used to find positions from
         *      \p x within the given cutoff.
         * \throws    std::bad_alloc if out of memory.
         */
        AnalysisNeighborhoodSearch
        initSearch(const t_pbc *pbc, int n, const rvec x[]);
        /*! \brief
         * Initializes neighborhood search for a set of positions.
         *
         * \param[in] pbc PBC information for the frame.
         * \param[in] p   Reference positions for the frame.
         * \returns   Search object that can be used to find positions from
         *      \p p within the given cutoff.
         * \throws    std::bad_alloc if out of memory.
         */
        AnalysisNeighborhoodSearch
        initSearch(const t_pbc *pbc, const gmx_ana_pos_t *p);

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
        AnalysisNeighborhoodPair() : refIndex_(-1), testIndex_(0) {}
        //! Initializes a pair object with the given data.
        AnalysisNeighborhoodPair(int refIndex, int testIndex)
            : refIndex_(refIndex), testIndex_(testIndex)
        {
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

    private:
        int                     refIndex_;
        int                     testIndex_;
};

/*! \brief
 * Initialized neighborhood search with a fixed set of reference positions.
 *
 * An instance of this class is obtained through
 * AnalysisNeighborhood::initSearch(), and can be used to do multiple searches
 * against the provided set of reference positions.
 * Currently, it is not possible to call startPairSearch() while a previous
 * AnalysisNeighborhoodPairSearch object obtained from startPairSearch() of the
 * same instance still exists.
 *
 * This class works like a pointer: copies of it point to the same search.
 * In general, avoid creating copies, and only use the copy/assignment support
 * for moving the variable around.  With C++11, this class would best be
 * movable.
 *
 * Methods in this class do not throw.
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
        //! Internal short-hand type for a pointer to the implementation class.
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
        //! Creates a shallow copy (copy references the same implementation).
        AnalysisNeighborhoodSearch(const AnalysisNeighborhoodSearch &other);
        ~AnalysisNeighborhoodSearch()
        {
            reset();
        }

        //! Assign a shallow copy (copy references the same implementation).
        AnalysisNeighborhoodSearch &operator=(const AnalysisNeighborhoodSearch &other);

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
         * \param[in] x  Test position.
         * \returns   true if the test position is within the cutoff of any
         *     reference position.
         */
        bool isWithin(const rvec x) const;
        /*! \brief
         * Check whether a point is within a neighborhood.
         *
         * \param[in] p  Test positions.
         * \param[in] i  Use the i'th position in \p p for testing.
         * \returns   true if the test position is within the cutoff of any
         *     reference position.
         */
        bool isWithin(const gmx_ana_pos_t *p, int i) const;
        /*! \brief
         * Calculates the minimum distance from the reference points.
         *
         * \param[in] x  Test position.
         * \returns   The distance to the nearest reference position, or the
         *     cutoff value if there are no reference positions within the
         *     cutoff.
         */
        real minimumDistance(const rvec x) const;
        /*! \brief
         * Calculates the minimum distance from the reference points.
         *
         * \param[in] p  Test positions.
         * \param[in] i  Use the i'th position in \p p for testing.
         * \returns   The distance to the nearest reference position, or the
         *     cutoff value if there are no reference positions within the
         *     cutoff.
         */
        real minimumDistance(const gmx_ana_pos_t *p, int i) const;
        /*! \brief
         * Finds the closest reference point.
         *
         * \param[in] x  Test position.
         * \returns   The reference index identifies the reference position
         *     that is closest to the test position.
         *     The test index is always zero.  The returned pair is invalid if
         *     no reference position is within the cutoff.
         */
        AnalysisNeighborhoodPair nearestPoint(const rvec x) const;
        /*! \brief
         * Finds the closest reference point.
         *
         * \param[in] p  Test positions.
         * \param[in] i  Use the i'th position in \p p for testing.
         * \returns   The reference index identifies the reference position
         *     that is closest to the test position.
         *     The test index is always \p i.  The returned pair is invalid if
         *     no reference position is within the cutoff.
         */
        AnalysisNeighborhoodPair nearestPoint(const gmx_ana_pos_t *p, int i) const;

        /*! \brief
         * Start a search to find reference positions within a cutoff.
         *
         * \param[in] x  Test position to search the neighbors for.
         * \returns   Initialized search object to loop through all reference
         *     positions within the configured cutoff.
         *
         * In the AnalysisNeighborhoodPair objects returned by the search, the
         * test index is always zero.
         */
        AnalysisNeighborhoodPairSearch startPairSearch(const rvec x);
        /*! \brief
         * Start a search to find reference positions within a cutoff.
         *
         * \param[in] p  Test positions.
         * \param[in] i  Use the i'th position in \p p for testing.
         * \returns   Initialized search object to loop through all reference
         *     positions within the configured cutoff.
         *
         * In the AnalysisNeighborhoodPair objects returned by the search, the
         * test index is always \p i.
         */
        AnalysisNeighborhoodPairSearch startPairSearch(const gmx_ana_pos_t *p, int i);

    private:
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
   gmx::AnalysisNeighborhoodSearch search = nb.initSearch(pbc, nref, xref);
   gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(x);
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
        //! Internal short-hand type for a pointer to the implementation class.
        typedef boost::shared_ptr<internal::AnalysisNeighborhoodPairSearchImpl>
            ImplPointer;

        /*! \brief
         * Internally initialize the search.
         *
         * Used to implement AnalysisNeighborhoodSearch::startPairSearch().
         * Cannot be called from user code.
         */
        explicit AnalysisNeighborhoodPairSearch(const ImplPointer &impl);
        //! Creates a shallow copy (copy references the same implementation).
        AnalysisNeighborhoodPairSearch(const AnalysisNeighborhoodPairSearch &other);
        ~AnalysisNeighborhoodPairSearch();

        //! Assign a shallow copy (copy references the same implementation).
        AnalysisNeighborhoodPairSearch &operator=(const AnalysisNeighborhoodPairSearch &other);

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

    private:
        ImplPointer             impl_;
};

} // namespace gmx

#endif
