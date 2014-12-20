/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Tests selection neighborhood searching.
 *
 * \todo
 * Increase coverage of these tests for different corner cases: other PBC cases
 * than full 3D, large cutoffs (larger than half the box size), etc.
 * At least some of these probably don't work correctly.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/nbsearch.h"

#include <cmath>

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * NeighborhoodSearchTestData
 */

class NeighborhoodSearchTestData
{
    public:
        struct RefPair
        {
            RefPair(int refIndex, real distance)
                : refIndex(refIndex), distance(distance), bFound(false),
                  bExcluded(false), bIndexed(true)
            {
            }

            bool operator<(const RefPair &other) const
            {
                return refIndex < other.refIndex;
            }

            int                 refIndex;
            real                distance;
            // The variables below are state variables that are only used
            // during the actual testing after creating a copy of the reference
            // pair list, not as part of the reference data.
            // Simpler to have just a single structure for both purposes.
            bool                bFound;
            bool                bExcluded;
            bool                bIndexed;
        };

        struct TestPosition
        {
            explicit TestPosition(const rvec x)
                : refMinDist(0.0), refNearestPoint(-1)
            {
                copy_rvec(x, this->x);
            }

            rvec                 x;
            real                 refMinDist;
            int                  refNearestPoint;
            std::vector<RefPair> refPairs;
        };

        typedef std::vector<TestPosition> TestPositionList;

        NeighborhoodSearchTestData(int seed, real cutoff);
        ~NeighborhoodSearchTestData();

        gmx::AnalysisNeighborhoodPositions refPositions() const
        {
            return gmx::AnalysisNeighborhoodPositions(refPos_);
        }
        gmx::AnalysisNeighborhoodPositions testPositions() const
        {
            if (testPos_.empty())
            {
                testPos_.reserve(testPositions_.size());
                for (size_t i = 0; i < testPositions_.size(); ++i)
                {
                    testPos_.push_back(testPositions_[i].x);
                }
            }
            return gmx::AnalysisNeighborhoodPositions(testPos_);
        }
        gmx::AnalysisNeighborhoodPositions testPosition(int index) const
        {
            return testPositions().selectSingleFromArray(index);
        }

        void addTestPosition(const rvec x)
        {
            GMX_RELEASE_ASSERT(testPos_.empty(),
                               "Cannot add positions after testPositions() call");
            testPositions_.push_back(TestPosition(x));
        }
        gmx::RVec generateRandomPosition();
        std::vector<int> generateIndex(int count) const;
        void generateRandomRefPositions(int count);
        void generateRandomTestPositions(int count);
        void computeReferences(t_pbc *pbc)
        {
            computeReferencesInternal(pbc, false);
        }
        void computeReferencesXY(t_pbc *pbc)
        {
            computeReferencesInternal(pbc, true);
        }

        bool containsPair(int testIndex, const RefPair &pair) const
        {
            const std::vector<RefPair>          &refPairs = testPositions_[testIndex].refPairs;
            std::vector<RefPair>::const_iterator foundRefPair
                = std::lower_bound(refPairs.begin(), refPairs.end(), pair);
            if (foundRefPair == refPairs.end() || foundRefPair->refIndex != pair.refIndex)
            {
                return false;
            }
            return true;
        }

        gmx_rng_t                        rng_;
        real                             cutoff_;
        matrix                           box_;
        t_pbc                            pbc_;
        int                              refPosCount_;
        std::vector<gmx::RVec>           refPos_;
        TestPositionList                 testPositions_;

    private:
        void computeReferencesInternal(t_pbc *pbc, bool bXY);

        mutable std::vector<gmx::RVec>   testPos_;
};

//! Shorthand for a collection of reference pairs.
typedef std::vector<NeighborhoodSearchTestData::RefPair> RefPairList;

NeighborhoodSearchTestData::NeighborhoodSearchTestData(int seed, real cutoff)
    : rng_(NULL), cutoff_(cutoff), refPosCount_(0)
{
    // TODO: Handle errors.
    rng_ = gmx_rng_init(seed);
    clear_mat(box_);
    set_pbc(&pbc_, epbcNONE, box_);
}

NeighborhoodSearchTestData::~NeighborhoodSearchTestData()
{
    if (rng_ != NULL)
    {
        gmx_rng_destroy(rng_);
    }
}

gmx::RVec NeighborhoodSearchTestData::generateRandomPosition()
{
    rvec fx, x;
    fx[XX] = gmx_rng_uniform_real(rng_);
    fx[YY] = gmx_rng_uniform_real(rng_);
    fx[ZZ] = gmx_rng_uniform_real(rng_);
    mvmul(box_, fx, x);
    // Add a small displacement to allow positions outside the box
    x[XX] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
    x[YY] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
    x[ZZ] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
    return x;
}

std::vector<int> NeighborhoodSearchTestData::generateIndex(int count) const
{
    std::vector<int> result;
    for (int i = 0; i < count; ++i)
    {
        if (gmx_rng_uniform_real(rng_) > 0.5)
        {
            result.push_back(i);
        }
    }
    return result;
}

void NeighborhoodSearchTestData::generateRandomRefPositions(int count)
{
    refPosCount_ = count;
    refPos_.reserve(count);
    for (int i = 0; i < count; ++i)
    {
        refPos_.push_back(generateRandomPosition());
    }
}

void NeighborhoodSearchTestData::generateRandomTestPositions(int count)
{
    testPositions_.reserve(count);
    for (int i = 0; i < count; ++i)
    {
        addTestPosition(generateRandomPosition());
    }
}

void NeighborhoodSearchTestData::computeReferencesInternal(t_pbc *pbc, bool bXY)
{
    real cutoff = cutoff_;
    if (cutoff <= 0)
    {
        cutoff = std::numeric_limits<real>::max();
    }
    TestPositionList::iterator i;
    for (i = testPositions_.begin(); i != testPositions_.end(); ++i)
    {
        i->refMinDist      = cutoff;
        i->refNearestPoint = -1;
        i->refPairs.clear();
        for (int j = 0; j < refPosCount_; ++j)
        {
            rvec dx;
            if (pbc != NULL)
            {
                pbc_dx(pbc, i->x, refPos_[j], dx);
            }
            else
            {
                rvec_sub(i->x, refPos_[j], dx);
            }
            // TODO: This may not work intuitively for 2D with the third box
            // vector not parallel to the Z axis, but neither does the actual
            // neighborhood search.
            const real dist =
                !bXY ? norm(dx) : sqrt(sqr(dx[XX]) + sqr(dx[YY]));
            if (dist < i->refMinDist)
            {
                i->refMinDist      = dist;
                i->refNearestPoint = j;
            }
            if (dist <= cutoff)
            {
                RefPair pair(j, dist);
                GMX_RELEASE_ASSERT(i->refPairs.empty() || i->refPairs.back() < pair,
                                   "Reference pairs should be generated in sorted order");
                i->refPairs.push_back(pair);
            }
        }
    }
}

/********************************************************************
 * ExclusionsHelper
 */

class ExclusionsHelper
{
    public:
        static void markExcludedPairs(RefPairList *refPairs, int testIndex,
                                      const t_blocka *excls);

        ExclusionsHelper(int refPosCount, int testPosCount);

        void generateExclusions();

        const t_blocka *exclusions() const { return &excls_; }

        gmx::ConstArrayRef<int> refPosIds() const
        {
            return gmx::constArrayRefFromVector<int>(exclusionIds_.begin(),
                                                     exclusionIds_.begin() + refPosCount_);
        }
        gmx::ConstArrayRef<int> testPosIds() const
        {
            return gmx::constArrayRefFromVector<int>(exclusionIds_.begin(),
                                                     exclusionIds_.begin() + testPosCount_);
        }

    private:
        int              refPosCount_;
        int              testPosCount_;
        std::vector<int> exclusionIds_;
        std::vector<int> exclsIndex_;
        std::vector<int> exclsAtoms_;
        t_blocka         excls_;
};

// static
void ExclusionsHelper::markExcludedPairs(RefPairList *refPairs, int testIndex,
                                         const t_blocka *excls)
{
    int count = 0;
    for (int i = excls->index[testIndex]; i < excls->index[testIndex + 1]; ++i)
    {
        const int                           excludedIndex = excls->a[i];
        NeighborhoodSearchTestData::RefPair searchPair(excludedIndex, 0.0);
        RefPairList::iterator               excludedRefPair
            = std::lower_bound(refPairs->begin(), refPairs->end(), searchPair);
        if (excludedRefPair != refPairs->end()
            && excludedRefPair->refIndex == excludedIndex)
        {
            excludedRefPair->bFound    = true;
            excludedRefPair->bExcluded = true;
            ++count;
        }
    }
}

ExclusionsHelper::ExclusionsHelper(int refPosCount, int testPosCount)
    : refPosCount_(refPosCount), testPosCount_(testPosCount)
{
    // Generate an array of 0, 1, 2, ...
    // TODO: Make the tests work also with non-trivial exclusion IDs,
    // and test that.
    exclusionIds_.resize(std::max(refPosCount, testPosCount), 1);
    exclusionIds_[0] = 0;
    std::partial_sum(exclusionIds_.begin(), exclusionIds_.end(),
                     exclusionIds_.begin());

    excls_.nr           = 0;
    excls_.index        = NULL;
    excls_.nra          = 0;
    excls_.a            = NULL;
    excls_.nalloc_index = 0;
    excls_.nalloc_a     = 0;
}

void ExclusionsHelper::generateExclusions()
{
    // TODO: Consider a better set of test data, where the density of the
    // particles would be higher, or where the exclusions would not be random,
    // to make a higher percentage of the exclusions to actually be within the
    // cutoff.
    exclsIndex_.reserve(testPosCount_ + 1);
    exclsAtoms_.reserve(testPosCount_ * 20);
    exclsIndex_.push_back(0);
    for (int i = 0; i < testPosCount_; ++i)
    {
        for (int j = 0; j < 20; ++j)
        {
            exclsAtoms_.push_back(i + j*3);
        }
        exclsIndex_.push_back(exclsAtoms_.size());
    }
    excls_.nr    = exclsIndex_.size();
    excls_.index = &exclsIndex_[0];
    excls_.nra   = exclsAtoms_.size();
    excls_.a     = &exclsAtoms_[0];
}

/********************************************************************
 * NeighborhoodSearchTest
 */

class NeighborhoodSearchTest : public ::testing::Test
{
    public:
        void testIsWithin(gmx::AnalysisNeighborhoodSearch  *search,
                          const NeighborhoodSearchTestData &data);
        void testMinimumDistance(gmx::AnalysisNeighborhoodSearch  *search,
                                 const NeighborhoodSearchTestData &data);
        void testNearestPoint(gmx::AnalysisNeighborhoodSearch  *search,
                              const NeighborhoodSearchTestData &data);
        void testPairSearch(gmx::AnalysisNeighborhoodSearch  *search,
                            const NeighborhoodSearchTestData &data);
        void testPairSearchIndexed(gmx::AnalysisNeighborhood        *nb,
                                   const NeighborhoodSearchTestData &data);
        void testPairSearchFull(gmx::AnalysisNeighborhoodSearch          *search,
                                const NeighborhoodSearchTestData         &data,
                                const gmx::AnalysisNeighborhoodPositions &pos,
                                const t_blocka                           *excls,
                                const gmx::ConstArrayRef<int>            &refIndices,
                                const gmx::ConstArrayRef<int>            &testIndices);

        gmx::AnalysisNeighborhood        nb_;
};

void NeighborhoodSearchTest::testIsWithin(
        gmx::AnalysisNeighborhoodSearch  *search,
        const NeighborhoodSearchTestData &data)
{
    NeighborhoodSearchTestData::TestPositionList::const_iterator i;
    for (i = data.testPositions_.begin(); i != data.testPositions_.end(); ++i)
    {
        const bool bWithin = (i->refMinDist <= data.cutoff_);
        EXPECT_EQ(bWithin, search->isWithin(i->x))
        << "Distance is " << i->refMinDist;
    }
}

void NeighborhoodSearchTest::testMinimumDistance(
        gmx::AnalysisNeighborhoodSearch  *search,
        const NeighborhoodSearchTestData &data)
{
    NeighborhoodSearchTestData::TestPositionList::const_iterator i;
    for (i = data.testPositions_.begin(); i != data.testPositions_.end(); ++i)
    {
        const real refDist = i->refMinDist;
        EXPECT_REAL_EQ_TOL(refDist, search->minimumDistance(i->x),
                           gmx::test::ulpTolerance(20));
    }
}

void NeighborhoodSearchTest::testNearestPoint(
        gmx::AnalysisNeighborhoodSearch  *search,
        const NeighborhoodSearchTestData &data)
{
    NeighborhoodSearchTestData::TestPositionList::const_iterator i;
    for (i = data.testPositions_.begin(); i != data.testPositions_.end(); ++i)
    {
        const gmx::AnalysisNeighborhoodPair pair = search->nearestPoint(i->x);
        if (pair.isValid())
        {
            EXPECT_EQ(i->refNearestPoint, pair.refIndex());
            EXPECT_EQ(0, pair.testIndex());
            EXPECT_REAL_EQ_TOL(i->refMinDist, sqrt(pair.distance2()),
                               gmx::test::ulpTolerance(64));
        }
        else
        {
            EXPECT_EQ(i->refNearestPoint, -1);
        }
    }
}

//! Helper function for formatting test failure messages.
std::string formatVector(const rvec x)
{
    return gmx::formatString("[%.3f, %.3f, %.3f]", x[XX], x[YY], x[ZZ]);
}

/*! \brief
 * Helper function to check that all expected pairs were found.
 */
void checkAllPairsFound(const RefPairList &refPairs,
                        const std::vector<gmx::RVec> &refPos,
                        int testPosIndex, const rvec testPos)
{
    // This could be elegantly expressed with Google Mock matchers, but that
    // has a significant effect on the runtime of the tests...
    int                         count = 0;
    RefPairList::const_iterator first;
    for (RefPairList::const_iterator i = refPairs.begin(); i != refPairs.end(); ++i)
    {
        if (!i->bFound)
        {
            ++count;
            first = i;
        }
    }
    if (count > 0)
    {
        ADD_FAILURE()
        << "Some pairs (" << count << "/" << refPairs.size() << ") "
        << "within the cutoff were not found. First pair:\n"
        << " Ref: " << first->refIndex << " at "
        << formatVector(refPos[first->refIndex]) << "\n"
        << "Test: " << testPosIndex << " at " << formatVector(testPos) << "\n"
        << "Dist: " << first->distance;
    }
}

void NeighborhoodSearchTest::testPairSearch(
        gmx::AnalysisNeighborhoodSearch  *search,
        const NeighborhoodSearchTestData &data)
{
    testPairSearchFull(search, data, data.testPositions(), NULL,
                       gmx::EmptyArrayRef(), gmx::EmptyArrayRef());
}

void NeighborhoodSearchTest::testPairSearchIndexed(
        gmx::AnalysisNeighborhood        *nb,
        const NeighborhoodSearchTestData &data)
{
    std::vector<int>                refIndices(data.generateIndex(data.refPos_.size()));
    std::vector<int>                testIndices(data.generateIndex(data.testPositions_.size()));
    gmx::AnalysisNeighborhoodSearch search =
        nb->initSearch(&data.pbc_,
                       data.refPositions().indexed(refIndices));
    testPairSearchFull(&search, data, data.testPositions(), NULL,
                       refIndices, testIndices);
}

void NeighborhoodSearchTest::testPairSearchFull(
        gmx::AnalysisNeighborhoodSearch          *search,
        const NeighborhoodSearchTestData         &data,
        const gmx::AnalysisNeighborhoodPositions &pos,
        const t_blocka                           *excls,
        const gmx::ConstArrayRef<int>            &refIndices,
        const gmx::ConstArrayRef<int>            &testIndices)
{
    // TODO: Some parts of this code do not work properly if pos does not
    // initially contain all the test positions.
    std::set<int> remainingTestPositions;
    gmx::AnalysisNeighborhoodPositions  posCopy(pos);
    if (testIndices.empty())
    {
        for (size_t i = 0; i < data.testPositions_.size(); ++i)
        {
            remainingTestPositions.insert(i);
        }
    }
    else
    {
        remainingTestPositions.insert(testIndices.begin(), testIndices.end());
        posCopy.indexed(testIndices);
    }
    gmx::AnalysisNeighborhoodPairSearch pairSearch
        = search->startPairSearch(posCopy);
    gmx::AnalysisNeighborhoodPair       pair;
    // There is an ordering assumption here that all pairs for a test position
    // are returned consencutively; with the current optimizations in the
    // search code, this is reasoable, as the set of grid cell pairs searched
    // depends on the test position.
    RefPairList refPairs;
    int         prevTestPos = -1;
    while (pairSearch.findNextPair(&pair))
    {
        const int testIndex =
            (testIndices.empty() ? pair.testIndex() : testIndices[pair.testIndex()]);
        const int refIndex =
            (refIndices.empty() ? pair.refIndex() : refIndices[pair.refIndex()]);
        if (testIndex != prevTestPos)
        {
            if (prevTestPos != -1)
            {
                checkAllPairsFound(refPairs, data.refPos_, prevTestPos,
                                   data.testPositions_[prevTestPos].x);
            }
            if (remainingTestPositions.count(testIndex) == 0)
            {
                ADD_FAILURE()
                << "Pairs for test position " << testIndex
                << " are returned more than once.";
            }
            remainingTestPositions.erase(testIndex);
            refPairs = data.testPositions_[testIndex].refPairs;
            if (excls != NULL)
            {
                ExclusionsHelper::markExcludedPairs(&refPairs, testIndex, excls);
            }
            if (!refIndices.empty())
            {
                RefPairList::iterator refPair;
                for (refPair = refPairs.begin(); refPair != refPairs.end(); ++refPair)
                {
                    refPair->bIndexed = false;
                }
                for (size_t i = 0; i < refIndices.size(); ++i)
                {
                    NeighborhoodSearchTestData::RefPair searchPair(refIndices[i], 0.0);
                    refPair = std::lower_bound(refPairs.begin(), refPairs.end(), searchPair);
                    if (refPair != refPairs.end() && refPair->refIndex == refIndices[i])
                    {
                        refPair->bIndexed = true;
                    }
                }
                for (refPair = refPairs.begin(); refPair != refPairs.end(); ++refPair)
                {
                    if (!refPair->bIndexed)
                    {
                        refPair->bFound = true;
                    }
                }
            }
            prevTestPos = testIndex;
        }

        NeighborhoodSearchTestData::RefPair searchPair(refIndex,
                                                       sqrt(pair.distance2()));
        RefPairList::iterator               foundRefPair
            = std::lower_bound(refPairs.begin(), refPairs.end(), searchPair);
        if (foundRefPair == refPairs.end() || foundRefPair->refIndex != refIndex)
        {
            ADD_FAILURE()
            << "Expected: Pair (ref: " << refIndex << ", test: " << testIndex
            << ") is not within the cutoff.\n"
            << "  Actual: It is returned.";
        }
        else if (foundRefPair->bExcluded)
        {
            ADD_FAILURE()
            << "Expected: Pair (ref: " << refIndex << ", test: " << testIndex
            << ") is excluded from the search.\n"
            << "  Actual: It is returned.";
        }
        else if (!foundRefPair->bIndexed)
        {
            ADD_FAILURE()
            << "Expected: Pair (ref: " << refIndex << ", test: " << testIndex
            << ") is not part of the indexed set.\n"
            << "  Actual: It is returned.";
        }
        else if (foundRefPair->bFound)
        {
            ADD_FAILURE()
            << "Expected: Pair (ref: " << refIndex << ", test: " << testIndex
            << ") is returned only once.\n"
            << "  Actual: It is returned multiple times.";
            return;
        }
        else
        {
            foundRefPair->bFound = true;
            EXPECT_REAL_EQ_TOL(foundRefPair->distance, searchPair.distance,
                               gmx::test::ulpTolerance(64))
            << "Distance computed by the neighborhood search does not match.";
        }
    }
    checkAllPairsFound(refPairs, data.refPos_, prevTestPos,
                       data.testPositions_[prevTestPos].x);
    for (std::set<int>::const_iterator i = remainingTestPositions.begin();
         i != remainingTestPositions.end(); ++i)
    {
        if (!data.testPositions_[*i].refPairs.empty())
        {
            ADD_FAILURE()
            << "Expected: Pairs would be returned for test position " << *i << ".\n"
            << "  Actual: None were returned.";
            break;
        }
    }
}

/********************************************************************
 * Test data generation
 */

class TrivialTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static TrivialTestData singleton;
            return singleton.data_;
        }

        TrivialTestData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 5.0;
            data_.box_[YY][YY] = 5.0;
            data_.box_[ZZ][ZZ] = 5.0;
            data_.generateRandomRefPositions(10);
            data_.generateRandomTestPositions(5);
            set_pbc(&data_.pbc_, epbcXYZ, data_.box_);
            data_.computeReferences(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class TrivialNoPBCTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static TrivialNoPBCTestData singleton;
            return singleton.data_;
        }

        TrivialNoPBCTestData() : data_(12345, 1.0)
        {
            data_.generateRandomRefPositions(10);
            data_.generateRandomTestPositions(5);
            data_.computeReferences(NULL);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class RandomBoxFullPBCData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBoxFullPBCData singleton;
            return singleton.data_;
        }

        RandomBoxFullPBCData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 10.0;
            data_.box_[YY][YY] = 5.0;
            data_.box_[ZZ][ZZ] = 7.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            data_.generateRandomRefPositions(1000);
            data_.generateRandomTestPositions(100);
            set_pbc(&data_.pbc_, epbcXYZ, data_.box_);
            data_.computeReferences(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class RandomBoxXYFullPBCData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBoxXYFullPBCData singleton;
            return singleton.data_;
        }

        RandomBoxXYFullPBCData() : data_(54321, 1.0)
        {
            data_.box_[XX][XX] = 10.0;
            data_.box_[YY][YY] = 5.0;
            data_.box_[ZZ][ZZ] = 7.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            data_.generateRandomRefPositions(1000);
            data_.generateRandomTestPositions(100);
            set_pbc(&data_.pbc_, epbcXYZ, data_.box_);
            data_.computeReferencesXY(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class RandomTriclinicFullPBCData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomTriclinicFullPBCData singleton;
            return singleton.data_;
        }

        RandomTriclinicFullPBCData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 5.0;
            data_.box_[YY][XX] = 2.5;
            data_.box_[YY][YY] = 2.5*sqrt(3.0);
            data_.box_[ZZ][XX] = 2.5;
            data_.box_[ZZ][YY] = 2.5*sqrt(1.0/3.0);
            data_.box_[ZZ][ZZ] = 5.0*sqrt(2.0/3.0);
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            data_.generateRandomRefPositions(1000);
            data_.generateRandomTestPositions(100);
            set_pbc(&data_.pbc_, epbcXYZ, data_.box_);
            data_.computeReferences(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class RandomBox2DPBCData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBox2DPBCData singleton;
            return singleton.data_;
        }

        RandomBox2DPBCData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 10.0;
            data_.box_[YY][YY] = 7.0;
            data_.box_[ZZ][ZZ] = 5.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            data_.generateRandomRefPositions(1000);
            data_.generateRandomTestPositions(100);
            set_pbc(&data_.pbc_, epbcXY, data_.box_);
            data_.computeReferences(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class RandomBoxNoPBCData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBoxNoPBCData singleton;
            return singleton.data_;
        }

        RandomBoxNoPBCData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 10.0;
            data_.box_[YY][YY] = 5.0;
            data_.box_[ZZ][ZZ] = 7.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            data_.generateRandomRefPositions(1000);
            data_.generateRandomTestPositions(100);
            set_pbc(&data_.pbc_, epbcNONE, data_.box_);
            data_.computeReferences(NULL);
        }

    private:
        NeighborhoodSearchTestData data_;
};

/********************************************************************
 * Actual tests
 */

TEST_F(NeighborhoodSearchTest, SimpleSearch)
{
    const NeighborhoodSearchTestData &data = RandomBoxFullPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Simple);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);

    search.reset();
    testPairSearchIndexed(&nb_, data);
}

TEST_F(NeighborhoodSearchTest, SimpleSearchXY)
{
    const NeighborhoodSearchTestData &data = RandomBoxXYFullPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Simple);
    nb_.setXYMode(true);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, GridSearchBox)
{
    const NeighborhoodSearchTestData &data = RandomBoxFullPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);

    search.reset();
    testPairSearchIndexed(&nb_, data);
}

TEST_F(NeighborhoodSearchTest, GridSearchTriclinic)
{
    const NeighborhoodSearchTestData &data = RandomTriclinicFullPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, GridSearch2DPBC)
{
    const NeighborhoodSearchTestData &data = RandomBox2DPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, GridSearchNoPBC)
{
    const NeighborhoodSearchTestData &data = RandomBoxNoPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, GridSearchXYBox)
{
    const NeighborhoodSearchTestData &data = RandomBoxXYFullPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    nb_.setXYMode(true);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, HandlesConcurrentSearches)
{
    const NeighborhoodSearchTestData &data = TrivialTestData::get();

    nb_.setCutoff(data.cutoff_);
    gmx::AnalysisNeighborhoodSearch search1 =
        nb_.initSearch(&data.pbc_, data.refPositions());
    gmx::AnalysisNeighborhoodSearch search2 =
        nb_.initSearch(&data.pbc_, data.refPositions());

    gmx::AnalysisNeighborhoodPairSearch pairSearch1 =
        search1.startPairSearch(data.testPosition(0));
    gmx::AnalysisNeighborhoodPairSearch pairSearch2 =
        search1.startPairSearch(data.testPosition(1));

    testPairSearch(&search2, data);

    gmx::AnalysisNeighborhoodPair pair;
    ASSERT_TRUE(pairSearch1.findNextPair(&pair))
    << "Test data did not contain any pairs for position 0 (problem in the test).";
    EXPECT_EQ(0, pair.testIndex());
    {
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), sqrt(pair.distance2()));
        EXPECT_TRUE(data.containsPair(0, searchPair));
    }

    ASSERT_TRUE(pairSearch2.findNextPair(&pair))
    << "Test data did not contain any pairs for position 1 (problem in the test).";
    EXPECT_EQ(1, pair.testIndex());
    {
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), sqrt(pair.distance2()));
        EXPECT_TRUE(data.containsPair(1, searchPair));
    }
}

TEST_F(NeighborhoodSearchTest, HandlesNoPBC)
{
    const NeighborhoodSearchTestData &data = TrivialNoPBCTestData::get();

    nb_.setCutoff(data.cutoff_);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, HandlesNullPBC)
{
    const NeighborhoodSearchTestData &data = TrivialNoPBCTestData::get();

    nb_.setCutoff(data.cutoff_);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(NULL, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, HandlesSkippingPairs)
{
    const NeighborhoodSearchTestData &data = TrivialTestData::get();

    nb_.setCutoff(data.cutoff_);
    gmx::AnalysisNeighborhoodSearch     search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    gmx::AnalysisNeighborhoodPairSearch pairSearch =
        search.startPairSearch(data.testPositions());
    gmx::AnalysisNeighborhoodPair       pair;
    // TODO: This test needs to be adjusted if the grid search gets optimized
    // to loop over the test positions in cell order (first, the ordering
    // assumption here breaks, and second, it then needs to be tested
    // separately for simple and grid searches).
    int currentIndex = 0;
    while (pairSearch.findNextPair(&pair))
    {
        while (currentIndex < pair.testIndex())
        {
            ++currentIndex;
        }
        EXPECT_EQ(currentIndex, pair.testIndex());
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), sqrt(pair.distance2()));
        EXPECT_TRUE(data.containsPair(currentIndex, searchPair));
        pairSearch.skipRemainingPairsForTestPosition();
        ++currentIndex;
    }
}

TEST_F(NeighborhoodSearchTest, SimpleSearchExclusions)
{
    const NeighborhoodSearchTestData &data = RandomBoxFullPBCData::get();

    ExclusionsHelper                  helper(data.refPosCount_, data.testPositions_.size());
    helper.generateExclusions();

    nb_.setCutoff(data.cutoff_);
    nb_.setTopologyExclusions(helper.exclusions());
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Simple);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_,
                       data.refPositions().exclusionIds(helper.refPosIds()));
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testPairSearchFull(&search, data,
                       data.testPositions().exclusionIds(helper.testPosIds()),
                       helper.exclusions(), gmx::EmptyArrayRef(), gmx::EmptyArrayRef());
}

TEST_F(NeighborhoodSearchTest, GridSearchExclusions)
{
    const NeighborhoodSearchTestData &data = RandomBoxFullPBCData::get();

    ExclusionsHelper                  helper(data.refPosCount_, data.testPositions_.size());
    helper.generateExclusions();

    nb_.setCutoff(data.cutoff_);
    nb_.setTopologyExclusions(helper.exclusions());
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_,
                       data.refPositions().exclusionIds(helper.refPosIds()));
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testPairSearchFull(&search, data,
                       data.testPositions().exclusionIds(helper.testPosIds()),
                       helper.exclusions(), gmx::EmptyArrayRef(), gmx::EmptyArrayRef());
}

} // namespace
