/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include <map>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
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

        NeighborhoodSearchTestData(gmx_uint64_t seed, real cutoff);

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
                    testPos_.emplace_back(testPositions_[i].x);
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
            testPositions_.emplace_back(x);
        }
        gmx::RVec generateRandomPosition();
        std::vector<int> generateIndex(int count, gmx_uint64_t seed) const;
        void generateRandomRefPositions(int count);
        void generateRandomTestPositions(int count);
        void useRefPositionsAsTestPositions();
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

        // Return a tolerance that accounts for the magnitudes of the coordinates
        // when doing subtractions, so that we set the ULP tolerance relative to the
        // coordinate values rather than their difference.
        // i.e., 10.0-9.9999999 will achieve a few ULP accuracy relative
        // to 10.0, but a much larger error relative to the difference.
        gmx::test::FloatingPointTolerance relativeTolerance() const
        {
            real magnitude = std::max(box_[XX][XX], std::max(box_[YY][YY], box_[ZZ][ZZ]));
            return gmx::test::relativeToleranceAsUlp(magnitude, 4);
        }

        gmx::DefaultRandomEngine         rng_;
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

NeighborhoodSearchTestData::NeighborhoodSearchTestData(gmx_uint64_t seed, real cutoff)
    : rng_(seed), cutoff_(cutoff), refPosCount_(0)
{
    clear_mat(box_);
    set_pbc(&pbc_, epbcNONE, box_);
}

gmx::RVec NeighborhoodSearchTestData::generateRandomPosition()
{
    gmx::UniformRealDistribution<real>  dist;
    rvec fx, x;
    fx[XX] = dist(rng_);
    fx[YY] = dist(rng_);
    fx[ZZ] = dist(rng_);
    mvmul(box_, fx, x);
    // Add a small displacement to allow positions outside the box
    x[XX] += 0.2 * dist(rng_) - 0.1;
    x[YY] += 0.2 * dist(rng_) - 0.1;
    x[ZZ] += 0.2 * dist(rng_) - 0.1;
    return x;
}

std::vector<int> NeighborhoodSearchTestData::generateIndex(int count, gmx_uint64_t seed) const
{
    gmx::DefaultRandomEngine             rngIndex(seed);
    gmx::UniformRealDistribution<real>   dist;
    std::vector<int>                     result;

    for (int i = 0; i < count; ++i)
    {
        if (dist(rngIndex) > 0.5)
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

void NeighborhoodSearchTestData::useRefPositionsAsTestPositions()
{
    testPositions_.reserve(refPosCount_);
    for (const auto &refPos : refPos_)
    {
        addTestPosition(refPos);
    }
}

void NeighborhoodSearchTestData::computeReferencesInternal(t_pbc *pbc, bool bXY)
{
    real cutoff = cutoff_;
    if (cutoff <= 0)
    {
        cutoff = std::numeric_limits<real>::max();
    }
    for (TestPosition &testPos : testPositions_)
    {
        testPos.refMinDist      = cutoff;
        testPos.refNearestPoint = -1;
        testPos.refPairs.clear();
        for (int j = 0; j < refPosCount_; ++j)
        {
            rvec dx;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, testPos.x, refPos_[j], dx);
            }
            else
            {
                rvec_sub(testPos.x, refPos_[j], dx);
            }
            // TODO: This may not work intuitively for 2D with the third box
            // vector not parallel to the Z axis, but neither does the actual
            // neighborhood search.
            const real dist =
                !bXY ? norm(dx) : std::hypot(dx[XX], dx[YY]);
            if (dist < testPos.refMinDist)
            {
                testPos.refMinDist      = dist;
                testPos.refNearestPoint = j;
            }
            if (dist > 0 && dist <= cutoff)
            {
                RefPair pair(j, dist);
                GMX_RELEASE_ASSERT(testPos.refPairs.empty() || testPos.refPairs.back() < pair,
                                   "Reference pairs should be generated in sorted order");
                testPos.refPairs.push_back(pair);
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

        gmx::ArrayRef<const int> refPosIds() const
        {
            return gmx::makeArrayRef(exclusionIds_).subArray(0, refPosCount_);
        }
        gmx::ArrayRef<const int> testPosIds() const
        {
            return gmx::makeArrayRef(exclusionIds_).subArray(0, testPosCount_);
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
    excls_.index        = nullptr;
    excls_.nra          = 0;
    excls_.a            = nullptr;
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
    excls_.index = exclsIndex_.data();
    excls_.nra   = exclsAtoms_.size();
    excls_.a     = exclsAtoms_.data();
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
                                   const NeighborhoodSearchTestData &data,
                                   gmx_uint64_t                      seed);
        void testPairSearchFull(gmx::AnalysisNeighborhoodSearch          *search,
                                const NeighborhoodSearchTestData         &data,
                                const gmx::AnalysisNeighborhoodPositions &pos,
                                const t_blocka                           *excls,
                                const gmx::ArrayRef<const int>           &refIndices,
                                const gmx::ArrayRef<const int>           &testIndices,
                                bool                                      selfPairs);

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
        EXPECT_REAL_EQ_TOL(refDist, search->minimumDistance(i->x), data.relativeTolerance());
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
            EXPECT_REAL_EQ_TOL(i->refMinDist, std::sqrt(pair.distance2()), data.relativeTolerance());
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
    testPairSearchFull(search, data, data.testPositions(), nullptr,
                       gmx::EmptyArrayRef(), gmx::EmptyArrayRef(), false);
}

void NeighborhoodSearchTest::testPairSearchIndexed(
        gmx::AnalysisNeighborhood        *nb,
        const NeighborhoodSearchTestData &data,
        gmx_uint64_t                      seed)
{
    std::vector<int>                refIndices(data.generateIndex(data.refPos_.size(), seed++));
    std::vector<int>                testIndices(data.generateIndex(data.testPositions_.size(), seed++));
    gmx::AnalysisNeighborhoodSearch search =
        nb->initSearch(&data.pbc_,
                       data.refPositions().indexed(refIndices));
    testPairSearchFull(&search, data, data.testPositions(), nullptr,
                       refIndices, testIndices, false);
}

void NeighborhoodSearchTest::testPairSearchFull(
        gmx::AnalysisNeighborhoodSearch          *search,
        const NeighborhoodSearchTestData         &data,
        const gmx::AnalysisNeighborhoodPositions &pos,
        const t_blocka                           *excls,
        const gmx::ArrayRef<const int>           &refIndices,
        const gmx::ArrayRef<const int>           &testIndices,
        bool                                      selfPairs)
{
    std::map<int, RefPairList> refPairs;
    // TODO: Some parts of this code do not work properly if pos does not
    // initially contain all the test positions.
    if (testIndices.empty())
    {
        for (size_t i = 0; i < data.testPositions_.size(); ++i)
        {
            refPairs[i] = data.testPositions_[i].refPairs;
        }
    }
    else
    {
        for (int index : testIndices)
        {
            refPairs[index] = data.testPositions_[index].refPairs;
        }
    }
    if (excls != nullptr)
    {
        GMX_RELEASE_ASSERT(!selfPairs, "Self-pairs testing not implemented with exclusions");
        for (auto &entry : refPairs)
        {
            const int testIndex = entry.first;
            ExclusionsHelper::markExcludedPairs(&entry.second, testIndex, excls);
        }
    }
    if (!refIndices.empty())
    {
        GMX_RELEASE_ASSERT(!selfPairs, "Self-pairs testing not implemented with indexing");
        for (auto &entry : refPairs)
        {
            for (auto &refPair : entry.second)
            {
                refPair.bIndexed = false;
            }
            for (int index : refIndices)
            {
                NeighborhoodSearchTestData::RefPair searchPair(index, 0.0);
                auto refPair = std::lower_bound(entry.second.begin(), entry.second.end(), searchPair);
                if (refPair != entry.second.end() && refPair->refIndex == index)
                {
                    refPair->bIndexed = true;
                }
            }
            for (auto &refPair : entry.second)
            {
                if (!refPair.bIndexed)
                {
                    refPair.bFound = true;
                }
            }
        }
    }

    gmx::AnalysisNeighborhoodPositions  posCopy(pos);
    if (!testIndices.empty())
    {
        posCopy.indexed(testIndices);
    }
    gmx::AnalysisNeighborhoodPairSearch pairSearch
        = selfPairs
            ? search->startSelfPairSearch()
            : search->startPairSearch(posCopy);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const int testIndex =
            (testIndices.empty() ? pair.testIndex() : testIndices[pair.testIndex()]);
        const int refIndex =
            (refIndices.empty() ? pair.refIndex() : refIndices[pair.refIndex()]);

        if (refPairs.count(testIndex) == 0)
        {
            ADD_FAILURE()
            << "Expected: No pairs are returned for test position " << testIndex << ".\n"
            << "  Actual: Pair with ref " << refIndex << " is returned.";
            continue;
        }
        NeighborhoodSearchTestData::RefPair searchPair(refIndex,
                                                       std::sqrt(pair.distance2()));
        const auto foundRefPair
            = std::lower_bound(refPairs[testIndex].begin(), refPairs[testIndex].end(),
                               searchPair);
        if (foundRefPair == refPairs[testIndex].end() || foundRefPair->refIndex != refIndex)
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

            EXPECT_REAL_EQ_TOL(foundRefPair->distance, searchPair.distance, data.relativeTolerance())
            << "Distance computed by the neighborhood search does not match.";
            if (selfPairs)
            {
                searchPair = NeighborhoodSearchTestData::RefPair(testIndex, 0.0);
                const auto otherRefPair
                    = std::lower_bound(refPairs[refIndex].begin(), refPairs[refIndex].end(),
                                       searchPair);
                GMX_RELEASE_ASSERT(otherRefPair != refPairs[refIndex].end(),
                                   "Precomputed reference data is not symmetric");
                otherRefPair->bFound = true;
            }
        }
    }

    for (auto &entry : refPairs)
    {
        const int testIndex = entry.first;
        checkAllPairsFound(entry.second, data.refPos_, testIndex,
                           data.testPositions_[testIndex].x);
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
            // Make the box so small we are virtually guaranteed to have
            // several neighbors for the five test positions
            data_.box_[XX][XX] = 3.0;
            data_.box_[YY][YY] = 3.0;
            data_.box_[ZZ][ZZ] = 3.0;
            data_.generateRandomRefPositions(10);
            data_.generateRandomTestPositions(5);
            set_pbc(&data_.pbc_, epbcXYZ, data_.box_);
            data_.computeReferences(&data_.pbc_);
        }

    private:
        NeighborhoodSearchTestData data_;
};

class TrivialSelfPairsTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static TrivialSelfPairsTestData singleton;
            return singleton.data_;
        }

        TrivialSelfPairsTestData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 3.0;
            data_.box_[YY][YY] = 3.0;
            data_.box_[ZZ][ZZ] = 3.0;
            data_.generateRandomRefPositions(20);
            data_.useRefPositionsAsTestPositions();
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
            data_.computeReferences(nullptr);
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

class RandomBoxSelfPairsData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBoxSelfPairsData singleton;
            return singleton.data_;
        }

        RandomBoxSelfPairsData() : data_(12345, 1.0)
        {
            data_.box_[XX][XX] = 10.0;
            data_.box_[YY][YY] = 5.0;
            data_.box_[ZZ][ZZ] = 7.0;
            data_.generateRandomRefPositions(1000);
            data_.useRefPositionsAsTestPositions();
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
            data_.box_[YY][YY] = 2.5*std::sqrt(3.0);
            data_.box_[ZZ][XX] = 2.5;
            data_.box_[ZZ][YY] = 2.5*std::sqrt(1.0/3.0);
            data_.box_[ZZ][ZZ] = 5.0*std::sqrt(2.0/3.0);
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
            data_.computeReferences(nullptr);
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
    testPairSearchIndexed(&nb_, data, 123);
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
    testPairSearchIndexed(&nb_, data, 456);
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

TEST_F(NeighborhoodSearchTest, SimpleSelfPairsSearch)
{
    const NeighborhoodSearchTestData &data = TrivialSelfPairsTestData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Simple);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Simple, search.mode());

    testPairSearchFull(&search, data, data.testPositions(), nullptr,
                       gmx::EmptyArrayRef(), gmx::EmptyArrayRef(), true);
}

TEST_F(NeighborhoodSearchTest, GridSelfPairsSearch)
{
    const NeighborhoodSearchTestData &data = RandomBoxSelfPairsData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPositions());
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testPairSearchFull(&search, data, data.testPositions(), nullptr,
                       gmx::EmptyArrayRef(), gmx::EmptyArrayRef(), true);
}

TEST_F(NeighborhoodSearchTest, HandlesConcurrentSearches)
{
    const NeighborhoodSearchTestData &data = TrivialTestData::get();

    nb_.setCutoff(data.cutoff_);
    gmx::AnalysisNeighborhoodSearch search1 =
        nb_.initSearch(&data.pbc_, data.refPositions());
    gmx::AnalysisNeighborhoodSearch search2 =
        nb_.initSearch(&data.pbc_, data.refPositions());

    // These checks are fragile, and unfortunately depend on the random
    // engine used to create the test positions. There is no particular reason
    // why exactly particles 0 & 2 should have neighbors, but in this case they do.
    gmx::AnalysisNeighborhoodPairSearch pairSearch1 =
        search1.startPairSearch(data.testPosition(0));
    gmx::AnalysisNeighborhoodPairSearch pairSearch2 =
        search1.startPairSearch(data.testPosition(2));

    testPairSearch(&search2, data);

    gmx::AnalysisNeighborhoodPair pair;
    ASSERT_TRUE(pairSearch1.findNextPair(&pair))
    << "Test data did not contain any pairs for position 0 (problem in the test).";
    EXPECT_EQ(0, pair.testIndex());
    {
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), std::sqrt(pair.distance2()));
        EXPECT_TRUE(data.containsPair(0, searchPair));
    }

    ASSERT_TRUE(pairSearch2.findNextPair(&pair))
    << "Test data did not contain any pairs for position 2 (problem in the test).";
    EXPECT_EQ(2, pair.testIndex());
    {
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), std::sqrt(pair.distance2()));
        EXPECT_TRUE(data.containsPair(2, searchPair));
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
        nb_.initSearch(nullptr, data.refPositions());
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
        NeighborhoodSearchTestData::RefPair searchPair(pair.refIndex(), std::sqrt(pair.distance2()));
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
                       helper.exclusions(), gmx::EmptyArrayRef(),
                       gmx::EmptyArrayRef(), false);
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
                       helper.exclusions(), gmx::EmptyArrayRef(),
                       gmx::EmptyArrayRef(), false);
}

} // namespace
