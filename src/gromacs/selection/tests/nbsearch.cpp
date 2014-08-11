/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include <gtest/gtest.h>

#include <cmath>

#include <limits>
#include <set>
#include <vector>

#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"

#include "gromacs/selection/nbsearch.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * NeighborhoodSearchTestData
 */

class NeighborhoodSearchTestData
{
    public:
        struct TestPosition
        {
            explicit TestPosition(const rvec x)
                : refMinDist(0.0), refNearestPoint(-1)
            {
                copy_rvec(x, this->x);
            }

            rvec                x;
            real                refMinDist;
            int                 refNearestPoint;
            std::set<int>       refPairs;
        };
        typedef std::vector<TestPosition> TestPositionList;

        NeighborhoodSearchTestData(int seed, real cutoff);
        ~NeighborhoodSearchTestData();

        gmx::AnalysisNeighborhoodPositions refPositions() const
        {
            return gmx::AnalysisNeighborhoodPositions(refPos_, refPosCount_);
        }
        gmx::AnalysisNeighborhoodPositions testPositions() const
        {
            if (testPos_ == NULL)
            {
                snew(testPos_, testPositions_.size());
                for (size_t i = 0; i < testPositions_.size(); ++i)
                {
                    copy_rvec(testPositions_[i].x, testPos_[i]);
                }
            }
            return gmx::AnalysisNeighborhoodPositions(testPos_,
                                                      testPositions_.size());
        }
        gmx::AnalysisNeighborhoodPositions testPosition(int index) const
        {
            return testPositions().selectSingleFromArray(index);
        }

        void addTestPosition(const rvec x)
        {
            GMX_RELEASE_ASSERT(testPos_ == NULL,
                               "Cannot add positions after testPositions() call");
            testPositions_.push_back(TestPosition(x));
        }
        void generateRandomPosition(rvec x);
        void generateRandomRefPositions(int count);
        void generateRandomTestPositions(int count);
        void computeReferences(t_pbc *pbc);

        gmx_rng_t                        rng_;
        real                             cutoff_;
        matrix                           box_;
        t_pbc                            pbc_;
        int                              refPosCount_;
        rvec                            *refPos_;
        TestPositionList                 testPositions_;

    private:
        mutable rvec                    *testPos_;
};

NeighborhoodSearchTestData::NeighborhoodSearchTestData(int seed, real cutoff)
    : rng_(NULL), cutoff_(cutoff), refPosCount_(0), refPos_(NULL), testPos_(NULL)
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
    sfree(refPos_);
    sfree(testPos_);
}

void NeighborhoodSearchTestData::generateRandomPosition(rvec x)
{
    rvec fx;
    fx[XX] = gmx_rng_uniform_real(rng_);
    fx[YY] = gmx_rng_uniform_real(rng_);
    fx[ZZ] = gmx_rng_uniform_real(rng_);
    mvmul(box_, fx, x);
    // Add a small displacement to allow positions outside the box
    x[XX] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
    x[YY] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
    x[ZZ] += 0.2 * gmx_rng_uniform_real(rng_) - 0.1;
}

void NeighborhoodSearchTestData::generateRandomRefPositions(int count)
{
    refPosCount_ = count;
    snew(refPos_, refPosCount_);
    for (int i = 0; i < refPosCount_; ++i)
    {
        generateRandomPosition(refPos_[i]);
    }
}

void NeighborhoodSearchTestData::generateRandomTestPositions(int count)
{
    testPositions_.reserve(count);
    for (int i = 0; i < count; ++i)
    {
        rvec x;
        generateRandomPosition(x);
        addTestPosition(x);
    }
}

void NeighborhoodSearchTestData::computeReferences(t_pbc *pbc)
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
            const real dist = norm(dx);
            if (dist < i->refMinDist)
            {
                i->refMinDist      = dist;
                i->refNearestPoint = j;
            }
            if (dist <= cutoff)
            {
                i->refPairs.insert(j);
            }
        }
    }
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
        }
        else
        {
            EXPECT_EQ(i->refNearestPoint, -1);
        }
    }
}

void NeighborhoodSearchTest::testPairSearch(
        gmx::AnalysisNeighborhoodSearch  *search,
        const NeighborhoodSearchTestData &data)
{
    NeighborhoodSearchTestData::TestPositionList::const_iterator i;
    for (i = data.testPositions_.begin(); i != data.testPositions_.end(); ++i)
    {
        std::set<int> checkSet                         = i->refPairs;
        gmx::AnalysisNeighborhoodPairSearch pairSearch =
            search->startPairSearch(i->x);
        gmx::AnalysisNeighborhoodPair       pair;
        while (pairSearch.findNextPair(&pair))
        {
            EXPECT_EQ(0, pair.testIndex());
            if (checkSet.erase(pair.refIndex()) == 0)
            {
                // TODO: Check whether the same pair was returned more than
                // once and give a better error message if so.
                ADD_FAILURE()
                << "Expected: Position " << pair.refIndex()
                << " is within cutoff.\n"
                << "  Actual: It is not.";
            }
        }
        EXPECT_TRUE(checkSet.empty()) << "Some positions were not returned by the pair search.";
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
    // Currently, grid searching not supported with 2D PBC.
    //ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

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
    pairSearch1.findNextPair(&pair);
    EXPECT_EQ(0, pair.testIndex());
    EXPECT_TRUE(data.testPositions_[0].refPairs.count(pair.refIndex()) == 1);

    pairSearch2.findNextPair(&pair);
    EXPECT_EQ(1, pair.testIndex());
    EXPECT_TRUE(data.testPositions_[1].refPairs.count(pair.refIndex()) == 1);
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
        EXPECT_TRUE(data.testPositions_[currentIndex].refPairs.count(pair.refIndex()) == 1);
        pairSearch.skipRemainingPairsForTestPosition();
        ++currentIndex;
    }
}

} // namespace
