/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/gmx_random.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/vec.h"

#include "gromacs/selection/nbsearch.h"

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
            TestPosition() : refMinDist(0.0), refNearestPoint(-1)
            {
                clear_rvec(x);
            }
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

        void addTestPosition(const rvec x)
        {
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
};

NeighborhoodSearchTestData::NeighborhoodSearchTestData(int seed, real cutoff)
    : rng_(NULL), cutoff_(cutoff), refPosCount_(0), refPos_(NULL)
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
        EXPECT_NEAR_REL(refDist, search->minimumDistance(i->x), 20*GMX_REAL_EPS);
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
            EXPECT_EQ(i->refNearestPoint, pair.firstIndex());
            EXPECT_EQ(0, pair.secondIndex());
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
            EXPECT_EQ(0, pair.secondIndex());
            if (checkSet.erase(pair.firstIndex()) == 0)
            {
                // TODO: Check whether the same pair was returned more than
                // once and give a better error message if so.
                ADD_FAILURE()
                << "Expected: Position " << pair.firstIndex()
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

class RandomBoxFullPBCData : public NeighborhoodSearchTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBoxFullPBCData singleton;
            return singleton;
        }

        RandomBoxFullPBCData()
            : NeighborhoodSearchTestData(12345, 1.0)
        {
            box_[XX][XX] = 10.0;
            box_[YY][YY] = 5.0;
            box_[ZZ][ZZ] = 7.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            generateRandomRefPositions(1000);
            generateRandomTestPositions(100);
            set_pbc(&pbc_, epbcXYZ, box_);
            computeReferences(&pbc_);
        }
};

class RandomTriclinicFullPBCData : public NeighborhoodSearchTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomTriclinicFullPBCData singleton;
            return singleton;
        }

        RandomTriclinicFullPBCData()
            : NeighborhoodSearchTestData(12345, 1.0)
        {
            box_[XX][XX] = 5.0;
            box_[YY][XX] = 2.5;
            box_[YY][YY] = 2.5*sqrt(3.0);
            box_[ZZ][XX] = 2.5;
            box_[ZZ][YY] = 2.5*sqrt(1.0/3.0);
            box_[ZZ][ZZ] = 5.0*sqrt(2.0/3.0);
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            generateRandomRefPositions(1000);
            generateRandomTestPositions(100);
            set_pbc(&pbc_, epbcXYZ, box_);
            computeReferences(&pbc_);
        }
};

class RandomBox2DPBCData : public NeighborhoodSearchTestData
{
    public:
        static const NeighborhoodSearchTestData &get()
        {
            static RandomBox2DPBCData singleton;
            return singleton;
        }

        RandomBox2DPBCData()
            : NeighborhoodSearchTestData(12345, 1.0)
        {
            box_[XX][XX] = 10.0;
            box_[YY][YY] = 7.0;
            box_[ZZ][ZZ] = 5.0;
            // TODO: Consider whether manually picking some positions would give better
            // test coverage.
            generateRandomRefPositions(1000);
            generateRandomTestPositions(100);
            set_pbc(&pbc_, epbcXY, box_);
            computeReferences(&pbc_);
        }
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
        nb_.initSearch(&data.pbc_, data.refPosCount_, data.refPos_);
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
        nb_.initSearch(&data.pbc_, data.refPosCount_, data.refPos_);
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
        nb_.initSearch(&data.pbc_, data.refPosCount_, data.refPos_);
    ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testPairSearch(&search, data);
}

TEST_F(NeighborhoodSearchTest, GridSearch2DPBC)
{
    const NeighborhoodSearchTestData &data = RandomBox2DPBCData::get();

    nb_.setCutoff(data.cutoff_);
    nb_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
    gmx::AnalysisNeighborhoodSearch search =
        nb_.initSearch(&data.pbc_, data.refPosCount_, data.refPos_);
    // Currently, grid searching not supported with 2D PBC.
    //ASSERT_EQ(gmx::AnalysisNeighborhood::eSearchMode_Grid, search.mode());

    testIsWithin(&search, data);
    testMinimumDistance(&search, data);
    testNearestPoint(&search, data);
    testPairSearch(&search, data);
}

} // namespace
