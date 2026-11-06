/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Tests for H5MD particle block class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_particleblock.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_timedatablock.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

// Helper variables
constexpr int64_t oneParticle    = 1;
constexpr int64_t sevenParticles = 7;

using H5mdParticleBlockTest = H5mdTestBase;

TEST_F(H5mdParticleBlockTest, DefaultBlocksAreEmpty)
{
    H5mdParticleBlock block(sevenParticles);
    EXPECT_FALSE(block.position().has_value());
    EXPECT_FALSE(block.velocity().has_value());
    EXPECT_FALSE(block.force().has_value());
    EXPECT_FALSE(block.box().has_value());
}

TEST_F(H5mdParticleBlockTest, NumParticlesWorks)
{
    H5mdParticleBlock block(sevenParticles);
    ASSERT_EQ(block.numParticles(), sevenParticles);
}

TEST_F(H5mdParticleBlockTest, SettersWork)
{
    const DataSetDims frameDims = { sevenParticles };
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock").withFrameDimension(frameDims).build();
    H5mdFrameDataSetBuilder<real>(fileid(), "matrixDataSetBlock").withFrameDimension({ DIM, DIM }).build();

    {
        SCOPED_TRACE("Set position data block only");
        H5mdParticleBlock block(sevenParticles);
        block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));

        ASSERT_TRUE(block.position().has_value());
        EXPECT_EQ(block.position()->frameDims(), frameDims);

        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.force().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set velocity data block only");
        H5mdParticleBlock block(sevenParticles);
        block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));

        ASSERT_TRUE(block.velocity().has_value());
        EXPECT_EQ(block.velocity()->frameDims(), frameDims);

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.force().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set force data block only");
        H5mdParticleBlock block(sevenParticles);
        block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));

        ASSERT_TRUE(block.force().has_value());
        EXPECT_EQ(block.force()->frameDims(), frameDims);

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set simulation box data block only");
        H5mdParticleBlock block(sevenParticles);
        block.addBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock"));

        ASSERT_TRUE(block.box().has_value());
        EXPECT_EQ(block.box()->frameDims(), (DataSetDims{ DIM, DIM }));

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.force().has_value());
    }
}

TEST_F(H5mdParticleBlockTest, HasBlockMethodsWork)
{
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock").withFrameDimension({ oneParticle }).build();
    H5mdFrameDataSetBuilder<real>(fileid(), "matrixDataSetBlock").withFrameDimension({ DIM, DIM }).build();

    {
        SCOPED_TRACE("Empty block");
        H5mdParticleBlock block(oneParticle);
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with positions");
        H5mdParticleBlock block(oneParticle);
        block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
        EXPECT_TRUE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with velocities");
        H5mdParticleBlock block(oneParticle);
        block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
        EXPECT_FALSE(block.hasPosition());
        EXPECT_TRUE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with forces");
        H5mdParticleBlock block(oneParticle);
        block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_TRUE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with boxes");
        H5mdParticleBlock block(oneParticle);
        block.addBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock"));
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_TRUE(block.hasBox());
    }
}

TEST_F(H5mdParticleBlockTest, AddBoxThrowsForBadDims)
{
    H5mdParticleBlock block(sevenParticles);

    ASSERT_NO_THROW_GMX(block.addBox(
            H5mdFrameDataSetBuilder<real>(fileid(), "3x3").withFrameDimension({ DIM, DIM }).build()))
            << "Matrices are { 3, 3 }, so this should work and all other cases should throw";

    EXPECT_THROW_GMX(block.addBox(H5mdFrameDataSetBuilder<real>(fileid(), "2x3")
                                          .withFrameDimension({ DIM - 1, DIM })
                                          .build()),
                     H5mdError);
    EXPECT_THROW_GMX(block.addBox(H5mdFrameDataSetBuilder<real>(fileid(), "3x2")
                                          .withFrameDimension({ DIM, DIM - 1 })
                                          .build()),
                     H5mdError);
    EXPECT_THROW_GMX(block.addBox(H5mdFrameDataSetBuilder<real>(fileid(), "4x3")
                                          .withFrameDimension({ DIM + 1, DIM })
                                          .build()),
                     H5mdError);
    EXPECT_THROW_GMX(block.addBox(H5mdFrameDataSetBuilder<real>(fileid(), "3x4")
                                          .withFrameDimension({ DIM, DIM + 1 })
                                          .build()),
                     H5mdError);
    EXPECT_THROW_GMX(
            block.addBox(
                    H5mdFrameDataSetBuilder<real>(fileid(), "3x1").withFrameDimension({ DIM, 1 }).build()),
            H5mdError);
    EXPECT_THROW_GMX(
            block.addBox(
                    H5mdFrameDataSetBuilder<real>(fileid(), "1x3").withFrameDimension({ 1, DIM }).build()),
            H5mdError);
}

TEST_F(H5mdParticleBlockTest, AddBlockThrowsForNon1DFrameDims)
{
    {
        SCOPED_TRACE("1-dimensional frames must work (the value is number of particles)");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "1d").withFrameDimension({ oneParticle }).build();

        H5mdParticleBlock block(oneParticle);
        EXPECT_NO_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
        EXPECT_NO_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
        EXPECT_NO_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
    }
    {
        SCOPED_TRACE("Zero-dimensional (empty) frames must throw");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "Empty").build();

        H5mdParticleBlock block(oneParticle);
        EXPECT_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), H5mdError);
        EXPECT_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), H5mdError);
        EXPECT_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), H5mdError);
    }
    {
        SCOPED_TRACE("Two-dimensional frames must throw");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "2d")
                .withFrameDimension({ oneParticle, oneParticle })
                .build();

        H5mdParticleBlock block(oneParticle);
        EXPECT_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "2d")), H5mdError);
        EXPECT_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "2d")), H5mdError);
        EXPECT_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "2d")), H5mdError);
    }
}

TEST_F(H5mdParticleBlockTest, AddBlockThrowsForFrameDimsNotMatchingNumParticles)
{
    H5mdParticleBlock block(sevenParticles);

    H5mdTimeDataBlockBuilder<RVec>(fileid(), "goodNumParticles")
            .withFrameDimension({ sevenParticles })
            .build();
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "lowerNumParticles")
            .withFrameDimension({ sevenParticles - 1 })
            .build();
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "higherNumParticles")
            .withFrameDimension({ sevenParticles + 1 })
            .build();

    EXPECT_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "lowerNumParticles")), H5mdError);
    EXPECT_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "higherNumParticles")), H5mdError);
    EXPECT_NO_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "goodNumParticles")))
            << "Sanity check failed: must be able to add block with correct number of particles";

    EXPECT_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "lowerNumParticles")), H5mdError);
    EXPECT_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "higherNumParticles")), H5mdError);
    EXPECT_NO_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "goodNumParticles")))
            << "Sanity check failed: must be able to add block with correct number of particles";

    EXPECT_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "lowerNumParticles")), H5mdError);
    EXPECT_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "higherNumParticles")), H5mdError);
    EXPECT_NO_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "goodNumParticles")))
            << "Sanity check failed: must be able to add block with correct number of particles";
}

TEST_F(H5mdParticleBlockTest, AddBlockMethodsThrowIfCalledMoreThanOnce)
{
    H5mdParticleBlock block(sevenParticles);

    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock")
            .withFrameDimension({ sevenParticles })
            .build();
    H5mdFrameDataSetBuilder<real>(fileid(), "matrixDataSetBlock").withFrameDimension({ DIM, DIM }).build();

    block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    EXPECT_THROW_GMX(block.addPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock")), H5mdError);

    block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    EXPECT_THROW_GMX(block.addVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock")), H5mdError);

    block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    EXPECT_THROW_GMX(block.addForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock")), H5mdError);

    block.addBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock"));
    EXPECT_THROW_GMX(block.addBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock")), H5mdError);
}

} // namespace
} // namespace test
} // namespace gmx
