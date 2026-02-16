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

#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_timedatablock.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

using H5mdParticleBlockTest = H5mdTestBase;

TEST_F(H5mdParticleBlockTest, DefaultBlockGettersAreEmpty)
{
    H5mdParticleBlock block = H5mdParticleBlockBuilder().build();
    EXPECT_FALSE(block.position().has_value());
    EXPECT_FALSE(block.velocity().has_value());
    EXPECT_FALSE(block.force().has_value());
    EXPECT_FALSE(block.box().has_value());
}

TEST_F(H5mdParticleBlockTest, BuilderSettersWorkTogether)
{
    const DataSetDims frameDims = { 7 };
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock").withFrameDimension(frameDims).build();

    H5mdParticleBlockBuilder builder;
    builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"));
    H5mdParticleBlock block = builder.build();

    EXPECT_TRUE(block.position().has_value());
    EXPECT_EQ(block.position()->frameDims(), frameDims);
    EXPECT_TRUE(block.velocity().has_value());
    EXPECT_EQ(block.velocity()->frameDims(), frameDims);
    EXPECT_TRUE(block.force().has_value());
    EXPECT_EQ(block.force()->frameDims(), frameDims);
}

TEST_F(H5mdParticleBlockTest, BuilderSettersWorkIndividually)
{
    const DataSetDims frameDims = { 7 };
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock").withFrameDimension(frameDims).build();
    H5mdFrameDataSetBuilder<real>(fileid(), "matrixDataSetBlock").withFrameDimension({ DIM, DIM }).build();

    {
        SCOPED_TRACE("Set position data block only");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();

        EXPECT_TRUE(block.position().has_value());
        EXPECT_EQ(block.position()->frameDims(), frameDims);

        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.force().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set velocity data block only");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();

        EXPECT_TRUE(block.velocity().has_value());
        EXPECT_EQ(block.velocity()->frameDims(), frameDims);

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.force().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set force data block only");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();

        EXPECT_TRUE(block.force().has_value());
        EXPECT_EQ(block.force()->frameDims(), frameDims);

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.box().has_value());
    }
    {
        SCOPED_TRACE("Set simulation box data block only");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock"))
                        .build();
        EXPECT_TRUE(block.box().has_value());
        EXPECT_EQ(block.box()->frameDims(), (DataSetDims{ DIM, DIM }));

        EXPECT_FALSE(block.position().has_value());
        EXPECT_FALSE(block.velocity().has_value());
        EXPECT_FALSE(block.force().has_value());
    }
}

TEST_F(H5mdParticleBlockTest, HasBlockMethodsWork)
{
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "rvecDataSetBlock").withFrameDimension({ 1 }).build();
    H5mdFrameDataSetBuilder<real>(fileid(), "matrixDataSetBlock").withFrameDimension({ DIM, DIM }).build();

    {
        SCOPED_TRACE("Empty block");
        H5mdParticleBlock block = H5mdParticleBlockBuilder().build();
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with positions");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setPosition(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();
        EXPECT_TRUE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with velocities");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();
        EXPECT_FALSE(block.hasPosition());
        EXPECT_TRUE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with forces");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setForce(H5mdTimeDataBlock<RVec>(fileid(), "rvecDataSetBlock"))
                        .build();
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_TRUE(block.hasForce());
        EXPECT_FALSE(block.hasBox());
    }
    {
        SCOPED_TRACE("Block with boxes");
        H5mdParticleBlock block =
                H5mdParticleBlockBuilder()
                        .setBox(H5mdFrameDataSet<real>(fileid(), "matrixDataSetBlock"))
                        .build();
        EXPECT_FALSE(block.hasPosition());
        EXPECT_FALSE(block.hasVelocity());
        EXPECT_FALSE(block.hasForce());
        EXPECT_TRUE(block.hasBox());
    }
}

TEST_F(H5mdParticleBlockTest, SetBoxThrowsForBadDims)
{
    H5mdParticleBlockBuilder builder;

    ASSERT_NO_THROW(H5mdFrameDataSetBuilder<RVec>(fileid(), "3x3").withFrameDimension({ DIM, DIM }).build())
            << "Matrices are { 3, 3 }, so this should work and all other cases should throw";

    EXPECT_THROW(builder.setBox(H5mdFrameDataSetBuilder<real>(fileid(), "2x3")
                                        .withFrameDimension({ DIM - 1, DIM })
                                        .build()),
                 gmx::FileIOError);
    EXPECT_THROW(builder.setBox(H5mdFrameDataSetBuilder<real>(fileid(), "3x2")
                                        .withFrameDimension({ DIM, DIM - 1 })
                                        .build()),
                 gmx::FileIOError);
    EXPECT_THROW(builder.setBox(H5mdFrameDataSetBuilder<real>(fileid(), "4x3")
                                        .withFrameDimension({ DIM + 1, DIM })
                                        .build()),
                 gmx::FileIOError);
    EXPECT_THROW(builder.setBox(H5mdFrameDataSetBuilder<real>(fileid(), "3x4")
                                        .withFrameDimension({ DIM, DIM + 1 })
                                        .build()),
                 gmx::FileIOError);
    EXPECT_THROW(
            builder.setBox(
                    H5mdFrameDataSetBuilder<real>(fileid(), "3x1").withFrameDimension({ DIM, 1 }).build()),
            gmx::FileIOError);
    EXPECT_THROW(
            builder.setBox(
                    H5mdFrameDataSetBuilder<real>(fileid(), "1x3").withFrameDimension({ 1, DIM }).build()),
            gmx::FileIOError);
}

TEST_F(H5mdParticleBlockTest, SetBlockThrowsForBadDims)
{
    {
        SCOPED_TRACE("1-dimensional frames must work (the value is number of particles)");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "1d").withFrameDimension({ 1 }).build();

        H5mdParticleBlockBuilder builder;
        EXPECT_NO_THROW(builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
        EXPECT_NO_THROW(builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
        EXPECT_NO_THROW(builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "1d")));
    }
    {
        SCOPED_TRACE("Zero-dimensional (empty) frames must throw");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "Empty").build();

        H5mdParticleBlockBuilder builder;
        EXPECT_THROW(builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), gmx::FileIOError);
        EXPECT_THROW(builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), gmx::FileIOError);
        EXPECT_THROW(builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "Empty")), gmx::FileIOError);
    }
    {
        SCOPED_TRACE("Two-dimensional frames must throw");
        H5mdTimeDataBlockBuilder<RVec>(fileid(), "2d").withFrameDimension({ 1, 1 }).build();

        H5mdParticleBlockBuilder builder;
        EXPECT_THROW(builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "2d")), gmx::FileIOError);
        EXPECT_THROW(builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "2d")), gmx::FileIOError);
        EXPECT_THROW(builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "2d")), gmx::FileIOError);
    }
}

TEST_F(H5mdParticleBlockTest, SetBlockThrowsForInconsistentFrameDims)
{
    const DataSetDims frameDims3 = { 3 };
    const DataSetDims frameDims4 = { 4 };
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "frameDims3").withFrameDimension(frameDims3).build();
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "frameDims4").withFrameDimension(frameDims4).build();

    {
        SCOPED_TRACE("Position before velocity and force");
        H5mdParticleBlockBuilder builder;
        builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "frameDims3")); // write
        EXPECT_THROW(builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "frameDims4")), gmx::FileIOError);
        EXPECT_THROW(builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "frameDims4")), gmx::FileIOError);
        EXPECT_EQ(builder.build().numParticles(), 3)
                << "numParticles should be as for the position block";
    }
    {
        SCOPED_TRACE("Velocity before position and force");
        H5mdParticleBlockBuilder builder;
        builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "frameDims4"));
        EXPECT_THROW(builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "frameDims3")), gmx::FileIOError);
        EXPECT_THROW(builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "frameDims3")), gmx::FileIOError);
        EXPECT_EQ(builder.build().numParticles(), 4)
                << "numParticles should be as for the velocity block";
        ;
    }
    {
        SCOPED_TRACE("Force before position and velocity");
        H5mdParticleBlockBuilder builder;
        builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "frameDims3"));
        EXPECT_THROW(builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "frameDims4")), gmx::FileIOError);
        EXPECT_THROW(builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "frameDims4")), gmx::FileIOError);
        EXPECT_EQ(builder.build().numParticles(), 3)
                << "numParticles should be as for the force block";
        ;
    }
}

TEST_F(H5mdParticleBlockTest, NumParticles)
{
    const int64_t numParticles = 13;
    H5mdTimeDataBlockBuilder<RVec>(fileid(), "testBlock").withFrameDimension({ numParticles }).build();

    {
        SCOPED_TRACE("numParticles is 0 for empty block");
        H5mdParticleBlock block = H5mdParticleBlockBuilder().build();
        EXPECT_EQ(block.numParticles(), 0);
    }
    {
        SCOPED_TRACE("numParticles can be set by the position block");
        H5mdParticleBlockBuilder builder;
        builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        H5mdParticleBlock block = builder.build();
        EXPECT_EQ(block.numParticles(), numParticles);
    }
    {
        SCOPED_TRACE("numParticles can be set by the velocity block");
        H5mdParticleBlockBuilder builder;
        builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        H5mdParticleBlock block = builder.build();
        EXPECT_EQ(block.numParticles(), numParticles);
    }
    {
        SCOPED_TRACE("numParticles can be set by the force block");
        H5mdParticleBlockBuilder builder;
        builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        H5mdParticleBlock block = builder.build();
        EXPECT_EQ(block.numParticles(), numParticles);
    }
    {
        SCOPED_TRACE("Unaffected by setting blocks multiple times");
        H5mdParticleBlockBuilder builder;
        builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), "testBlock"));
        H5mdParticleBlock block = builder.build();
        EXPECT_EQ(block.numParticles(), numParticles);
    }
}

TEST_F(H5mdParticleBlockTest, HasTime)
{
    // Create the data structure for two H5mdTimeDataBlocks:
    // - One with a time data set
    // - One without a time data set
    // By pointing our H5mdTimeDataBlock builders to the corresponding
    // names when constructing our H5mdParticleBlocks we can select
    // whether they have a time data set or not in the tests below.

    constexpr char withTimeName[]        = "withTime";
    const auto [withTime, withTimeGuard] = makeH5mdGroupGuard(createGroup(fileid(), withTimeName));
    H5mdFrameDataSetBuilder<RVec>(withTime, "value").withFrameDimension({ 1 }).build();
    H5mdFrameDataSetBuilder<int64_t>(withTime, "step").build();
    H5mdFrameDataSetBuilder<double>(withTime, "time").build();

    constexpr char withoutTimeName[] = "withoutTime";
    const auto [withoutTime, withoutTimeGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), withoutTimeName));
    H5mdFrameDataSetBuilder<RVec>(withoutTime, "value").withFrameDimension({ 1 }).build();
    H5mdFrameDataSetBuilder<int64_t>(withoutTime, "step").build();

    {
        SCOPED_TRACE("Empty blocks do not have time");
        H5mdParticleBlock block = H5mdParticleBlockBuilder().build();
        EXPECT_FALSE(block.hasTime());
    }
    {
        SCOPED_TRACE("Adding blocks without time does not affect hasTime");
        H5mdParticleBlockBuilder builder;
        builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), withoutTimeName));
        builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), withoutTimeName));
        builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), withoutTimeName));
        EXPECT_FALSE(builder.build().hasTime());
    }
    {
        SCOPED_TRACE("Adding any block with time sets it to true");
        H5mdParticleBlockBuilder builder;
        builder.setPosition(H5mdTimeDataBlock<RVec>(fileid(), withoutTimeName));
        // Velocity data set has a time data set
        builder.setVelocity(H5mdTimeDataBlock<RVec>(fileid(), withTimeName));
        builder.setForce(H5mdTimeDataBlock<RVec>(fileid(), withoutTimeName));
        EXPECT_TRUE(builder.build().hasTime());
    }
}

} // namespace
} // namespace test
} // namespace gmx
