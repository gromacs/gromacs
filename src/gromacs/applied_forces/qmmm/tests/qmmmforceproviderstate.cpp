/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for QMMMForceProviderState class used for checkpointing.
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/qmmm/qmmmforceproviderstate.h"

#include <gtest/gtest.h>

#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

namespace gmx
{

class QMMMForceProviderStateTest : public ::testing::Test
{
protected:
    RVec defaultQMTrans_ = { 1.0, -2.5, 3.25 };
};

TEST_F(QMMMForceProviderStateTest, WritesTranslationVectorToCheckpoint)
{
    QMMMForceProviderState state;
    state.setQMTrans(defaultQMTrans_);

    // Try writing to KVT
    KeyValueTreeBuilder builder;
    state.writeState(builder.rootObject(), "qmmm-cp2k");
    const KeyValueTreeObject kvt = builder.build();

    // Check that data is there
    const auto values = kvt["qmmm-cp2k-qmtrans"].asArray().values();
    ASSERT_EQ(DIM, values.size());
    for (int i = 0; i < DIM; i++)
    {
        EXPECT_DOUBLE_EQ(static_cast<double>(defaultQMTrans_[i]), values[i].cast<double>());
    }
}

TEST_F(QMMMForceProviderStateTest, ReadsTranslationVectorWhenPresent)
{
    // Prepare KVT with data
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder objectBuilder = builder.rootObject();
    auto arrayAdder = objectBuilder.addUniformArray<double>("qmmm-cp2k-qmtrans");
    for (auto i = 0; i < DIM; i++)
    {
        arrayAdder.addValue(static_cast<double>(defaultQMTrans_[i]));
    }

    // Try reading from KVT
    QMMMForceProviderState state;
    state.readState(builder.build(), "qmmm-cp2k");

    // Check that state was read correctly
    EXPECT_TRUE(state.isStateRead());
    const RVec& qmTrans = state.qmTrans();
    for (int i = 0; i < DIM; i++)
    {
        EXPECT_DOUBLE_EQ(static_cast<double>(defaultQMTrans_[i]), static_cast<double>(qmTrans[i]));
    }
}

TEST_F(QMMMForceProviderStateTest, MissingStateWhenNoData)
{
    // Prepare empty KVT
    KeyValueTreeBuilder    builder;
    QMMMForceProviderState state;

    // Try reading from KVT
    state.readState(builder.build(), "qmmm-cp2k");

    // Check that state was not read
    EXPECT_FALSE(state.isStateRead());
    const RVec& qmTrans = state.qmTrans();
    for (int i = 0; i < DIM; i++)
    {
        EXPECT_DOUBLE_EQ(0.0, static_cast<double>(qmTrans[i]));
    }
}

} // namespace gmx
