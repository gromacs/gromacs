/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Tests handling of time units with gmx::TimeUnitManager.
 *
 * Also related functionality in gmx::DoubleOptionStorage is tested.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/timeunitmanager.h"

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"

#include "testutils/testasserts.h"

namespace
{

TEST(TimeUnitManagerTest, BasicOperations)
{
    gmx::TimeUnitManager manager;
    EXPECT_EQ(gmx::eTimeUnit_ps, manager.timeUnit());
    EXPECT_DOUBLE_EQ(1.0, manager.timeScaleFactor());
    manager.setTimeUnit(gmx::eTimeUnit_ns);
    EXPECT_EQ(gmx::eTimeUnit_ns, manager.timeUnit());
    EXPECT_DOUBLE_EQ(1e3, manager.timeScaleFactor());
    EXPECT_DOUBLE_EQ(1e-3, manager.inverseTimeScaleFactor());
}

TEST(TimeUnitManagerTest, ScalesAssignedOptionValue)
{
    gmx::TimeUnitManager manager;

    gmx::Options         options(NULL, NULL);
    double               value = 0.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW_GMX(options.addOption(DoubleOption("p").store(&value).timeValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("1.5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());

    EXPECT_DOUBLE_EQ(1.5, value);
    manager.setTimeUnit(gmx::eTimeUnit_ns);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500, value);

    EXPECT_NO_THROW_GMX(options.finish());

    manager.setTimeUnit(gmx::eTimeUnit_us);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500000, value);

    manager.setTimeUnit(gmx::eTimeUnit_fs);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(0.0015, value);

    manager.setTimeUnit(gmx::eTimeUnit_ps);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1.5, value);
}

TEST(TimeUnitManagerTest, DoesNotScaleDefaultValues)
{
    gmx::TimeUnitManager manager;

    gmx::Options         options(NULL, NULL);
    double               value = 1.5, value2 = 0.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW_GMX(options.addOption(DoubleOption("p").store(&value).timeValue()));
    ASSERT_NO_THROW_GMX(options.addOption(DoubleOption("q").store(&value2).timeValue()
                                              .defaultValueIfSet(2.5)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("q"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_DOUBLE_EQ(2.5, value2);
    manager.setTimeUnit(gmx::eTimeUnit_ns);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1.5, value);
    EXPECT_DOUBLE_EQ(2.5, value2);
}

TEST(TimeUnitManagerTest, ScalesUserInputWithMultipleSources)
{
    gmx::TimeUnitManager manager;

    gmx::Options         options(NULL, NULL);
    double               value = 0.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW_GMX(options.addOption(DoubleOption("p").store(&value).timeValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("1.5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    gmx::OptionsAssigner assigner2(&options);
    EXPECT_NO_THROW_GMX(assigner2.start());
    EXPECT_NO_THROW_GMX(assigner2.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_DOUBLE_EQ(1.5, value);
    manager.setTimeUnit(gmx::eTimeUnit_ns);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500, value);
}

TEST(TimeUnitManagerTest, TimeUnitOptionWorks)
{
    gmx::TimeUnitManager manager;

    gmx::Options         options(NULL, NULL);
    double               value = 0.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW_GMX(options.addOption(DoubleOption("p").store(&value).timeValue()));
    ASSERT_NO_THROW_GMX(manager.addTimeUnitOption(&options, "tu"));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("1.5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("tu"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("ns"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());

    EXPECT_DOUBLE_EQ(1.5, value);
    EXPECT_EQ(gmx::eTimeUnit_ns, manager.timeUnit());
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500, value);

    EXPECT_NO_THROW_GMX(options.finish());
}

} // namespace
