/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
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
#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/options/timeunitmanager.h"

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
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value).timeValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("1.5"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());

    EXPECT_DOUBLE_EQ(1.5, value);
    manager.setTimeUnit(gmx::eTimeUnit_ns);
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500, value);

    EXPECT_NO_THROW(options.finish());

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
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value).timeValue()));
    ASSERT_NO_THROW(options.addOption(DoubleOption("q").store(&value2).timeValue()
                                          .defaultValueIfSet(2.5)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("q"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

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
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value).timeValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("1.5"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    gmx::OptionsAssigner assigner2(&options);
    EXPECT_NO_THROW(assigner2.start());
    EXPECT_NO_THROW(assigner2.finish());
    EXPECT_NO_THROW(options.finish());

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
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value).timeValue()));
    ASSERT_NO_THROW(manager.addTimeUnitOption(&options, "tu"));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("1.5"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("tu"));
    ASSERT_NO_THROW(assigner.appendValue("ns"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());

    EXPECT_DOUBLE_EQ(1.5, value);
    EXPECT_EQ(gmx::eTimeUnit_ns, manager.timeUnit());
    manager.scaleTimeOptions(&options);
    EXPECT_DOUBLE_EQ(1500, value);

    EXPECT_NO_THROW(options.finish());
}

} // namespace
