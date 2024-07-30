/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Tests functionality related to repeating sections.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/repeatingsection.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::RepeatingOptionSection;

struct SectionData
{
    int value;
};

TEST(RepeatingOptionSectionTest, HandlesNoInstance)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ(0U, values.size());
}

TEST(RepeatingOptionSectionTest, HandlesNoInstanceWithRequiredOption)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value).required()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ(0U, values.size());
}

TEST(RepeatingOptionSectionTest, HandlesSingleInstance)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("4"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(1U, values.size());
    EXPECT_EQ(4, values[0].value);
}

TEST(RepeatingOptionSectionTest, HandlesDefaultValue)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value).defaultValue(3)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(1U, values.size());
    EXPECT_EQ(3, values[0].value);
}

TEST(RepeatingOptionSectionTest, HandlesTwoInstances)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("4"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(4, values[0].value);
    EXPECT_EQ(5, values[1].value);
}

TEST(RepeatingOptionSectionTest, HandlesUnsetOptionWithImplicitDefault)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("4"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(4, values[0].value);
    EXPECT_EQ(0, values[1].value);
}

TEST(RepeatingOptionSectionTest, HandlesUnsetOptionWithExplicitDefault)
{
    std::vector<SectionData> values;
    gmx::Options             options;
    auto sec = options.addSection(RepeatingOptionSection<SectionData>("section").storeVector(&values));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value).defaultValue(1)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("4"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(4, values[0].value);
    EXPECT_EQ(1, values[1].value);
}

struct NestedSectionData
{
    int                      value{};
    std::vector<SectionData> subsec;
};

TEST(RepeatingOptionSectionTest, HandlesNestedSections)
{
    std::vector<NestedSectionData> values;
    gmx::Options                   options;
    auto sec = options.addSection(RepeatingOptionSection<NestedSectionData>("section").storeVector(&values));
    auto subsec =
            sec.addSection(RepeatingOptionSection<SectionData>("subsec").storeVector(&sec.bind().subsec));
    using gmx::IntegerOption;
    ASSERT_NO_THROW_GMX(sec.addOption(IntegerOption("p").store(&sec.bind().value)));
    ASSERT_NO_THROW_GMX(subsec.addOption(IntegerOption("p").store(&subsec.bind().value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startSection("section"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("4"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startSection("subsec"));
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finishSection());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(1U, values.size());
    EXPECT_EQ(4, values[0].value);
    ASSERT_EQ(1U, values[0].subsec.size());
    EXPECT_EQ(5, values[0].subsec[0].value);
}

} // namespace
} // namespace test
} // namespace gmx
