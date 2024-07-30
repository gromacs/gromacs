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
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreetransform.h"

#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

class TreeValueTransformTest : public ::testing::Test
{
public:
    void testTransform(const gmx::KeyValueTreeObject& input, const gmx::KeyValueTreeTransformer& transform)
    {
        gmx::KeyValueTreeTransformResult result = transform.transform(input, nullptr);
        gmx::KeyValueTreeObject          object = result.object();

        gmx::test::TestReferenceData    data;
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(input, "Input");
        auto mappedPaths = transform.mappedPaths();
        checker.checkSequence(
                mappedPaths.begin(), mappedPaths.end(), "MappedPaths", &TreeValueTransformTest::checkMappedPath);
        checker.checkKeyValueTreeObject(object, "Tree");
        checkBackMapping(&checker, object, result.backMapping());
    }

private:
    static void checkMappedPath(gmx::test::TestReferenceChecker* checker, const gmx::KeyValueTreePath& path)
    {
        checker->checkString(path.toString(), nullptr);
    }
    void checkBackMapping(gmx::test::TestReferenceChecker*     checker,
                          const gmx::KeyValueTreeObject&       object,
                          const gmx::IKeyValueTreeBackMapping& mapping)
    {
        auto compound(checker->checkCompound("BackMapping", "Mapping"));
        checkBackMappingImpl(&compound, object, mapping, gmx::KeyValueTreePath());
    }

    void checkBackMappingImpl(gmx::test::TestReferenceChecker*     checker,
                              const gmx::KeyValueTreeObject&       object,
                              const gmx::IKeyValueTreeBackMapping& mapping,
                              const gmx::KeyValueTreePath&         prefix)
    {
        for (const auto& prop : object.properties())
        {
            gmx::KeyValueTreePath path = prefix;
            path.append(prop.key());
            if (prop.value().isObject())
            {
                checkBackMappingImpl(checker, prop.value().asObject(), mapping, path);
            }
            else
            {
                gmx::KeyValueTreePath orgPath = mapping.originalPath(path);
                checker->checkString(orgPath.toString(), path.toString().c_str());
            }
        }
    }
};

TEST_F(TreeValueTransformTest, SimpleTransforms)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "1");
    builder.rootObject().addValue<std::string>("b", "2");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().from<std::string>("/a").to<int>("/i").transformWith(
            &gmx::fromStdString<int>);
    transform.rules()->addRule().from<std::string>("/b").to<int>("/j").transformWith(
            &gmx::fromStdString<int>);

    testTransform(input, transform);
}

TEST_F(TreeValueTransformTest, SimpleTransformsCaseAndDashInsensitive)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a-x", "1");
    builder.rootObject().addValue<std::string>("by", "2");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
    transform.rules()->addRule().from<std::string>("/Ax").to<int>("/i").transformWith(
            &gmx::fromStdString<int>);
    transform.rules()->addRule().from<std::string>("/B-Y").to<int>("/j").transformWith(
            &gmx::fromStdString<int>);

    testTransform(input, transform);
}

TEST_F(TreeValueTransformTest, SimpleTransformsToObject)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "1");
    builder.rootObject().addValue<std::string>("b", "2");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().from<std::string>("/a").to<int>("/foo/i").transformWith(
            &gmx::fromStdString<int>);
    transform.rules()->addRule().from<std::string>("/b").to<int>("/foo/j").transformWith(
            &gmx::fromStdString<int>);

    testTransform(input, transform);
}


TEST_F(TreeValueTransformTest, ObjectFromString)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "1 2");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().from<std::string>("/a").toObject("/foo").transformWith(
            [](gmx::KeyValueTreeObjectBuilder* builder, const std::string& value) {
                std::vector<std::string> values = gmx::splitString(value);
                builder->addValue<int>("a", gmx::fromString<int>(values[0]));
                builder->addValue<int>("b", gmx::fromString<int>(values[1]));
            });

    testTransform(input, transform);
}

TEST_F(TreeValueTransformTest, ObjectFromMultipleStrings)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "1");
    builder.rootObject().addValue<std::string>("b", "2 3");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().from<std::string>("/a").to<int>("/foo/a").transformWith(
            &gmx::fromStdString<int>);
    transform.rules()->addRule().from<std::string>("/b").toObject("/foo").transformWith(
            [](gmx::KeyValueTreeObjectBuilder* builder, const std::string& value) {
                std::vector<std::string> values = gmx::splitString(value);
                builder->addValue<int>("b", gmx::fromString<int>(values[0]));
                builder->addValue<int>("c", gmx::fromString<int>(values[1]));
            });

    testTransform(input, transform);
}

TEST_F(TreeValueTransformTest, ScopedTransformRules)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "1");
    builder.rootObject().addValue<std::string>("b", "2");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    auto                         scope = transform.rules()->scopedTransform("/foo");
    scope.rules()->addRule().from<std::string>("/a").to<int>("/i").transformWith(&gmx::fromStdString<int>);
    scope.rules()->addRule().from<std::string>("/b").to<int>("/j").transformWith(&gmx::fromStdString<int>);

    testTransform(input, transform);
}

TEST_F(TreeValueTransformTest, CanAssignUserMultiValue)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("numValues", "3");
    auto arrayBuilder = builder.rootObject().addObjectArray("values");
    for (int i = 0; i < 3; ++i)
    {
        auto valueABuilder = arrayBuilder.addObject();
        auto valueBBuilder = arrayBuilder.addObject();
        auto keyA          = gmx::formatString("TestA-%d", i);
        auto keyB          = gmx::formatString("TestB-%d", i);
        valueABuilder.addValue<std::string>(keyA, std::to_string(i + 1));
        valueBBuilder.addValue<std::string>(keyB, std::to_string(i + 2));
    }
    auto input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
    transform.rules()->addRule().from<std::string>("/numValues").to<int>("/i").transformWith(&gmx::fromStdString<int>);
    for (int i = 0; i < 3; ++i)
    {
        auto keyA    = gmx::formatString("/TestA-%d", i);
        auto keyB    = gmx::formatString("/TestB-%d", i);
        auto resultA = gmx::formatString("/j-%d", i);
        auto resultB = gmx::formatString("/j-%d", i + 2);
        transform.rules()->addRule().from<std::string>(keyA).to<int>(resultA).transformWith(
                &gmx::fromStdString<int>);
        transform.rules()->addRule().from<std::string>(keyB).to<int>(resultB).transformWith(
                &gmx::fromStdString<int>);
    }
    testTransform(input, transform);
}

/********************************************************************
 * Tests for errors
 */

TEST(TreeValueTransformErrorTest, ConversionError)
{
    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<std::string>("a", "foo");
    gmx::KeyValueTreeObject input = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().from<std::string>("/a").to<int>("/i").transformWith(
            &gmx::fromStdString<int>);

    EXPECT_THROW_GMX(transform.transform(input, nullptr), gmx::InvalidInputError);
}

} // namespace
} // namespace test
} // namespace gmx
