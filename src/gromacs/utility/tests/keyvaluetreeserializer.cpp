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

#include "gromacs/utility/keyvaluetreeserializer.h"

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{
namespace
{

//! Dummy function to raise assert for use of unimplemented function.
void raiseAssert()
{
    GMX_RELEASE_ASSERT(false, "Unimplemented function");
}

class RefDataSerializer : public gmx::ISerializer
{
public:
    RefDataSerializer(gmx::test::TestReferenceChecker* parentChecker, const char* id) :
        checker_(parentChecker->checkCompound("SerializedData", id))
    {
    }

    bool reading() const override { return false; }

    void doBool(bool* value) override { checker_.checkBoolean(*value, nullptr); }
    void doUChar(unsigned char* value) override { checker_.checkUChar(*value, nullptr); }
    void doChar(char* value) override { checker_.checkChar(*value, nullptr); }
    void doUShort(unsigned short* /* value */) override { raiseAssert(); }
    void doInt(int* value) override { checker_.checkInteger(*value, nullptr); }
    void doInt32(int32_t* value) override { checker_.checkInt32(*value, nullptr); }
    void doInt64(int64_t* value) override { checker_.checkInt64(*value, nullptr); }
    void doFloat(float* value) override { checker_.checkFloat(*value, nullptr); }
    void doDouble(double* value) override { checker_.checkDouble(*value, nullptr); }
    void doString(std::string* value) override { checker_.checkString(*value, nullptr); }
    void doOpaque(char* /* value */, std::size_t /* size */) override { raiseAssert(); }
    void doReal(real* /* value */) override { raiseAssert(); }
    void doIvec(ivec* /* value */) override { raiseAssert(); }
    void doRvec(rvec* /* value */) override { raiseAssert(); }

private:
    gmx::test::TestReferenceChecker checker_;
};

class KeyValueTreeSerializerTest : public ::testing::Test
{
public:
    void runTest()
    {
        gmx::KeyValueTreeObject         input(builder_.build());
        gmx::test::TestReferenceData    data;
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(input, "Input");
        {
            RefDataSerializer serializer(&checker, "Stream");
            gmx::serializeKeyValueTree(input, &serializer);
        }
        std::vector<char> buffer = serializeTree(input);
        {
            gmx::InMemoryDeserializer deserializer(buffer, false);
            gmx::KeyValueTreeObject   output = gmx::deserializeKeyValueTree(&deserializer);
            checker.checkKeyValueTreeObject(output, "Input");
        }
    }

    gmx::KeyValueTreeBuilder builder_;

private:
    static std::vector<char> serializeTree(const gmx::KeyValueTreeObject& tree)
    {
        gmx::InMemorySerializer serializer;
        gmx::serializeKeyValueTree(tree, &serializer);
        return serializer.finishAndGetBuffer();
    }
};

TEST_F(KeyValueTreeSerializerTest, EmptyTree)
{
    runTest();
}

TEST_F(KeyValueTreeSerializerTest, SimpleObject)
{
    builder_.rootObject().addValue<char>("c", 'x');
    builder_.rootObject().addValue<unsigned char>("u", 0170);
    builder_.rootObject().addValue<int>("foo", 1);
    builder_.rootObject().addValue<int32_t>("foo32", 1);
    builder_.rootObject().addValue<int64_t>("foo64", 1);
    builder_.rootObject().addValue<std::string>("bar", "a");
    builder_.rootObject().addValue<float>("f", 1.5);
    builder_.rootObject().addValue<double>("d", 2.5);
    runTest();
}

TEST_F(KeyValueTreeSerializerTest, ObjectWithArrays)
{
    auto arr1 = builder_.rootObject().addUniformArray<int>("a");
    arr1.addValue(1);
    arr1.addValue(2);
    auto arr2 = builder_.rootObject().addUniformArray<std::string>("b");
    arr2.addValue("foo");
    arr2.addValue("bar");
    runTest();
}

TEST_F(KeyValueTreeSerializerTest, ObjectWithObjects)
{
    auto obj1 = builder_.rootObject().addObject("obj");
    obj1.addValue<int>("a", 1);
    obj1.addValue<std::string>("b", "foo");
    auto obj2 = builder_.rootObject().addObject("obj2");
    obj2.addValue<int>("c", 2);
    obj2.addValue<std::string>("d", "bar");
    builder_.rootObject().addValue<int>("foo", 3);
    runTest();
}

} // namespace
} // namespace test
} // namespace gmx
