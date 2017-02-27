/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreeserializer.h"

#include <gtest/gtest.h>

#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "testutils/refdata.h"

namespace
{

class RefDataWriteSerializer : public gmx::ISerializer
{
    public:
        RefDataWriteSerializer(gmx::test::TestReferenceChecker *parentChecker,
                               const char                      *id)
            : checker_(parentChecker->checkCompound("SerializedData", id))
        {
        }

        virtual bool reading() const { return false; }

        virtual void doUChar(unsigned char *value)
        {
            checker_.checkUChar(*value, nullptr);
        }
        virtual void doInt(int *value)
        {
            checker_.checkInteger(*value, nullptr);
        }
        virtual void doInt64(gmx_int64_t *value)
        {
            checker_.checkInt64(*value, nullptr);
        }
        virtual void doFloat(float *value)
        {
            checker_.checkFloat(*value, nullptr);
        }
        virtual void doDouble(double *value)
        {
            checker_.checkDouble(*value, nullptr);
        }
        virtual void doString(std::string *value)
        {
            checker_.checkString(*value, nullptr);
        }

    private:
        gmx::test::TestReferenceChecker checker_;
};

class RefDataReadSerializer : public gmx::ISerializer
{
    public:
        RefDataReadSerializer(gmx::test::TestReferenceChecker *parentChecker,
                              const char                      *id)
            : checker_(parentChecker->checkCompound("SerializedData", id))
        {
        }

        virtual bool reading() const { return true; }

        virtual void doUChar(unsigned char *value)
        {
            *value = checker_.readUChar(nullptr);
        }
        virtual void doInt(int *value)
        {
            *value = checker_.readInteger(nullptr);
        }
        virtual void doInt64(gmx_int64_t *value)
        {
            *value = checker_.readInt64(nullptr);
        }
        virtual void doFloat(float *value)
        {
            *value = checker_.readFloat(nullptr);
        }
        virtual void doDouble(double *value)
        {
            *value = checker_.readDouble(nullptr);
        }
        virtual void doString(std::string *value)
        {
            *value = checker_.readString(nullptr);
        }

    private:
        gmx::test::TestReferenceChecker checker_;
};

class KeyValueTreeSerializerTest : public ::testing::Test
{
    public:
        void runTest()
        {
            gmx::KeyValueTreeObject           input(builder_.build());
            gmx::test::TestReferenceData      data;
            gmx::test::TestReferenceChecker   checker(data.rootChecker());
            checker.checkKeyValueTreeObject(input, "Input");
            {
                RefDataWriteSerializer        serializer(&checker, "Stream");
                gmx::serializeKeyValueTree(input, &serializer);
            }
            {
                RefDataReadSerializer         serializer(&checker, "Stream");
                gmx::KeyValueTreeObject       output
                    = gmx::deserializeKeyValueTree(&serializer);
                checker.checkKeyValueTreeObject(output, "Input");
            }
        }

        gmx::KeyValueTreeBuilder builder_;
};

TEST_F(KeyValueTreeSerializerTest, EmptyTree)
{
    runTest();
}

TEST_F(KeyValueTreeSerializerTest, SimpleObject)
{
    builder_.rootObject().addValue<int>("foo", 1);
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
