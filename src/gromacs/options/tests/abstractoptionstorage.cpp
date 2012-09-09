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
 * Tests proper handling of option storage.
 *
 * These tests check that methods in storage objects are called properly in all
 * situations, and also that the OptionStorageTemplate class behaves properly.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include <vector>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionstoragetemplate.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/testexceptions.h"

namespace
{

class MockOption;
class MockOptionStorage;

class MockOptionInfo : public gmx::OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit MockOptionInfo(MockOptionStorage *option);

        MockOptionStorage &option();
};

/*! \internal \brief
 * Mock implementation of an option storage class for unit testing.
 *
 * Provides facilities for checking that correct methods are called, and for
 * controlling how they add values using the base class methods.
 *
 * \ingroup module_options
 */
class MockOptionStorage : public gmx::OptionStorageTemplate<std::string>
{
    public:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings   Storage settings.
         */
        MockOptionStorage(const MockOption &settings);

        /*! \brief
         * Calls addValue("dummy") in the base class.
         */
        void addDummyValue()
        {
            addValue("dummy");
        }
        // using MyBase::markAsSet;
        // using MyBase::addValue;
        // using MyBase::commitValues;
        // "using" is correct but MSVC gives error C2248. Workaround:
        void markAsSet()
        {
            MyBase::markAsSet();
        }
        void addValue(const std::string &value)
        {
            MyBase::addValue(value);
        }
        void commitValues()
        {
            MyBase::commitValues();
        }

        virtual gmx::OptionInfo &optionInfo() { return info_; }
        // These are not used.
        virtual const char *typeString() const { return "mock"; }
        virtual std::string formatSingleValue(const std::string &/*value*/) const
        {
            return "";
        }

        MOCK_METHOD1(convertValue, void(const std::string &value));
        MOCK_METHOD1(processSetValues, void(ValueList *values));
        MOCK_METHOD0(processAll, void());

    private:
        MockOptionInfo          info_;
};

/*! \internal \brief
 * Specifies an option that has a mock storage object for unit testing.
 *
 * \ingroup module_options
 */
class MockOption : public gmx::OptionTemplate<std::string, MockOption>
{
    public:
        //! OptionInfo subclass corresponding to this option type.
        typedef MockOptionInfo InfoType;

        //! Initializes an option with the given name.
        explicit MockOption(const char *name)
            : MyBase(name)
        {
        }

    private:
        virtual gmx::AbstractOptionStoragePointer createStorage() const
        {
            return gmx::AbstractOptionStoragePointer(new MockOptionStorage(*this));
        }
};

MockOptionStorage::MockOptionStorage(const MockOption &settings)
    : MyBase(settings), info_(this)
{
    using ::testing::_;
    using ::testing::Invoke;
    using ::testing::WithArg;
    ON_CALL(*this, convertValue(_))
        .WillByDefault(WithArg<0>(Invoke(this, &MockOptionStorage::addValue)));
}

MockOptionInfo::MockOptionInfo(MockOptionStorage *option)
    : gmx::OptionInfo(option)
{
}

MockOptionStorage &MockOptionInfo::option()
{
    return static_cast<MockOptionStorage &>(gmx::OptionInfo::option());
}

/*
 * Tests that finish() can set a required option even if the user has not
 * provided it.
 */
TEST(AbstractOptionStorageTest, HandlesSetInFinish)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionInfo *info = options.addOption(
            MockOption("name").required().storeVector(&values));
    MockOptionStorage *mock = &info->option();
    {
        ::testing::InSequence dummy;
        using ::testing::DoAll;
        using ::testing::InvokeWithoutArgs;
        EXPECT_CALL(*mock, processAll())
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::markAsSet),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            InvokeWithoutArgs(mock, &MockOptionStorage::commitValues)));
    }

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(1U, values.size());
    EXPECT_EQ("dummy", values[0]);
}

/*
 * Tests that storage works if the storage object does not add a value in a
 * call to appendValue().
 */
TEST(AbstractOptionStorageTest, HandlesValueRemoval)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionInfo *info = options.addOption(
            MockOption("name").storeVector(&values).multiValue());
    MockOptionStorage *mock = &info->option();
    {
        ::testing::InSequence dummy;
        using ::testing::ElementsAre;
        using ::testing::Pointee;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a"));
        EXPECT_CALL(*mock, convertValue("b"))
            .WillOnce(Return());
        EXPECT_CALL(*mock, convertValue("c"));
        EXPECT_CALL(*mock, processSetValues(Pointee(ElementsAre("a", "c"))));
        EXPECT_CALL(*mock, processAll());
    }

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("name"));
    EXPECT_NO_THROW(assigner.appendValue("a"));
    EXPECT_NO_THROW(assigner.appendValue("b"));
    EXPECT_NO_THROW(assigner.appendValue("c"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ("a", values[0]);
    EXPECT_EQ("c", values[1]);
}

/*
 * Tests that storage works if the storage object adds more than one value in
 * one call to appendValue().
 */
TEST(AbstractOptionStorageTest, HandlesValueAddition)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionInfo *info = options.addOption(
            MockOption("name").storeVector(&values).multiValue());
    MockOptionStorage *mock = &info->option();
    {
        ::testing::InSequence dummy;
        using ::testing::DoAll;
        using ::testing::ElementsAre;
        using ::testing::InvokeWithoutArgs;
        using ::testing::Pointee;
        EXPECT_CALL(*mock, convertValue("a"));
        EXPECT_CALL(*mock, convertValue("b"))
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue)));
        EXPECT_CALL(*mock, processSetValues(Pointee(ElementsAre("a", "dummy", "dummy"))));
        EXPECT_CALL(*mock, processAll());
    }

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("name"));
    EXPECT_NO_THROW(assigner.appendValue("a"));
    EXPECT_NO_THROW(assigner.appendValue("b"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(3U, values.size());
    EXPECT_EQ("a", values[0]);
    EXPECT_EQ("dummy", values[1]);
    EXPECT_EQ("dummy", values[2]);
}

/*
 * Tests that storage works if the storage object adds more than one value in
 * one call to appendValue(), and this results in too many values.
 */
TEST(AbstractOptionStorageTest, HandlesTooManyValueAddition)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionInfo *info = options.addOption(
            MockOption("name").storeVector(&values).valueCount(2));
    MockOptionStorage *mock = &info->option();
    {
        ::testing::InSequence dummy;
        using ::testing::DoAll;
        using ::testing::ElementsAre;
        using ::testing::InvokeWithoutArgs;
        using ::testing::Pointee;
        EXPECT_CALL(*mock, convertValue("a"));
        EXPECT_CALL(*mock, convertValue("b"))
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue)));
        EXPECT_CALL(*mock, processAll());
    }

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("name"));
    EXPECT_NO_THROW(assigner.appendValue("a"));
    EXPECT_THROW(assigner.appendValue("b"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_TRUE(values.empty());
}

/*
 * Tests that the storage object is properly invoked even if no output values
 * should be produced.
 */
TEST(AbstractOptionStorageTest, AllowsEmptyValues)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionInfo *info = options.addOption(
            MockOption("name").storeVector(&values).valueCount(0));
    MockOptionStorage *mock = &info->option();
    {
        ::testing::InSequence dummy;
        using ::testing::DoAll;
        using ::testing::ElementsAre;
        using ::testing::Pointee;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a"))
            .WillOnce(Return());
        EXPECT_CALL(*mock, processSetValues(Pointee(ElementsAre())));
        EXPECT_CALL(*mock, processAll());
    }

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.startOption("name"));
    EXPECT_NO_THROW(assigner.appendValue("a"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(0U, values.size());
}

} // namespace
