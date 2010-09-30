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

#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/errorreporting/emptyerrorreporter.h"
#include "gromacs/options/abstractoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionstoragetemplate.h"
#include "gromacs/options/optionsassigner.h"

class MockOption;

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
         * \param[in] options    Options object.
         * \retval 0 on success.
         */
        int init(const MockOption &settings, gmx::Options *options);

        /*! \brief
         * Calls addValue() in the base class and expects it to succeed.
         */
        void addValue(const std::string &value)
        {
            EXPECT_EQ(0, MyBase::addValue(value));
        }
        /*! \brief
         * Calls addValue() in the base class and expects it to fail.
         */
        void addValueExpectFail(const std::string &value)
        {
            EXPECT_NE(0, MyBase::addValue(value));
        }
        /*! \brief
         * Calls processAll() in the base class.
         */
        int processAllBase(gmx::AbstractErrorReporter *errors)
        {
            return MyBase::processAll(errors);
        }
        /*! \brief
         * Calls addValue("dummy") in the base class and expects it to succeed.
         */
        void addDummyValue()
        {
            addValue("dummy");
        }
        /*! \brief
         * Calls addValue("dummy") in the base class and expects it to fail.
         */
        void addDummyValueExpectFail()
        {
            addValueExpectFail("dummy");
        }
        /*! \brief
         * Calls setFlag(efSet).
         */
        void setOption()
        {
            setFlag(gmx::efSet);
        }

        virtual const char *typeString() const { return "mock"; }
        virtual std::string formatValue(int /*i*/) const { return ""; }

        MOCK_METHOD2(convertValue, int(const std::string &value,
                                       gmx::AbstractErrorReporter *errors));
        MOCK_METHOD2(processSet, int(int nvalues,
                                     gmx::AbstractErrorReporter *errors));
        MOCK_METHOD1(processAll, int(gmx::AbstractErrorReporter *errors));
};

/*! \internal \brief
 * Specifies an option that has a mock storage object for unit testing.
 *
 * \ingroup module_options
 */
class MockOption : public gmx::OptionTemplate<std::string, MockOption>
{
    public:
        //! Initializes an option with the given name.
        explicit MockOption(const char *name)
            : MyBase(name), _storagePtr(NULL)
        {
        }

        //! Sets the required flags to support storage that may not add values.
        MyClass &mayNotAddValues()
        { setFlag(gmx::efConversionMayNotAddValues); return me(); }
        //! Sets an output pointer to give access to the created storage object.
        MyClass &storageObject(MockOptionStorage **storagePtr)
        { _storagePtr = storagePtr; return me(); }

    private:
        virtual int createDefaultStorage(gmx::Options *options,
                                         gmx::AbstractOptionStorage **storage) const
        {
            int rc = gmx::createOptionStorage<MockOption, MockOptionStorage>(this, options, storage);
            if (_storagePtr != NULL)
            {
                *_storagePtr = static_cast<MockOptionStorage *>(*storage);
            }
            return rc;
        }

        MockOptionStorage     **_storagePtr;
};

int MockOptionStorage::init(const MockOption &settings, gmx::Options *options)
{
    using ::testing::_;
    using ::testing::DoAll;
    using ::testing::Invoke;
    using ::testing::Return;
    using ::testing::WithArg;
    ON_CALL(*this, convertValue(_, _))
        .WillByDefault(DoAll(WithArg<0>(Invoke(this, &MockOptionStorage::addValue)),
                             Return(0)));
    ON_CALL(*this, processAll(_))
        .WillByDefault(Invoke(this, &MockOptionStorage::processAllBase));
    return MyBase::init(settings, options);
}

namespace
{

/*
 * Tests that finish() can set a required option even if the user has not
 * provided it.
 */
TEST(AbstractOptionStorageTest, HandlesSetInFinish)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionStorage          *mock;
    options.addOption(MockOption("name").storageObject(&mock).required()
                          .storeVector(&values));

    {
        ::testing::InSequence dummy;
        using ::testing::_;
        using ::testing::DoAll;
        using ::testing::Return;
        using ::testing::InvokeWithoutArgs;
        EXPECT_CALL(*mock, processAll(_))
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::setOption),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            Return(0)));
    }

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

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
    MockOptionStorage          *mock;
    options.addOption(MockOption("name").storageObject(&mock).mayNotAddValues()
                          .storeVector(&values).multiValue());

    {
        ::testing::InSequence dummy;
        using ::testing::_;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a", _));
        EXPECT_CALL(*mock, convertValue("b", _))
            .WillOnce(Return(0));
        EXPECT_CALL(*mock, convertValue("c", _));
        EXPECT_CALL(*mock, processSet(2, _));
        EXPECT_CALL(*mock, processAll(_));
    }

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("name"));
    EXPECT_EQ(0, assigner.appendValue("a"));
    EXPECT_EQ(0, assigner.appendValue("b"));
    EXPECT_EQ(0, assigner.appendValue("c"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

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
    MockOptionStorage          *mock;
    options.addOption(MockOption("name").storageObject(&mock)
                          .storeVector(&values).multiValue());

    {
        ::testing::InSequence dummy;
        using ::testing::_;
        using ::testing::DoAll;
        using ::testing::InvokeWithoutArgs;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a", _));
        EXPECT_CALL(*mock, convertValue("b", _))
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            Return(0)));
        EXPECT_CALL(*mock, processSet(3, _));
        EXPECT_CALL(*mock, processAll(_));
    }

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("name"));
    EXPECT_EQ(0, assigner.appendValue("a"));
    EXPECT_EQ(0, assigner.appendValue("b"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

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
    MockOptionStorage          *mock;
    options.addOption(MockOption("name").storageObject(&mock)
                          .storeVector(&values).valueCount(2));

    {
        ::testing::InSequence dummy;
        using ::testing::_;
        using ::testing::DoAll;
        using ::testing::InvokeWithoutArgs;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a", _));
        EXPECT_CALL(*mock, convertValue("b", _))
            .WillOnce(DoAll(InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValue),
                            InvokeWithoutArgs(mock, &MockOptionStorage::addDummyValueExpectFail),
                            Return(gmx::eeInvalidInput)));
        EXPECT_CALL(*mock, processSet(2, _));
        EXPECT_CALL(*mock, processAll(_));
    }

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("name"));
    EXPECT_EQ(0, assigner.appendValue("a"));
    EXPECT_NE(0, assigner.appendValue("b"));
    EXPECT_NE(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ("a", values[0]);
    EXPECT_EQ("dummy", values[1]);
}

/*
 * Tests that the storage object is properly invoked even if no output values
 * should be produced.
 */
TEST(AbstractOptionStorageTest, AllowsEmptyValues)
{
    gmx::Options                options(NULL, NULL);
    std::vector<std::string>    values;
    MockOptionStorage          *mock;
    options.addOption(MockOption("name").storageObject(&mock).mayNotAddValues()
                          .storeVector(&values).valueCount(0));

    {
        ::testing::InSequence dummy;
        using ::testing::_;
        using ::testing::DoAll;
        using ::testing::InvokeWithoutArgs;
        using ::testing::Return;
        EXPECT_CALL(*mock, convertValue("a", _))
            .WillOnce(Return(0));
        EXPECT_CALL(*mock, processSet(0, _));
        EXPECT_CALL(*mock, processAll(_));
    }

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("name"));
    EXPECT_EQ(0, assigner.appendValue("a"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    ASSERT_EQ(0U, values.size());
}

} // namespace
