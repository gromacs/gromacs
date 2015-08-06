/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements classes and functions from refdata.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "refdata.h"

#include <cstdlib>

#include <limits>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata-checkers.h"
#include "testutils/refdata-impl.h"
#include "testutils/refdata-xml.h"
#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * TestReferenceData::Impl declaration
 */

namespace internal
{

/*! \internal \brief
 * Private implementation class for TestReferenceData.
 *
 * \ingroup module_testutils
 */
class TestReferenceDataImpl
{
    public:
        //! Initializes a checker in the given mode.
        TestReferenceDataImpl(ReferenceDataMode mode, bool bSelfTestMode);

        //! Performs final reference data processing when test ends.
        void onTestEnd();

        //! Full path of the reference data file.
        std::string             fullFilename_;
        /*! \brief
         * Root entry for the reference data.
         *
         * If null after construction, the reference data is not present.
         */
        ReferenceDataEntry::EntryPointer  rootEntry_;
        /*! \brief
         * Whether the reference data is being written (true) or compared
         * (false).
         */
        bool                    bWrite_;
        //! `true` if self-testing (enables extra failure messages).
        bool                    bSelfTestMode_;
        /*! \brief
         * Whether any reference checkers have been created for this data.
         */
        bool                    bInUse_;
};

}       // namespace internal

/********************************************************************
 * Internal helpers
 */

namespace
{

//! Convenience typedef for a smart pointer to TestReferenceDataImpl.
typedef boost::shared_ptr<internal::TestReferenceDataImpl>
    TestReferenceDataImplPointer;

/*! \brief
 * Global reference data instance.
 *
 * The object is destructed between tests using an event listener.
 */
TestReferenceDataImplPointer g_referenceData;
//! Global reference data mode set with setReferenceDataMode().
// TODO: Make this a real enum; requires solving a TODO in StringOption.
int                          g_referenceDataMode = erefdataCompare;

//! Returns the global reference data mode.
ReferenceDataMode getReferenceDataMode()
{
    return static_cast<ReferenceDataMode>(g_referenceDataMode);
}

//! Returns a reference to the global reference data object.
TestReferenceDataImplPointer getReferenceDataInstance()
{
    GMX_RELEASE_ASSERT(!g_referenceData,
                       "Test cannot create multiple TestReferenceData instances");
    g_referenceData.reset(new internal::TestReferenceDataImpl(getReferenceDataMode(), false));
    return g_referenceData;
}

//! Handles reference data creation for self-tests.
TestReferenceDataImplPointer initReferenceDataInstanceForSelfTest(ReferenceDataMode mode)
{
    if (g_referenceData)
    {
        g_referenceData->onTestEnd();
        g_referenceData.reset();
    }
    g_referenceData.reset(new internal::TestReferenceDataImpl(mode, true));
    return g_referenceData;
}

class ReferenceDataTestEventListener : public ::testing::EmptyTestEventListener
{
    public:
        virtual void OnTestEnd(const ::testing::TestInfo &)
        {
            if (g_referenceData)
            {
                GMX_RELEASE_ASSERT(g_referenceData.unique(),
                                   "Test leaked TestRefeferenceData objects");
                g_referenceData->onTestEnd();
                g_referenceData.reset();
            }
        }

        // Frees internal buffers allocated by libxml2.
        virtual void OnTestProgramEnd(const ::testing::UnitTest &)
        {
            cleanupReferenceData();
        }
};

}       // namespace

void initReferenceData(IOptionsContainer *options)
{
    // Needs to correspond to the enum order in refdata.h.
    const char *const refDataEnum[] = { "check", "create", "update" };
    options->addOption(
            StringOption("ref-data").enumValue(refDataEnum)
                .defaultEnumIndex(0)
                .storeEnumIndex(&g_referenceDataMode)
                .description("Operation mode for tests that use reference data"));
    ::testing::UnitTest::GetInstance()->listeners().Append(
            new ReferenceDataTestEventListener);
}

/********************************************************************
 * TestReferenceDataImpl definition
 */

namespace internal
{

TestReferenceDataImpl::TestReferenceDataImpl(
        ReferenceDataMode mode, bool bSelfTestMode)
    : bWrite_(false), bSelfTestMode_(bSelfTestMode), bInUse_(false)
{
    const std::string dirname =
        bSelfTestMode
        ? TestFileManager::getGlobalOutputTempDirectory()
        : TestFileManager::getInputDataDirectory();
    const std::string filename = TestFileManager::getTestSpecificFileName(".xml");
    fullFilename_ = Path::join(dirname, "refdata", filename);

    bWrite_ = true;
    if (mode != erefdataUpdateAll)
    {
        if (File::exists(fullFilename_))
        {
            bWrite_ = false;
        }
        else if (mode == erefdataCompare)
        {
            bWrite_ = false;
            return;
        }
    }
    if (bWrite_)
    {
        rootEntry_ = ReferenceDataEntry::createRoot();
    }
    else
    {
        rootEntry_ = readReferenceDataFile(fullFilename_);
    }
}

void TestReferenceDataImpl::onTestEnd()
{
    if (bWrite_ && bInUse_ && rootEntry_)
    {
        std::string dirname = Path::getParentPath(fullFilename_);
        if (!Directory::exists(dirname))
        {
            if (Directory::create(dirname) != 0)
            {
                GMX_THROW(TestException("Creation of reference data directory failed: " + dirname));
            }
        }
        writeReferenceDataFile(fullFilename_, *rootEntry_);
    }
}

}       // namespace internal


/********************************************************************
 * TestReferenceChecker::Impl
 */

/*! \internal \brief
 * Private implementation class for TestReferenceChecker.
 *
 * \ingroup module_testutils
 */
class TestReferenceChecker::Impl
{
    public:
        //! String constant for naming XML elements for boolean values.
        static const char * const    cBooleanNodeName;
        //! String constant for naming XML elements for string values.
        static const char * const    cStringNodeName;
        //! String constant for naming XML elements for integer values.
        static const char * const    cIntegerNodeName;
        //! String constant for naming XML elements for int64 values.
        static const char * const    cInt64NodeName;
        //! String constant for naming XML elements for unsigned int64 values.
        static const char * const    cUInt64NodeName;
        //! String constant for naming XML elements for floating-point values.
        static const char * const    cRealNodeName;
        //! String constant for naming XML attribute for value identifiers.
        static const char * const    cIdAttrName;
        //! String constant for naming compounds for vectors.
        static const char * const    cVectorType;
        //! String constant for naming compounds for sequences.
        static const char * const    cSequenceType;
        //! String constant for value identifier for sequence length.
        static const char * const    cSequenceLengthName;

        //! Creates a checker that does nothing.
        explicit Impl(bool bWrite);
        //! Creates a checker with a given root entry.
        Impl(const std::string &path, ReferenceDataEntry *rootEntry, bool bWrite,
             bool bSelfTestMode, const FloatingPointTolerance &defaultTolerance);

        //! Returns the path of this checker with \p id appended.
        std::string appendPath(const char *id) const;

        //! Returns whether an iterator returned by findEntry() is valid.
        bool isValidEntry(const ReferenceDataEntry::ChildIterator &iter) const
        {
            if (rootEntry_ == NULL)
            {
                return false;
            }
            return iter != rootEntry_->children().end();
        }

        /*! \brief
         * Finds a reference data entry.
         *
         * \param[in]  type   Type of entry to find (can be NULL, in which case
         *      any type is matched).
         * \param[in]  id     Unique identifier of the entry (can be NULL, in
         *      which case the next entry without an id is matched).
         * \returns    Matching entry, or an invalid iterator (see
         *      isValidEntry()) if no matching entry found.
         *
         * Searches for an entry in the reference data that matches the given
         * \p name and \p id.  Searching starts from the entry that follows the
         * previously matched entry (relevant for performance, and if there are
         * nodes without ids).  Note that the match pointer is not updated by
         * this method.
         */
        ReferenceDataEntry::ChildIterator
        findEntry(const char *type, const char *id) const;
        /*! \brief
         * Finds/creates a reference data entry to match against.
         *
         * \param[in]  type   Type of entry to find.
         * \param[in]  id     Unique identifier of the entry (can be NULL, in
         *      which case the next entry without an id is matched).
         * \param[out] bFound Whether the entry was found (false if the entry
         *      was created in write mode).
         * \returns Matching entry, or NULL if no matching entry found
         *      (NULL is never returned in write mode).
         *
         * Finds an entry using findEntry() and updates the match pointer if a
         * match is found.  If a match is not found, the method returns NULL in
         * read mode and creates a new entry in write mode.
         */
        ReferenceDataEntry *findOrCreateEntry(const char *type, const char *id,
                                              bool *bFound);
        /*! \brief
         * Helper method for checking a reference data value.
         *
         * \param[in]  name   Type of entry to find.
         * \param[in]  id     Unique identifier of the entry (can be NULL, in
         *     which case the next entry without an id is matched).
         * \param[in]  checker  Checker that provides logic specific to the
         *     type of the entry.
         * \returns    Whether the reference data matched, including details
         *     of the mismatch if the comparison failed.
         * \throws     TestException if there is a problem parsing the
         *     reference data.
         *
         * Performs common tasks in checking a reference value, such as
         * finding or creating the correct entry.
         * Caller needs to provide a checker object that provides the string
         * value for a newly created entry and performs the actual comparison
         * against a found entry.
         */
        ::testing::AssertionResult
        processItem(const char *name, const char *id,
                    const IReferenceDataEntryChecker &checker);
        /*! \brief
         * Whether the checker should ignore all validation calls.
         *
         * This is used to ignore any calls within compounds for which
         * reference data could not be found, such that only one error is
         * issued for the missing compound, instead of every individual value.
         */
        bool shouldIgnore() const { return rootEntry_ == NULL; }

        //! Default floating-point comparison tolerance.
        FloatingPointTolerance  defaultTolerance_;
        /*! \brief
         * Human-readable path to the root node of this checker.
         *
         * For the root checker, this will be "/", and for each compound, the
         * id of the compound is added.  Used for reporting comparison
         * mismatches.
         */
        std::string             path_;
        /*! \brief
         * Current node under which reference data is searched.
         *
         * Points to either the TestReferenceDataImpl::rootEntry_, or to
         * a compound entry in the tree rooted at that entry.
         *
         * Can be NULL, in which case this checker does nothing (doesn't even
         * report errors, see shouldIgnore()).
         */
        ReferenceDataEntry     *rootEntry_;
        /*! \brief
         * Iterator to a child of \a rootEntry_ that was last found.
         *
         * If isValidEntry() returns false, no entry has been found yet.
         * After every check, is updated to point to the entry that was used
         * for the check.
         * Subsequent checks start the search for the matching node on this
         * node.
         */
        ReferenceDataEntry::ChildIterator prevFoundNode_;
        /*! \brief
         * Whether the reference data is being written (true) or compared
         * (false).
         */
        bool                    bWrite_;
        //! `true` if self-testing (enables extra failure messages).
        bool                    bSelfTestMode_;
        /*! \brief
         * Current number of unnamed elements in a sequence.
         *
         * It is the index of the next added unnamed element.
         */
        int                     seqIndex_;
};

const char *const TestReferenceChecker::Impl::cBooleanNodeName    = "Bool";
const char *const TestReferenceChecker::Impl::cStringNodeName     = "String";
const char *const TestReferenceChecker::Impl::cIntegerNodeName    = "Int";
const char *const TestReferenceChecker::Impl::cInt64NodeName      = "Int64";
const char *const TestReferenceChecker::Impl::cUInt64NodeName     = "UInt64";
const char *const TestReferenceChecker::Impl::cRealNodeName       = "Real";
const char *const TestReferenceChecker::Impl::cIdAttrName         = "Name";
const char *const TestReferenceChecker::Impl::cVectorType         = "Vector";
const char *const TestReferenceChecker::Impl::cSequenceType       = "Sequence";
const char *const TestReferenceChecker::Impl::cSequenceLengthName = "Length";


TestReferenceChecker::Impl::Impl(bool bWrite)
    : defaultTolerance_(defaultRealTolerance()),
      rootEntry_(NULL), bWrite_(bWrite),
      bSelfTestMode_(false), seqIndex_(0)
{
}


TestReferenceChecker::Impl::Impl(const std::string &path, ReferenceDataEntry *rootEntry,
                                 bool bWrite, bool bSelfTestMode,
                                 const FloatingPointTolerance &defaultTolerance)
    : defaultTolerance_(defaultTolerance), path_(path + "/"),
      rootEntry_(rootEntry), prevFoundNode_(rootEntry->children().end()),
      bWrite_(bWrite), bSelfTestMode_(bSelfTestMode), seqIndex_(0)
{
}


std::string
TestReferenceChecker::Impl::appendPath(const char *id) const
{
    std::string printId = (id != NULL) ? id : formatString("[%d]", seqIndex_);
    return path_ + printId;
}


ReferenceDataEntry::ChildIterator
TestReferenceChecker::Impl::findEntry(const char *type, const char *id) const
{
    const ReferenceDataEntry::ChildList &children = rootEntry_->children();
    if (children.empty())
    {
        return children.end();
    }
    ReferenceDataEntry::ChildIterator  node  = prevFoundNode_;
    bool                               bWrap = true;
    if (node != children.end())
    {
        if (id == NULL)
        {
            if ((*node)->id().empty())
            {
                if (type == NULL || (*node)->type() == type)
                {
                    bWrap = false;
                    ++node;
                    if (node == children.end())
                    {
                        return children.end();
                    }
                }
            }
        }
    }
    else
    {
        node  = children.begin();
        bWrap = false;
    }
    do
    {
        if (type == NULL || (*node)->type() == type)
        {
            if (id == NULL && (*node)->id().empty())
            {
                return node;
            }
            if (!(*node)->id().empty())
            {
                if (id != NULL && (*node)->id() == id)
                {
                    return node;
                }
            }
        }
        ++node;
        if (bWrap && node == children.end())
        {
            node = children.begin();
        }
    }
    while (node != children.end() && node != prevFoundNode_);
    return children.end();
}


ReferenceDataEntry *
TestReferenceChecker::Impl::findOrCreateEntry(const char *type,
                                              const char *id,
                                              bool       *bFound)
{
    *bFound = false;
    if (rootEntry_ == NULL)
    {
        return NULL;
    }
    ReferenceDataEntry::ChildIterator node = findEntry(type, id);
    if (isValidEntry(node))
    {
        *bFound        = true;
        prevFoundNode_ = node;
    }
    else
    {
        if (bWrite_)
        {
            ReferenceDataEntry::EntryPointer newEntry(
                    new ReferenceDataEntry(type, id));
            node           = rootEntry_->addChild(move(newEntry));
            prevFoundNode_ = node;
        }
    }
    seqIndex_ = (id == NULL) ? seqIndex_+1 : 0;

    return isValidEntry(node) ? node->get() : NULL;
}

::testing::AssertionResult
TestReferenceChecker::Impl::processItem(const char *type, const char *id,
                                        const IReferenceDataEntryChecker &checker)
{
    if (shouldIgnore())
    {
        return ::testing::AssertionSuccess();
    }
    std::string         fullId = appendPath(id);
    bool                bFound;
    ReferenceDataEntry *entry = findOrCreateEntry(type, id, &bFound);
    if (entry == NULL)
    {
        return ::testing::AssertionFailure()
               << "Reference data item " << fullId << " not found";
    }
    if (bWrite_ && !bFound)
    {
        checker.fillEntry(entry);
        return ::testing::AssertionSuccess();
    }
    else
    {
        ::testing::AssertionResult result(checker.checkEntry(*entry, fullId));
        if (bSelfTestMode_ && !result)
        {
            ReferenceDataEntry expected(type, id);
            checker.fillEntry(&expected);
            result << std::endl
            << "String value: " << expected.value() << std::endl
            << " Ref. string: " << entry->value();
        }
        return result;
    }
}


/********************************************************************
 * TestReferenceData
 */

TestReferenceData::TestReferenceData()
    : impl_(getReferenceDataInstance())
{
}


TestReferenceData::TestReferenceData(ReferenceDataMode mode)
    : impl_(initReferenceDataInstanceForSelfTest(mode))
{
}


TestReferenceData::~TestReferenceData()
{
}


bool TestReferenceData::isWriteMode() const
{
    return impl_->bWrite_;
}


TestReferenceChecker TestReferenceData::rootChecker()
{
    if (!isWriteMode() && !impl_->bInUse_ && !impl_->rootEntry_)
    {
        ADD_FAILURE() << "Reference data file not found: "
        << impl_->fullFilename_;
    }
    impl_->bInUse_ = true;
    if (!impl_->rootEntry_)
    {
        return TestReferenceChecker(new TestReferenceChecker::Impl(isWriteMode()));
    }
    return TestReferenceChecker(
            new TestReferenceChecker::Impl("", impl_->rootEntry_.get(),
                                           isWriteMode(), impl_->bSelfTestMode_,
                                           defaultRealTolerance()));
}


/********************************************************************
 * TestReferenceChecker
 */

TestReferenceChecker::TestReferenceChecker(Impl *impl)
    : impl_(impl)
{
}


TestReferenceChecker::TestReferenceChecker(const TestReferenceChecker &other)
    : impl_(new Impl(*other.impl_))
{
}


TestReferenceChecker &
TestReferenceChecker::operator=(const TestReferenceChecker &other)
{
    impl_.reset(new Impl(*other.impl_));
    return *this;
}


TestReferenceChecker::~TestReferenceChecker()
{
}


bool TestReferenceChecker::isWriteMode() const
{
    return impl_->bWrite_;
}


void TestReferenceChecker::setDefaultTolerance(
        const FloatingPointTolerance &tolerance)
{
    impl_->defaultTolerance_ = tolerance;
}


bool TestReferenceChecker::checkPresent(bool bPresent, const char *id)
{
    if (isWriteMode() || impl_->shouldIgnore())
    {
        return bPresent;
    }
    ReferenceDataEntry::ChildIterator  node   = impl_->findEntry(NULL, id);
    bool                               bFound = impl_->isValidEntry(node);
    if (bFound != bPresent)
    {
        ADD_FAILURE() << "Mismatch while checking reference data item '"
        << impl_->appendPath(id) << "'\n"
        << "Expected: " << (bPresent ? "it is present.\n" : "it is absent.\n")
        << "  Actual: " << (bFound ? "it is present." : "it is absent.");
    }
    if (bFound && bPresent)
    {
        impl_->prevFoundNode_ = node;
        return true;
    }
    return false;
}


TestReferenceChecker TestReferenceChecker::checkCompound(const char *type, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    std::string         fullId = impl_->appendPath(id);
    bool                bFound;
    ReferenceDataEntry *newNode = impl_->findOrCreateEntry(type, id, &bFound);
    if (newNode == NULL)
    {
        ADD_FAILURE() << "Reference data item " << fullId << " not found";
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    return TestReferenceChecker(
            new Impl(fullId, newNode, isWriteMode(),
                     impl_->bSelfTestMode_, impl_->defaultTolerance_));
}


void TestReferenceChecker::checkBoolean(bool value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cBooleanNodeName, id,
                                   ExactStringChecker(value ? "true" : "false")));
}


void TestReferenceChecker::checkString(const char *value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cStringNodeName, id,
                                   ExactStringChecker(value)));
}


void TestReferenceChecker::checkString(const std::string &value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cStringNodeName, id,
                                   ExactStringChecker(value)));
}


void TestReferenceChecker::checkStringBlock(const std::string &value,
                                            const char        *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cStringNodeName, id,
                                   ExactStringBlockChecker(value)));
}


void TestReferenceChecker::checkInteger(int value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cIntegerNodeName, id,
                                   ExactStringChecker(formatString("%d", value))));
}

void TestReferenceChecker::checkInt64(gmx_int64_t value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cInt64NodeName, id,
                                   ExactStringChecker(formatString("%" GMX_PRId64, value))));
}

void TestReferenceChecker::checkUInt64(gmx_uint64_t value, const char *id)
{
    EXPECT_TRUE(impl_->processItem(Impl::cUInt64NodeName, id,
                                   ExactStringChecker(formatString("%" GMX_PRIu64, value))));
}

void TestReferenceChecker::checkDouble(double value, const char *id)
{
    FloatingPointChecker<double> checker(value, impl_->defaultTolerance_);
    EXPECT_TRUE(impl_->processItem(Impl::cRealNodeName, id, checker));
}


void TestReferenceChecker::checkFloat(float value, const char *id)
{
    FloatingPointChecker<float> checker(value, impl_->defaultTolerance_);
    EXPECT_TRUE(impl_->processItem(Impl::cRealNodeName, id, checker));
}


void TestReferenceChecker::checkReal(float value, const char *id)
{
    checkFloat(value, id);
}


void TestReferenceChecker::checkReal(double value, const char *id)
{
    checkDouble(value, id);
}


void TestReferenceChecker::checkVector(const int value[3], const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkInteger(value[0], "X");
    compound.checkInteger(value[1], "Y");
    compound.checkInteger(value[2], "Z");
}


void TestReferenceChecker::checkVector(const float value[3], const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkReal(value[0], "X");
    compound.checkReal(value[1], "Y");
    compound.checkReal(value[2], "Z");
}


void TestReferenceChecker::checkVector(const double value[3], const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkReal(value[0], "X");
    compound.checkReal(value[1], "Y");
    compound.checkReal(value[2], "Z");
}


TestReferenceChecker
TestReferenceChecker::checkSequenceCompound(const char *id, size_t length)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    return compound;
}

} // namespace test
} // namespace gmx
