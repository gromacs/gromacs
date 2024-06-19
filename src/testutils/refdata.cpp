/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
 * Implements classes and functions from refdata.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/refdata.h"

#include <cctype>
#include <cinttypes>
#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"
#include "testutils/testfilemanager.h"

#include "refdata_checkers.h"
#include "refdata_impl.h"
#include "refdata_xml.h"

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
    TestReferenceDataImpl(ReferenceDataMode                    mode,
                          bool                                 bSelfTestMode,
                          std::optional<std::filesystem::path> testNameOverride);

    //! Performs final reference data processing when test ends.
    void onTestEnd(bool testPassed) const;

    //! Full path of the reference data file.
    std::filesystem::path fullFilename_;
    /*! \brief
     * Root entry for comparing the reference data.
     *
     * Null after construction iff in compare mode and reference data was
     * not loaded successfully.
     * In all write modes, copies are present for nodes added to
     * \a outputRootEntry_, and ReferenceDataEntry::correspondingOutputEntry()
     * points to the copy in the output tree.
     */
    ReferenceDataEntry::EntryPointer compareRootEntry_;
    /*! \brief
     * Root entry for writing new reference data.
     *
     * Null if only comparing against existing data.  Otherwise, starts
     * always as empty.
     * When creating new reference data, this is maintained as a copy of
     * \a compareRootEntry_.
     * When updating existing data, entries are added either by copying
     * from \a compareRootEntry_ (if they exist and comparison passes), or
     * by creating new ones.
     */
    ReferenceDataEntry::EntryPointer outputRootEntry_;
    /*! \brief
     * Whether updating existing reference data.
     */
    bool updateMismatchingEntries_;
    //! `true` if self-testing (enables extra failure messages).
    bool bSelfTestMode_;
    /*! \brief
     * Whether any reference checkers have been created for this data.
     */
    bool bInUse_;
};

} // namespace internal

/********************************************************************
 * Internal helpers
 */

namespace
{

//! Convenience typedef for a smart pointer to TestReferenceDataImpl.
typedef std::shared_ptr<internal::TestReferenceDataImpl> TestReferenceDataImplPointer;

/*! \brief
 * Global reference data instance.
 *
 * The object is created when the test creates a TestReferenceData, and the
 * object is destructed (and other post-processing is done) at the end of each
 * test by ReferenceDataTestEventListener (which is installed as a Google Test
 * test listener).
 */
TestReferenceDataImplPointer g_referenceData;
//! Global reference data mode set by the `-ref-data` command-line option
ReferenceDataMode g_referenceDataMode = ReferenceDataMode::Compare;

} // namespace

ReferenceDataMode referenceDataMode()
{
    return g_referenceDataMode;
}

namespace
{

//! Returns a reference to the global reference data object.
TestReferenceDataImplPointer initReferenceDataInstance(std::optional<std::filesystem::path> testNameOverride)
{
    GMX_RELEASE_ASSERT(!g_referenceData, "Test cannot create multiple TestReferenceData instances");
    g_referenceData.reset(new internal::TestReferenceDataImpl(
            referenceDataMode(), false, std::move(testNameOverride)));
    return g_referenceData;
}

//! Handles reference data creation for self-tests.
TestReferenceDataImplPointer initReferenceDataInstanceForSelfTest(ReferenceDataMode mode)
{
    if (g_referenceData)
    {
        GMX_RELEASE_ASSERT(g_referenceData.use_count() == 1,
                           "Test cannot create multiple TestReferenceData instances");
        g_referenceData->onTestEnd(true);
        g_referenceData.reset();
    }
    g_referenceData.reset(new internal::TestReferenceDataImpl(mode, true, std::nullopt));
    return g_referenceData;
}

class ReferenceDataTestEventListener : public ::testing::EmptyTestEventListener
{
public:
    void OnTestEnd(const ::testing::TestInfo& test_info) override
    {
        if (g_referenceData)
        {
            GMX_RELEASE_ASSERT(g_referenceData.use_count() == 1,
                               "Test leaked TestRefeferenceData objects");
            g_referenceData->onTestEnd(test_info.result()->Passed());
            g_referenceData.reset();
        }
    }

    void OnTestProgramEnd(const ::testing::UnitTest& /*unused*/) override
    {
        // Could be used e.g. to free internal buffers allocated by an XML parsing library
    }
};

//! Formats a path to a reference data entry with a non-null id.
std::string formatEntryPath(const std::string& prefix, const std::string& id)
{
    return prefix + "/" + id;
}

//! Formats a path to a reference data entry with a null id.
std::string formatSequenceEntryPath(const std::string& prefix, int seqIndex)
{
    return formatString("%s/[%d]", prefix.c_str(), seqIndex + 1);
}

//! Finds all entries that have not been checked under a given root.
void gatherUnusedEntries(const ReferenceDataEntry& root,
                         const std::string&        rootPath,
                         std::vector<std::string>* unusedPaths)
{
    if (!root.hasBeenChecked())
    {
        unusedPaths->push_back(rootPath);
        return;
    }
    int seqIndex = 0;
    for (const auto& child : root.children())
    {
        std::string path;
        if (child->id().empty())
        {
            path = formatSequenceEntryPath(rootPath, seqIndex);
            ++seqIndex;
        }
        else
        {
            path = formatEntryPath(rootPath, child->id());
        }
        gatherUnusedEntries(*child, path, unusedPaths);
    }
}

//! Produces a GTest assertion of any entries under given root have not been checked.
void checkUnusedEntries(const ReferenceDataEntry& root, const std::string& rootPath)
{
    std::vector<std::string> unusedPaths;
    gatherUnusedEntries(root, rootPath, &unusedPaths);
    if (!unusedPaths.empty())
    {
        std::string paths;
        if (unusedPaths.size() > 5)
        {
            paths = joinStrings(unusedPaths.begin(), unusedPaths.begin() + 5, "\n  ");
            paths = "  " + paths + "\n  ...";
        }
        else
        {
            paths = joinStrings(unusedPaths.begin(), unusedPaths.end(), "\n  ");
            paths = "  " + paths;
        }
        ADD_FAILURE() << "Reference data items not used in test:" << std::endl << paths;
    }
}

} // namespace

void initReferenceData(IOptionsContainer* options)
{
    static const gmx::EnumerationArray<ReferenceDataMode, const char*> s_refDataNames = {
        { "check", "create", "update-changed", "update-all" }
    };
    options->addOption(EnumOption<ReferenceDataMode>("ref-data")
                               .enumValue(s_refDataNames)
                               .store(&g_referenceDataMode)
                               .description("Operation mode for tests that use reference data"));
    ::testing::UnitTest::GetInstance()->listeners().Append(new ReferenceDataTestEventListener);
}

/********************************************************************
 * TestReferenceDataImpl definition
 */

namespace internal
{

TestReferenceDataImpl::TestReferenceDataImpl(ReferenceDataMode                    mode,
                                             bool                                 bSelfTestMode,
                                             std::optional<std::filesystem::path> testNameOverride) :
    updateMismatchingEntries_(false), bSelfTestMode_(bSelfTestMode), bInUse_(false)
{
    const std::filesystem::path dirname = bSelfTestMode
                                                  ? TestFileManager::getGlobalOutputTempDirectory()
                                                  : TestFileManager::getInputDataDirectory();
    const std::filesystem::path filename =
            testNameOverride.has_value() ? testNameOverride.value()
                                         : TestFileManager::getTestSpecificFileName(".xml");
    fullFilename_ = dirname / "refdata" / filename;

    switch (mode)
    {
        case ReferenceDataMode::Compare:
            if (File::exists(fullFilename_, File::throwOnError))
            {
                compareRootEntry_ = readReferenceDataFile(fullFilename_.string());
            }
            break;
        case ReferenceDataMode::CreateMissing:
            if (File::exists(fullFilename_, File::throwOnError))
            {
                compareRootEntry_ = readReferenceDataFile(fullFilename_.string());
            }
            else
            {
                compareRootEntry_ = ReferenceDataEntry::createRoot();
                outputRootEntry_  = ReferenceDataEntry::createRoot();
            }
            break;
        case ReferenceDataMode::UpdateChanged:
            if (File::exists(fullFilename_, File::throwOnError))
            {
                compareRootEntry_ = readReferenceDataFile(fullFilename_.string());
            }
            else
            {
                compareRootEntry_ = ReferenceDataEntry::createRoot();
            }
            outputRootEntry_          = ReferenceDataEntry::createRoot();
            updateMismatchingEntries_ = true;
            break;
        case ReferenceDataMode::UpdateAll:
            compareRootEntry_ = ReferenceDataEntry::createRoot();
            outputRootEntry_  = ReferenceDataEntry::createRoot();
            break;
        case ReferenceDataMode::Count: GMX_THROW(InternalError("Invalid reference data mode"));
    }
}

void TestReferenceDataImpl::onTestEnd(bool testPassed) const
{
    if (!bInUse_)
    {
        return;
    }
    // TODO: Only write the file with update-changed if there were actual changes.
    if (outputRootEntry_)
    {
        if (testPassed)
        {
            auto dirname = fullFilename_.parent_path();
            if (!std::filesystem::exists(dirname))
            {
                if (!std::filesystem::create_directory(dirname))
                {
                    GMX_THROW(TestException(gmx::formatString(
                            "Creation of reference data directory failed: %s", dirname.string().c_str())));
                }
            }
            writeReferenceDataFile(fullFilename_.string(), *outputRootEntry_);
        }
    }
    else if (compareRootEntry_)
    {
        checkUnusedEntries(*compareRootEntry_, "");
    }
}

} // namespace internal


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
    static const char* const cBooleanNodeName;
    //! String constant for naming XML elements for string values.
    static const char* const cStringNodeName;
    //! String constant for naming XML elements for single char values.
    static const char* const cCharNodeName;
    //! String constant for naming XML elements for unsigned char values.
    static const char* const cUCharNodeName;
    //! String constant for naming XML elements for integer values.
    static const char* const cIntegerNodeName;
    //! String constant for naming XML elements for int32 values.
    static const char* const cInt32NodeName;
    //! String constant for naming XML elements for unsigned int32 values.
    static const char* const cUInt32NodeName;
    //! String constant for naming XML elements for int32 values.
    static const char* const cInt64NodeName;
    //! String constant for naming XML elements for unsigned int64 values.
    static const char* const cUInt64NodeName;
    //! String constant for naming XML elements for floating-point values.
    static const char* const cRealNodeName;
    //! String constant for naming XML attribute for value identifiers.
    static const char* const cIdAttrName;
    //! String constant for naming compounds for vectors.
    static const char* const cVectorType;
    //! String constant for naming compounds for key-value tree objects.
    static const char* const cObjectType;
    //! String constant for naming compounds for sequences.
    static const char* const cSequenceType;
    //! String constant for value identifier for sequence length.
    static const char* const cSequenceLengthName;

    //! Creates a checker that does nothing.
    explicit Impl(bool initialized);
    //! Creates a checker with a given root entry.
    Impl(const std::string&            path,
         ReferenceDataEntry*           compareRootEntry,
         ReferenceDataEntry*           outputRootEntry,
         bool                          updateMismatchingEntries,
         bool                          bSelfTestMode,
         const FloatingPointTolerance& defaultTolerance);

    //! Returns the path of this checker with \p id appended.
    std::string appendPath(const char* id) const;

    //! Creates an entry with given parameters and fills it with \p checker.
    static ReferenceDataEntry::EntryPointer createEntry(const char*                       type,
                                                        const char*                       id,
                                                        const IReferenceDataEntryChecker& checker)
    {
        ReferenceDataEntry::EntryPointer entry(new ReferenceDataEntry(type, id));
        checker.fillEntry(entry.get());
        return entry;
    }
    //! Checks an entry for correct type and using \p checker.
    static ::testing::AssertionResult checkEntry(const ReferenceDataEntry&         entry,
                                                 const std::string&                fullId,
                                                 const char*                       type,
                                                 const IReferenceDataEntryChecker& checker)
    {
        if (entry.type() != type)
        {
            return ::testing::AssertionFailure() << "Mismatching reference data item type" << std::endl
                                                 << "  In item: " << fullId << std::endl
                                                 << "   Actual: " << type << std::endl
                                                 << "Reference: " << entry.type();
        }
        return checker.checkEntry(entry, fullId);
    }
    //! Finds an entry by id and updates the last found entry pointer.
    ReferenceDataEntry* findEntry(const char* id);
    /*! \brief
     * Finds/creates a reference data entry to match against.
     *
     * \param[in]  type   Type of entry to create.
     * \param[in]  id     Unique identifier of the entry (can be NULL, in
     *      which case the next entry without an id is matched).
     * \param[out] checker  Checker to use for filling out created entries.
     * \returns    Matching entry, or NULL if no matching entry found
     *      (NULL is never returned in write mode; new entries are created
     *      instead).
     */
    ReferenceDataEntry* findOrCreateEntry(const char*                       type,
                                          const char*                       id,
                                          const IReferenceDataEntryChecker& checker);
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
    ::testing::AssertionResult processItem(const char*                       name,
                                           const char*                       id,
                                           const IReferenceDataEntryChecker& checker);
    /*! \brief
     * Whether the checker is initialized.
     */
    bool initialized() const { return initialized_; }
    /*! \brief
     * Whether the checker should ignore all validation calls.
     *
     * This is used to ignore any calls within compounds for which
     * reference data could not be found, such that only one error is
     * issued for the missing compound, instead of every individual value.
     */
    bool shouldIgnore() const
    {
        GMX_RELEASE_ASSERT(initialized(), "Accessing uninitialized reference data checker.");
        return compareRootEntry_ == nullptr;
    }

    //! Whether initialized with other means than the default constructor.
    bool initialized_;
    //! Default floating-point comparison tolerance.
    FloatingPointTolerance defaultTolerance_;
    /*! \brief
     * Human-readable path to the root node of this checker.
     *
     * For the root checker, this will be "/", and for each compound, the
     * id of the compound is added.  Used for reporting comparison
     * mismatches.
     */
    std::string path_;
    /*! \brief
     * Current entry under which reference data is searched for comparison.
     *
     * Points to either the TestReferenceDataImpl::compareRootEntry_, or to
     * a compound entry in the tree rooted at that entry.
     *
     * Can be NULL, in which case this checker does nothing (doesn't even
     * report errors, see shouldIgnore()).
     */
    ReferenceDataEntry* compareRootEntry_;
    /*! \brief
     * Current entry under which entries for writing are created.
     *
     * Points to either the TestReferenceDataImpl::outputRootEntry_, or to
     * a compound entry in the tree rooted at that entry.  NULL if only
     * comparing, or if shouldIgnore() returns `false`.
     */
    ReferenceDataEntry* outputRootEntry_;
    /*! \brief
     * Iterator to a child of \a compareRootEntry_ that was last found.
     *
     * If `compareRootEntry_->isValidChild()` returns false, no entry has
     * been found yet.
     * After every check, is updated to point to the entry that was used
     * for the check.
     * Subsequent checks start the search for the matching node on this
     * node.
     */
    ReferenceDataEntry::ChildIterator lastFoundEntry_;
    /*! \brief
     * Whether the reference data is being written (true) or compared
     * (false).
     */
    bool updateMismatchingEntries_;
    //! `true` if self-testing (enables extra failure messages).
    bool bSelfTestMode_;
    /*! \brief
     * Current number of unnamed elements in a sequence.
     *
     * It is the index of the current unnamed element.
     */
    int seqIndex_;
};

const char* const TestReferenceChecker::Impl::cBooleanNodeName    = "Bool";
const char* const TestReferenceChecker::Impl::cStringNodeName     = "String";
const char* const TestReferenceChecker::Impl::cCharNodeName       = "Char";
const char* const TestReferenceChecker::Impl::cUCharNodeName      = "UChar";
const char* const TestReferenceChecker::Impl::cIntegerNodeName    = "Int";
const char* const TestReferenceChecker::Impl::cInt32NodeName      = "Int32";
const char* const TestReferenceChecker::Impl::cUInt32NodeName     = "UInt32";
const char* const TestReferenceChecker::Impl::cInt64NodeName      = "Int64";
const char* const TestReferenceChecker::Impl::cUInt64NodeName     = "UInt64";
const char* const TestReferenceChecker::Impl::cRealNodeName       = "Real";
const char* const TestReferenceChecker::Impl::cIdAttrName         = "Name";
const char* const TestReferenceChecker::Impl::cVectorType         = "Vector";
const char* const TestReferenceChecker::Impl::cObjectType         = "Object";
const char* const TestReferenceChecker::Impl::cSequenceType       = "Sequence";
const char* const TestReferenceChecker::Impl::cSequenceLengthName = "Length";


TestReferenceChecker::Impl::Impl(bool initialized) :
    initialized_(initialized),
    defaultTolerance_(defaultRealTolerance()),
    compareRootEntry_(nullptr),
    outputRootEntry_(nullptr),
    updateMismatchingEntries_(false),
    bSelfTestMode_(false),
    seqIndex_(-1)
{
}


TestReferenceChecker::Impl::Impl(const std::string&            path,
                                 ReferenceDataEntry*           compareRootEntry,
                                 ReferenceDataEntry*           outputRootEntry,
                                 bool                          updateMismatchingEntries,
                                 bool                          bSelfTestMode,
                                 const FloatingPointTolerance& defaultTolerance) :
    initialized_(true),
    defaultTolerance_(defaultTolerance),
    path_(path),
    compareRootEntry_(compareRootEntry),
    outputRootEntry_(outputRootEntry),
    lastFoundEntry_(compareRootEntry->children().end()),
    updateMismatchingEntries_(updateMismatchingEntries),
    bSelfTestMode_(bSelfTestMode),
    seqIndex_(-1)
{
}


std::string TestReferenceChecker::Impl::appendPath(const char* id) const
{
    return id != nullptr ? formatEntryPath(path_, id) : formatSequenceEntryPath(path_, seqIndex_);
}


ReferenceDataEntry* TestReferenceChecker::Impl::findEntry(const char* id)
{
    ReferenceDataEntry::ChildIterator entry = compareRootEntry_->findChild(id, lastFoundEntry_);
    seqIndex_                               = (id == nullptr) ? seqIndex_ + 1 : -1;
    if (compareRootEntry_->isValidChild(entry))
    {
        lastFoundEntry_ = entry;
        return entry->get();
    }
    return nullptr;
}

ReferenceDataEntry* TestReferenceChecker::Impl::findOrCreateEntry(const char* type,
                                                                  const char* id,
                                                                  const IReferenceDataEntryChecker& checker)
{
    ReferenceDataEntry* entry = findEntry(id);
    if (entry == nullptr && outputRootEntry_ != nullptr)
    {
        lastFoundEntry_ = compareRootEntry_->addChild(createEntry(type, id, checker));
        entry           = lastFoundEntry_->get();
    }
    return entry;
}

::testing::AssertionResult TestReferenceChecker::Impl::processItem(const char* type,
                                                                   const char* id,
                                                                   const IReferenceDataEntryChecker& checker)
{
    if (shouldIgnore())
    {
        return ::testing::AssertionSuccess();
    }
    std::string         fullId = appendPath(id);
    ReferenceDataEntry* entry  = findOrCreateEntry(type, id, checker);
    if (entry == nullptr)
    {
        return ::testing::AssertionFailure() << "Reference data item " << fullId << " not found";
    }
    entry->setChecked();
    ::testing::AssertionResult result(checkEntry(*entry, fullId, type, checker));
    if (outputRootEntry_ != nullptr && entry->correspondingOutputEntry() == nullptr)
    {
        if (!updateMismatchingEntries_ || result)
        {
            outputRootEntry_->addChild(entry->cloneToOutputEntry());
        }
        else
        {
            ReferenceDataEntry::EntryPointer outputEntry(createEntry(type, id, checker));
            entry->setCorrespondingOutputEntry(outputEntry.get());
            outputRootEntry_->addChild(std::move(outputEntry));
            return ::testing::AssertionSuccess();
        }
    }
    if (bSelfTestMode_ && !result)
    {
        ReferenceDataEntry expected(type, id);
        checker.fillEntry(&expected);
        result << std::endl
               << "String value: '" << expected.value() << "'" << std::endl
               << " Ref. string: '" << entry->value() << "'";
    }
    return result;
}


/********************************************************************
 * TestReferenceData
 */

TestReferenceData::TestReferenceData() : impl_(initReferenceDataInstance(std::nullopt)) {}


TestReferenceData::TestReferenceData(std::string testNameOverride) :
    impl_(initReferenceDataInstance(std::move(testNameOverride)))
{
}

TestReferenceData::TestReferenceData(ReferenceDataMode mode) :
    impl_(initReferenceDataInstanceForSelfTest(mode))
{
}


TestReferenceData::~TestReferenceData() {}


TestReferenceChecker TestReferenceData::rootChecker()
{
    if (!impl_->bInUse_ && !impl_->compareRootEntry_)
    {
        ADD_FAILURE() << "Reference data file not found: " << impl_->fullFilename_;
    }
    impl_->bInUse_ = true;
    if (!impl_->compareRootEntry_)
    {
        return TestReferenceChecker(new TestReferenceChecker::Impl(true));
    }
    impl_->compareRootEntry_->setChecked();
    return TestReferenceChecker(new TestReferenceChecker::Impl("",
                                                               impl_->compareRootEntry_.get(),
                                                               impl_->outputRootEntry_.get(),
                                                               impl_->updateMismatchingEntries_,
                                                               impl_->bSelfTestMode_,
                                                               defaultRealTolerance()));
}


/********************************************************************
 * TestReferenceChecker
 */

TestReferenceChecker::TestReferenceChecker() : impl_(new Impl(false)) {}

TestReferenceChecker::TestReferenceChecker(Impl* impl) : impl_(impl) {}

TestReferenceChecker::TestReferenceChecker(const TestReferenceChecker& other) :
    impl_(new Impl(*other.impl_))
{
}

TestReferenceChecker::TestReferenceChecker(TestReferenceChecker&& other) noexcept :
    impl_(std::move(other.impl_))
{
}

TestReferenceChecker& TestReferenceChecker::operator=(TestReferenceChecker&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

TestReferenceChecker::~TestReferenceChecker() {}

bool TestReferenceChecker::isValid() const
{
    return impl_->initialized();
}


void TestReferenceChecker::setDefaultTolerance(const FloatingPointTolerance& tolerance)
{
    impl_->defaultTolerance_ = tolerance;
}


void TestReferenceChecker::checkUnusedEntries()
{
    if (impl_->compareRootEntry_)
    {
        gmx::test::checkUnusedEntries(*impl_->compareRootEntry_, impl_->path_);
        // Mark them checked so that they are reported only once.
        impl_->compareRootEntry_->setCheckedIncludingChildren();
    }
}

void TestReferenceChecker::disableUnusedEntriesCheck()
{
    if (impl_->compareRootEntry_)
    {
        impl_->compareRootEntry_->setCheckedIncludingChildren();
    }
}


bool TestReferenceChecker::checkPresent(bool bPresent, const char* id)
{
    if (impl_->shouldIgnore() || impl_->outputRootEntry_ != nullptr)
    {
        return bPresent;
    }
    ReferenceDataEntry::ChildIterator entry =
            impl_->compareRootEntry_->findChild(id, impl_->lastFoundEntry_);
    const bool bFound = impl_->compareRootEntry_->isValidChild(entry);
    if (bFound != bPresent)
    {
        ADD_FAILURE() << "Mismatch while checking reference data item '" << impl_->appendPath(id) << "'\n"
                      << "Expected: " << (bPresent ? "it is present.\n" : "it is absent.\n")
                      << "  Actual: " << (bFound ? "it is present." : "it is absent.");
    }
    if (bFound && bPresent)
    {
        impl_->lastFoundEntry_ = entry;
        return true;
    }
    return false;
}


TestReferenceChecker TestReferenceChecker::checkCompound(const char* type, const char* id)
{
    if (impl_->shouldIgnore())
    {
        return TestReferenceChecker(new Impl(true));
    }
    std::string         fullId = impl_->appendPath(id);
    NullChecker         checker;
    ReferenceDataEntry* entry = impl_->findOrCreateEntry(type, id, checker);
    if (entry == nullptr)
    {
        ADD_FAILURE() << "Reference data item " << fullId << " not found";
        return TestReferenceChecker(new Impl(true));
    }
    entry->setChecked();
    if (impl_->updateMismatchingEntries_)
    {
        entry->makeCompound(type);
    }
    else
    {
        ::testing::AssertionResult result(impl_->checkEntry(*entry, fullId, type, checker));
        EXPECT_PLAIN(result);
        if (!result)
        {
            return TestReferenceChecker(new Impl(true));
        }
    }
    if (impl_->outputRootEntry_ != nullptr && entry->correspondingOutputEntry() == nullptr)
    {
        impl_->outputRootEntry_->addChild(entry->cloneToOutputEntry());
    }
    return TestReferenceChecker(new Impl(fullId,
                                         entry,
                                         entry->correspondingOutputEntry(),
                                         impl_->updateMismatchingEntries_,
                                         impl_->bSelfTestMode_,
                                         impl_->defaultTolerance_));
}

TestReferenceChecker TestReferenceChecker::checkCompound(const char* type, const std::string& id)
{
    return checkCompound(type, id.c_str());
}

/*! \brief Throw a TestException if the caller tries to write particular refdata that can't work.
 *
 * If the string to write is non-empty and has only whitespace,
 * TinyXML2 can't read it correctly, so throw an exception for this
 * case, so that we can't accidentally use it and run into mysterious
 * problems.
 *
 * \todo Eliminate this limitation of TinyXML2. See
 * e.g. https://github.com/leethomason/tinyxml2/issues/432
 */
static void throwIfNonEmptyAndOnlyWhitespace(const std::string& s, const char* id)
{
    if (!s.empty() && std::all_of(s.cbegin(), s.cend(), [](const char& c) { return std::isspace(c); }))
    {
        std::string message("String '" + s + "' with ");
        message += (id != nullptr) ? "null " : "";
        message += "ID ";
        message += (id != nullptr) ? "" : id;
        message +=
                " cannot be handled. We must refuse to write a refdata String"
                "field for a non-empty string that contains only whitespace, "
                "because it will not be read correctly by TinyXML2.";
        GMX_THROW(TestException(message));
    }
}

void TestReferenceChecker::checkBoolean(bool value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cBooleanNodeName, id, ExactStringChecker(value ? "true" : "false")));
}


void TestReferenceChecker::checkString(const char* value, const char* id)
{
    throwIfNonEmptyAndOnlyWhitespace(value, id);
    EXPECT_PLAIN(impl_->processItem(Impl::cStringNodeName, id, ExactStringChecker(value)));
}


void TestReferenceChecker::checkString(const std::string& value, const char* id)
{
    throwIfNonEmptyAndOnlyWhitespace(value, id);
    EXPECT_PLAIN(impl_->processItem(Impl::cStringNodeName, id, ExactStringChecker(value)));
}


void TestReferenceChecker::checkTextBlock(const std::string& value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(Impl::cStringNodeName, id, ExactStringBlockChecker(value)));
}


void TestReferenceChecker::checkChar(char value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cCharNodeName, id, ExactStringChecker(formatString("%c", value))));
}


void TestReferenceChecker::checkUChar(unsigned char value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cUCharNodeName, id, ExactStringChecker(formatString("%d", value))));
}

void TestReferenceChecker::checkInteger(int value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cIntegerNodeName, id, ExactStringChecker(formatString("%d", value))));
}

void TestReferenceChecker::checkInt32(int32_t value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cInt32NodeName, id, ExactStringChecker(formatString("%" PRId32, value))));
}

void TestReferenceChecker::checkUInt32(uint32_t value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cUInt32NodeName, id, ExactStringChecker(formatString("%" PRIu32, value))));
}

void TestReferenceChecker::checkInt64(int64_t value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cInt64NodeName, id, ExactStringChecker(formatString("%" PRId64, value))));
}

void TestReferenceChecker::checkUInt64(uint64_t value, const char* id)
{
    EXPECT_PLAIN(impl_->processItem(
            Impl::cUInt64NodeName, id, ExactStringChecker(formatString("%" PRIu64, value))));
}

void TestReferenceChecker::checkDouble(double value, const char* id)
{
    FloatingPointChecker<double> checker(value, impl_->defaultTolerance_);
    EXPECT_PLAIN(impl_->processItem(Impl::cRealNodeName, id, checker));
}


void TestReferenceChecker::checkFloat(float value, const char* id)
{
    FloatingPointChecker<float> checker(value, impl_->defaultTolerance_);
    EXPECT_PLAIN(impl_->processItem(Impl::cRealNodeName, id, checker));
}


void TestReferenceChecker::checkReal(float value, const char* id)
{
    checkFloat(value, id);
}


void TestReferenceChecker::checkReal(double value, const char* id)
{
    checkDouble(value, id);
}


void TestReferenceChecker::checkRealFromString(const std::string& value, const char* id)
{
    FloatingPointFromStringChecker<real> checker(value, impl_->defaultTolerance_);
    EXPECT_PLAIN(impl_->processItem(Impl::cRealNodeName, id, checker));
}


void TestReferenceChecker::checkVector(const BasicVector<int>& value, const char* id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkInteger(value[0], "X");
    compound.checkInteger(value[1], "Y");
    compound.checkInteger(value[2], "Z");
}


void TestReferenceChecker::checkVector(const BasicVector<float>& value, const char* id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkReal(value[0], "X");
    compound.checkReal(value[1], "Y");
    compound.checkReal(value[2], "Z");
}


void TestReferenceChecker::checkVector(const BasicVector<double>& value, const char* id)
{
    TestReferenceChecker compound(checkCompound(Impl::cVectorType, id));
    compound.checkReal(value[0], "X");
    compound.checkReal(value[1], "Y");
    compound.checkReal(value[2], "Z");
}


void TestReferenceChecker::checkVector(const int value[3], const char* id)
{
    checkVector(BasicVector<int>(value), id);
}


void TestReferenceChecker::checkVector(const float value[3], const char* id)
{
    checkVector(BasicVector<float>(value), id);
}


void TestReferenceChecker::checkVector(const double value[3], const char* id)
{
    checkVector(BasicVector<double>(value), id);
}


void TestReferenceChecker::checkAny(const Any& any, const char* id)
{
    if (any.isType<bool>())
    {
        checkBoolean(any.cast<bool>(), id);
    }
    else if (any.isType<char>())
    {
        checkChar(any.cast<char>(), id);
    }
    else if (any.isType<unsigned char>())
    {
        checkUChar(any.cast<unsigned char>(), id);
    }
    else if (any.isType<int>())
    {
        checkInteger(any.cast<int>(), id);
    }
    else if (any.isType<int32_t>())
    {
        checkInt32(any.cast<int32_t>(), id);
    }
    else if (any.isType<uint32_t>())
    {
        checkInt32(any.cast<uint32_t>(), id);
    }
    else if (any.isType<int64_t>())
    {
        checkInt64(any.cast<int64_t>(), id);
    }
    else if (any.isType<uint64_t>())
    {
        checkInt64(any.cast<uint64_t>(), id);
    }
    else if (any.isType<float>())
    {
        checkFloat(any.cast<float>(), id);
    }
    else if (any.isType<double>())
    {
        checkDouble(any.cast<double>(), id);
    }
    else if (any.isType<std::string>())
    {
        checkString(any.cast<std::string>(), id);
    }
    else
    {
        GMX_THROW(TestException("Unsupported any type"));
    }
}


void TestReferenceChecker::checkKeyValueTreeObject(const KeyValueTreeObject& tree, const char* id)
{
    TestReferenceChecker compound(checkCompound(Impl::cObjectType, id));
    for (const auto& prop : tree.properties())
    {
        compound.checkKeyValueTreeValue(prop.value(), prop.key().c_str());
    }
    compound.checkUnusedEntries();
}


void TestReferenceChecker::checkKeyValueTreeValue(const KeyValueTreeValue& value, const char* id)
{
    if (value.isObject())
    {
        checkKeyValueTreeObject(value.asObject(), id);
    }
    else if (value.isArray())
    {
        const auto& values = value.asArray().values();
        checkSequence(values.begin(), values.end(), id);
    }
    else
    {
        checkAny(value.asAny(), id);
    }
}


TestReferenceChecker TestReferenceChecker::checkSequenceCompound(const char* id, size_t length)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    return compound;
}


unsigned char TestReferenceChecker::readUChar(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    int value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cUCharNodeName, id, ValueExtractor<int>(&value)));
    return value;
}


int TestReferenceChecker::readInteger(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    int value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cIntegerNodeName, id, ValueExtractor<int>(&value)));
    return value;
}


int32_t TestReferenceChecker::readInt32(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    int32_t value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cInt32NodeName, id, ValueExtractor<int32_t>(&value)));
    return value;
}


int64_t TestReferenceChecker::readInt64(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    int64_t value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cInt64NodeName, id, ValueExtractor<int64_t>(&value)));
    return value;
}


float TestReferenceChecker::readFloat(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    float value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cRealNodeName, id, ValueExtractor<float>(&value)));
    return value;
}


double TestReferenceChecker::readDouble(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    double value = 0;
    EXPECT_PLAIN(impl_->processItem(Impl::cRealNodeName, id, ValueExtractor<double>(&value)));
    return value;
}


std::string TestReferenceChecker::readString(const char* id)
{
    if (impl_->shouldIgnore())
    {
        GMX_THROW(TestException("Trying to read from non-existent reference data value"));
    }
    std::string value;
    EXPECT_PLAIN(impl_->processItem(Impl::cStringNodeName, id, ValueExtractor<std::string>(&value)));
    return value;
}

} // namespace test
} // namespace gmx
