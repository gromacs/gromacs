/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
#include "refdata.h"

#include <cstdio>
#include <cstdlib>

#include <limits>
#include <string>
#include <sstream>

#include <gtest/gtest.h>
#include <rapidxml.hpp>
#include <rapidxml_print.hpp>

#include "gromacs/fileio/path.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/uniqueptr.h"

#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"
#include "testutils/testfilemanager.h"

namespace
{

//@{
/*! \brief Helper typedef
 *
 * A previous version of this implementation used libxml2 and we might
 * one day use it again. This implementation uses types that are named
 * similarly to libxml2, but are based on RapidXml to do the actual
 * XML parsing. */
typedef char xmlChar;
typedef rapidxml::xml_document<xmlChar> xmlDoc;
typedef gmx::gmx_unique_ptr<xmlDoc>::type xmlDocPtr;
typedef rapidxml::xml_node<xmlChar> xmlNode;
typedef xmlNode *xmlNodePtr;
typedef rapidxml::xml_attribute<xmlChar> xmlAttr;
typedef xmlAttr *xmlAttrPtr;
//@}

//! Global reference data mode set with gmx::test::setReferenceDataMode().
// TODO: Make this a real enum; requires solving a TODO in StringOption.
int g_referenceDataMode = gmx::test::erefdataCompare;

} // namespace

namespace gmx
{
namespace test
{

ReferenceDataMode getReferenceDataMode()
{
    return static_cast<ReferenceDataMode>(g_referenceDataMode);
}

void setReferenceDataMode(ReferenceDataMode mode)
{
    g_referenceDataMode = mode;
}

std::string getReferenceDataPath()
{
    return TestFileManager::getInputFilePath("refdata");
}

void initReferenceData(Options *options)
{
    // Needs to correspond to the enum order in refdata.h.
    const char *const refDataEnum[] = { "check", "create", "update" };
    options->addOption(
            StringOption("ref-data").enumValue(refDataEnum)
                .defaultEnumIndex(0)
                .storeEnumIndex(&g_referenceDataMode)
                .description("Operation mode for tests that use reference data"));
}

/********************************************************************
 * TestReferenceData::Impl
 */

/*! \internal \brief
 * Private implementation class for TestReferenceData.
 *
 * \ingroup module_testutils
 */
class TestReferenceData::Impl
{
    public:
        //! String constant for output XML version string.
        static const xmlChar * const cXmlVersion;
        //! String constant for XML stylesheet processing instruction name.
        static const xmlChar * const cXmlStyleSheetNodeName;
        //! String constant for XML stylesheet reference.
        static const xmlChar * const cXmlStyleSheetContent;
        //! String constant for naming the root XML element.
        static const xmlChar * const cRootNodeName;

        //! Initializes a checker in the given mode.
        Impl(ReferenceDataMode mode, bool bSelfTestMode);
        ~Impl();

        //! Full path of the reference data file.
        std::string             fullFilename_;
        /*! \brief
         * XML document for the reference data.
         *
         * May be NULL if there was an I/O error in initialization.
         */
        xmlDocPtr               refDoc_;
        /*! \brief String used to hold the text read from the
         * reference data file. Since the XML parsing modifies the
         * string, and the lifetime of the string must match that of
         * the document parsed from it, this needs to be a data
         * member. */
        std::string             fileText_;
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

const xmlChar * const TestReferenceData::Impl::cXmlVersion =
    (const xmlChar *)"1.0";
const xmlChar * const TestReferenceData::Impl::cXmlStyleSheetNodeName =
    (const xmlChar *)"xml-stylesheet";
const xmlChar * const TestReferenceData::Impl::cXmlStyleSheetContent =
    (const xmlChar *)"type=\"text/xsl\" href=\"referencedata.xsl\"";
const xmlChar * const TestReferenceData::Impl::cRootNodeName =
    (const xmlChar *)"ReferenceData";


TestReferenceData::Impl::Impl(ReferenceDataMode mode, bool bSelfTestMode)
    : refDoc_(new xmlDoc), fileText_(), bWrite_(false), bSelfTestMode_(bSelfTestMode), bInUse_(false)
{
    std::string dirname  = getReferenceDataPath();
    std::string filename = TestFileManager::getTestSpecificFileName(".xml");
    fullFilename_ = Path::join(dirname, filename);

    bWrite_ = true;
    if (mode != erefdataUpdateAll)
    {
        FILE *fp = std::fopen(fullFilename_.c_str(), "r");
        if (fp != NULL)
        {
            bWrite_ = false;
            fclose(fp);
        }
        else if (mode == erefdataCompare)
        {
            bWrite_ = false;
            return;
        }
    }
    if (bWrite_)
    {
        try
        {
            xmlNodePtr rootNode = refDoc_->allocate_node(rapidxml::node_element, cRootNodeName);
            refDoc_->append_node(rootNode);
            // TODO Somehow, writing the style-sheet information leads
            // RapidXml to struggle to find the cRootNodeName element
            // later on, but supporting the XSL is not all that
            // important either.
            // xmlNodePtr xslNode = refDoc_->allocate_node(rapidxml::node_pi, cXmlStyleSheetNodeName, cXmlStyleSheetContent);
            // refDoc_->prepend_node(xslNode);
        }
        catch (const std::bad_alloc &)
        {
            GMX_THROW(TestException("XML node creation failed"));
        }
    }
    else
    {
        try
        {
            fileText_ = gmx::File::readToString(fullFilename_);
        }
        catch (const FileIOError &)
        {
            ADD_FAILURE() << "Reference data file could not be read: " + fullFilename_;
        }
        if (fileText_.empty())
        {
            ADD_FAILURE() << "Reference data file produced no content: " + fullFilename_;
        }
        try
        {
            // Parsing needs to modify the string in place, but
            // apparently does so without needing to reallocate, so
            // casting away the const-ness of the pointer returned by
            // std::string::c_str works OK.
            refDoc_->parse<rapidxml::parse_default>(const_cast<xmlChar *>(fileText_.c_str()));
        }
        catch (const rapidxml::parse_error &parseError)
        {
            GMX_THROW(TestException("Reference data not parsed successfully: " + fullFilename_ + " because " + parseError.what()));
        }
        xmlNodePtr rootNode = refDoc_->first_node();
        if (rootNode == NULL)
        {
            GMX_THROW(TestException("Reference data is empty: " + fullFilename_));
        }
        if (0 != std::strcmp(rootNode->name(), cRootNodeName))
        {
            GMX_THROW(TestException("Invalid root node type in " + fullFilename_));
        }
    }
}


TestReferenceData::Impl::~Impl()
{
    if (bWrite_ && bInUse_ && refDoc_ != NULL)
    {
        std::string dirname = getReferenceDataPath();
        if (!Directory::exists(dirname))
        {
            if (Directory::create(dirname) != 0)
            {
                ADD_FAILURE() << "Creation of reference data directory failed for " << dirname;
            }
        }
        std::ostringstream fileText;
        xmlNodePtr         rootNode = refDoc_->first_node();
        rapidxml::print<xmlChar>(fileText, *rootNode);
        try
        {
            gmx::File::writeFileFromString(fullFilename_, fileText.str());
        }
        catch (const gmx::FileIOError &)
        {
            ADD_FAILURE() << "Saving reference data failed for " + fullFilename_;
        }
    }
}


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
        static const xmlChar * const cBooleanNodeName;
        //! String constant for naming XML elements for string values.
        static const xmlChar * const cStringNodeName;
        //! String constant for naming XML elements for integer values.
        static const xmlChar * const cIntegerNodeName;
        //! String constant for naming XML elements for int64 values.
        static const xmlChar * const cInt64NodeName;
        //! String constant for naming XML elements for unsigned int64 values.
        static const xmlChar * const cUInt64NodeName;
        //! String constant for naming XML elements for floating-point values.
        static const xmlChar * const cRealNodeName;
        //! String constant for naming XML attribute for value identifiers.
        static const xmlChar * const cIdAttrName;
        //! String constant for naming compounds for vectors.
        static const char * const    cVectorType;
        //! String constant for naming compounds for sequences.
        static const char * const    cSequenceType;
        //! String constant for value identifier for sequence length.
        static const char * const    cSequenceLengthName;

        //! Creates a checker that does nothing.
        explicit Impl(bool bWrite);
        //! Creates a checker with a given root node.
        Impl(const std::string &path, xmlNodePtr rootNode, bool bWrite,
             bool bSelfTestMode, const FloatingPointTolerance &defaultTolerance);

        //! Returns a string for SCOPED_TRACE() for checking element \p id.
        std::string traceString(const char *id) const;
        //! Returns the path of this checker with \p id appended.
        std::string appendPath(const char *id) const;

        /*! \brief
         * Finds a reference data node.
         *
         * \param[in]  name   Type of node to find (can be NULL, in which case
         *      any type is matched).
         * \param[in]  id     Unique identifier of the node (can be NULL, in
         *      which case the next node without an id is matched).
         * \returns    Matching node, or NULL if no matching node found.
         *
         * Searches for a node in the reference data that matches the given
         * \p name and \p id.  Searching starts from the node that follows the
         * previously matched node (relevant for performance, and if there are
         * duplicate ids or nodes without ids).  Note that the match pointer is
         * not updated by this method.
         */
        xmlNodePtr findNode(const xmlChar *name, const char *id) const;
        /*! \brief
         * Finds/creates a reference data node to match against.
         *
         * \param[in]  name   Type of node to find.
         * \param[in]  id     Unique identifier of the node (can be NULL, in
         *      which case the next node without an id is matched).
         * \param[out] bFound Whether the node was found (false if the node was
         *      created in write mode).
         * \returns Matching node, or NULL if no matching node found
         *      (NULL is never returned in write mode).
         * \throws  TestException if node creation fails in write mode.
         *
         * Finds a node using findNode() and updates the match pointer is a
         * match is found.  If a match is not found, the method returns NULL in
         * read mode and creates a new node in write mode.  If the creation
         * fails in write mode, throws.
         */
        xmlNodePtr findOrCreateNode(const xmlChar *name, const char *id,
                                    bool *bFound);
        /*! \brief
         * Helper method for checking a reference data value.
         *
         * \param[in]  name   Type of node to find.
         * \param[in]  id     Unique identifier of the node (can be NULL, in
         *      which case the next node without an id is matched).
         * \param[in]  value  String value of the value to be compared.
         * \param[out] bFound true if a matching value was found.
         * \returns String value for the reference value.
         * \throws  TestException if node creation fails in write mode.
         *
         * Performs common tasks in checking a reference value:
         * finding/creating the correct XML node and reading/writing its string
         * value.  Caller is responsible for converting the value to and from
         * string where necessary and performing the actual comparison.
         *
         * In read mode, if a value is not found, adds a Google Test failure
         * and returns an empty string.  If the reference value is found,
         * returns it (\p value is not used in this case).
         *
         * In write mode, creates the node if it is not found, sets its value
         * as \p value and returns \p value.
         */
        std::string processItem(const xmlChar *name, const char *id,
                                const char *value, bool *bFound);
        //! Convenience wrapper that takes a std::string.
        std::string processItem(const xmlChar *name, const char *id,
                                const std::string &value, bool *bFound);
        /*! \brief
         * Whether the checker should ignore all validation calls.
         *
         * This is used to ignore any calls within compounds for which
         * reference data could not be found, such that only one error is
         * issued for the missing compound, instead of every individual value.
         */
        bool shouldIgnore() const;

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
         * Points to either the root of TestReferenceData::Impl::refDoc_, or to
         * a compound node.
         *
         * Can be NULL, in which case this checker does nothing (doesn't even
         * report errors, see shouldIgnore()).
         */
        xmlNodePtr              currNode_;
        /*! \brief
         * Points to a child of \a currNode_ that was last found.
         *
         * On initialization, is initialized to NULL.  After every check, is
         * updated to point to the node that was used for the check.
         * Subsequent checks start the search for the matching node on this
         * node.
         *
         * Is NULL if \a currNode_ contains no children or if no checks have
         * yet been made.
         * Otherwise, always points to a direct child of \a currNode_.
         */
        xmlNodePtr              prevFoundNode_;
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

const xmlChar * const TestReferenceChecker::Impl::cBooleanNodeName =
    (const xmlChar *)"Bool";
const xmlChar * const TestReferenceChecker::Impl::cStringNodeName =
    (const xmlChar *)"String";
const xmlChar * const TestReferenceChecker::Impl::cIntegerNodeName  =
    (const xmlChar *)"Int";
const xmlChar * const TestReferenceChecker::Impl::cInt64NodeName  =
    (const xmlChar *)"Int64";
const xmlChar * const TestReferenceChecker::Impl::cUInt64NodeName  =
    (const xmlChar *)"UInt64";
const xmlChar * const TestReferenceChecker::Impl::cRealNodeName =
    (const xmlChar *)"Real";
const xmlChar * const TestReferenceChecker::Impl::cIdAttrName =
    (const xmlChar *)"Name";
const char * const    TestReferenceChecker::Impl::cVectorType =
    "Vector";
const char * const    TestReferenceChecker::Impl::cSequenceType =
    "Sequence";
const char * const    TestReferenceChecker::Impl::cSequenceLengthName =
    "Length";


TestReferenceChecker::Impl::Impl(bool bWrite)
    : defaultTolerance_(defaultRealTolerance()),
      currNode_(NULL), prevFoundNode_(NULL), bWrite_(bWrite),
      bSelfTestMode_(false), seqIndex_(0)
{
}


TestReferenceChecker::Impl::Impl(const std::string &path, xmlNodePtr rootNode,
                                 bool bWrite, bool bSelfTestMode,
                                 const FloatingPointTolerance &defaultTolerance)
    : defaultTolerance_(defaultTolerance), path_(path + "/"),
      currNode_(rootNode), prevFoundNode_(NULL), bWrite_(bWrite),
      bSelfTestMode_(bSelfTestMode), seqIndex_(0)
{
}


std::string
TestReferenceChecker::Impl::traceString(const char *id) const
{
    return "Checking '" + appendPath(id) + "'";
}


std::string
TestReferenceChecker::Impl::appendPath(const char *id) const
{
    std::string printId = (id != NULL) ? id : formatString("[%d]", seqIndex_);
    return path_ + printId;
}


xmlNodePtr
TestReferenceChecker::Impl::findNode(const xmlChar *name, const char *id) const
{
    if (currNode_ == NULL || currNode_->first_node() == NULL)
    {
        return NULL;
    }
    const xmlChar *xmlId = reinterpret_cast<const xmlChar *>(id);
    xmlNodePtr     node  = prevFoundNode_;
    bool           bWrap = true;
    if (node != NULL)
    {
        if (id == NULL)
        {
            xmlAttrPtr refId = node->first_attribute(cIdAttrName);
            if (refId == NULL)
            {
                if (name == NULL || 0 == std::strcmp(node->name(), name))
                {
                    bWrap = false;
                    node  = node->next_sibling();
                    if (node == NULL)
                    {
                        return NULL;
                    }
                }
            }
        }
    }
    else
    {
        node  = currNode_->first_node();
        bWrap = false;
    }
    do
    {
        if (name == NULL || 0 == std::strcmp(node->name(), name))
        {
            xmlAttrPtr refId = node->first_attribute(cIdAttrName);
            if (xmlId == NULL && refId == NULL)
            {
                return node;
            }
            if (refId != NULL)
            {
                if (xmlId != NULL && 0 == std::strcmp(refId->value(), xmlId))
                {
                    return node;
                }
            }
        }
        node = node->next_sibling();
        if (bWrap && node == NULL)
        {
            node = currNode_->first_node();
        }
    }
    while (node != NULL && node != prevFoundNode_);
    return NULL;
}


xmlNodePtr
TestReferenceChecker::Impl::findOrCreateNode(const xmlChar *name,
                                             const char    *id,
                                             bool          *bFound)
{
    *bFound = false;
    xmlNodePtr node = findNode(name, id);
    if (node != NULL)
    {
        *bFound        = true;
        prevFoundNode_ = node;
    }
    else
    {
        if (bWrite_)
        {
            try
            {
                node = currNode_->document()->allocate_node(rapidxml::node_element, name);
                currNode_->append_node(node);
                if (id != NULL)
                {
                    const xmlChar *xmlId = reinterpret_cast<const xmlChar *>(id);
                    xmlAttrPtr     prop  = currNode_->document()->allocate_attribute(cIdAttrName, xmlId);
                    node->append_attribute(prop);
                }
            }
            catch (const std::bad_alloc &)
            {
                GMX_THROW(TestException("XML node creation failed"));
            }
            prevFoundNode_ = node;
        }
        else
        {
            ADD_FAILURE() << "Reference data item not found";
        }
    }
    seqIndex_ = (id == NULL) ? seqIndex_+1 : 0;

    return node;
}


std::string
TestReferenceChecker::Impl::processItem(const xmlChar *name, const char *id,
                                        const char *value, bool *bFound)
{
    xmlNodePtr node = findOrCreateNode(name, id, bFound);
    if (node == NULL)
    {
        return std::string();
    }
    if (bWrite_ && !*bFound)
    {
        node->value(node->document()->allocate_string(value));
        *bFound = true;
        return std::string(value);
    }
    else
    {
        xmlChar    *refXmlValue = node->value();
        std::string refValue(reinterpret_cast<const char *>(refXmlValue));
        return refValue;
    }
}


std::string
TestReferenceChecker::Impl::processItem(const xmlChar *name, const char *id,
                                        const std::string &value, bool *bFound)
{
    return processItem(name, id, value.c_str(), bFound);
}


bool
TestReferenceChecker::Impl::shouldIgnore() const
{
    return currNode_ == NULL;
}


/********************************************************************
 * TestReferenceData
 */

TestReferenceData::TestReferenceData()
    : impl_(new Impl(getReferenceDataMode(), false))
{
}


TestReferenceData::TestReferenceData(ReferenceDataMode mode)
    : impl_(new Impl(mode, true))
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
    if (!isWriteMode() && !impl_->bInUse_ && impl_->refDoc_ == NULL)
    {
        ADD_FAILURE() << "Reference data file not found: "
        << impl_->fullFilename_;
    }
    impl_->bInUse_ = true;
    if (impl_->refDoc_ == NULL)
    {
        return TestReferenceChecker(new TestReferenceChecker::Impl(isWriteMode()));
    }
    xmlNodePtr rootNode = impl_->refDoc_->first_node();
    // TODO: The default tolerance for double-precision builds that explicitly
    // call checkFloat() may not be ideal.
    return TestReferenceChecker(
            new TestReferenceChecker::Impl("", rootNode, isWriteMode(),
                                           impl_->bSelfTestMode_,
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
    xmlNodePtr node   = impl_->findNode(NULL, id);
    bool       bFound = (node != NULL);
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
    SCOPED_TRACE(impl_->traceString(id));
    if (impl_->shouldIgnore())
    {
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    const xmlChar *xmlNodeName = reinterpret_cast<const xmlChar *>(type);
    bool           bFound;
    xmlNodePtr     newNode     = impl_->findOrCreateNode(xmlNodeName, id, &bFound);
    if (newNode == NULL)
    {
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    return TestReferenceChecker(
            new Impl(impl_->appendPath(id), newNode, isWriteMode(),
                     impl_->bSelfTestMode_, impl_->defaultTolerance_));
}


void TestReferenceChecker::checkBoolean(bool value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    const char *strValue    = value ? "true" : "false";
    std::string refStrValue =
        impl_->processItem(Impl::cBooleanNodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}


void TestReferenceChecker::checkString(const char *value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    std::string refStrValue =
        impl_->processItem(Impl::cStringNodeName, id, value, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, value);
    }
}


void TestReferenceChecker::checkString(const std::string &value, const char *id)
{
    checkString(value.c_str(), id);
}


void TestReferenceChecker::checkStringBlock(const std::string &value,
                                            const char        *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool       bFound;
    xmlNodePtr node = impl_->findOrCreateNode(Impl::cStringNodeName, id, &bFound);
    if (node == NULL)
    {
        return;
    }
    // An extra newline is written in the beginning to make lines align
    // in the output xml (otherwise, the first line would be off by the length
    // of the starting CDATA tag).
    if (isWriteMode() && !bFound)
    {
        std::string    adjustedValue = "\n" + value;
        const xmlChar *xmlValue
            = node->document()->allocate_string(reinterpret_cast<const xmlChar *>(adjustedValue.c_str()));

        // TODO: Figure out if \r and \r\n can be handled without them changing
        // to \n in the roundtrip
        xmlNodePtr cdata
            = node->document()->allocate_node(rapidxml::node_cdata, NULL, xmlValue);
        node->append_node(cdata);
    }
    else
    {
        xmlNodePtr cdata = node->first_node();
        while (cdata != NULL && rapidxml::node_cdata != cdata->type())
        {
            cdata = cdata->next_sibling();
        }
        if (cdata == NULL)
        {
            ADD_FAILURE() << "Invalid string block element";
            return;
        }
        xmlChar *refXmlValue = cdata->value();
        if (refXmlValue[0] != '\n')
        {
            ADD_FAILURE() << "Invalid string block element";
            return;
        }
        std::string refValue(reinterpret_cast<const char *>(refXmlValue + 1));
        EXPECT_EQ(refValue, value);
    }
}


void TestReferenceChecker::checkInteger(int value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    std::string strValue    = formatString("%d", value);
    std::string refStrValue =
        impl_->processItem(Impl::cIntegerNodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}

void TestReferenceChecker::checkInt64(gmx_int64_t value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    std::string strValue    = formatString("%" GMX_PRId64, value);
    std::string refStrValue =
        impl_->processItem(Impl::cInt64NodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}

void TestReferenceChecker::checkUInt64(gmx_uint64_t value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    std::string strValue    = formatString("%" GMX_PRIu64, value);
    std::string refStrValue =
        impl_->processItem(Impl::cUInt64NodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}

void TestReferenceChecker::checkDouble(double value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    const int   prec        = std::numeric_limits<double>::digits10 + 2;
    std::string strValue    = formatString("%.*g", prec, value);
    std::string refStrValue =
        impl_->processItem(Impl::cRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        char  *endptr;
        double refValue = std::strtod(refStrValue.c_str(), &endptr);
        EXPECT_EQ('\0', *endptr);
        if (impl_->bSelfTestMode_)
        {
            EXPECT_DOUBLE_EQ_TOL(refValue, value, impl_->defaultTolerance_)
            << "String value: " << strValue << std::endl
            << " Ref. string: " << refStrValue;
        }
        else
        {
            EXPECT_DOUBLE_EQ_TOL(refValue, value, impl_->defaultTolerance_);
        }
    }
}


void TestReferenceChecker::checkFloat(float value, const char *id)
{
    if (impl_->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(impl_->traceString(id));
    bool        bFound      = false;
    const int   prec        = std::numeric_limits<float>::digits10 + 2;
    std::string strValue    = formatString("%.*g", prec, value);
    std::string refStrValue =
        impl_->processItem(Impl::cRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        char  *endptr;
        float  refValue = static_cast<float>(std::strtod(refStrValue.c_str(), &endptr));
        EXPECT_EQ('\0', *endptr);
        if (impl_->bSelfTestMode_)
        {
            EXPECT_FLOAT_EQ_TOL(refValue, value, impl_->defaultTolerance_)
            << "String value: " << strValue << std::endl
            << " Ref. string: " << refStrValue;
        }
        else
        {
            EXPECT_FLOAT_EQ_TOL(refValue, value, impl_->defaultTolerance_);
        }
    }
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
