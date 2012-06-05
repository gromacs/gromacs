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
 * Implements classes and functions from refdata.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "refdata.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <new>
#include <string>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/format.h"
#include "testutils/datapath.h"
#include "testutils/testexceptions.h"

#include "refdata-impl.h"

namespace
{

class TestReferenceDataEnvironment : public ::testing::Environment
{
    public:
        virtual void TearDown()
        {
            xmlCleanupParser();
        }
};

} // namespace

namespace gmx
{
namespace test
{

static ReferenceDataMode g_referenceDataMode = erefdataCompare;

ReferenceDataMode getReferenceDataMode()
{
    return g_referenceDataMode;
}

void setReferenceDataMode(ReferenceDataMode mode)
{
    g_referenceDataMode = mode;
}

std::string getReferenceDataPath()
{
    return TestFileManager::getTestFilePath("refdata");
}

void initReferenceData(int *argc, char **argv)
{
    int i, newi;

    for (i = newi = 1; i < *argc; ++i, ++newi)
    {
        argv[newi] = argv[i];
        if (!std::strcmp(argv[i], "--create-ref-data"))
        {
            setReferenceDataMode(erefdataCreateMissing);
            --newi;
        }
        else if (!std::strcmp(argv[i], "--update-ref-data"))
        {
            setReferenceDataMode(erefdataUpdateAll);
            --newi;
        }
    }
    *argc = newi;
    try
    {
        ::testing::AddGlobalTestEnvironment(new TestReferenceDataEnvironment);
    }
    catch (const std::bad_alloc &)
    {
        std::fprintf(stderr, "Out of memory\n");
        std::exit(1);
    }
}

/********************************************************************
 * TestReferenceData::Impl
 */

const xmlChar * const TestReferenceData::Impl::cXmlVersion =
    (const xmlChar *)"1.0";
const xmlChar * const TestReferenceData::Impl::cXmlStyleSheetNodeName =
    (const xmlChar *)"xml-stylesheet";
const xmlChar * const TestReferenceData::Impl::cXmlStyleSheetContent =
    (const xmlChar *)"type=\"text/xsl\" href=\"referencedata.xsl\"";
const xmlChar * const TestReferenceData::Impl::cRootNodeName =
    (const xmlChar *)"ReferenceData";
const xmlChar * const TestReferenceChecker::Impl::cBooleanNodeName =
    (const xmlChar *)"Bool";
const xmlChar * const TestReferenceChecker::Impl::cStringNodeName =
    (const xmlChar *)"String";
const xmlChar * const TestReferenceChecker::Impl::cIntegerNodeName  =
    (const xmlChar *)"Int";
const xmlChar * const TestReferenceChecker::Impl::cRealNodeName =
    (const xmlChar *)"Real";
const xmlChar * const TestReferenceChecker::Impl::cIdAttrName =
    (const xmlChar *)"Name";
const char * const TestReferenceChecker::Impl::cVectorType =
    "Vector";
const char * const TestReferenceChecker::Impl::cSequenceType =
    "Sequence";
const char * const TestReferenceChecker::Impl::cSequenceLengthName =
    "Length";


TestReferenceData::Impl::Impl(ReferenceDataMode mode)
    : _refDoc(NULL), _bWrite(false), _bInUse(false)
{
    std::string dirname = getReferenceDataPath();
    std::string filename = TestFileManager::getTestSpecificFileName(".xml");
    _fullFilename = Path::join(dirname, filename);

    _bWrite = true;
    if (mode != erefdataUpdateAll)
    {
        FILE *fp = std::fopen(_fullFilename.c_str(), "r");
        if (fp != NULL)
        {
            _bWrite = false;
            fclose(fp);
        }
        else if (mode == erefdataCompare)
        {
            _bWrite = false;
            return;
        }
    }
    if (_bWrite)
    {
        // TODO: Error checking
        _refDoc = xmlNewDoc(cXmlVersion);
        xmlNodePtr rootNode = xmlNewDocNode(_refDoc, NULL, cRootNodeName, NULL);
        xmlDocSetRootElement(_refDoc, rootNode);
        xmlNodePtr xslNode = xmlNewDocPI(_refDoc, cXmlStyleSheetNodeName,
                                         cXmlStyleSheetContent);
        xmlAddPrevSibling(rootNode, xslNode);
    }
    else
    {
        _refDoc = xmlParseFile(_fullFilename.c_str());
        if (_refDoc == NULL)
        {
            GMX_THROW(TestException("Reference data not parsed successfully: " + _fullFilename));
        }
        xmlNodePtr rootNode = xmlDocGetRootElement(_refDoc);
        if (rootNode == NULL)
        {
            xmlFreeDoc(_refDoc);
            GMX_THROW(TestException("Reference data is empty: " + _fullFilename));
        }
        if (xmlStrcmp(rootNode->name, cRootNodeName) != 0)
        {
            xmlFreeDoc(_refDoc);
            GMX_THROW(TestException("Invalid root node type in " + _fullFilename));
        }
    }
}


TestReferenceData::Impl::~Impl()
{
    if (_bWrite && _bInUse && _refDoc != NULL)
    {
        std::string dirname = getReferenceDataPath();
        if (!Directory::exists(dirname))
        {
            if (Directory::create(dirname) != 0)
            {
                ADD_FAILURE() << "Creation of reference data directory failed for " << dirname;
            }
        }
        if (xmlSaveFormatFile(_fullFilename.c_str(), _refDoc, 1) == -1)
        {
            ADD_FAILURE() << "Saving reference data failed for " + _fullFilename;
        }
    }
    if (_refDoc != NULL)
    {
        xmlFreeDoc(_refDoc);
    }
}


/********************************************************************
 * TestReferenceChecker::Impl
 */

TestReferenceChecker::Impl::Impl(bool bWrite)
    : _currNode(NULL), _nextSearchNode(NULL), _bWrite(bWrite), _seqIndex(0)
{
}


TestReferenceChecker::Impl::Impl(const std::string &path, xmlNodePtr rootNode,
                                 bool bWrite)
    : _path(path + "/"), _currNode(rootNode),
      _nextSearchNode(rootNode->xmlChildrenNode),
      _bWrite(bWrite), _seqIndex(0)
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
    std::string printId = (id != NULL) ? id : formatString("[%d]", _seqIndex);
    return _path + printId;
}


xmlNodePtr
TestReferenceChecker::Impl::findNode(const xmlChar *name, const char *id) const
{
    const xmlChar *xmlId = reinterpret_cast<const xmlChar *>(id);
    xmlNodePtr node = _nextSearchNode;
    if (node == NULL)
    {
        return NULL;
    }
    do
    {
        if (name == NULL || xmlStrcmp(node->name, name) == 0)
        {
            xmlChar *refId = xmlGetProp(node, cIdAttrName);
            if (xmlId == NULL && refId == NULL)
            {
                return node;
            }
            if (refId != NULL)
            {
                if (xmlId != NULL && xmlStrcmp(refId, xmlId) == 0)
                {
                    xmlFree(refId);
                    return node;
                }
                xmlFree(refId);
            }
        }
        node = node->next;
        if (node == NULL && _nextSearchNode != _currNode->xmlChildrenNode)
        {
            node = _currNode->xmlChildrenNode;
        }
    }
    while (node != NULL && node != _nextSearchNode);
    return NULL;
}


xmlNodePtr
TestReferenceChecker::Impl::findOrCreateNode(const xmlChar *name, const char *id)
{
    xmlNodePtr node = findNode(name, id);
    if (node == NULL)
    {
        if (_bWrite)
        {
            node = xmlNewTextChild(_currNode, NULL, name, NULL);
            if (node != NULL && id != NULL)
            {
                const xmlChar *xmlId = reinterpret_cast<const xmlChar *>(id);
                xmlAttrPtr prop = xmlNewProp(node, cIdAttrName, xmlId);
                if (prop == NULL)
                {
                    xmlFreeNode(node);
                    node = NULL;
                }
            }
            if (node == NULL)
            {
                GMX_THROW(TestException("XML node creation failed"));
            }
        }
        else
        {
            node = NULL;
        }
    }
    else
    {
        _nextSearchNode = node->next;
        if (_nextSearchNode == NULL)
        {
            _nextSearchNode = _currNode->xmlChildrenNode;
        }

    }
    if (node == NULL)
    {
        GMX_RELEASE_ASSERT(!_bWrite, "Node creation failed without exception");
        ADD_FAILURE() << "Reference data item not found";
    }
    _seqIndex = (id == NULL) ? _seqIndex+1 : 0;

    return node;
}


std::string
TestReferenceChecker::Impl::processItem(const xmlChar *name, const char *id,
                                        const char *value, bool *bFound)
{
    *bFound = false;
    xmlNodePtr node = findOrCreateNode(name, id);
    if (node == NULL)
    {
        return std::string();
    }
    *bFound = true;
    if (_bWrite)
    {
        xmlNodeAddContent(node, reinterpret_cast<const xmlChar *>(value));
        return std::string(value);
    }
    else
    {
        xmlChar *refXmlValue = xmlNodeGetContent(node);
        std::string refValue(reinterpret_cast<const char *>(refXmlValue));
        xmlFree(refXmlValue);
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
    return _currNode == NULL;
}


/********************************************************************
 * TestReferenceData
 */

TestReferenceData::TestReferenceData()
    : _impl(new Impl(getReferenceDataMode()))
{
}


TestReferenceData::TestReferenceData(ReferenceDataMode mode)
    : _impl(new Impl(mode))
{
}


TestReferenceData::~TestReferenceData()
{
}


bool TestReferenceData::isWriteMode() const
{
    return _impl->_bWrite;
}


TestReferenceChecker TestReferenceData::rootChecker()
{
    if (!isWriteMode() && !_impl->_bInUse && _impl->_refDoc == NULL)
    {
        ADD_FAILURE() << "Reference data file not found: "
                      << _impl->_fullFilename;
    }
    _impl->_bInUse = true;
    if (_impl->_refDoc == NULL)
    {
        return TestReferenceChecker(new TestReferenceChecker::Impl(isWriteMode()));
    }
    xmlNodePtr rootNode = xmlDocGetRootElement(_impl->_refDoc);
    return TestReferenceChecker(
            new TestReferenceChecker::Impl("", rootNode, isWriteMode()));
}


/********************************************************************
 * TestReferenceChecker
 */

TestReferenceChecker::TestReferenceChecker(Impl *impl)
    : _impl(impl)
{
}


TestReferenceChecker::TestReferenceChecker(const TestReferenceChecker &other)
    : _impl(new Impl(*other._impl))
{
}


TestReferenceChecker &
TestReferenceChecker::operator =(const TestReferenceChecker &other)
{
    _impl.reset(new Impl(*other._impl));
    return *this;
}


TestReferenceChecker::~TestReferenceChecker()
{
}


bool TestReferenceChecker::isWriteMode() const
{
    return _impl->_bWrite;
}


bool TestReferenceChecker::checkPresent(bool bPresent, const char *id)
{
    if (isWriteMode())
    {
        return bPresent;
    }
    xmlNodePtr node = _impl->findNode(NULL, id);
    bool bFound = (node != NULL);
    if (bFound != bPresent)
    {
        ADD_FAILURE() << "Mismatch while checking reference data item'"
                          << _impl->appendPath(id) << "'\n"
                      << "Expected: " << (bPresent ? "it is present.\n" : "it is absent.\n")
                      << "  Actual: " << (bFound ? "it is present." : "it is absent.");
    }
    if (bFound && bPresent)
    {
        _impl->_nextSearchNode = node;
        return true;
    }
    return false;
}


TestReferenceChecker TestReferenceChecker::checkCompound(const char *type, const char *id)
{
    SCOPED_TRACE(_impl->traceString(id));
    if (_impl->shouldIgnore())
    {
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    const xmlChar *xmlNodeName = reinterpret_cast<const xmlChar *>(type);
    xmlNodePtr newNode = _impl->findOrCreateNode(xmlNodeName, id);
    if (newNode == NULL)
    {
        return TestReferenceChecker(new Impl(isWriteMode()));
    }
    return TestReferenceChecker(
            new Impl(_impl->appendPath(id), newNode, isWriteMode()));
}


void TestReferenceChecker::checkBoolean(bool value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(_impl->traceString(id));
    bool bFound = false;
    const char *strValue = value ? "true" : "false";
    std::string refStrValue =
        _impl->processItem(Impl::cBooleanNodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}


void TestReferenceChecker::checkString(const char *value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(_impl->traceString(id));
    bool bFound = false;
    std::string refStrValue =
        _impl->processItem(Impl::cStringNodeName, id, value, &bFound);
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
                                            const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(_impl->traceString(id));
    xmlNodePtr node = _impl->findOrCreateNode(Impl::cStringNodeName, id);
    if (node == NULL)
    {
        return;
    }
    // An extra newline is written in the beginning to make lines align
    // in the output xml (otherwise, the first line would be off by the length
    // of the starting CDATA tag).
    if (isWriteMode())
    {
        std::string adjustedValue = "\n" + value;
        const xmlChar *xmlValue
            = reinterpret_cast<const xmlChar *>(adjustedValue.c_str());
        // TODO: Figure out if \r and \r\n can be handled without them changing
        // to \n in the roundtrip
        xmlNodePtr cdata
            = xmlNewCDataBlock(node->doc, xmlValue,
                               static_cast<int>(adjustedValue.length()));
        xmlAddChild(node, cdata);
    }
    else
    {
        xmlNodePtr cdata = node->children;
        while (cdata != NULL && cdata->type != XML_CDATA_SECTION_NODE)
        {
            cdata = cdata->next;
        }
        if (cdata == NULL)
        {
            ADD_FAILURE() << "Invalid string block element";
            return;
        }
        xmlChar *refXmlValue = xmlNodeGetContent(cdata);
        if (refXmlValue[0] != '\n')
        {
            ADD_FAILURE() << "Invalid string block element";
            xmlFree(refXmlValue);
            return;
        }
        std::string refValue(reinterpret_cast<const char *>(refXmlValue + 1));
        xmlFree(refXmlValue);
        EXPECT_EQ(refValue, value);
    }
}


void TestReferenceChecker::checkInteger(int value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(_impl->traceString(id));
    bool bFound = false;
    std::string strValue = formatString("%d", value);
    std::string refStrValue =
        _impl->processItem(Impl::cIntegerNodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}


void TestReferenceChecker::checkDouble(double value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    SCOPED_TRACE(_impl->traceString(id));
    bool bFound = false;
    std::string strValue = formatString("%f", value);
    std::string refStrValue =
        _impl->processItem(Impl::cRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        char *endptr;
        double refValue = std::strtod(refStrValue.c_str(), &endptr);
        EXPECT_EQ('\0', *endptr);
        EXPECT_NEAR(refValue, value, 0.0001);
    }
}


void TestReferenceChecker::checkFloat(float value, const char *id)
{
    checkDouble(value, id);
}


void TestReferenceChecker::checkReal(float value, const char *id)
{
    checkDouble(value, id);
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
