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
 * Implements gmx::test::TestReferenceData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "refdata.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/utility/path.h"
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

namespace internal
{

void addGlobalReferenceDataEnvironment()
{
    ::testing::AddGlobalTestEnvironment(new TestReferenceDataEnvironment);
}

} // namespace internal

/********************************************************************
 * TestReferenceData::Impl
 */

const xmlChar * const TestReferenceData::Impl::cXmlVersion =
    (const xmlChar *)"1.0";
const xmlChar * const TestReferenceData::Impl::cRootNodeName =
    (const xmlChar *)"ReferenceData";
const xmlChar * const TestReferenceChecker::Impl::cCompoundNodeName =
    (const xmlChar *)"Compound";
const xmlChar * const TestReferenceChecker::Impl::cBooleanNodeName =
    (const xmlChar *)"Bool";
const xmlChar * const TestReferenceChecker::Impl::cStringNodeName =
    (const xmlChar *)"String";
const xmlChar * const TestReferenceChecker::Impl::cIntegerNodeName  =
    (const xmlChar *)"Int";
const xmlChar * const TestReferenceChecker::Impl::cRealNodeName =
    (const xmlChar *)"Real";
const xmlChar * const TestReferenceChecker::Impl::cVectorIntegerNodeName =
    (const xmlChar *)"Vector";
const xmlChar * const TestReferenceChecker::Impl::cVectorRealNodeName =
    (const xmlChar *)"Vector";
const xmlChar * const TestReferenceChecker::Impl::cIdAttrName =
    (const xmlChar *)"Name";
const xmlChar * const TestReferenceChecker::Impl::cCompoundTypeAttrName =
    (const xmlChar *)"Subtype";
const char * const TestReferenceChecker::Impl::cSequenceIntegerType =
    "SequenceInteger";
const char * const TestReferenceChecker::Impl::cSequenceRealType =
    "SequenceReal";
const char * const TestReferenceChecker::Impl::cSequenceVectorType =
    "SequenceVector";
const char * const TestReferenceChecker::Impl::cSequenceLengthName =
    "Length";


TestReferenceData::Impl::Impl(ReferenceDataMode mode)
    : _refDoc(NULL), _bWrite(false), _bInUse(false)
{
    const ::testing::TestInfo *test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    std::string dirname = getReferenceDataPath();
    std::string filename = std::string(test_info->test_case_name())
        + "_" + test_info->name() + ".xml";
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

TestReferenceChecker::Impl::Impl(xmlNodePtr rootNode, bool bWrite)
    : _currNode(rootNode),
      _nextSearchNode(rootNode != NULL ? rootNode->xmlChildrenNode : NULL),
      _bWrite(bWrite)
{
}


xmlNodePtr
TestReferenceChecker::Impl::findOrCreateNode(const xmlChar *name, const char *id)
{
    const xmlChar *xmlId = reinterpret_cast<const xmlChar *>(id);
    xmlNodePtr node = _nextSearchNode;
    while (node != NULL && node->next != _nextSearchNode)
    {
        if (xmlStrcmp(node->name, name) == 0)
        {
            xmlChar *refId = xmlGetProp(node, cIdAttrName);
            if (xmlId == NULL && refId == NULL)
            {
                break;
            }
            if (refId != NULL)
            {
                if (xmlId != NULL && xmlStrcmp(refId, xmlId) == 0)
                {
                    xmlFree(refId);
                    break;
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
    if (node == NULL || node->next == _nextSearchNode)
    {
        if (_bWrite)
        {
            node = xmlNewTextChild(_currNode, NULL, name, NULL);
            if (node != NULL && xmlId != NULL)
            {
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
        GMX_RELEASE_ASSERT(!_bWrite, "Node creation failed without exception");
        ADD_FAILURE() << "Reference data item not found";
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
    delete _impl;
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
    xmlNodePtr rootNode =
        _impl->_refDoc != NULL ? xmlDocGetRootElement(_impl->_refDoc) : NULL;
    return TestReferenceChecker(new TestReferenceChecker::Impl(rootNode, isWriteMode()));
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
    Impl *newImpl = new Impl(*other._impl);
    std::swap(_impl, newImpl);
    delete newImpl;
    return *this;
}


TestReferenceChecker::~TestReferenceChecker()
{
    delete _impl;
}


bool TestReferenceChecker::isWriteMode() const
{
    return _impl->_bWrite;
}


TestReferenceChecker TestReferenceChecker::checkCompound(const char *type, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return TestReferenceChecker(new Impl(NULL, isWriteMode()));
    }
    xmlNodePtr newNode = _impl->findOrCreateNode(Impl::cCompoundNodeName, id);
    if (newNode == NULL)
    {
        GMX_RELEASE_ASSERT(!isWriteMode(), "Node creation failed without exception");
        ADD_FAILURE() << "Reference data item not found";
        return TestReferenceChecker(new Impl(NULL, isWriteMode()));
    }
    if (isWriteMode())
    {
        if (xmlNewProp(newNode, Impl::cCompoundTypeAttrName,
                       reinterpret_cast<const xmlChar *>(type)) == NULL)
        {
            GMX_THROW(TestException("XML property creation failed"));
        }
    }
    return TestReferenceChecker(new Impl(newNode, isWriteMode()));
}


void TestReferenceChecker::checkBoolean(bool value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
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


void TestReferenceChecker::checkInteger(int value, const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    bool bFound = false;
    char strValue[20];
    snprintf(strValue, 20, "%d", value);
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
    bool bFound = false;
    char strValue[20];
    snprintf(strValue, 20, "%f", value);
    std::string refStrValue =
        _impl->processItem(Impl::cRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        char *endptr = NULL;
        double refValue = std::strtof(refStrValue.c_str(), &endptr);
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
    if (_impl->shouldIgnore())
    {
        return;
    }
    bool bFound = false;
    char strValue[50];
    snprintf(strValue, 50, "%d %d %d", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorIntegerNodeName, id, strValue, &bFound);
    if (bFound)
    {
        EXPECT_EQ(refStrValue, strValue);
    }
}


void TestReferenceChecker::checkVector(const float value[3], const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    bool bFound = false;
    char strValue[50];
    snprintf(strValue, 50, "%f %f %f", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        float refX, refY, refZ;
        int count = std::sscanf(refStrValue.c_str(), " %g %g %g", &refX, &refY, &refZ);
        if (count != 3)
        {
            GMX_THROW(TestException("Corrupt reference vector data"));
        }
        EXPECT_NEAR(refX, value[0], 0.0001);
        EXPECT_NEAR(refY, value[1], 0.0001);
        EXPECT_NEAR(refZ, value[2], 0.0001);
    }
}


void TestReferenceChecker::checkVector(const double value[3], const char *id)
{
    if (_impl->shouldIgnore())
    {
        return;
    }
    bool bFound = false;
    char strValue[50];
    snprintf(strValue, 50, "%f %f %f", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorRealNodeName, id, strValue, &bFound);
    if (bFound)
    {
        double refX, refY, refZ;
        int count = std::sscanf(refStrValue.c_str(), " %lg %lg %lg", &refX, &refY, &refZ);
        if (count != 3)
        {
            GMX_THROW(TestException("Corrupt reference vector data"));
        }
        EXPECT_NEAR(refX, value[0], 0.0001);
        EXPECT_NEAR(refY, value[1], 0.0001);
        EXPECT_NEAR(refZ, value[2], 0.0001);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const int *values,
                                              const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceIntegerType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkInteger(values[i], NULL);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const float *values,
                                              const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceRealType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkFloat(values[i], NULL);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const double *values,
                                              const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceRealType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkDouble(values[i], NULL);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const int values[][3],
                                              const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceVectorType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkVector(values[i], NULL);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const float values[][3],
                                              const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceVectorType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkVector(values[i], NULL);
    }
}


void TestReferenceChecker::checkSequenceArray(size_t length, const double values[][3],
                                             const char *id)
{
    TestReferenceChecker compound(checkCompound(Impl::cSequenceVectorType, id));
    compound.checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        compound.checkVector(values[i], NULL);
    }
}

} // namespace test
} // namespace gmx
