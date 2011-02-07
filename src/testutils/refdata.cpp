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

#include "gromacs/utility/path.h"

#include "refdata-impl.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * TestReferenceData::Impl
 */

const xmlChar * const TestReferenceData::Impl::cXmlVersion =
    (const xmlChar *)"1.0";
const xmlChar * const TestReferenceData::Impl::cRootNodeName =
    (const xmlChar *)"ReferenceData";
const xmlChar * const TestReferenceData::Impl::cCompoundNodeName =
    (const xmlChar *)"Compound";
const xmlChar * const TestReferenceData::Impl::cBooleanNodeName =
    (const xmlChar *)"Bool";
const xmlChar * const TestReferenceData::Impl::cStringNodeName =
    (const xmlChar *)"String";
const xmlChar * const TestReferenceData::Impl::cIntegerNodeName  =
    (const xmlChar *)"Int";
const xmlChar * const TestReferenceData::Impl::cRealNodeName =
    (const xmlChar *)"Real";
const xmlChar * const TestReferenceData::Impl::cVectorIntegerNodeName =
    (const xmlChar *)"Vector";
const xmlChar * const TestReferenceData::Impl::cVectorRealNodeName =
    (const xmlChar *)"Vector";
const xmlChar * const TestReferenceData::Impl::cIdAttrName =
    (const xmlChar *)"Name";
const xmlChar * const TestReferenceData::Impl::cCompoundTypeAttrName =
    (const xmlChar *)"Subtype";
const char * const TestReferenceData::Impl::cSequenceIntegerType =
    "SequenceInteger";
const char * const TestReferenceData::Impl::cSequenceRealType =
    "SequenceReal";
const char * const TestReferenceData::Impl::cSequenceVectorType =
    "SequenceVector";
const char * const TestReferenceData::Impl::cSequenceLengthName =
    "Length";


TestReferenceData::Impl::Impl()
    : _bWrite(false), _refDoc(NULL), _currNode(NULL), _nextSearchNode(NULL)
{
    const ::testing::TestInfo *test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    std::string dirname = getReferenceDataPath();
    std::string filename = std::string(test_info->test_case_name())
        + "_" + test_info->name() + ".xml";
    _fullFilename = Path::join(dirname, filename);
}


void TestReferenceData::Impl::Initialize(ReferenceDataMode mode)
{
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
            FAIL() << "Could not find reference data file: " << _fullFilename;
        }
    }
    if (_bWrite)
    {
        _refDoc = xmlNewDoc(cXmlVersion);
        _currNode = xmlNewDocNode(_refDoc, NULL, cRootNodeName, NULL);
        _nextSearchNode = NULL;
        xmlDocSetRootElement(_refDoc, _currNode);
    }
    else
    {
        _refDoc = xmlParseFile(_fullFilename.c_str());
        if (_refDoc == NULL)
        {
            FAIL() << "Reference data not parsed successfully: " << _fullFilename;
        }
        _currNode = xmlDocGetRootElement(_refDoc);
        if (_currNode == NULL)
        {
            FAIL() << "Reference data is empty: " << _fullFilename;
        }
        if (xmlStrcmp(_currNode->name, cRootNodeName))
        {
            FAIL() << "Invalid root node type in " << _fullFilename;
        }
        _nextSearchNode = _currNode->xmlChildrenNode;
    }
}


TestReferenceData::Impl::~Impl()
{
    if (_bWrite)
    {
        EXPECT_NE(-1, xmlSaveFormatFile(_fullFilename.c_str(), _refDoc, 1));
    }
    if (_refDoc != NULL)
    {
        xmlFreeDoc(_refDoc);
    }
}


xmlNodePtr
TestReferenceData::Impl::findOrCreateNode(const xmlChar *name, const char *id)
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
TestReferenceData::Impl::processItem(const xmlChar *name, const char *id,
                                     const char *value)
{
    xmlNodePtr node = findOrCreateNode(name, id);
    if (node == NULL)
    {
        return std::string();
    }
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
TestReferenceData::Impl::processItem(const xmlChar *name, const char *id,
                                     const std::string &value)
{
    return processItem(name, id, value.c_str());
}


/********************************************************************
 * TestReferenceData
 */

TestReferenceData::TestReferenceData()
    : _impl(new Impl)
{
    _impl->Initialize(getReferenceDataMode());
}


TestReferenceData::TestReferenceData(ReferenceDataMode mode)
    : _impl(new Impl)
{
    _impl->Initialize(mode);
}


TestReferenceData::~TestReferenceData()
{
    delete _impl;
}


bool TestReferenceData::isWriteMode() const
{
    return _impl->_bWrite;
}


void TestReferenceData::startCompound(const char *type, const char *id)
{
    xmlNodePtr newNode = _impl->findOrCreateNode(Impl::cCompoundNodeName, id);
    ASSERT_TRUE(NULL != newNode);
    _impl->_currNode = newNode;
    _impl->_nextSearchNode = newNode->xmlChildrenNode;
    if (isWriteMode())
    {
        ASSERT_TRUE(NULL != xmlNewProp(newNode, Impl::cCompoundTypeAttrName,
                                       reinterpret_cast<const xmlChar *>(type)));
    }
}


void TestReferenceData::finishCompound()
{
    ASSERT_TRUE(NULL != _impl->_currNode);
    ASSERT_TRUE(NULL != _impl->_currNode->parent);
    _impl->_nextSearchNode = _impl->_currNode->next;
    _impl->_currNode = _impl->_currNode->parent;
    if (_impl->_nextSearchNode == NULL)
    {
        _impl->_nextSearchNode = _impl->_currNode->xmlChildrenNode;
    }
}


void TestReferenceData::checkBoolean(bool value, const char *id)
{
    const char *strValue = value ? "true" : "false";
    std::string refStrValue =
        _impl->processItem(Impl::cBooleanNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    EXPECT_EQ(refStrValue, strValue);
}


void TestReferenceData::checkString(const char *value, const char *id)
{
    std::string refStrValue =
        _impl->processItem(Impl::cStringNodeName, id, value);
    ASSERT_FALSE(refStrValue.empty());
    EXPECT_EQ(refStrValue, value);
}


void TestReferenceData::checkString(const std::string &value, const char *id)
{
    checkString(value.c_str(), id);
}


void TestReferenceData::checkInteger(int value, const char *id)
{
    char strValue[20];
    snprintf(strValue, 20, "%d", value);
    std::string refStrValue =
        _impl->processItem(Impl::cIntegerNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    EXPECT_EQ(refStrValue, strValue);
}


void TestReferenceData::checkDouble(double value, const char *id)
{
    char strValue[20];
    snprintf(strValue, 20, "%f", value);
    std::string refStrValue =
        _impl->processItem(Impl::cRealNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    char *endptr = NULL;
    double refValue = std::strtof(refStrValue.c_str(), &endptr);
    EXPECT_EQ('\0', *endptr);
    EXPECT_NEAR(refValue, value, 0.0001);
}


void TestReferenceData::checkFloat(float value, const char *id)
{
    checkDouble(value, id);
}


void TestReferenceData::checkReal(float value, const char *id)
{
    checkDouble(value, id);
}


void TestReferenceData::checkReal(double value, const char *id)
{
    checkDouble(value, id);
}


void TestReferenceData::checkVector(int value[3], const char *id)
{
    char strValue[50];
    snprintf(strValue, 50, "%d %d %d", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorIntegerNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    EXPECT_EQ(refStrValue, strValue);
}


void TestReferenceData::checkVector(float value[3], const char *id)
{
    char strValue[50];
    snprintf(strValue, 50, "%f %f %f", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorRealNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    double refX, refY, refZ;
    int count = std::sscanf(refStrValue.c_str(), " %lg %lg %lg", &refX, &refY, &refZ);
    ASSERT_EQ(3, count);
    EXPECT_NEAR(refX, value[0], 0.0001);
    EXPECT_NEAR(refY, value[1], 0.0001);
    EXPECT_NEAR(refZ, value[2], 0.0001);
}


void TestReferenceData::checkVector(double value[3], const char *id)
{
    char strValue[50];
    snprintf(strValue, 50, "%f %f %f", value[0], value[1], value[2]);
    std::string refStrValue =
        _impl->processItem(Impl::cVectorRealNodeName, id, strValue);
    ASSERT_FALSE(refStrValue.empty());
    float refX, refY, refZ;
    int count = std::sscanf(refStrValue.c_str(), " %g %g %g", &refX, &refY, &refZ);
    ASSERT_EQ(3, count);
    EXPECT_NEAR(refX, value[0], 0.0001);
    EXPECT_NEAR(refY, value[1], 0.0001);
    EXPECT_NEAR(refZ, value[2], 0.0001);
}


void TestReferenceData::checkSequenceInteger(size_t length, int *values,
                                             const char *id)
{
    startCompound(Impl::cSequenceIntegerType, id);
    checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        checkInteger(values[i], NULL);
    }
    finishCompound();
}


void TestReferenceData::checkSequenceDouble(size_t length, double *values,
                                            const char *id)
{
    startCompound(Impl::cSequenceRealType, id);
    checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        checkDouble(values[i], NULL);
    }
    finishCompound();
}


void TestReferenceData::checkSequenceVector(size_t length, int values[][3],
                                            const char *id)
{
    startCompound(Impl::cSequenceVectorType, id);
    checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        checkVector(values[i], NULL);
    }
    finishCompound();
}


void TestReferenceData::checkSequenceVector(size_t length, float values[][3],
                                            const char *id)
{
    startCompound(Impl::cSequenceVectorType, id);
    checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        checkVector(values[i], NULL);
    }
    finishCompound();
}


void TestReferenceData::checkSequenceVector(size_t length, double values[][3],
                                            const char *id)
{
    startCompound(Impl::cSequenceVectorType, id);
    checkInteger(static_cast<int>(length), Impl::cSequenceLengthName);
    for (size_t i = 0; i < length; ++i)
    {
        checkVector(values[i], NULL);
    }
    finishCompound();
}

} // namespace test
} // namespace gmx
