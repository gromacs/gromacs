/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Implements reference data XML persistence.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "refdata-xml.h"

#include <tinyxml2.h>

#include "gromacs/utility/exceptions.h"

#include "testutils/refdata-impl.h"
#include "testutils/testexceptions.h"

namespace gmx
{
namespace test
{

namespace
{

//! XML version declaration used when writing the reference data.
const char *const c_VersionDeclarationString = "xml version=\"1.0\"";
//! XML stylesheet declaration used for writing the reference data.
const char *const c_StyleSheetDeclarationString = "xml-stylesheet type=\"text/xsl\" href=\"referencedata.xsl\"";
//! Name of the root element in reference data XML files.
const char *const c_RootNodeName = "ReferenceData";
//! Name of the XML attribute used to store identifying strings for reference data elements.
const char *const c_IdAttrName = "Name";

}       // namespace

/********************************************************************
 * XML reading
 */

namespace
{

//! Convenience typedef
typedef tinyxml2::XMLDocument *XMLDocumentPtr;
//! Convenience typedef
typedef tinyxml2::XMLNode *XMLNodePtr;
//! Convenience typedef
typedef tinyxml2::XMLElement *XMLElementPtr;
//! Convenience typedef
typedef tinyxml2::XMLText *XMLTextPtr;

//! \name Helper functions for XML reading
//! \{

void readEntry(XMLNodePtr node, ReferenceDataEntry *entry);

XMLNodePtr getCDataChildNode(XMLNodePtr node)
{
    XMLNodePtr cdata = node->FirstChild();
    while (cdata != nullptr &&
           cdata->ToText() != nullptr &&
           !cdata->ToText()->CData())
    {
        cdata = cdata->NextSibling();
    }
    return cdata;
}

bool hasCDataContent(XMLNodePtr node)
{
    return getCDataChildNode(node) != nullptr;
}

//! Return a node convertible to text, either \c childNode or its first such sibling.
XMLNodePtr getNextTextChildNode(XMLNodePtr childNode)
{
    // Note that when reading, we don't have to care if it is in a
    // CDATA section, or not.
    while (childNode != nullptr)
    {
        if (childNode->ToText() != nullptr)
        {
            break;
        }
        childNode = childNode->NextSibling();
    }
    return childNode;
}

//! Return the concatenation of all the text children of \c node, including multiple CDATA children.
std::string getValueFromLeafElement(XMLNodePtr node)
{
    std::string value;

    XMLNodePtr  childNode = getNextTextChildNode(node->FirstChild());
    while (childNode != nullptr)
    {
        value += std::string(childNode->Value());

        childNode = getNextTextChildNode(childNode->NextSibling());
    }

    if (hasCDataContent(node))
    {
        // Prepare to strip the convenience newline added in
        // createElementContents, when writing CDATA content for
        // StringBlock data.
        if (value.empty() || value[0] != '\n')
        {
            GMX_THROW(TestException("Invalid string block in reference data"));
        }
        value.erase(0, 1);
    }

    return value;
}

//! Make a new entry from \c element.
ReferenceDataEntry::EntryPointer createEntry(XMLElementPtr element)
{
    const char *id = element->Attribute(c_IdAttrName);
    ReferenceDataEntry::EntryPointer entry(new ReferenceDataEntry(element->Value(), id));
    return entry;
}

//! Read the child entries of \c parentElement and transfer the contents to \c entry
void readChildEntries(XMLNodePtr parentElement, ReferenceDataEntry *entry)
{
    XMLElementPtr childElement = parentElement->FirstChildElement();
    while (childElement != nullptr)
    {
        ReferenceDataEntry::EntryPointer child(createEntry(childElement));
        readEntry(childElement, child.get());
        entry->addChild(move(child));
        childElement = childElement->NextSiblingElement();
    }
}

//! Return whether \c node has child XML elements (rather than text content).
bool isCompoundElement(XMLNodePtr node)
{
    return node->FirstChildElement() != nullptr;
}

//! Read \c element and transfer the contents to \c entry
void readEntry(XMLNodePtr element, ReferenceDataEntry *entry)
{
    if (isCompoundElement(element))
    {
        readChildEntries(element, entry);
    }
    else if (hasCDataContent(element))
    {
        entry->setTextBlockValue(getValueFromLeafElement(element));
    }
    else
    {
        entry->setValue(getValueFromLeafElement(element));
    }
}

//! \}

}       // namespace

//! \cond internal
ReferenceDataEntry::EntryPointer
readReferenceDataFile(const std::string &path)
{
    tinyxml2::XMLDocument document;
    document.LoadFile(path.c_str());
    if (document.Error())
    {
        const char *errorStr1 = document.GetErrorStr1();
        const char *errorStr2 = document.GetErrorStr2();
        std::string errorString("Error was ");
        if (errorStr1)
        {
            errorString += errorStr1;
        }
        if (errorStr2)
        {
            errorString += errorStr2;
        }
        if (!errorStr1 && !errorStr2)
        {
            errorString += "not specified.";
        }
        GMX_THROW(TestException("Reference data not parsed successfully: " + path + "\n." + errorString + "\n"));
    }
    XMLElementPtr rootNode = document.RootElement();
    if (rootNode == nullptr)
    {
        GMX_THROW(TestException("Reference data is empty: " + path));
    }
    if (std::strcmp(rootNode->Value(), c_RootNodeName) != 0)
    {
        GMX_THROW(TestException("Invalid root node type in " + path));
    }

    ReferenceDataEntry::EntryPointer rootEntry(ReferenceDataEntry::createRoot());
    readEntry(rootNode, rootEntry.get());
    return rootEntry;
}
//! \endcond

/********************************************************************
 * XML writing
 */

namespace
{

//! \name Helper functions for XML writing
//! \{

void createElementAndContents(XMLElementPtr             parentElement,
                              const ReferenceDataEntry &entry);

void setIdAttribute(XMLElementPtr element, const std::string &id)
{
    if (!id.empty())
    {
        element->SetAttribute(c_IdAttrName, id.c_str()); // If this fails, it throws std::bad_alloc
    }
}

XMLElementPtr createElement(XMLElementPtr parentElement, const ReferenceDataEntry &entry)
{
    XMLElementPtr element = parentElement->GetDocument()->NewElement(entry.type().c_str());
    parentElement->InsertEndChild(element);
    setIdAttribute(element, entry.id()); // If this fails, it throws std::bad_alloc
    return element;
}

void createChildElements(XMLElementPtr parentElement, const ReferenceDataEntry &entry)
{
    const ReferenceDataEntry::ChildList &children(entry.children());
    ReferenceDataEntry::ChildIterator    child;
    for (child = children.begin(); child != children.end(); ++child)
    {
        createElementAndContents(parentElement, **child);
    }
}

/*! \brief Handle \c input intended to be written as CDATA
 *
 * This method searches for any ']]>' sequences embedded in \c input,
 * because this must always end a CDATA field. If any are found, it
 * breaks the string so that instead multiple CDATA fields will be
 * written with that token sequence split across the fields. Note that
 * tinyxml2 does not handle such things itself.
 *
 * This is an edge case that is unimportant for GROMACS refdata, but
 * it is preferable to know that the infrastructure won't break.
 */
std::vector<std::string> breakUpAnyCdataEndTags(const std::string &input)
{
    std::vector<std::string> strings;
    std::size_t              startPos = 0;
    std::size_t              endPos;

    do
    {
        endPos = input.find("]]>", startPos);
        if (endPos != std::string::npos)
        {
            // We found an embedded CDATA end tag, so arrange to split it into multiple CDATA blocks
            endPos++;
        }
        strings.push_back(input.substr(startPos, endPos));
        startPos = endPos;
    }
    while (endPos != std::string::npos);

    return strings;
}

void createElementContents(XMLElementPtr element, const ReferenceDataEntry &entry)
{
    // TODO: Figure out if \r and \r\n can be handled without them
    // changing to \n in the roundtrip.
    if (entry.isCompound())
    {
        createChildElements(element, entry);
    }
    else if (entry.isTextBlock())
    {
        // An extra newline is written in the beginning to make lines align
        // in the output xml (otherwise, the first line would be off by the length
        // of the starting CDATA tag).
        const std::string        adjustedValue = "\n" + entry.value();
        std::vector<std::string> cdataStrings  = breakUpAnyCdataEndTags(adjustedValue);
        for (auto const &s : cdataStrings)
        {
            XMLTextPtr textNode = element->GetDocument()->NewText(s.c_str());
            textNode->SetCData(true);
            element->InsertEndChild(textNode);
        }
    }
    else
    {
        XMLTextPtr textNode = element->GetDocument()->NewText(entry.value().c_str());
        element->InsertEndChild(textNode);
    }
}

void createElementAndContents(XMLElementPtr parentElement, const ReferenceDataEntry &entry)
{
    XMLElementPtr element = createElement(parentElement, entry);
    createElementContents(element, entry);
}

XMLElementPtr createRootElement(XMLDocumentPtr document)
{
    XMLElementPtr rootElement = document->NewElement(c_RootNodeName);
    document->InsertEndChild(rootElement);
    return rootElement;
}

//! \}

}       // namespace

//! \cond internal
void writeReferenceDataFile(const std::string        &path,
                            const ReferenceDataEntry &rootEntry)
{
    // TODO: Error checking
    tinyxml2::XMLDocument     document;

    tinyxml2::XMLDeclaration *declaration = document.NewDeclaration(c_VersionDeclarationString);
    document.InsertEndChild(declaration);

    declaration = document.NewDeclaration(c_StyleSheetDeclarationString);
    document.InsertEndChild(declaration);

    XMLElementPtr rootElement = createRootElement(&document);
    createChildElements(rootElement, rootEntry);

    if (document.SaveFile(path.c_str()) != tinyxml2::XML_NO_ERROR)
    {
        GMX_THROW(TestException("Reference data saving failed in " + path));
    }
}
//! \endcond

} // namespace test
} // namespace gmx
