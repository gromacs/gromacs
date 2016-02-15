/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "external/tinyxml2/tinyxml2.h"
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

typedef tinyxml2::XMLDocument *XMLDocumentPtr;
typedef tinyxml2::XMLNode *XMLNodePtr;
typedef tinyxml2::XMLElement *XMLElementPtr;
typedef tinyxml2::XMLText *XMLTextPtr;

//! \name Helper functions for XML reading
//! \{

void readEntry(XMLNodePtr node, ReferenceDataEntry *entry);

XMLTextPtr getTextChildNode(XMLNodePtr node)
{
    // Note that when reading, we don't have to care if it is in a
    // CDATA section, or not.
    XMLNodePtr textNode = node->FirstChild();
    if (textNode == nullptr)
    {
        return nullptr;
    }
    while (textNode->ToText() == nullptr)
    {
        textNode = textNode->NextSibling();
        if (textNode == nullptr)
        {
            return nullptr;
        }
    }
    // TODO: Consider checking that there is no more than one CDATA section.
    return textNode->ToText();
}

bool hasTextContent(XMLNodePtr node)
{
    return getTextChildNode(node) != nullptr;
}

std::string getValueFromLeafElement(XMLNodePtr node)
{
    std::string value;
    if (hasTextContent(node))
    {
        XMLTextPtr contentNode = getTextChildNode(node);
        value = std::string(contentNode->Value());
        if (value.empty())
        {
            GMX_THROW(TestException("Invalid string block in reference data"));
        }
        if (value[0] == '\n')
        {
            value.erase(0, 1);
        }
    }
    return value;
}

ReferenceDataEntry::EntryPointer createEntry(XMLElementPtr element)
{
    const char *id = element->Attribute(c_IdAttrName);
    ReferenceDataEntry::EntryPointer entry(new ReferenceDataEntry(element->Value(), id));
    return entry;
}

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

bool isCompoundElement(XMLNodePtr node)
{
    return node->FirstChildElement() != nullptr;
}


void readEntry(XMLNodePtr element, ReferenceDataEntry *entry)
{
    if (isCompoundElement(element))
    {
        readChildEntries(element, entry);
    }
    else if (hasTextContent(element))
    {
        entry->setTextBlockValue(getValueFromLeafElement(element));
    }
    else
    {
        entry->setValue("");
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
        std::string errorString("Error was ");
        errorString += document.GetErrorStr1();
        errorString += document.GetErrorStr2();
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
    if (element == nullptr)
    {
        GMX_THROW(TestException("XML element creation failed"));
    }
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

void createElementContents(XMLElementPtr element, const ReferenceDataEntry &entry)
{
    if (entry.isCompound())
    {
        createChildElements(element, entry);
    }
    else if (entry.isTextBlock())
    {
        // An extra newline is written in the beginning to make lines align
        // in the output xml
        const std::string adjustedValue = "\n" + entry.value();
        // TODO: Figure out if \r and \r\n can be handled without them changing
        // to \n in the roundtrip
        XMLTextPtr textNode = element->GetDocument()->NewText(adjustedValue.c_str());
        if (textNode == nullptr)
        {
            GMX_THROW(TestException("Could not create text node for " + entry.value()));
        }
        textNode->SetCData(true);
        element->InsertEndChild(textNode);
    }
    else
    {
        XMLTextPtr textNode = element->GetDocument()->NewText(entry.value().c_str());
        if (textNode == nullptr)
        {
            GMX_THROW(TestException("Could not create text node for " + entry.value()));
        }
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
    if (rootElement == nullptr)
    {
        GMX_THROW(TestException("Could not create root element"));
    }
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
    if (declaration == nullptr)
    {
        GMX_THROW(TestException("Could not create XML declaration"));
    }
    document.InsertEndChild(declaration);

    declaration = document.NewDeclaration(c_StyleSheetDeclarationString);
    if (declaration == nullptr)
    {
        GMX_THROW(TestException("Could not create XML declaration"));
    }
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
