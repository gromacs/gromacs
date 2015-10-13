/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "refdata-xml.h"

#include <libxml/parser.h>
#include <libxml/xmlmemory.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/scoped_cptr.h"

#include "testutils/refdata-impl.h"
#include "testutils/testexceptions.h"

namespace gmx
{
namespace test
{

namespace
{

//! XML version used for writing the reference data.
const xmlChar *const cXmlVersion =
    reinterpret_cast<const xmlChar *>("1.0");
//! Name of the XML processing instruction used for XSLT reference.
const xmlChar *const cXmlStyleSheetNodeName =
    reinterpret_cast<const xmlChar *>("xml-stylesheet");
//! XSLT reference written to the reference data XML files.
const xmlChar *const cXmlStyleSheetContent =
    reinterpret_cast<const xmlChar *>("type=\"text/xsl\" href=\"referencedata.xsl\"");
//! Name of the root element in reference data XML files.
const xmlChar *const cRootNodeName =
    reinterpret_cast<const xmlChar *>("ReferenceData");
//! Name of the XML attribute used to store identifying strings for reference data elements.
const xmlChar *const cIdAttrName =
    reinterpret_cast<const xmlChar *>("Name");

/********************************************************************
 * Generic helper functions and classes
 */

//! Helper function to convert strings to xmlChars.
const xmlChar *toXmlString(const std::string &str)
{
    // TODO: Consider asserting that str is ASCII.
    return reinterpret_cast<const xmlChar *>(str.c_str());
}

//! Helper function to convert strings from xmlChars.
const char *fromXmlString(const xmlChar *str)
{
    // TODO: Consider asserting that str is ASCII.
    return reinterpret_cast<const char *>(str);
}

class XmlString
{
    public:
        explicit XmlString(xmlChar *str) : str_(str) {}
        ~XmlString()
        {
            if (str_ != NULL)
            {
                xmlFree(str_);
            }
        }

        std::string toString() const
        {
            if (str_ == NULL)
            {
                return std::string();
            }
            return std::string(fromXmlString(str_));
        }

    private:
        xmlChar *str_;
};

//! C++ wrapper for xmlDocPtr for exception safety.
typedef scoped_cptr<xmlDoc, xmlFreeDoc> XmlDocumentPointer;

}       // namespace

/********************************************************************
 * XML reading
 */

namespace
{

//! \name Helper functions for XML reading
//! \{

void readEntry(xmlNodePtr node, ReferenceDataEntry *entry);

xmlNodePtr getCDataChildNode(xmlNodePtr node)
{
    xmlNodePtr cdata = node->children;
    while (cdata != NULL && cdata->type != XML_CDATA_SECTION_NODE)
    {
        cdata = cdata->next;
    }
    // TODO: Consider checking that there is no more than one CDATA section.
    return cdata;
}

bool hasCDataContent(xmlNodePtr node)
{
    return getCDataChildNode(node) != NULL;
}

xmlNodePtr findContentNode(xmlNodePtr node)
{
    xmlNodePtr cdata = getCDataChildNode(node);
    return cdata != NULL ? cdata : node;
}

std::string getValueFromLeafElement(xmlNodePtr node)
{
    xmlNodePtr  contentNode = findContentNode(node);
    XmlString   content(xmlNodeGetContent(contentNode));
    std::string value(content.toString());
    if (hasCDataContent(node))
    {
        if (value.empty() || value[0] != '\n')
        {
            GMX_THROW(TestException("Invalid CDATA string block in reference data"));
        }
        value.erase(0, 1);
    }
    return value;
}

ReferenceDataEntry::EntryPointer createEntry(xmlNodePtr element)
{
    XmlString id(xmlGetProp(element, cIdAttrName));
    ReferenceDataEntry::EntryPointer entry(
            new ReferenceDataEntry(fromXmlString(element->name),
                                   id.toString().c_str()));
    return move(entry);
}

void readChildEntries(xmlNodePtr parentElement, ReferenceDataEntry *entry)
{
    xmlNodePtr childElement = xmlFirstElementChild(parentElement);
    while (childElement != NULL)
    {
        ReferenceDataEntry::EntryPointer child(createEntry(childElement));
        readEntry(childElement, child.get());
        entry->addChild(move(child));
        childElement = xmlNextElementSibling(childElement);
    }
}

bool isCompoundElement(xmlNodePtr node)
{
    return xmlFirstElementChild(node) != NULL;
}

void readEntry(xmlNodePtr element, ReferenceDataEntry *entry)
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
    XmlDocumentPointer document(xmlParseFile(path.c_str()));
    if (!document)
    {
        GMX_THROW(TestException("Reference data not parsed successfully: " + path));
    }
    xmlNodePtr rootNode = xmlDocGetRootElement(document.get());
    if (rootNode == NULL)
    {
        GMX_THROW(TestException("Reference data is empty: " + path));
    }
    if (xmlStrcmp(rootNode->name, cRootNodeName) != 0)
    {
        GMX_THROW(TestException("Invalid root node type in " + path));
    }

    ReferenceDataEntry::EntryPointer rootEntry(ReferenceDataEntry::createRoot());
    readEntry(rootNode, rootEntry.get());
    return move(rootEntry);
}
//! \endcond

/********************************************************************
 * XML writing
 */

namespace
{

//! \name Helper functions for XML writing
//! \{

void createElementAndContents(xmlNodePtr                parentNode,
                              const ReferenceDataEntry &entry);

void setIdAttribute(xmlNodePtr node, const std::string &id)
{
    if (!id.empty())
    {
        const xmlChar *xmlId = toXmlString(id);
        xmlAttrPtr     prop  = xmlNewProp(node, cIdAttrName, xmlId);
        if (prop == NULL)
        {
            GMX_THROW(TestException("XML attribute creation failed"));
        }
    }
}

xmlNodePtr createElement(xmlNodePtr parentNode, const ReferenceDataEntry &entry)
{
    xmlNodePtr node = xmlNewTextChild(parentNode, NULL, toXmlString(entry.type()), NULL);
    if (node == NULL)
    {
        GMX_THROW(TestException("XML element creation failed"));
    }
    setIdAttribute(node, entry.id());
    return node;
}

void createChildElements(xmlNodePtr parentNode, const ReferenceDataEntry &entry)
{
    const ReferenceDataEntry::ChildList &children(entry.children());
    ReferenceDataEntry::ChildIterator    child;
    for (child = children.begin(); child != children.end(); ++child)
    {
        createElementAndContents(parentNode, **child);
    }
}

void createElementContents(xmlNodePtr node, const ReferenceDataEntry &entry)
{
    if (entry.isCompound())
    {
        createChildElements(node, entry);
    }
    else if (entry.isTextBlock())
    {
        // An extra newline is written in the beginning to make lines align
        // in the output xml (otherwise, the first line would be off by the length
        // of the starting CDATA tag).
        const std::string adjustedValue = "\n" + entry.value();
        // TODO: Figure out if \r and \r\n can be handled without them changing
        // to \n in the roundtrip
        xmlNodePtr cdata
            = xmlNewCDataBlock(node->doc, toXmlString(adjustedValue),
                               static_cast<int>(adjustedValue.length()));
        xmlAddChild(node, cdata);
    }
    else
    {
        xmlNodeAddContent(node, toXmlString(entry.value()));
    }
}

void createElementAndContents(xmlNodePtr parentNode, const ReferenceDataEntry &entry)
{
    xmlNodePtr node = createElement(parentNode, entry);
    createElementContents(node, entry);
}

xmlNodePtr createRootElement(xmlDocPtr document)
{
    xmlNodePtr rootElement = xmlNewDocNode(document, NULL, cRootNodeName, NULL);
    xmlDocSetRootElement(document, rootElement);
    return rootElement;
}

void createXsltReference(xmlDocPtr document, xmlNodePtr rootElement)
{
    xmlNodePtr xslNode = xmlNewDocPI(document, cXmlStyleSheetNodeName,
                                     cXmlStyleSheetContent);
    xmlAddPrevSibling(rootElement, xslNode);
}

//! \}

}       // namespace

//! \cond internal
void writeReferenceDataFile(const std::string        &path,
                            const ReferenceDataEntry &rootEntry)
{
    // TODO: Error checking
    XmlDocumentPointer  document(xmlNewDoc(cXmlVersion));
    xmlNodePtr          rootElement = createRootElement(document.get());
    createXsltReference(document.get(), rootElement);
    createChildElements(rootElement, rootEntry);

    if (xmlSaveFormatFile(path.c_str(), document.get(), 1) == -1)
    {
        GMX_THROW(TestException("Reference data saving failed in " + path));
    }
}
//! \endcond

/********************************************************************
 * Cleanup
 */

//! \cond internal
void cleanupReferenceData()
{
    xmlCleanupParser();
}
//! \endcond

} // namespace test
} // namespace gmx
