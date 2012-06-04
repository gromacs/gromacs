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
 * Declares private implementation class for gmx::test::TestReferenceData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_IMPL_H
#define GMX_TESTUTILS_REFDATA_IMPL_H

#include <string>

#include <libxml/parser.h>

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{

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
        explicit Impl(ReferenceDataMode mode);
        ~Impl();

        //! Full path of the reference data file.
        std::string             _fullFilename;
        /*! \brief
         * XML document for the reference data.
         *
         * May be NULL if there was an I/O error in initialization.
         */
        xmlDocPtr               _refDoc;
        /*! \brief
         * Whether the reference data is being written (true) or compared
         * (false).
         */
        bool                    _bWrite;
        /*! \brief
         * Whether any reference checkers have been created for this data.
         */
        bool                    _bInUse;
};


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
        //! String constant for naming XML elements for floating-point values.
        static const xmlChar * const cRealNodeName;
        //! String constant for naming XML attribute for value identifiers.
        static const xmlChar * const cIdAttrName;
        //! String constant for naming compounds for vectors.
        static const char * const cVectorType;
        //! String constant for naming compounds for sequences.
        static const char * const cSequenceType;
        //! String constant for value identifier for sequence length.
        static const char * const cSequenceLengthName;

        //! Creates a checker that does nothing.
        explicit Impl(bool bWrite);
        //! Creates a checker with a given root node.
        Impl(const std::string &path, xmlNodePtr rootNode, bool bWrite);

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
         * \returns Matching node, or NULL if no matching node found
         *      (NULL is never returned in write mode).
         * \throws  TestException if node creation fails in write mode.
         *
         * Finds a node using findNode() and updates the match pointer is a
         * match is found.  If a match is not found, the method returns NULL in
         * read mode and creates a new node in write mode.  If the creation
         * fails in write mode, throws.
         */
        xmlNodePtr findOrCreateNode(const xmlChar *name, const char *id);
        /*! \brief
         * Helper method for checking a reference data value.
         *
         * \param[in]  name   Type of node to find.
         * \param[in]  id     Unique identifier of the node (can be NULL, in
         *      which case the next node without an id is matched).
         * \param[in]  value  String value of the value to be compared.
         * \param[out] bFound true if a matchin value was found.
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

        /*! \brief
         * Human-readable path to the root node of this checker.
         *
         * For the root checker, this will be "/", and for each compound, the
         * id of the compound is added.  Used for reporting comparison
         * mismatches.
         */
        std::string             _path;
        /*! \brief
         * Current node under which reference data is searched.
         *
         * Points to either the root of TestReferenceData::Impl::_refDoc, or to
         * a compound node.
         *
         * Can be NULL, in which case this checker does nothing (doesn't even
         * report errors, see shouldIgnore()).
         */
        xmlNodePtr              _currNode;
        /*! \brief
         * Points to a child of \a _currNode where the next search should start.
         *
         * On initialization, points to the first child of \a _currNode.  After
         * every check, is updated to point to the node following the one
         * found, with possible wrapping.
         *
         * Is NULL if and only if \a _currNode contains no children.
         * Otherwise, always points to a direct child of \a _currNode.
         */
        xmlNodePtr              _nextSearchNode;
        /*! \brief
         * Whether the reference data is being written (true) or compared
         * (false).
         */
        bool                    _bWrite;
        /*! \brief
         * Current number of unnamed elements in a sequence.
         *
         * It is the index of the next added unnamed element.
         */
        int                     _seqIndex;
};

} // namespace test
} // namespace gmx

#endif
