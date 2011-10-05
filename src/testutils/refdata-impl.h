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
        static const xmlChar * const cXmlVersion;
        static const xmlChar * const cRootNodeName;
        static const xmlChar * const cCompoundNodeName;
        static const xmlChar * const cBooleanNodeName;
        static const xmlChar * const cStringNodeName;
        static const xmlChar * const cIntegerNodeName;
        static const xmlChar * const cRealNodeName;
        static const xmlChar * const cVectorIntegerNodeName;
        static const xmlChar * const cVectorRealNodeName;
        static const xmlChar * const cIdAttrName;
        static const xmlChar * const cCompoundTypeAttrName;
        static const char * const cSequenceIntegerType;
        static const char * const cSequenceRealType;
        static const char * const cSequenceVectorType;
        static const char * const cSequenceLengthName;

        Impl(ReferenceDataMode mode);
        ~Impl();

        xmlNodePtr findOrCreateNode(const xmlChar *name, const char *id);
        std::string processItem(const xmlChar *name, const char *id,
                                const char *value, bool *bFound);
        std::string processItem(const xmlChar *name, const char *id,
                                const std::string &value, bool *bFound);
        bool shouldIgnore() const;

        //! Full path of the reference data file.
        std::string             _fullFilename;
        /*! \brief
         * Whether the reference data is being written (true) or compared
         * (false).
         */
        bool                    _bWrite;
        /*! \brief
         * XML document for the reference data.
         *
         * May be NULL if there was an I/O error in initialization.
         */
        xmlDocPtr               _refDoc;
        /*! \brief
         * Current nore in \a _refDoc under which reference data is searched.
         *
         * On initialization, points to the root of \a _refDoc.  Is changed only
         * when compound groups are started or finished.
         *
         * Is NULL if and only if \a _refDoc is NULL.
         */
        xmlNodePtr              _currNode;
        /*! \brief
         * Points to a child of \a _currNode where the next search shoudl start.
         *
         * On initialization or when entering a compound node, points to the
         * first child of \a _currNode.  After every check, is updated to point
         * to the node following the one found, with possible wrapping.
         * On leaving a compound, is updated to point to the node following the
         * compound.
         *
         * Is NULL if and only if \a _currNode contains no children.
         * Otherwise, always points to a direct child of \a _currNode.
         */
        xmlNodePtr              _nextSearchNode;
        /*! \brief
         * Count of currently failed startCompound() calls.
         *
         * If larger than 0, currently we are within a non-existent compound,
         * so all comparisons should be ignored.
         */
        int                     _failedCompounds;
};

} // namespace test
} // namespace gmx

#endif
