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
 * Declares classes for (partial) parsing of reStructuredText.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_RSTPARSER_H
#define GMX_ONLINEHELP_RSTPARSER_H

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class TextLineWrapperSettings;

/*! \internal
 * \brief
 * Iterator over reStructuredText paragraphs.
 *
 * After initialization, nextParagraph() needs to be called to access the first
 * paragraph.  Subsequence paragraphs can be accessed by repeated calls to
 * nextParagraph().  After the last paragraph, nextParagraph() returns `false`.
 *
 * After each call to nextParagraph(), other methods can be called to query
 * details of the current paragraph.
 *
 * \ingroup module_onlinehelp
 */
class RstParagraphIterator
{
    public:
        /*! \brief
         * Initializes an iterator for given input text.
         *
         * Does not throw.
         */
        explicit RstParagraphIterator(const std::string &text);

        /*! \brief
         * Advances the iterator to the next paragraph.
         *
         * \returns `false` if there were no more paragraphs.
         *
         * Does not throw (except std::bad_alloc if std::string::compare()
         * throws).
         */
        bool nextParagraph();

        //! Returns the indentation for first line of this paragraph.
        int firstLineIndent() const { return firstLineIndent_; }
        //! Returns the indentation for subsequent lines of this paragraph.
        int indent() const { return indent_; }
        /*! \brief
         * Returns the text
         *
         * \param[out] result  Variable to receive the paragraph text.
         * \throws std::bad_alloc if out of memory.
         *
         * Indentation and internal line breaks have been stripped from the
         * paragraph text (except for literal blocks etc.).  For literal
         * blocks, the common indentation has been stripped and is returned in
         * indent() instead.
         *
         * Leading newlines are returned to indicate necessary separation from
         * the preceding paragraph.
         */
        void getParagraphText(std::string *result) const;

    private:
        enum ParagraphType
        {
            eParagraphType_Normal,
            eParagraphType_Literal,
            eParagraphType_Title
        };

        //! The text to iterate over.
        const std::string &text_;

        //! Start of the current paragraph.
        size_t             begin_;
        //! End of the current paragraph (C++-style iterator).
        size_t             end_;
        //! Type of the current paragraph.
        ParagraphType      type_;
        //! Number of newlines to print before the current paragraph.
        int                breakSize_;
        //! Indentation of the first line of this paragraph.
        int                firstLineIndent_;
        //! (Minimum) indentation of other lines in this paragraph.
        int                indent_;

        //! Start of the next paragrah.
        size_t             nextBegin_;
        //! Number of newlines to print after the current paragraph.
        int                nextBreakSize_;
        /*! \brief
         * Indentation of the preceding paragraph that contained `::`.
         *
         * If the next paragraph is not a literal block, the value is `-1`.
         */
        int                literalIndent_;

        GMX_DISALLOW_COPY_AND_ASSIGN(RstParagraphIterator);
};

} // namespace gmx

#endif
