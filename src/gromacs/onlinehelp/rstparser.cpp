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
 * Implements classes from rstparser.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "rstparser.h"

#include <cctype>

#include <algorithm>

#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \brief
 * Counts the number of leading spaces in a text range.
 *
 * Does not throw.
 */
int countLeadingSpace(const std::string &text, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        if (!std::isspace(text[i]))
        {
            return i - start;
        }
    }
    return end - start;
}

/*! \brief
 * Returns `true` if a list item starts in \p text at \p index.
 *
 * Does not throw.
 */
bool startsListItem(const std::string &text, size_t index)
{
    if (text.length() <= index + 1)
    {
        return false;
    }
    if (text[index] == '*' && std::isspace(text[index+1]))
    {
        return true;
    }
    if (std::isdigit(text[index]))
    {
        while (index < text.length() && std::isdigit(text[index]))
        {
            ++index;
        }
        if (text.length() > index + 1 && text[index] == '.'
            && std::isspace(text[index+1]))
        {
            return true;
        }
    }
    return false;
}

/*! \brief
 * Returns `true` if a table starts in \p text at \p index.
 *
 * The function only inspects the first line for something that looks like a
 * reStructuredText table, and accepts also some malformed tables.
 * Any issues should be apparent when Sphinx parses the reStructuredText
 * export, so full validation is not done here.
 *
 * Does not throw.
 */
bool startsTable(const std::string &text, size_t index)
{
    if (text[index] == '=')
    {
        while (index < text.length() && text[index] != '\n')
        {
            if (text[index] != '=' && !std::isspace(text[index]))
            {
                return false;
            }
            ++index;
        }
        return true;
    }
    else if (text[index] == '+')
    {
        while (index < text.length() && text[index] != '\n')
        {
            if (text[index] != '-' && text[index] != '+')
            {
                return false;
            }
            ++index;
        }
        return true;
    }
    return false;
}

/*! \brief
 * Returns `true` if a line in \p text starting at \p index is a title underline.
 *
 * Does not throw.
 */
bool isTitleUnderline(const std::string &text, size_t index)
{
    const char firstChar = text[index];
    if (std::ispunct(firstChar))
    {
        while (index < text.length() && text[index] != '\n')
        {
            if (text[index] != firstChar)
            {
                return false;
            }
            ++index;
        }
        return true;
    }
    return false;
}

}    // namespace

/********************************************************************
 * RstParagraphIterator
 */

RstParagraphIterator::RstParagraphIterator(const std::string &text)
    : text_(text), begin_(0), end_(0), type_(eParagraphType_Normal),
      breakSize_(0), firstLineIndent_(0), indent_(0),
      nextBegin_(0), nextBreakSize_(0), literalIndent_(-1)
{
}

bool RstParagraphIterator::nextParagraph()
{
    begin_     = nextBegin_;
    type_      = eParagraphType_Normal;
    breakSize_ = nextBreakSize_;
    // Skip leading newlines (includes those separating paragraphs).
    while (begin_ < text_.length() && text_[begin_] == '\n')
    {
        ++begin_;
    }
    if (begin_ == text_.length())
    {
        end_       = begin_;
        breakSize_ = 0;
        nextBegin_ = begin_;
        return false;
    }
    if (literalIndent_ >= 0)
    {
        type_ = eParagraphType_Literal;
    }
    // Loop over lines in input until the end of the current paragraph.
    size_t i         = begin_;
    int    lineCount = 0;
    while (true)
    {
        const bool   bFirstLine = (lineCount == 0);
        const size_t lineStart  = i;
        const size_t lineEnd    = std::min(text_.find('\n', i), text_.length());
        const int    lineIndent = countLeadingSpace(text_, lineStart, lineEnd);
        const size_t textStart  = lineStart + lineIndent;
        const bool   bListItem  = startsListItem(text_, textStart);
        // Return each list item as a separate paragraph to make the behavior
        // the same always; the item text could even contain multiple
        // paragraphs, that would anyways produce breaks.
        if (bListItem && !bFirstLine)
        {
            // Since there was no empty line in input, do not produce one in
            // the output, either.
            nextBreakSize_ = 1;
            // end_ is not updated to break the paragraph before the current line.
            break;
        }
        // Now we will actually use this line as part of this paragraph.
        end_ = lineEnd;
        ++lineCount;
        // Update indentation.
        if (bFirstLine)
        {
            firstLineIndent_ = indent_ = lineIndent;
            if (bListItem)
            {
                // Find the indentation of the actual text after the
                // bullet/number.
                int prefixLength = 0;
                while (!std::isspace(text_[textStart + prefixLength]))
                {
                    ++prefixLength;
                }
                while (textStart + prefixLength < text_.length()
                       && std::isspace(text_[textStart + prefixLength]))
                {
                    ++prefixLength;
                }
                indent_ += prefixLength;
            }
        }
        else
        {
            indent_ = std::min(indent_, lineIndent);
        }
        // We need to check for the title underline before checking for the
        // paragraph break so that the title is correctly recognized.
        if (lineCount == 2 && isTitleUnderline(text_, lineStart))
        {
            type_ = eParagraphType_Title;
        }
        // Check for end-of-input or an empty line, i.e., a normal paragraph
        // break.
        if (lineEnd + 1 >= text_.length() || text_[lineEnd + 1] == '\n')
        {
            nextBreakSize_ = 2;
            break;
        }
        // Always return the title as a separate paragraph, as it requires
        // different processing.
        // TODO: This should allow nicer formatting that shares
        // implementation with writeTitle() and honors the nesting depths etc.,
        // but that is not implemented.
        if (type_ == eParagraphType_Title)
        {
            // If we are here, there was no actual paragraph break, so do not
            // produce one in the output either.
            nextBreakSize_ = 1;
            break;
        }
        // Next loop starts at the character after the newline.
        i = lineEnd + 1;
    }
    nextBegin_ = end_;
    // Check if the next paragraph should be treated as a literal paragraph,
    // and deal with transformations for the :: marker.
    if (end_ - begin_ >= 2 && text_.compare(end_ - 2, 2, "::") == 0)
    {
        literalIndent_ = indent_;
        // Return the actual literal block if the paragraph was just an "::".
        if (end_ - begin_ == 2)
        {
            // Avoid leading whitespace at the beginning; breakSize_ == 0
            // only for the first paragraph.
            if (breakSize_ == 0)
            {
                nextBreakSize_ = 0;
            }
            return nextParagraph();
        }
        // Remove one of the colons, or both if preceded by whitespace.
        const bool bRemoveDoubleColon = (text_[end_ - 3] == ' ');
        end_ -= (bRemoveDoubleColon ? 3 : 1);
    }
    else
    {
        literalIndent_ = -1;
    }
    // Treat a table like a literal block (preserve newlines).
    if (startsTable(text_, begin_ + firstLineIndent_))
    {
        type_ = eParagraphType_Literal;
    }
    return true;
}

void RstParagraphIterator::getParagraphText(std::string *result) const
{
    result->clear();
    result->reserve(end_ - begin_);
    result->append(breakSize_, '\n');
    const bool bPreserveNewlines = (type_ != eParagraphType_Normal);
    size_t     i                 = begin_;
    while (i < end_)
    {
        const bool   bFirstLine = (i == begin_);
        const size_t lineStart  = i + (bFirstLine ? firstLineIndent_ : indent_);
        const size_t lineEnd    = std::min(text_.find('\n', i), end_);
        if (!bFirstLine)
        {
            if (bPreserveNewlines)
            {
                result->push_back('\n');
            }
            else if (!std::isspace((*result)[result->length() - 1]))
            {
                result->push_back(' ');
            }
        }
        result->append(text_, lineStart, lineEnd - lineStart);
        i = lineEnd + 1;
    }
}

} // namespace gmx
