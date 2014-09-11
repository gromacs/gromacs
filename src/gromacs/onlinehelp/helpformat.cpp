/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Implements functions in helpformat.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "helpformat.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * TextTableFormatter::Impl
 */

/*! \internal \brief
 * Private implementation class for TextTableFormatter.
 *
 * \ingroup module_onlinehelp
 */
class TextTableFormatter::Impl
{
    public:
        /*! \internal \brief
         * Manages a single column for TextTableFormatter.
         *
         * \ingroup module_onlinehelp
         */
        struct ColumnData
        {
            //! Initializes a text table column with given values.
            ColumnData(const char *title, int width, bool bWrap)
                : title_(title != NULL ? title : ""),
                  width_(width), bWrap_(bWrap), firstLine_(0),
                  nextLineIndex_(0), nextLineOffset_(0)
            {
                GMX_ASSERT(width >= 0, "Negative width not possible");
                GMX_ASSERT(title_.length() <= static_cast<size_t>(width),
                           "Title too long for column width");
            }

            //! Returns the title of the column.
            const std::string &title() const { return title_; }
            //! Returns the width of the column.
            int                width() const { return width_; }
            /*! \brief
             * Returns the first line offset for the current row.
             *
             * Note that the return value may be outside the printed lines if
             * there is no text.
             */
            int firstLine() const { return firstLine_; }

            /*! \brief
             * Resets the formatting state.
             *
             * After this call, textForNextLine() and hasLinesRemaining() can
             * be used to format the lines for the column.
             */
            void startFormatting()
            {
                nextLineIndex_  = (!lines_.empty() ? -firstLine_ : 0);
                nextLineOffset_ = 0;
            }
            //! Whether there are lines remaining for textForNextLine().
            bool hasLinesRemaining() const
            {
                return nextLineIndex_ < static_cast<int>(lines_.size());
            }
            /*! \brief
             * Returns the text for the next line.
             *
             * \param[in] columnWidth  Width to wrap the text to.
             * \returns   Text for the next line, or empty string if there is
             *   no text for this column.
             */
            std::string textForNextLine(int columnWidth)
            {
                if (nextLineIndex_ < 0 || !hasLinesRemaining())
                {
                    ++nextLineIndex_;
                    return std::string();
                }
                if (bWrap_)
                {
                    TextLineWrapperSettings  settings;
                    settings.setLineLength(columnWidth);
                    TextLineWrapper          wrapper(settings);
                    const std::string       &currentLine = lines_[nextLineIndex_];
                    const size_t             prevOffset  = nextLineOffset_;
                    const size_t             nextOffset
                        = wrapper.findNextLine(currentLine, prevOffset);
                    if (nextOffset >= currentLine.size())
                    {
                        ++nextLineIndex_;
                        nextLineOffset_ = 0;
                    }
                    else
                    {
                        nextLineOffset_ = nextOffset;
                    }
                    return wrapper.formatLine(currentLine, prevOffset, nextOffset);
                }
                else
                {
                    return lines_[nextLineIndex_++];
                }
            }

            //! Statit data: title of the column.
            std::string                 title_;
            //! Static data: width of the column.
            int                         width_;
            //! Static data: whether to automatically wrap input text.
            bool                        bWrap_;
            //! First line offset for the current row.
            int                         firstLine_;
            //! Text lines for the current row.
            std::vector<std::string>    lines_;
            //! Formatting state: index in `lines_` for the next line.
            int                         nextLineIndex_;
            //! Formatting state: offset within line `nextLineIndex_` for the next line.
            size_t                      nextLineOffset_;
        };

        //! Container type for column data.
        typedef std::vector<ColumnData> ColumnList;

        //! Initializes data for an empty formatter.
        Impl();

        /*! \brief
         * Convenience method for checked access to data for a column.
         *
         * \param[in] index  Zero-based column index.
         * \returns   \c columns_[index]
         */
        ColumnData &columnData(int index)
        {
            GMX_ASSERT(index >= 0 && index < static_cast<int>(columns_.size()),
                       "Invalid column index");
            return columns_[index];
        }
        //! \copydoc columnData()
        const ColumnData &columnData(int index) const
        {
            return const_cast<Impl *>(this)->columnData(index);
        }

        //! Container for column data.
        ColumnList              columns_;
        //! Indentation before the first column.
        int                     firstColumnIndent_;
        //! Indentation before the last column if folded.
        int                     foldLastColumnToNextLineIndent_;
        //! If true, no output has yet been produced.
        bool                    bFirstRow_;
        //! If true, a header will be printed before the first row.
        bool                    bPrintHeader_;
};

TextTableFormatter::Impl::Impl()
    : firstColumnIndent_(0), foldLastColumnToNextLineIndent_(-1),
      bFirstRow_(true), bPrintHeader_(false)
{
}

/********************************************************************
 * TextTableFormatter
 */

TextTableFormatter::TextTableFormatter()
    : impl_(new Impl)
{
}

TextTableFormatter::~TextTableFormatter()
{
}

void TextTableFormatter::addColumn(const char *title, int width, bool bWrap)
{
    if (title != NULL && title[0] != '\0')
    {
        impl_->bPrintHeader_ = true;
    }
    impl_->columns_.push_back(Impl::ColumnData(title, width, bWrap));
}

void TextTableFormatter::setFirstColumnIndent(int indent)
{
    GMX_RELEASE_ASSERT(indent >= 0, "Negative indentation not allowed");
    impl_->firstColumnIndent_ = indent;
}

void TextTableFormatter::setFoldLastColumnToNextLine(int indent)
{
    impl_->foldLastColumnToNextLineIndent_ = indent;
}

bool TextTableFormatter::didOutput() const
{
    return !impl_->bFirstRow_;
}

void TextTableFormatter::clear()
{
    Impl::ColumnList::iterator i;
    for (i = impl_->columns_.begin(); i != impl_->columns_.end(); ++i)
    {
        i->firstLine_ = 0;
        i->lines_.clear();
    }
}

void TextTableFormatter::addColumnLine(int index, const std::string &text)
{
    Impl::ColumnData        &column = impl_->columnData(index);
    TextLineWrapper          wrapper;
    std::vector<std::string> lines(wrapper.wrapToVector(text));
    column.lines_.insert(column.lines_.end(), lines.begin(), lines.end());
}

void TextTableFormatter::addColumnHelpTextBlock(
        int index, const HelpWriterContext &context, const std::string &text)
{
    Impl::ColumnData       &column = impl_->columnData(index);
    TextLineWrapperSettings settings;
    // TODO: If in the future, there is actually a coupling between the markup
    // and the wrapping, this must be postponed into formatRow(), where we do
    // the actual line wrapping.
    std::vector<std::string> lines(
            context.substituteMarkupAndWrapToVector(settings, text));
    column.lines_.insert(column.lines_.end(), lines.begin(), lines.end());
}

void TextTableFormatter::setColumnFirstLineOffset(int index, int firstLine)
{
    GMX_ASSERT(firstLine >= 0, "Invalid first line");
    Impl::ColumnData &column = impl_->columnData(index);
    column.firstLine_ = firstLine;
}

std::string TextTableFormatter::formatRow()
{
    std::string                result;
    Impl::ColumnList::iterator column;
    // Print a header if this is the first line.
    if (impl_->bPrintHeader_ && impl_->bFirstRow_)
    {
        size_t totalWidth = 0;
        result.append(impl_->firstColumnIndent_, ' ');
        for (column  = impl_->columns_.begin();
             column != impl_->columns_.end();
             ++column)
        {
            std::string title(column->title());
            if (column != impl_->columns_.end() - 1)
            {
                title.resize(column->width() + 1, ' ');
                totalWidth += title.length();
            }
            else
            {
                totalWidth += std::min(column->width(),
                                       static_cast<int>(title.length() + 13));
            }
            result.append(title);
        }
        result.append("\n");
        result.append(impl_->firstColumnIndent_, ' ');
        result.append(totalWidth, '-');
        result.append("\n");
    }

    // Format all the lines, one column at a time.
    std::vector<std::string> lines;
    std::vector<std::string> columnLines;
    int                      currentWidth    = 0;
    bool                     bFoldLastColumn = false;
    for (column  = impl_->columns_.begin();
         column != impl_->columns_.end();
         ++column)
    {
        // Format the column into columnLines.
        column->startFormatting();
        columnLines.clear();
        columnLines.reserve(lines.size());
        for (size_t line = 0; column->hasLinesRemaining(); ++line)
        {
            int columnWidth = column->width();
            if (line < lines.size())
            {
                const int overflow = static_cast<int>(lines[line].length()) - currentWidth;
                if (overflow > 0)
                {
                    if (overflow > columnWidth && column->bWrap_)
                    {
                        columnLines.push_back(std::string());
                        continue;
                    }
                    columnWidth -= overflow;
                }
            }
            columnLines.push_back(column->textForNextLine(columnWidth));
        }
        if (column == impl_->columns_.end() - 1
            && impl_->foldLastColumnToNextLineIndent_ >= 0
            && columnLines.size() >= lines.size() + column->lines_.size())
        {
            bFoldLastColumn = true;
            currentWidth   += column->width();
            break;
        }
        // Add columnLines into lines.
        if (lines.size() < columnLines.size())
        {
            lines.resize(columnLines.size());
        }
        for (size_t line = 0; line < columnLines.size(); ++line)
        {
            if (column != impl_->columns_.begin() && !columnLines[line].empty())
            {
                lines[line].append(" ");
                if (static_cast<int>(lines[line].length()) < currentWidth)
                {
                    lines[line].resize(currentWidth, ' ');
                }
            }
            lines[line].append(columnLines[line]);
        }
        currentWidth += column->width() + 1;
    }

    // Construct the result by concatenating all the lines.
    std::vector<std::string>::const_iterator line;
    for (line = lines.begin(); line != lines.end(); ++line)
    {
        result.append(impl_->firstColumnIndent_, ' ');
        result.append(*line);
        result.append("\n");
    }

    if (bFoldLastColumn)
    {
        Impl::ColumnList::reference lastColumn = impl_->columns_.back();
        const int                   totalIndent
            = impl_->firstColumnIndent_ + impl_->foldLastColumnToNextLineIndent_;
        lastColumn.startFormatting();
        currentWidth -= impl_->foldLastColumnToNextLineIndent_;
        while (lastColumn.hasLinesRemaining())
        {
            result.append(totalIndent, ' ');
            result.append(lastColumn.textForNextLine(currentWidth));
            result.append("\n");
        }
    }

    impl_->bFirstRow_ = false;
    clear();
    return result;
}

} // namespace gmx
