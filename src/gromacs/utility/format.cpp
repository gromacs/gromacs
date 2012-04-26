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
 *
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
 * Implements functions in format.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "gromacs/utility/format.h"

#include <cstdio>
#include <cstdarg>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

std::string formatString(const char *fmt, ...)
{
    va_list ap;
    char staticBuf[1024];
    int length = 1024;
    std::vector<char> dynamicBuf;
    char *buf = staticBuf;

    // TODO: There may be a better way of doing this on Windows, Microsoft
    // provides their own way of doing things...
    while (1)
    {
        va_start(ap, fmt);
        int n = vsnprintf(buf, length, fmt, ap);
        va_end(ap);
        if (n > -1 && n < length)
        {
            std::string result(buf);
            return result;
        }
        if (n > -1)
        {
            length = n + 1;
        }
        else
        {
            length *= 2;
        }
        dynamicBuf.resize(length);
        buf = &dynamicBuf[0];
    }
}

/********************************************************************
 * TextLineWrapper::Impl
 */

/*! \internal \brief
 * Private implementation class for TextLineWrapper.
 *
 * \ingroup module_utility
 */
class TextLineWrapper::Impl
{
    public:
        //! Initialize default values for the wrapper.
        Impl() : maxLength_(0) {}

        /*! \brief
         * Helper method to find the next wrapped line.
         *
         * \param[in]     input         Full input string.
         * \param[in,out] lineStartPtr
         *      Index of first character for the line to wrap.
         *      On exit, index of the first character of the next line.
         * \returns   Next output line.
         */
        std::string wrapNextLine(const std::string &input,
                                 size_t *lineStartPtr) const;

        //! Maximum length of output lines, or <= 0 if no limit.
        int                     maxLength_;
};

std::string
TextLineWrapper::Impl::wrapNextLine(const std::string &input,
                                    size_t *lineStartPtr) const
{
    size_t lineStart = *lineStartPtr;
    size_t lineEnd = std::string::npos;
    size_t nextNewline
        = std::min(input.find('\n', lineStart), input.length());
    if (maxLength_ <= 0 || nextNewline <= lineStart + maxLength_)
    {
        lineEnd = nextNewline;
    }
    else
    {
        size_t bestSpace = input.rfind(' ', lineStart + maxLength_);
        if (bestSpace < lineStart || bestSpace == std::string::npos)
        {
            bestSpace = input.find(' ', lineStart);
        }
        lineEnd = std::min(bestSpace, nextNewline);
    }
    size_t lineLength = lineEnd - lineStart;
    *lineStartPtr
        = (lineEnd < std::string::npos) ? lineEnd + 1 : std::string::npos;
    return input.substr(lineStart, lineLength);
}

/********************************************************************
 * TextLineWrapper
 */

TextLineWrapper::TextLineWrapper()
    : impl_(new Impl)
{
}

TextLineWrapper::~TextLineWrapper()
{
}

TextLineWrapper &TextLineWrapper::setLineLength(int length)
{
    impl_->maxLength_ = length;
    return *this;
}

std::string
TextLineWrapper::wrapToString(const std::string &input) const
{
    std::string result;
    size_t lineStart = 0;
    size_t length = input.length();
    while (lineStart < length)
    {
        result.append(impl_->wrapNextLine(input, &lineStart));
        if (lineStart < length
            || (lineStart == length && input[length - 1] == '\n'))
        {
            result.append("\n");
        }
    }
    return result;
}

std::vector<std::string>
TextLineWrapper::wrapToVector(const std::string &input) const
{
    std::vector<std::string> result;
    size_t lineStart = 0;
    while (lineStart < input.length())
    {
        result.push_back(impl_->wrapNextLine(input, &lineStart));
    }
    return result;
}

/********************************************************************
 * TextTableFormatter::Impl
 */

/*! \internal \brief
 * Private implementation class for TextTableFormatter.
 *
 * \ingroup module_utility
 */
class TextTableFormatter::Impl
{
    public:
        struct ColumnData
        {
            ColumnData(const char *title, int width, bool bWrap)
                : title_(title != NULL ? title : ""),
                  width_(width), bWrap_(bWrap), firstLine_(0)
            {
                GMX_ASSERT(width >= 0, "Negative width not possible");
                GMX_ASSERT(title_.length() <= static_cast<size_t>(width),
                           "Title too long for column width");
            }

            //! Returns the title of the column.
            const std::string &title() const { return title_; }
            //! Returns the width of the column.
            int width() const { return width_; }
            /*! \brief
             * Returns the first line offset for the current row.
             *
             * Note that the return value may be outside the printed lines if
             * there is no text.
             */
            int firstLine() const { return firstLine_; }
            /*! \brief
             * Returns the index of the last line with text for the current row.
             *
             * If there is no text, returns -1.
             */
            int lastLine() const
            {
                if (lines_.empty())
                {
                    return -1;
                }
                return firstLine_ + static_cast<int>(lines_.size()) - 1;
            }
            /*! \brief
             * Returns the text for a line.
             *
             * \param[in] line  Zero-based line index.
             * \returns   Text for line \p line, or empty string if \p line has
             *     no text for this column.
             */
            std::string textForLine(int line) const
            {
                // The second conditional matches if there are no lines
                if (line < firstLine() || line > lastLine())
                {
                    return std::string();
                }
                return lines_[line - firstLine()];
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
        //! If true, no output has yet been produced.
        bool                    bFirstRow_;
};

TextTableFormatter::Impl::Impl()
    : bFirstRow_(true)
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
    impl_->columns_.push_back(Impl::ColumnData(title, width, bWrap));
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
    Impl::ColumnData &column = impl_->columnData(index);
    TextLineWrapper wrapper;
    if (column.bWrap_)
    {
        wrapper.setLineLength(column.width());
    }
    std::vector<std::string> lines(wrapper.wrapToVector(text));
    column.lines_.insert(column.lines_.end(), lines.begin(), lines.end());
}

void TextTableFormatter::setColumnFirstLineOffset(int index, int firstLine)
{
    GMX_ASSERT(firstLine >= 0, "Invalid first line");
    Impl::ColumnData &column = impl_->columnData(index);
    column.firstLine_ = firstLine;
}

int TextTableFormatter::lastColumnLine(int index) const
{
    return impl_->columnData(index).lastLine();
}

std::string TextTableFormatter::formatRow()
{
    std::string result;
    Impl::ColumnList::const_iterator column;
    // Print a header if this is the first line.
    if (impl_->bFirstRow_)
    {
        size_t totalWidth = 0;
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
        result.append(totalWidth, '-');
        result.append("\n");
    }

    // Compute the last applicable line.
    int lastLine = -1;
    for (column  = impl_->columns_.begin();
         column != impl_->columns_.end();
         ++column)
    {
        lastLine = std::max(lastLine, column->lastLine());
    }

    // Format the actual row data.
    for (int line = 0; line <= lastLine; ++line)
    {
        std::string lineResult;
        size_t currentWidth = 0;
        for (column  = impl_->columns_.begin();
             column != impl_->columns_.end();
             ++column)
        {
            std::string value(column->textForLine(line));
            if (column != impl_->columns_.begin())
            {
                ++currentWidth;
                if (!value.empty())
                {
                    lineResult.append(" ");
                    if (lineResult.length() < currentWidth)
                    {
                        lineResult.resize(currentWidth, ' ');
                    }
                }
            }
            // TODO: Rewrap the text if wrapping is on and the previous columns
            // overflow.
            lineResult.append(value);
            currentWidth += column->width();
        }
        result.append(lineResult);
        result.append("\n");
    }

    impl_->bFirstRow_ = false;
    clear();
    return result;
}

} // namespace gmx
