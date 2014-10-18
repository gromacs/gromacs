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
/*! \libinternal \file
 * \brief
 * Declares common string formatting routines for online help.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPFORMAT_H
#define GMX_ONLINEHELP_HELPFORMAT_H

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class File;
class HelpWriterContext;

/*! \libinternal \brief
 * Formats rows of a table for text output.
 *
 * This utility class formats tabular data, mainly for console output.
 * Each row in the table can take multiple lines, and automatic text wrapping
 * is supported.  If text overflows the allocated width, the remaining columns
 * on that line become shifted.  To avoid this, it is possible to start the
 * output for different columns from different lines (it is the caller's
 * responsibility to check that overflows are avoided or are acceptable).
 *
 * Column definitions are first set with addColumn().
 * To format a row, first call clear().  Then call addColumnLine() to add text
 * to each column (can be called multiple times on a single column to add
 * multiple lines), and possibly setColumnFirstLineOffset() to adjust the line
 * from which the column output should start.  Finally, call formatRow() to
 * obtain the formatted row.
 *
 * A header (column titles and a horizontal line) is printed before the first
 * line.
 *
 * Typical usage:
 * \code
   gmx::TextTableFormatter formatter;
   formatter.addColumn("Name", 10, false);
   formatter.addColumn("Type", 10, false);
   formatter.addColumn("Description", 50, true);

   formatter.clear();
   formatter.addColumnLine(0, "name");
   formatter.addColumnLine(1, "type");
   formatter.addColumnLine(2, "Description for name");
   printf("%s", formatter.formatRow().c_str());

   formatter.clear();
   formatter.addColumnLine(0, "averylongname");
   formatter.addColumnLine(1, "type");
   formatter.setColumnFirstLineOffset(1, 1);
   formatter.addColumnLine(2, "Description for name");
   printf("%s", formatter.formatRow().c_str());

   // format other rows by repeating the above code
 * \endcode
 *
 * Methods in this class may throw std::bad_alloc if out of memory.
 * Other exceptions are not thrown.
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class TextTableFormatter
{
    public:
        //! Constructs an empty formatter.
        TextTableFormatter();
        ~TextTableFormatter();

        /*! \brief
         * Adds a column to the table.
         *
         * \param[in]  title  Title string for the column (used for header).
         * \param[in]  width  Width of the column (must be > 0).
         * \param[in]  bWrap  Whether text that exceeds \p width is
         *      automatically wrapped.
         *
         * The length of \p title must not exceed \p width.
         */
        void addColumn(const char *title, int width, bool bWrap);
        /*! \brief
         * Sets amount on indentation before the first column.
         *
         * \param[in]  indent  Number of spaces to use for indenting.
         *
         * Does not throw.
         */
        void setFirstColumnIndent(int indent);
        /*! \brief
         * Enables folding the last column to separate lines if it overflows.
         *
         * \param[in]  indent  Number of spaces to use for indenting the lines.
         *
         * If called with `indent >= 0`, the last column for each row is
         * treated specially: if it contains more lines than the other columns,
         * and if the text would fit more compactly as separate lines after the
         * row, then the whole last column is written after the row with the
         * given \p indent.  The column text then spans the whole space
         * reserved for the table, making long text fit into a smaller amount
         * of vertical space.
         * If not called, the last column is not treates specially.
         *
         * Does not throw.
         */
        void setFoldLastColumnToNextLine(int indent);

        /*! \brief
         * Whether formatRow() has been successfully called.
         *
         * This method can be used to determine after-the-fact whether anything
         * was written in the table.
         *
         * Does not throw.
         */
        bool didOutput() const;

        /*! \brief
         * Removes all text from all columns and resets the line offsets.
         *
         * Removes all text added using addColumnLine() and resets line offsets
         * set with setColumnFirstLineOffset() to zero.
         * Should be called before starting to add data for a row.
         *
         * Does not throw.
         */
        void clear();
        /*! \brief
         * Adds text to be printed in a column.
         *
         * \param[in]  index     Zero-based column index.
         * \param[in]  text      Text to add.
         *
         * Can be called multiple times.  Additional calls append \p text as
         * additional lines.  Any calls with \p text empty have no effect.
         * To add an empty line, use "\n" as \p text.
         *
         * If \p text contains newlines, the text is automatically splitted to
         * multiple lines.  The same happens if automatic wrapping is on for
         * the column and the text contains lines that are longer than what
         * fits the column.
         */
        void addColumnLine(int index, const std::string &text);
        /*! \brief
         * Adds text containing help markup to be printed in a column.
         *
         * \param[in]  index     Zero-based column index.
         * \param[in]  context   Context to use for markup processing.
         * \param[in]  text      Text to add.
         *
         * Works as addColumnLine(), except that it uses
         * HelpWriterContext::substituteMarkupAndWrapToVector() to process
         * markup in the input text instead of just wrapping it as plain text.
         */
        void addColumnHelpTextBlock(int index, const HelpWriterContext &context,
                                    const std::string &text);
        /*! \brief
         * Sets the first line to which text is printed for a column.
         *
         * \param[in]  index     Zero-based column index.
         * \param[in]  firstLine Zero-based line index from which to start the
         *      output.
         *
         * Can be called if there is no text for column \p index.
         * Does not affect the output in this case.
         *
         * Does not throw.
         */
        void setColumnFirstLineOffset(int index, int firstLine);
        /*! \brief
         * Formats the lines for the current row.
         *
         * \returns  Current row formatted as a single string
         *      (contains newlines).
         *
         * Formats the data as set after the previous clear()/formatRow() using
         * addColumnLine() and setColumnFirstLineOffset().
         *
         * If this is the first line to be formatted, a header is also added to
         * the beginning of the returned string if any column has a title.
         *
         * The return value always terminates with a newline.
         *
         * Calls clear() on successful return.
         */
        std::string formatRow();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
