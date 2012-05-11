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
/*! \libinternal \file
 * \brief
 * Declares common string formatting routines for online help.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_ONLINEHELP_HELPFORMAT_H
#define GMX_ONLINEHELP_HELPFORMAT_H

#include <string>

#include "../utility/common.h"

namespace gmx
{

class File;

/*! \cond libapi */
/*! \libinternal \brief
 * Make the string uppercase.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all characters transformed to uppercase.
 * \throws    std::bad_alloc if out of memory.
 *
 * \inlibraryapi
 */
std::string toUpperCase(const std::string &text);

/*! \libinternal \brief
 * Substitute markup used in help text for console output.
 *
 * \param[in] text  Text to substitute.
 * \returns   \p text with markup substituted.
 * \throws    std::bad_alloc if out of memory.
 *
 * \inlibraryapi
 */
std::string substituteMarkupForConsole(const std::string &text);
/*! \libinternal \brief
 * Format a help text block for console output.
 *
 * \param     file  File to write the formatted text to.
 * \param[in] text  Text to format.
 * \throws    std::bad_alloc if out of memory.
 * \throws    FileIOError on any I/O error.
 *
 * Calls substituteMarkupForConsole(), and also wraps the lines to 78
 * characters.
 *
 * \inlibraryapi
 */
void writeHelpTextForConsole(File *file, const std::string &text);
//! \endcond

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
 * To format a fow, first call clear().  Then call addColumnLine() to add text
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
         * Does nto throw.
         */
        void setFirstColumnIndent(int indent);

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

        /*! \brief
         * Returns the last line on which column \p index has text.
         *
         * \param[in] index  Zero-based column index.
         * \returns   Last line index (zero-based) on which \p index has text.
         *
         * The return value is the sum of the number of lines added with
         * addColumnLine() (taking into account possible wrapping) and the line
         * offset set with setColumnFirstLineOffset().
         *
         * Does not throw.
         */
        int lastColumnLine(int index) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
