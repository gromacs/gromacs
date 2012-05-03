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
/*! \file
 * \brief
 * Declares common string formatting routines.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FORMAT_H
#define GMX_UTILITY_FORMAT_H

#include <string>
#include <vector>

#include "common.h"

namespace gmx
{

/*! \brief
 * Format a string (snprintf() wrapper).
 *
 * \throws  std::bad_alloc if out of memory.
 *
 * This function works like sprintf(), except that it returns an std::string
 * instead of requiring a preallocated buffer.  Arbitrary length output is
 * supported.
 *
 * \inpublicapi
 */
std::string formatString(const char *fmt, ...);

/*! \brief
 * Joins strings in an array to a single string.
 *
 * \param[in] sarray  Array of strings to concatenate.
 * \param[in] count   Number of elements in \p sarray to concatenate.
 * \returns   All strings in \p sarray joined, ensuring at least one space
 *      between the strings.
 * \throws    std::bad_alloc if out of memory.
 *
 * The strings in the \p sarray array are concatenated, adding a single space
 * between the strings if there is no whitespace in the end of a string.
 * Terminal whitespace is removed.
 *
 * \inpublicapi
 */
std::string concatenateStrings(const char * const *sarray, size_t count);
/*! \brief
 * Convenience overload for joining strings in a C array (static data).
 *
 * \param[in] sarray  Array of strings to concatenate.
 * \tparam    count   Deduced number of elements in \p sarray.
 * \returns   All strings in \p sarray joined, ensuring at least one space
 *      between the strings.
 * \throws    std::bad_alloc if out of memory.
 *
 * \see concatenateStrings(const char * const *, size_t)
 *
 * \inpublicapi
 */
template <size_t count>
std::string concatenateStrings(const char * const (&sarray)[count])
{
    return concatenateStrings(sarray, count);
}

/*! \brief
 * Replace all occurrences of a string with another string.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \returns   \p input with all occurrences of \p from replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 *
 * The replacement is greedy and not recursive: starting from the beginning of
 * \p input, each match of \p from is replaced with \p to, and the search for
 * the next match begins after the end of the previous match.
 *
 * Compexity is O(N), where N is length of output.
 *
 * \see replaceAllWords()
 *
 * \inpublicapi
 */
std::string replaceAll(const std::string &input,
                       const char *from, const char *to);
/*! \brief
 * Replace whole words with others.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \returns   \p input with all \p from words replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 *
 * Works as replaceAll(), but a match is only considered if it is delimited by
 * non-alphanumeric characters.
 *
 * \see replaceAll()
 *
 * \inpublicapi
 */
std::string replaceAllWords(const std::string &input,
                            const char *from, const char *to);

/*! \brief
 * Wraps lines to a predefined length.
 *
 * This utility class wraps lines at word breaks to produce lines that are not
 * longer than a predefined length.  Explicit newlines ('\\n') are preserved.
 * Only space is considered a word separator.  If a single word exceeds the
 * maximum line length, it is still printed on a single line.
 * Extra whitespace is stripped from the start and end of produced lines.
 * If maximum line length is not set using setLineLength(), only wraps at
 * explicit newlines.
 *
 * Two output formats are possible: wrapToString() produces a single string
 * with embedded newlines, and wrapToVector() produces a vector of strings,
 * where each element is one line.
 *
 * Typical usage:
 * \code
gmx::TextLineWrapper wrapper;
wrapper.setLineLength(78);
printf("%s\n", wrapper.wrapToString(textToWrap).c_str());
 * \endcode
 *
 * Methods in this class may throw std::bad_alloc if out of memory.
 * Other exceptions are not thrown.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class TextLineWrapper
{
    public:
        //! Constructs a new line wrapper with no initial wrapping length.
        TextLineWrapper();
        ~TextLineWrapper();

        /*! \brief
         * Sets the maximum length for output lines.
         *
         * \param[in] length  Maximum length for the lines after wrapping.
         * \returns   *this
         *
         * If this method is not called, the wrapper has no maximum length
         * (only wraps at explicit line breaks).
         *
         * Does not throw.
         */
        TextLineWrapper &setLineLength(int length);

        /*! \brief
         * Formats a string, producing a single string with all the lines.
         *
         * \param[in] input  String to wrap.
         * \returns   \p input with added newlines such that maximum line
         *      length is not exceeded.
         *
         * Newlines in the input are preserved, including terminal newlines.
         * Note that if the input does not contain a terminal newline, the
         * output does not either.
         */
        std::string wrapToString(const std::string &input) const;
        /*! \brief
         * Formats a string, producing a vector with all the lines.
         *
         * \param[in] input  String to wrap.
         * \returns   \p input split into lines such that maximum line length
         *      is not exceeded.
         *
         * The strings in the returned vector do not contain newlines at the
         * end.
         * Note that a single terminal newline does not affect the output:
         * "line\\n" and "line" both produce the same output (but "line\\n\\n"
         * produces two lines, the second of which is empty).
         */
        std::vector<std::string> wrapToVector(const std::string &input) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \brief
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
 * \inpublicapi
 * \ingroup module_utility
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
         * the beginning of the returned string.
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
