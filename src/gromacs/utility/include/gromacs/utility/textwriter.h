/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Declares gmx::TextWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_TEXTWRITER_H
#define GMX_UTILITY_TEXTWRITER_H

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/utility/textstream.h"

namespace gmx
{

class TextLineWrapperSettings;

/*! \libinternal \brief
 * Writes text into a TextOutputStream.
 *
 * This class provides more formatting and line-oriented writing capabilities
 * than writing raw strings into the stream.
 *
 * All methods that write to the stream can throw any exceptions that the
 * underlying stream throws.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class TextWriter
{
public:
    /*! \brief
     * Convenience method for writing a file from a string in a single call.
     *
     * \param[in] filename  Name of the file to read.
     * \param[in] text      String to write to \p filename.
     * \throws    std::bad_alloc if out of memory.
     * \throws    FileIOError on any I/O error.
     *
     * If \p filename exists, it is overwritten.
     */
    static void writeFileFromString(const std::filesystem::path& filename, const std::string& text);

    /*! \brief
     * Creates a writer that writes to specified file.
     *
     * \param[in]  filename  Path to the file to open.
     * \throws     std::bad_alloc if out of memory.
     * \throws     FileIOError on any I/O error.
     *
     * This constructor is provided for convenience for writing directly to
     * a file, without the need to construct multiple objects.
     */
    explicit TextWriter(const std::filesystem::path& filename);
    /*! \brief
     * Creates a writer that writes to specified file.
     *
     * \param[in]  fp  File handle to write to.
     * \throws     std::bad_alloc if out of memory.
     * \throws     FileIOError on any I/O error.
     *
     * This constructor is provided for interoperability with C-like code
     * for writing directly to an already opened file, without the need to
     * construct multiple objects.
     *
     * The caller is responsible of closing \p fp; it is not allowed to
     * call close() on the writer.
     */
    explicit TextWriter(FILE* fp);
    /*! \brief
     * Creates a writer that writes to specified stream.
     *
     * \param[in]  stream  Stream to write to.
     * \throws     std::bad_alloc if out of memory.
     *
     * The caller is responsible of the lifetime of the stream (should
     * remain in existence as long as the writer exists).
     *
     * This constructor is provided for convenience for cases where the
     * stream is not allocated with `new` and/or not managed by a
     * std::shared_ptr (e.g., if the stream is an object on the stack).
     */
    explicit TextWriter(TextOutputStream* stream);
    /*! \brief
     * Creates a writer that writes to specified stream.
     *
     * \param[in]  stream  Stream to write to.
     * \throws     std::bad_alloc if out of memory.
     *
     * The writer keeps a reference to the stream, so the caller can pass
     * in a temporary if necessary.
     */
    explicit TextWriter(const TextOutputStreamPointer& stream);
    ~TextWriter();

    /*! \brief
     * Allows adjusting wrapping settings for the writer.
     *
     * \todo
     * Wrapping is not currently implemented for code that writes partial
     * lines with writeString().
     */
    TextLineWrapperSettings& wrapperSettings();

    /*! \brief
     * Writes a string to the stream.
     *
     * \param[in]  str  String to write.
     */
    void writeString(const char* str);
    //! \copydoc writeString(const char *)
    void writeString(const std::string& str);
    //! Writes a string to the stream, with printf-style formatting.
    void writeStringFormatted(const char* fmt, ...);
    /*! \brief
     * Writes a line to the stream.
     *
     * \param[in]  line  Line to write.
     *
     * If \p line does not end in a newline, one newline is appended.
     * Otherwise, works as writeString().
     */
    void writeLine(const char* line);
    //! \copydoc writeLine(const char *)
    void writeLine(const std::string& line);
    //! Writes a line to the stream, with printf-style formatting.
    void writeLineFormatted(const char* fmt, ...);
    //! Writes a newline to the stream.
    void writeLine();

    /*! \brief
     * Writes a newline if previous output did not end in one.
     *
     * If nothing has been written using the writer, this method does
     * nothing.
     */
    void ensureLineBreak();
    /*! \brief
     * Ensures that the next string written starts after an empty line.
     *
     * Always terminates the current line (as with ensureLineBreak()), but
     * the empty line is only written out when the next line is written,
     * so that trailing newlines after final output can be avoided.
     *
     * If nothing has been written using the writer, this method does
     * nothing.
     */
    void ensureEmptyLine();

    /*! \brief
     * Closes the underlying stream.
     */
    void close();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
