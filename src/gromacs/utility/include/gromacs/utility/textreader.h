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
 * Declares gmx::TextReader.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_TEXTREADER_H
#define GMX_UTILITY_TEXTREADER_H

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/utility/textstream.h"

namespace gmx
{

/*! \libinternal \brief
 * Reads text from a TextInputStream.
 *
 * This class provides more formatted reading capabilities than reading raw
 * lines from the stream (and a natural place to implement more such
 * capabilities).
 *
 * All methods that read from the stream can throw any exceptions that the
 * underlying stream throws.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class TextReader
{
public:
    /*! \brief
     * Reads contents of a file to a std::string.
     *
     * \param[in] filename  Name of the file to read.
     * \returns   The contents of \p filename.
     * \throws    std::bad_alloc if out of memory.
     * \throws    FileIOError on any I/O error.
     */
    static std::string readFileToString(const std::string& filename);
    //! \copydoc readFileToString(const std::string&)
    static std::string readFileToString(const std::filesystem::path& filename);

    /*! \brief
     * Creates a reader that reads from specified file.
     *
     * \param[in]  filename  Path to the file to open.
     * \throws     std::bad_alloc if out of memory.
     * \throws     FileIOError on any I/O error.
     *
     * This constructor is provided for convenience for reading directly
     * from a file, without the need to construct multiple objects.
     */
    explicit TextReader(const std::filesystem::path& filename);
    /*! \brief
     * Creates a reader that reads from specified stream.
     *
     * \param[in]  stream  Stream to read from.
     * \throws     std::bad_alloc if out of memory.
     *
     * The caller is responsible of the lifetime of the stream (should
     * remain in existence as long as the reader exists).
     *
     * This constructor is provided for convenience for cases where the
     * stream is not allocated with `new` and/or not managed by a
     * std::shared_ptr (e.g., if the stream is an object on the stack).
     */
    explicit TextReader(TextInputStream* stream);
    /*! \brief
     * Creates a reader that reads from specified stream.
     *
     * \param[in]  stream  Stream to read from.
     * \throws     std::bad_alloc if out of memory.
     *
     * The reader keeps a reference to the stream, so the caller can pass
     * in a temporary if necessary.
     */
    explicit TextReader(const TextInputStreamPointer& stream);
    ~TextReader();

    /*! \brief
     * Reads a single line (including newline) from the stream.
     *
     * \param[out] line    String to receive the line.
     * \returns    `false` if nothing was read because the file ended.
     *
     * On error or when false is returned, \p line will be empty.
     * Newlines will be returned as part of \p line if it was present in
     * the stream.
     * To loop over all lines in the stream, use:
     * \code
       std::string line;
       while (reader.readLine(&line))
       {
           // ...
       }
       \endcode
     *
     * Behaviours such as trimming whitespace or comments can be
     * configured by calling other methods before this one.
     */
    bool readLine(std::string* line);
    /*! \brief Sets whether the reader should trim leading whitespace
     * from a line before returning it.
     *
     * \param[in] doTrimming  Whether trimming should be active.
     */
    void setTrimLeadingWhiteSpace(bool doTrimming);
    /*! \brief Sets whether the reader should trim trailing whitespace
     * from a line before returning it.
     *
     * Note that comment trimming will precede whitespace trimming
     * when both are active.
     *
     * \param[in] doTrimming  Whether trimming should be active.
     */
    void setTrimTrailingWhiteSpace(bool doTrimming);
    /*! \brief Sets whether the reader should trim at trailing
     * comment from a line before returning it.
     *
     * Note that comment trimming will precede whitespace trimming
     * when both are active.
     *
     * \param[in]  commentChar  The character that begins a comment.
     *
     * \param[in] doTrimming  Whether trimming should be active.
     */
    void setTrimTrailingComment(bool doTrimming, char commentChar);
    /*! \brief
     * Reads all remaining lines from the stream as a single string.
     *
     * \returns   Full contents of the stream (from the current point to
     *     the end).
     */
    std::string readAll();

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
