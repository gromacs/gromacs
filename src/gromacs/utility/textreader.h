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

#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
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
        static std::string readFileToString(const char *filename);
        //! \copydoc readFileToString(const char *)
        static std::string readFileToString(const std::string &filename);

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
        explicit TextReader(const std::string &filename);
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
         * boost::shared_ptr (e.g., if the stream is an object on the stack).
         */
        explicit TextReader(TextInputStream *stream);
        /*! \brief
         * Creates a reader that reads from specified stream.
         *
         * \param[in]  stream  Stream to read from.
         * \throws     std::bad_alloc if out of memory.
         *
         * The reader keeps a reference to the stream, so the caller can pass
         * in a temporary if necessary.
         */
        explicit TextReader(const TextInputStreamPointer &stream);
        ~TextReader();

        /*! \brief
         * Reads a single line (including newline) from the stream.
         *
         * \returns    `false` if nothing was read because the file ended.
         *
         * To loop over all lines in the stream, use:
         * \code
           TextReader reader("file.txt");
           while (reader.readLine())
           {
               const std::string &line = reader.currentLine();
               // ...
           }
           \endcode
         */
        bool readLine();
        /*! \brief
         * Returns from the stream the last line read (including newline)
         *
         * Newlines will be returned as part of \p line if it was present in
         * the stream.
         *
         * \see readLine()
         */
        const std::string &currentLine() const;
        /*! \brief
         * Returns from the stream the last line read (including newline)
         *
         * Works as currentLine(), except that a copy of the string
         * will be returned, where trailing whitespace has been
         * removed.
         *
         * \see readLine()
         */
        std::string currentLineTrimmed() const;
        /*! \brief
         * Reads all remaining lines from the stream as a single string.
         *
         * \returns   Full contents of the stream (from the current point to
         *     the end).
         */
        std::string readAll();

        /*! \brief Constructs an InvalidInputError with an error
         * message, and adds custom context information useful for
         * trouble-shooting issues when reading line-based input from
         * files.
         *
         * This will produce an error message that reads
         * \verbatim
           <message>
            on line <number>, which was
            '<the actual line text>'
           \endverbatim
         *
         * \param[in]  message        String describing to the user the nature of the error.
         * \throws     std::bad_alloc if out of memory.
         */
        InvalidInputError makeError(const std::string &message) const;

        /*! \brief
         * Closes the underlying stream.
         */
        void close();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
