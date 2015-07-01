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
 * Declares gmx::TextWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_TEXTWRITER_H
#define GMX_UTILITY_TEXTWRITER_H

#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

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
         * Creates a writer that writes to specified file.
         *
         * \param[in]  filename  Path to the file to open.
         * \throws     std::bad_alloc if out of memory.
         * \throws     FileIOError on any I/O error.
         *
         * This constructor is provided for convenience for writing directly to
         * a file, without the need to construct multiple objects.
         */
        explicit TextWriter(const std::string &filename);
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
         * boost::shared_ptr (e.g., if the stream is an object on the stack).
         */
        explicit TextWriter(TextOutputStream *stream);
        /*! \brief
         * Creates a writer that writes to specified stream.
         *
         * \param[in]  stream  Stream to write to.
         * \throws     std::bad_alloc if out of memory.
         *
         * The writer keeps a reference to the stream, so the caller can pass
         * in a temporary if necessary.
         */
        explicit TextWriter(const TextOutputStreamPointer &stream);
        ~TextWriter();

        //! Returns the underlying stream for this writer.
        TextOutputStream &stream();

        /*! \brief
         * Writes a string to the stream.
         *
         * \param[in]  str  String to write.
         */
        void writeString(const char *str);
        //! \copydoc writeString(const char *)
        void writeString(const std::string &str);
        /*! \brief
         * Writes a line to the stream.
         *
         * \param[in]  line  Line to write.
         *
         * If \p line does not end in a newline, one newline is appended.
         * Otherwise, works as writeString().
         */
        void writeLine(const char *line);
        //! \copydoc writeLine(const char *)
        void writeLine(const std::string &line);
        //! Writes a newline to the stream.
        void writeLine();

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
