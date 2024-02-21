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
 * Declares implementations for textstream.h interfaces for file input/output.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FILESTREAM_H
#define GMX_UTILITY_FILESTREAM_H

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

namespace internal
{
class FileStreamImpl;
}

/*! \libinternal \brief
 * Text input stream implementation for reading from `stdin`.
 *
 * Implementations for the TextInputStream methods throw FileIOError on any
 * I/O error.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class StandardInputStream : public TextInputStream
{
public:
    /*! \brief
     * Returns whether `stdin` is an interactive terminal.
     *
     * Only works on Unix, otherwise always returns true.
     *
     * Does not throw.
     */
    static bool isInteractive();

    // From TextInputStream
    bool readLine(std::string* line) override;
    void close() override {}
};

/*! \libinternal \brief
 * Text input stream implementation for reading from a file.
 *
 * Implementations for the TextInputStream methods throw FileIOError on any
 * I/O error.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class TextInputFile : public TextInputStream
{
public:
    /*! \brief
     * Opens a file and returns an RAII-style `FILE` handle.
     *
     * \param[in] filename  Path of the file to open.
     * \throws    FileIOError on any I/O error.
     *
     * Instead of returning `NULL` on errors, throws an exception with
     * additional details (including the file name and `errno`).
     */
    static FilePtr openRawHandle(const std::filesystem::path& filename);

    /*! \brief
     * Opens a text file as a stream.
     *
     * \param[in]  filename  Path to the file to open.
     * \throws     std::bad_alloc if out of memory.
     * \throws     FileIOError on any I/O error.
     */
    explicit TextInputFile(const std::filesystem::path& filename);
    /*! \brief
     * Initializes file object from an existing file handle.
     *
     * \param[in]  fp     File handle to use.
     * \throws     std::bad_alloc if out of memory.
     *
     * The caller is responsible of closing the file; close() does nothing
     * for an object constructed this way.
     */
    explicit TextInputFile(FILE* fp);
    ~TextInputFile() override;

    /*! \brief
     * Returns a raw handle to the input file.
     *
     * This is provided for interoperability with older C-like code.
     */
    FILE* handle();

    // From TextInputStream
    bool readLine(std::string* line) override;
    void close() override;

private:
    std::unique_ptr<internal::FileStreamImpl> impl_;
};

/*! \libinternal \brief
 * Text output stream implementation for writing to a file.
 *
 * Implementations for the TextOutputStream methods throw FileIOError on any
 * I/O error.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class TextOutputFile : public TextOutputStream
{
public:
    //! \copydoc TextInputFile::TextInputFile(const std::filesystem::path&)
    explicit TextOutputFile(const std::filesystem::path& filename);
    //! \copydoc TextInputFile::TextInputFile(FILE *)
    explicit TextOutputFile(FILE* fp);
    ~TextOutputFile() override;

    // From TextOutputStream
    void write(const char* text) override;
    void close() override;

    /*! \brief
     * Returns a stream for accessing `stdout`.
     *
     * \throws    std::bad_alloc if out of memory (only on first call).
     */
    static TextOutputFile& standardOutput();
    /*! \brief
     * Returns a stream for accessing `stderr`.
     *
     * \throws    std::bad_alloc if out of memory (only on first call).
     */
    static TextOutputFile& standardError();

private:
    std::unique_ptr<internal::FileStreamImpl> impl_;
};

} // namespace gmx

#endif
