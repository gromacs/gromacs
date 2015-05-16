/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::File.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FILE_H
#define GMX_UTILITY_FILE_H

#include <cstdio>

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class File;

/*! \brief
 * Parameters for creating a File object.
 *
 * This class (mostly) replaces the ability to return a File object from a
 * function (since File is not copyable): returning a FileInitializer instead
 * allows the caller to construct the File object.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class FileInitializer
{
    public:
        /*! \brief
         * Creates the initializer with given parameters.
         *
         * The passed strings must remain valid until the initializer is used
         * to construct a File object.
         */
        FileInitializer(const char *filename, const char *mode)
            : filename_(filename), mode_(mode)
        {
        }

    private:
        const char *filename_;
        const char *mode_;

        /*! \brief
         * Needed to allow access to the parameters without otherwise
         * unnecessary accessors.
         */
        friend class File;
};

/*! \brief
 * Basic file object.
 *
 * This class provides basic file I/O functionality and uses exceptions
 * (FileIOError) for error reporting.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class File
{
    public:
        /*! \brief
         * Opens a file and returns a `FILE` handle.
         *
         * \param[in] filename  Path of the file to open.
         * \param[in] mode      Mode to open the file in (for fopen()).
         * \throws    FileIOError on any I/O error.
         *
         * Instead of returning `NULL` on errors, throws an exception with
         * additional details (including the file name and `errno`).
         */
        static FILE *openRawHandle(const char *filename, const char *mode);
        //! \copydoc openRawHandle(const char *, const char *)
        static FILE *openRawHandle(const std::string &filename, const char *mode);
        /*! \brief
         * Creates a file object and opens a file.
         *
         * \param[in] filename  Path of the file to open.
         * \param[in] mode      Mode to open the file in (for fopen()).
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * \see open(const char *, const char *)
         */
        File(const char *filename, const char *mode);
        //! \copydoc File(const char *, const char *)
        File(const std::string &filename, const char *mode);
        /*! \brief
         * Creates a file object and opens a file.
         *
         * \param[in] initializer  Parameters to open the file.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        File(const FileInitializer &initializer);
        /*! \brief
         * Destroys the file object.
         *
         * If the file is still open, it is closed.
         * Any error conditions will be ignored.
         */
        ~File();

        /*! \brief
         * Opens a file.
         *
         * \param[in] filename  Path of the file to open.
         * \param[in] mode      Mode to open the file in (for fopen()).
         * \throws    FileIOError on any I/O error.
         *
         * The file object must not be open.
         */
        void open(const char *filename, const char *mode);
        //! \copydoc open(const char *, const char *)
        void open(const std::string &filename, const char *mode);
        /*! \brief
         * Closes the file object.
         *
         * \throws  FileIOError on any I/O error.
         *
         * The file must be open.
         */
        void close();

        /*! \brief
         * Returns whether the file is an interactive terminal.
         *
         * Only works on Unix, otherwise always returns true.
         * It only makes sense to call this for File::standardInput() and
         * friends.
         *
         * Thie file must be open.
         * Does not throw.
         */
        bool isInteractive() const;
        /*! \brief
         * Returns a file handle for interfacing with C functions.
         *
         * The file must be open.
         * Does not throw.
         */
        FILE *handle();

        /*! \brief
         * Reads given number of bytes from the file.
         *
         * \param[out] buffer  Pointer to buffer that receives the bytes.
         * \param[in]  bytes   Number of bytes to read.
         * \throws     FileIOError on any I/O error.
         *
         * The file must be open.
         */
        void readBytes(void *buffer, size_t bytes);
        /*! \brief
         * Reads a single line from the file.
         *
         * \param[out] line    String to receive the line.
         * \returns    false if nothing was read because the file ended.
         * \throws     std::bad_alloc if out of memory.
         * \throws     FileIOError on any I/O error.
         *
         * On error or when false is returned, \p line will be empty.
         * Trailing space will be removed from the line.
         * To loop over all lines in the file, use:
         * \code
           std::string line;
           while (file.readLine(&line))
           {
               // ...
           }
           \endcode
         */
        bool readLine(std::string *line);
        /*! \brief
         * Reads a single line from the file.
         *
         * \param[out] line    String to receive the line.
         * \returns    false if nothing was read because the file ended.
         * \throws     std::bad_alloc if out of memory.
         * \throws     FileIOError on any I/O error.
         *
         * On error or when false is returned, \p line will be empty.
         * Works as readLine(), except that terminating newline will be present
         * in \p line if it was present in the file.
         *
         * \see readLine()
         */
        bool readLineWithTrailingSpace(std::string *line);

        /*! \brief
         * Writes a string to the file.
         *
         * \param[in]  str  String to write.
         * \throws     FileIOError on any I/O error.
         *
         * The file must be open.
         */
        void writeString(const char *str);
        //! \copydoc writeString(const char *)
        void writeString(const std::string &str) { writeString(str.c_str()); }
        /*! \brief
         * Writes a line to the file.
         *
         * \param[in]  line  Line to write.
         * \throws     FileIOError on any I/O error.
         *
         * If \p line does not end in a newline, one newline is appended.
         * Otherwise, works as writeString().
         *
         * The file must be open.
         */
        void writeLine(const char *line);
        //! \copydoc writeLine(const char *)
        void writeLine(const std::string &line) { writeLine(line.c_str()); }
        /*! \brief
         * Writes a newline to the file.
         *
         * \throws     FileIOError on any I/O error.
         */
        void writeLine();

        /*! \brief
         * Checks whether a file exists and is a regular file.
         *
         * \param[in] filename  Path to the file to check.
         * \returns   true if \p filename exists and is accessible.
         *
         * Does not throw.
         */
        static bool exists(const char *filename);
        //! \copydoc exists(const char *)
        static bool exists(const std::string &filename);

        /*! \brief
         * Returns a File object for accessing stdin.
         *
         * \throws    std::bad_alloc if out of memory (only on first call).
         */
        static File &standardInput();
        /*! \brief
         * Returns a File object for accessing stdout.
         *
         * \throws    std::bad_alloc if out of memory (only on first call).
         */
        static File &standardOutput();
        /*! \brief
         * Returns a File object for accessing stderr.
         *
         * \throws    std::bad_alloc if out of memory (only on first call).
         */
        static File &standardError();

        /*! \brief
         * Reads contents of a file to a std::string.
         *
         * \param[in] filename  Name of the file to read.
         * \returns   The contents of \p filename.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        static std::string readToString(const char *filename);
        //! \copydoc readToString(const char *)
        static std::string readToString(const std::string &filename);
        /*! \brief
         * Convenience method for writing a file from a string in a single call.
         *
         * \param[in] filename  Name of the file to read.
         * \param[in] text      String to write to \p filename.
         * \throws    FileIOError on any I/O error.
         *
         * If \p filename exists, it is overwritten.
         */
        static void writeFileFromString(const std::string &filename,
                                        const std::string &text);

    private:
        /*! \brief
         * Initialize file object from an existing file handle.
         *
         * \param[in]  fp     %File handle to use (may be NULL).
         * \param[in]  bClose Whether this object should close its file handle.
         * \throws     std::bad_alloc if out of memory.
         *
         * Used internally to implement standardOutput() and standardError().
         */
        File(FILE *fp, bool bClose);

        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
