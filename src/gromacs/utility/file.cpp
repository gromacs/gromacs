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
/*! \internal \file
 * \brief
 * Implements gmx::File.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "file.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \internal \brief
 * Private implementation class for File.
 *
 * \ingroup module_utility
 */
class File::Impl
{
    public:
        /*! \brief
         * Initialize a file object with the given handle.
         *
         * \param[in]  fp     %File handle to use (may be NULL).
         * \param[in]  bClose Whether this object should close its file handle.
         */
        Impl(FILE *fp, bool bClose);
        ~Impl();

        //! File handle for this object (may be NULL).
        FILE                   *fp_;
        /*! \brief
         * Whether \p fp_ should be closed by this object.
         *
         * Can be true if \p fp_ is NULL.
         */
        bool                    bClose_;
};

File::Impl::Impl(FILE *fp, bool bClose)
    : fp_(fp), bClose_(bClose)
{
}

File::Impl::~Impl()
{
    if (fp_ != NULL && bClose_)
    {
        if (fclose(fp_) != 0)
        {
            // TODO: Log the error somewhere
        }
    }
}

// static
FILE *File::openRawHandle(const char *filename, const char *mode)
{
    FILE *fp = fopen(filename, mode);
    if (fp == NULL)
    {
        GMX_THROW_WITH_ERRNO(
                FileIOError(formatString("Could not open file '%s'", filename)),
                "fopen", errno);
    }
    return fp;
}

// static
FILE *File::openRawHandle(const std::string &filename, const char *mode)
{
    return openRawHandle(filename.c_str(), mode);
}

File::File(const char *filename, const char *mode)
    : impl_(new Impl(NULL, true))
{
    open(filename, mode);
}

File::File(const std::string &filename, const char *mode)
    : impl_(new Impl(NULL, true))
{
    open(filename, mode);
}

File::File(const FileInitializer &initializer)
    : impl_(new Impl(NULL, true))
{
    open(initializer.filename_, initializer.mode_);
}

File::File(FILE *fp, bool bClose)
    : impl_(new Impl(fp, bClose))
{
}

File::~File()
{
}

void File::open(const char *filename, const char *mode)
{
    GMX_RELEASE_ASSERT(impl_->fp_ == NULL,
                       "Attempted to open the same file object twice");
    // TODO: Port all necessary functionality from gmx_ffopen() here.
    impl_->fp_ = openRawHandle(filename, mode);
}

void File::open(const std::string &filename, const char *mode)
{
    open(filename.c_str(), mode);
}

void File::close()
{
    GMX_RELEASE_ASSERT(impl_->fp_ != NULL,
                       "Attempted to close a file object that is not open");
    GMX_RELEASE_ASSERT(impl_->bClose_,
                       "Attempted to close a file object that should not be");
    bool bOk = (fclose(impl_->fp_) == 0);
    impl_->fp_ = NULL;
    if (!bOk)
    {
        GMX_THROW_WITH_ERRNO(
                FileIOError("Error while closing file"), "fclose", errno);
    }
}

bool File::isInteractive() const
{
    GMX_RELEASE_ASSERT(impl_->fp_ != NULL,
                       "Attempted to access a file object that is not open");
#ifdef HAVE_UNISTD_H
    return isatty(fileno(impl_->fp_));
#else
    return true;
#endif
}

FILE *File::handle()
{
    GMX_RELEASE_ASSERT(impl_->fp_ != NULL,
                       "Attempted to access a file object that is not open");
    return impl_->fp_;
}

void File::readBytes(void *buffer, size_t bytes)
{
    errno = 0;
    FILE  *fp = handle();
    // TODO: Retry based on errno or something else?
    size_t bytesRead = std::fread(buffer, 1, bytes, fp);
    if (bytesRead != bytes)
    {
        if (feof(fp))
        {
            GMX_THROW(FileIOError(
                              formatString("Premature end of file\n"
                                           "Attempted to read: %d bytes\n"
                                           "Successfully read: %d bytes",
                                           static_cast<int>(bytes),
                                           static_cast<int>(bytesRead))));
        }
        else
        {
            GMX_THROW_WITH_ERRNO(FileIOError("Error while reading file"),
                                 "fread", errno);
        }
    }
}

bool File::readLine(std::string *line)
{
    if (!readLineWithTrailingSpace(line))
    {
        return false;
    }
    size_t endPos = line->find_last_not_of(" \t\r\n");
    if (endPos != std::string::npos)
    {
        line->resize(endPos + 1);
    }
    return true;
}

bool File::readLineWithTrailingSpace(std::string *line)
{
    line->clear();
    const size_t bufsize = 256;
    std::string  result;
    char         buf[bufsize];
    buf[0] = '\0';
    FILE        *fp = handle();
    while (fgets(buf, bufsize, fp) != NULL)
    {
        size_t length = std::strlen(buf);
        result.append(buf, length);
        if (length < bufsize - 1 || buf[length - 1] == '\n')
        {
            break;
        }
    }
    if (ferror(fp))
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Error while reading file"),
                             "fgets", errno);
    }
    *line = result;
    return !result.empty() || !feof(fp);
}

void File::writeString(const char *str)
{
    if (fprintf(handle(), "%s", str) < 0)
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Writing to file failed"),
                             "fprintf", errno);
    }
}

void File::writeLine(const char *line)
{
    size_t length = std::strlen(line);

    writeString(line);
    if (length == 0 || line[length-1] != '\n')
    {
        writeString("\n");
    }
}

void File::writeLine()
{
    writeString("\n");
}

// static
bool File::exists(const char *filename)
{
    if (filename == NULL)
    {
        return false;
    }
    FILE *test = fopen(filename, "r");
    if (test == NULL)
    {
        return false;
    }
    else
    {
        fclose(test);
        // Windows doesn't allow fopen of directory, so we don't need to check
        // this separately.
#ifndef GMX_NATIVE_WINDOWS
        struct stat st_buf;
        int         status = stat(filename, &st_buf);
        if (status != 0 || !S_ISREG(st_buf.st_mode))
        {
            return false;
        }
#endif
        return true;
    }
}

// static
bool File::exists(const std::string &filename)
{
    return exists(filename.c_str());
}

// static
File &File::standardInput()
{
    static File stdinObject(stdin, false);
    return stdinObject;
}

// static
File &File::standardOutput()
{
    static File stdoutObject(stdout, false);
    return stdoutObject;
}

// static
File &File::standardError()
{
    static File stderrObject(stderr, false);
    return stderrObject;
}

// static
std::string File::readToString(const char *filename)
{
    // Binary mode is required on Windows to be able to determine a size
    // that can be passed to fread().
    File  file(filename, "rb");
    FILE *fp = file.handle();

    if (std::fseek(fp, 0L, SEEK_END) != 0)
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Seeking to end of file failed"),
                             "fseek", errno);
    }
    long len = std::ftell(fp);
    if (len == -1)
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Reading file length failed"),
                             "ftell", errno);
    }
    if (std::fseek(fp, 0L, SEEK_SET) != 0)
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Seeking to start of file failed"),
                             "fseek", errno);
    }

    std::vector<char> data(len);
    file.readBytes(&data[0], len);
    file.close();

    std::string result(&data[0], len);
    // The below is necessary on Windows to make newlines stay as '\n' on a
    // roundtrip.
    result = replaceAll(result, "\r\n", "\n");

    return result;
}

// static
std::string File::readToString(const std::string &filename)
{
    return readToString(filename.c_str());
}

// static
void File::writeFileFromString(const std::string &filename,
                               const std::string &text)
{
    File file(filename, "w");
    file.writeString(text);
    file.close();
}

} // namespace gmx
