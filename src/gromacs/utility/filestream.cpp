/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017,2018,2019, by the GROMACS development team, led by
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
 * Implements classes from filestream.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "filestream.h"

#include "config.h"

#include <cerrno>
#include <cstdio>

#ifdef HAVE_UNISTD_H
#    include <unistd.h>
#endif

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Helper function for implementing readLine() for input streams.
bool readLineImpl(FILE* fp, std::string* line)
{
    line->clear();
    const size_t bufsize = 256;
    std::string  result;
    char         buf[bufsize];
    buf[0] = '\0';
    while (std::fgets(buf, bufsize, fp) != nullptr)
    {
        const size_t length = std::strlen(buf);
        result.append(buf, length);
        if (length < bufsize - 1 || buf[length - 1] == '\n')
        {
            break;
        }
    }
    if (std::ferror(fp))
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Error while reading file"), "fgets", errno);
    }
    *line = result;
    return !result.empty() || (std::feof(fp) == 0);
}

} // namespace

namespace internal
{

/********************************************************************
 * FileStreamImpl
 */

class FileStreamImpl
{
public:
    explicit FileStreamImpl(FILE* fp) : fp_(fp), bClose_(false) {}
    FileStreamImpl(const char* filename, const char* mode) : fp_(nullptr), bClose_(true)
    {
        fp_ = std::fopen(filename, mode);
        if (fp_ == nullptr)
        {
            GMX_THROW_WITH_ERRNO(FileIOError(formatString("Could not open file '%s'", filename)),
                                 "fopen", errno);
        }
    }
    ~FileStreamImpl()
    {
        if (fp_ != nullptr && bClose_)
        {
            if (std::fclose(fp_) != 0)
            {
                // TODO: Log the error somewhere
            }
        }
    }

    FILE* handle()
    {
        GMX_RELEASE_ASSERT(fp_ != nullptr, "Attempted to access a file object that is not open");
        return fp_;
    }

    void close()
    {
        GMX_RELEASE_ASSERT(fp_ != nullptr, "Attempted to close a file object that is not open");
        GMX_RELEASE_ASSERT(bClose_, "Attempted to close a file object that should not be");
        const bool bOk = (std::fclose(fp_) == 0);
        fp_            = nullptr;
        if (!bOk)
        {
            GMX_THROW_WITH_ERRNO(FileIOError("Error while closing file"), "fclose", errno);
        }
    }

private:
    //! File handle for this object (NULL if the stream has been closed).
    FILE* fp_;
    //! Whether \p fp_ should be closed by this object.
    bool bClose_;
};

} // namespace internal

using internal::FileStreamImpl;

/********************************************************************
 * StandardInputStream
 */

bool StandardInputStream::isInteractive() const
{
#ifdef HAVE_UNISTD_H
    return isatty(fileno(stdin)) != 0;
#else
    return true;
#endif
}

bool StandardInputStream::readLine(std::string* line)
{
    return readLineImpl(stdin, line);
}

// static
StandardInputStream& StandardInputStream::instance()
{
    static StandardInputStream stdinObject;
    return stdinObject;
}

/********************************************************************
 * TextInputFile
 */

// static
FilePtr TextInputFile::openRawHandle(const char* filename)
{
    FilePtr fp(fopen(filename, "r"));
    if (fp == nullptr)
    {
        GMX_THROW_WITH_ERRNO(FileIOError(formatString("Could not open file '%s'", filename)),
                             "fopen", errno);
    }
    return fp;
}

// static
FilePtr TextInputFile::openRawHandle(const std::string& filename)
{
    return openRawHandle(filename.c_str());
}

TextInputFile::TextInputFile(const std::string& filename) :
    impl_(new FileStreamImpl(filename.c_str(), "r"))
{
}

TextInputFile::TextInputFile(FILE* fp) : impl_(new FileStreamImpl(fp)) {}

TextInputFile::~TextInputFile() {}

FILE* TextInputFile::handle()
{
    return impl_->handle();
}

bool TextInputFile::readLine(std::string* line)
{
    return readLineImpl(impl_->handle(), line);
}

void TextInputFile::close()
{
    impl_->close();
}

/********************************************************************
 * TextOutputFile
 */

TextOutputFile::TextOutputFile(const std::string& filename) :
    impl_(new FileStreamImpl(filename.c_str(), "w"))
{
}

TextOutputFile::TextOutputFile(FILE* fp) : impl_(new FileStreamImpl(fp)) {}

TextOutputFile::~TextOutputFile() {}

void TextOutputFile::write(const char* str)
{
    if (std::fprintf(impl_->handle(), "%s", str) < 0)
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Writing to file failed"), "fprintf", errno);
    }
}

void TextOutputFile::close()
{
    impl_->close();
}

// static
TextOutputFile& TextOutputFile::standardOutput()
{
    static TextOutputFile stdoutObject(stdout);
    return stdoutObject;
}

// static
TextOutputFile& TextOutputFile::standardError()
{
    static TextOutputFile stderrObject(stderr);
    return stderrObject;
}

} // namespace gmx
