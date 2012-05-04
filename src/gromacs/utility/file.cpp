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
 *
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
/*! \internal \file
 * \brief
 * Implements gmx::File.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "file.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/format.h"

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
         * \param[in]  fp     File handle to use (may be NULL).
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
    // TODO: Port all necessary functionality from ffopen() here.
    impl_->fp_ = fopen(filename, mode);
    if (impl_->fp_ == NULL)
    {
        GMX_THROW_WITH_ERRNO(
                FileIOError(formatString("Could not open file '%s'", filename)),
                "fopen", errno);
    }
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

FILE *File::handle()
{
    GMX_RELEASE_ASSERT(impl_->fp_ != NULL,
                       "Attempted to access a file object that is not open");
    return impl_->fp_;
}

void File::readBytes(void *buffer, size_t bytes)
{
    errno = 0;
    FILE *fp = handle();
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
    File file(filename, "rb");
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

} // namespace gmx
