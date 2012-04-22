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

#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/format.h"

namespace gmx
{

File::File(const char *filename, const char *mode)
    : fp_(NULL)
{
    open(filename, mode);
}

File::File(const std::string &filename, const char *mode)
    : fp_(NULL)
{
    open(filename, mode);
}

File::~File()
{
    if (fp_ != NULL)
    {
        if (fclose(fp_) != 0)
        {
            // TODO: Log the error somewhere
        }
    }
}

void File::open(const char *filename, const char *mode)
{
    GMX_RELEASE_ASSERT(fp_ == NULL,
                       "Attempted to open the same file object twice");
    // TODO: Port all necessary functionality from ffopen() here.
    fp_ = fopen(filename, mode);
    if (fp_ == NULL)
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
    GMX_RELEASE_ASSERT(fp_ != NULL,
                       "Attempted to close a file object that is not open");
    if (fclose(fp_) != 0)
    {
        GMX_THROW_WITH_ERRNO(
                FileIOError("Error while closing file"), "fclose", errno);
    }
}

FILE *File::handle()
{
    GMX_RELEASE_ASSERT(fp_ != NULL,
                       "Attempted to access a file object that is not open");
    return fp_;
}

void File::readBytes(void *buffer, size_t bytes)
{
    GMX_RELEASE_ASSERT(fp_ != NULL,
                       "Attempted to access a file object that is not open");
    errno = 0;
    // TODO: Retry based on errno or something else?
    size_t bytesRead = std::fread(buffer, 1, bytes, fp_);
    if (bytesRead != bytes)
    {
        if (feof(fp_))
        {
            GMX_THROW(FileIOError("Premature end of file"));
        }
        else
        {
            GMX_THROW_WITH_ERRNO(FileIOError("Error while reading file"),
                                 "fread", errno);
        }
    }
}

// static
std::string File::readToString(const char *filename)
{
    File file(filename, "r");
    FILE *fp = file.handle();

    // TODO: Full error checking.
    std::fseek(fp, 0L, SEEK_END);
    long len = std::ftell(fp);
    std::fseek(fp, 0L, SEEK_SET);

    std::vector<char> data(len);
    file.readBytes(&data[0], len);
    std::string result(&data[0], len);

    file.close();
    return result;
}

// static
std::string File::readToString(const std::string &filename)
{
    return readToString(filename.c_str());
}

} // namespace gmx
