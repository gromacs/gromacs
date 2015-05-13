/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::DirectoryEnumerator.
 *
 * \author Erik Lindahl (original C implementation)
 * \author Teemu Murtola <teemu.murtola@gmail.com> (C++ wrapper + errno handling)
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "directoryenumerator.h"

#include "config.h"

#include <cerrno>
#include <cstdio>

#include <algorithm>
#include <string>
#include <vector>

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif
#ifdef GMX_NATIVE_WINDOWS
#include <io.h>
#endif

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * DirectoryEnumerator::Impl
 */

// TODO: Consider whether checking the return value of closing would be useful,
// and what could we do if it fails?
#if defined GMX_NATIVE_WINDOWS
// TODO: Consider if Windows provides more error details through other APIs.
class DirectoryEnumerator::Impl
{
    public:
        static Impl *init(const char *dirname, bool bThrow)
        {
            std::string tmpname(dirname);
            // Remove possible trailing directory separator.
            // TODO: Use a method in gmx::Path instead.
            if (tmpname.back() == '/' || tmpname.back() == '\\')
            {
                tmpname.pop_back();
            }

            // Add wildcard.
            tmpname.append("/*");

            errno = 0;
            _finddata_t finddata;
            intptr_t    handle = _findfirst(tmpname.c_str(), &finddata);
            if (handle < 0L)
            {
                if (errno != ENOENT && bThrow)
                {
                    const int         code    = errno;
                    const std::string message =
                        formatString("Failed to list files in directory '%s'",
                                     dirname);
                    GMX_THROW_WITH_ERRNO(FileIOError(message), "_findfirst", code);
                }
                return NULL;
            }
            return new Impl(handle, finddata);
        }
        Impl(intptr_t handle, _finddata_t finddata)
            : windows_handle(handle), finddata(finddata), bFirst_(true)
        {
        }
        ~Impl()
        {
            _findclose(windows_handle);
        }

        bool nextFile(std::string *filename)
        {
            if (bFirst_)
            {
                *filename = finddata.name;
                bFirst_   = false;
                return true;
            }
            else
            {
                errno = 0;
                if (_findnext(windows_handle, &finddata) != 0)
                {
                    if (errno == 0 || errno == ENOENT)
                    {
                        filename->clear();
                        return false;
                    }
                    else
                    {
                        GMX_THROW_WITH_ERRNO(
                                FileIOError("Failed to list files in a directory"),
                                "_findnext", errno);
                    }
                }
                *filename = finddata.name;
                return true;
            }
        }

    private:
        intptr_t     windows_handle;
        _finddata_t  finddata;
        bool         bFirst_;
};
#elif defined HAVE_DIRENT_H
class DirectoryEnumerator::Impl
{
    public:
        static Impl *init(const char *dirname, bool bThrow)
        {
            errno       = 0;
            DIR *handle = opendir(dirname);
            if (handle == NULL)
            {
                if (bThrow)
                {
                    const int         code    = errno;
                    const std::string message =
                        formatString("Failed to list files in directory '%s'",
                                     dirname);
                    GMX_THROW_WITH_ERRNO(FileIOError(message), "opendir", code);
                }
                return NULL;
            }
            return new Impl(handle);
        }
        explicit Impl(DIR *handle) : dirent_handle(handle)
        {
            // TODO: Use memory allocation that throws, and handle
            // exception safety (close handle) in such a case.
            /* On some platforms no space is present for d_name in dirent.
             * Since d_name is guaranteed to be the last entry, allocating
             * extra space for dirent will allow more size for d_name.
             * GMX_MAX_PATH should always be >= the max possible d_name.
             */
            smalloc(direntp_large, sizeof(*direntp_large) + GMX_PATH_MAX);
        }
        ~Impl()
        {
            sfree(direntp_large);
            closedir(dirent_handle);
        }

        bool nextFile(std::string *filename)
        {
            errno = 0;
            dirent *p;
            int     rc = readdir_r(dirent_handle, direntp_large, &p);
            if (p == NULL && rc == 0)
            {
                filename->clear();
                return false;
            }
            else if (rc != 0)
            {
                GMX_THROW_WITH_ERRNO(
                        FileIOError("Failed to list files in a directory"),
                        "readdir_r", errno);
            }
            *filename = direntp_large->d_name;
            return true;
        }

    private:
        DIR    *dirent_handle;
        dirent *direntp_large;
};
#else
class DirectoryEnumerator::Impl
{
    public:
        static Impl *init(const char * /*dirname*/, bool /*bThrow*/)
        {
            std::string message(
                    "Source compiled without POSIX dirent or Windows support "
                    "- cannot scan directories. In the very unlikely event "
                    "this is not a compile-time mistake you could consider "
                    "implementing support for your platform in "
                    "directoryenumerator.cpp, but contact the developers "
                    "to make sure it's really necessary!");
            GMX_THROW(NotImplementedError(message));
        }

        bool nextFile(std::string * /*filename*/)
        {
            return false;
        }
};
#endif

/********************************************************************
 * DirectoryEnumerator
 */

// static
std::vector<std::string>
DirectoryEnumerator::enumerateFilesWithExtension(
        const char *dirname, const char *extension, bool bThrow)
{
    std::vector<std::string> result;
    DirectoryEnumerator      dir(dirname, bThrow);
    std::string              nextName;
    while (dir.nextFile(&nextName))
    {
        if (debug)
        {
            std::fprintf(debug, "dir '%s' file '%s'\n",
                         dirname, nextName.c_str());
        }
        // TODO: What about case sensitivity?
        if (endsWith(nextName, extension))
        {
            result.push_back(nextName);
        }
    }

    std::sort(result.begin(), result.end());
    return result;
}


DirectoryEnumerator::DirectoryEnumerator(const char *dirname, bool bThrow)
    : impl_(NULL)
{
    GMX_RELEASE_ASSERT(dirname != NULL && dirname[0] != '\0',
                       "Attempted to open empty/null directory path");
    impl_.reset(Impl::init(dirname, bThrow));
}

DirectoryEnumerator::DirectoryEnumerator(const std::string &dirname, bool bThrow)
    : impl_(NULL)
{
    GMX_RELEASE_ASSERT(!dirname.empty(),
                       "Attempted to open empty/null directory path");
    impl_.reset(Impl::init(dirname.c_str(), bThrow));
}

DirectoryEnumerator::~DirectoryEnumerator()
{
}

bool DirectoryEnumerator::nextFile(std::string *filename)
{
    if (!impl_.get())
    {
        filename->clear();
        return false;
    }
    return impl_->nextFile(filename);
}

} // namespace gmx
