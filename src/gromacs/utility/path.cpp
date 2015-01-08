/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements functions in path.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "path.h"

#include "config.h"

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include <sys/stat.h>

#ifdef GMX_NATIVE_WINDOWS
#include <Windows.h>
#include <direct.h>
#else
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#endif

#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace
{

//! Directory separator to use when joining paths.
const char cDirSeparator = '/';
//! Directory separators to use when parsing paths.
const char cDirSeparators[] = "/\\";
/*! \var cPathSeparator
 * \brief
 * Separator to use to split the PATH environment variable.
 *
 * When reading the PATH environment variable, Unix separates entries
 * with colon, while windows uses semicolon.
 */
#ifdef GMX_NATIVE_WINDOWS
const char cPathSeparator = ';';
#else
const char cPathSeparator = ':';
#endif

//! Check whether a given character is a directory separator.
bool isDirSeparator(char chr)
{
    return std::strchr(cDirSeparators, chr);
}

} // namespace

namespace gmx
{

/********************************************************************
 * Path
 */

bool Path::containsDirectory(const std::string &path)
{
    return path.find_first_of(cDirSeparators) != std::string::npos;
}

/* Check if the program name begins with "/" on unix/cygwin, or
 * with "\" or "X:\" on windows. If not, the program name
 * is relative to the current directory.
 */
bool Path::isAbsolute(const char *path)
{
    if (isDirSeparator(path[0]))
    {
        return true;
    }
#ifdef GMX_NATIVE_WINDOWS
    return path[0] != '\0' && path[1] == ':' && isDirSeparator(path[2]);
#else
    return false;
#endif
}

bool Path::isAbsolute(const std::string &path)
{
    return isAbsolute(path.c_str());
}

#ifdef GMX_NATIVE_WINDOWS
namespace
{
struct handle_wrapper
{
    HANDLE handle;
    handle_wrapper(HANDLE h)
        : handle(h){}
    ~handle_wrapper()
    {
        if (handle != INVALID_HANDLE_VALUE)
        {
            ::CloseHandle(handle);
        }
    }
};
}
#endif

bool Path::isEquivalent(const std::string &path1, const std::string &path2)
{
    //based on boost_1_56_0/libs/filesystem/src/operations.cpp under BSL
#ifdef GMX_NATIVE_WINDOWS
    // Note well: Physical location on external media is part of the
    // equivalence criteria. If there are no open handles, physical location
    // can change due to defragmentation or other relocations. Thus handles
    // must be held open until location information for both paths has
    // been retrieved.

    // p2 is done first, so any error reported is for p1
    // FixME: #1635
    handle_wrapper h2(
            CreateFile(
                    path2.c_str(),
                    0,
                    FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                    0,
                    OPEN_EXISTING,
                    FILE_FLAG_BACKUP_SEMANTICS,
                    0));

    handle_wrapper h1(
            CreateFile(
                    path1.c_str(),
                    0,
                    FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                    0,
                    OPEN_EXISTING,
                    FILE_FLAG_BACKUP_SEMANTICS,
                    0));

    if (h1.handle == INVALID_HANDLE_VALUE
        || h2.handle == INVALID_HANDLE_VALUE)
    {
        // if one is invalid and the other isn't, then they aren't equivalent,
        // but if both are invalid then it is an error
        if (h1.handle == INVALID_HANDLE_VALUE
            && h2.handle == INVALID_HANDLE_VALUE)
        {
            GMX_THROW(FileIOError("Path::isEquivalent called with two invalid files"));
        }

        return false;
    }

    // at this point, both handles are known to be valid

    BY_HANDLE_FILE_INFORMATION info1, info2;

    if (!GetFileInformationByHandle(h1.handle, &info1))
    {
        GMX_THROW(FileIOError("Path::isEquivalent: GetFileInformationByHandle failed"));
    }

    if (!GetFileInformationByHandle(h2.handle, &info2))
    {
        GMX_THROW(FileIOError("Path::isEquivalent: GetFileInformationByHandle failed"));
    }

    // In theory, volume serial numbers are sufficient to distinguish between
    // devices, but in practice VSN's are sometimes duplicated, so last write
    // time and file size are also checked.
    return
        info1.dwVolumeSerialNumber == info2.dwVolumeSerialNumber
        && info1.nFileIndexHigh == info2.nFileIndexHigh
        && info1.nFileIndexLow == info2.nFileIndexLow
        && info1.nFileSizeHigh == info2.nFileSizeHigh
        && info1.nFileSizeLow == info2.nFileSizeLow
        && info1.ftLastWriteTime.dwLowDateTime
        == info2.ftLastWriteTime.dwLowDateTime
        && info1.ftLastWriteTime.dwHighDateTime
        == info2.ftLastWriteTime.dwHighDateTime;
#else
    struct stat s1, s2;
    int         e2 = stat(path2.c_str(), &s2);
    int         e1 = stat(path1.c_str(), &s1);

    if (e1 != 0 || e2 != 0)
    {
        // if one is invalid and the other isn't then they aren't equivalent,
        // but if both are invalid then it is an error.
        if (e1 != 0 && e2 != 0)
        {
            GMX_THROW_WITH_ERRNO(
                    FileIOError("Path::isEquivalent called with two invalid files"),
                    "stat", errno);
        }
        return false;
    }

    // both stats now known to be valid
    return s1.st_dev == s2.st_dev && s1.st_ino == s2.st_ino
           // According to the POSIX stat specs, "The st_ino and st_dev fields
           // taken together uniquely identify the file within the system."
           // Just to be sure, size and mod time are also checked.
           && s1.st_size == s2.st_size && s1.st_mtime == s2.st_mtime;
#endif
}

std::string Path::join(const std::string &path1,
                       const std::string &path2)
{
    // TODO: Remove extra separators if they are present in the input paths.
    return path1 + cDirSeparator + path2;
}


std::string Path::join(const std::string &path1,
                       const std::string &path2,
                       const std::string &path3)
{
    // TODO: Remove extra separators if they are present in the input paths.
    return path1 + cDirSeparator + path2 + cDirSeparator + path3;
}

std::string Path::getParentPath(const std::string &path)
{
    /* Expects that the path doesn't contain "." or "..". If used on a path for
     * which this isn't guaranteed realpath needs to be called first. */
    size_t pos = path.find_last_of(cDirSeparators);
    if (pos == std::string::npos)
    {
        return std::string();
    }
    return path.substr(0, pos);
}

std::string Path::getFilename(const std::string &path)
{
    size_t pos = path.find_last_of(cDirSeparators);
    if (pos == std::string::npos)
    {
        return path;
    }
    return path.substr(pos+1);
}

bool Path::hasExtension(const std::string &path)
{
    return getFilename(path).find('.') != std::string::npos;
}

std::string Path::stripExtension(const std::string &path)
{
    size_t dirSeparatorPos = path.find_last_of(cDirSeparators);
    size_t extPos          = path.find_last_of('.');
    if (extPos == std::string::npos
        || (dirSeparatorPos != std::string::npos && extPos < dirSeparatorPos))
    {
        return path;
    }
    return path.substr(0, extPos);
}

std::string Path::normalize(const std::string &path)
{
    std::string result(path);
    if (DIR_SEPARATOR != '/')
    {
        std::replace(result.begin(), result.end(), '/', DIR_SEPARATOR);
    }
    return result;
}

bool Path::exists(const char *path)
{
    return gmx_fexist(path);
}

bool Path::exists(const std::string &path)
{
    return exists(path.c_str());
}

std::string Path::getWorkingDirectory()
{
    // TODO: Use exceptions instead of gmx_fatal().
    char cwd[GMX_PATH_MAX];
    gmx_getcwd(cwd, sizeof(cwd));
    return cwd;
}

void Path::splitPathEnvironment(const std::string        &pathEnv,
                                std::vector<std::string> *result)
{
    size_t prevPos = 0;
    size_t separator;
    do
    {
        separator = pathEnv.find(cPathSeparator, prevPos);
        result->push_back(pathEnv.substr(prevPos, separator - prevPos));
        prevPos = separator + 1;
    }
    while (separator != std::string::npos);
}

std::vector<std::string> Path::getExecutablePaths()
{
    std::vector<std::string> result;
#ifdef GMX_NATIVE_WINDOWS
    // Add the local dir since it is not in the path on Windows.
    result.push_back("");
#endif
    const char *path = std::getenv("PATH");
    if (path != NULL)
    {
        splitPathEnvironment(path, &result);
    }
    return result;
}

std::string Path::resolveSymlinks(const std::string &path)
{
    /* Does not fully resolve the path like realpath/boost::canonical would.
     * It doesn't resolve path elements (including "." or ".."), but only
     * resolves the entire path (it does that recursively). */
    std::string result(path);
#ifndef GMX_NATIVE_WINDOWS
    char        buf[GMX_PATH_MAX];
    int         length;
    while ((length = readlink(result.c_str(), buf, sizeof(buf)-1)) > 0)
    {
        buf[length] = '\0';
        if (isAbsolute(buf))
        {
            result = buf;
        }
        else
        {
            result = join(getParentPath(result), buf);
        }
    }
#endif
    return result;
}


/********************************************************************
 * Directory
 */

int Directory::create(const char *path)
{
    if (Directory::exists(path))
    {
        return 0;
    }
#ifdef GMX_NATIVE_WINDOWS
    if (_mkdir(path))
#else
    if (mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH) != 0)
#endif
    {
        // TODO: Proper error handling.
        return -1;
    }
    return 0;
}


int Directory::create(const std::string &path)
{
    return create(path.c_str());
}


bool Directory::exists(const char *path)
{
    struct stat info;
    if (stat(path, &info) != 0)
    {
        if (errno != ENOENT && errno != ENOTDIR)
        {
            // TODO: Proper error handling.
        }
        return false;
    }
#ifdef GMX_NATIVE_WINDOWS
    return ((_S_IFDIR & info.st_mode) != 0);
#else
    return S_ISDIR(info.st_mode);
#endif
}


bool Directory::exists(const std::string &path)
{
    return exists(path.c_str());
}

} // namespace gmx
