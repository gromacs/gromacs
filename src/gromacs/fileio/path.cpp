/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
 * \ingroup module_fileio
 */
#include "path.h"

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "config.h"

#include <sys/stat.h>
#ifdef GMX_NATIVE_WINDOWS
#include <direct.h>
#else
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#endif

#include "gromacs/fileio/futil.h"
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

bool Path::startsWith(const std::string &path, const std::string &prefix)
{
    return gmx::startsWith(normalize(path), normalize(prefix));
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

std::string Path::normalize(const std::string &path)
{
    std::string result(path);
    // TODO: Remove . and .. entries.
    if (DIR_SEPARATOR != '/')
    {
        std::replace(result.begin(), result.end(), '/', DIR_SEPARATOR);
    }
#ifdef GMX_NATIVE_WINDOWS
    if (std::isalpha(result[0]) && result[1] == ':')
    {
        result[0] = std::toupper(result[0]);
    }
#endif
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
