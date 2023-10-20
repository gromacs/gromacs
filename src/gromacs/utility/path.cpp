/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements functions in path.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/path.h"

#include "config.h"

#include <cctype>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sys/stat.h>

#if GMX_NATIVE_WINDOWS
#    include <Windows.h>
#    include <direct.h>
#else
#    ifdef HAVE_UNISTD_H
#        include <unistd.h>
#    endif
#endif

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace
{

//! Directory separators to use when parsing paths.
const char cDirSeparators[] = "/\\";
/*! \var cPathSeparator
 * \brief
 * Separator to use to split the PATH environment variable.
 *
 * When reading the PATH environment variable, Unix separates entries
 * with colon, while windows uses semicolon.
 */
#if GMX_NATIVE_WINDOWS
const char cPathSeparator = ';';
#else
const char cPathSeparator = ':';
#endif

//! Check whether a given character is a directory separator.
bool isDirSeparator(char chr)
{
    return std::strchr(cDirSeparators, chr) != nullptr;
}

} // namespace

namespace gmx
{

bool extensionMatches(const std::filesystem::path& path, const std::string_view extension)
{
    auto extensionWithSeparator = path.extension();
    return (!extensionWithSeparator.empty() && extensionWithSeparator.string().substr(1) == extension);
}

std::filesystem::path stripExtension(const std::filesystem::path& path)
{
    auto              parentPath = path.parent_path();
    const std::string stem       = path.stem().string();
    return parentPath.append(stem);
}

std::filesystem::path concatenateBeforeExtension(const std::filesystem::path& path, const std::string& addition)
{
    auto extension = path.extension();
    auto stripped  = path.parent_path().append(path.stem().string());
    return stripped.concat(addition).concat(extension.string());
}

std::string stripSourcePrefix(const char* path)
{
    const char* fallback           = path;
    const char* sep                = path + std::strlen(path);
    bool        gromacsSubdirFound = false;
    while (sep > path)
    {
        const char* prevSep = sep - 1;
        while (prevSep >= path && !isDirSeparator(*prevSep))
        {
            --prevSep;
        }
        const std::ptrdiff_t length = sep - prevSep - 1;
        if (gromacsSubdirFound)
        {
            if (std::strncmp(prevSep + 1, "src", length) == 0)
            {
                return prevSep + 1;
            }
            return fallback;
        }
        if (std::strncmp(prevSep + 1, "gromacs", length) == 0
            || std::strncmp(prevSep + 1, "programs", length) == 0
            || std::strncmp(prevSep + 1, "testutils", length) == 0)
        {
            gromacsSubdirFound = true;
        }
        if (fallback == path)
        {
            fallback = prevSep + 1;
        }
        sep = prevSep;
    }
    return fallback;
}

std::vector<std::filesystem::path> splitPathEnvironment(const std::string& pathEnv)
{
    std::vector<std::filesystem::path> result;
    size_t                             prevPos   = 0;
    size_t                             separator = 0;
    do
    {
        separator = pathEnv.find(cPathSeparator, prevPos);
        result.emplace_back(pathEnv.substr(prevPos, separator - prevPos));
        prevPos = separator + 1;
    } while (separator != std::string::npos);
    return result;
}

std::vector<std::filesystem::path> getSystemExecutablePaths()
{
    std::vector<std::filesystem::path> result;
#if GMX_NATIVE_WINDOWS
    // Add the local dir since it is not in the path on Windows.
    result.push_back("");
#endif
    const char* path = std::getenv("PATH");
    if (path != nullptr)
    {
        result = splitPathEnvironment(path);
    }
    return result;
}

/********************************************************************
 * File
 */

void File::returnFalseOnError(const NotFoundInfo& /*info*/) {}

void File::throwOnError(const NotFoundInfo& info)
{
    if (info.wasError_)
    {
        const std::string message = formatString(
                "Failed to access file '%s'.\n%s", info.filename_.string().c_str(), info.message_);
        GMX_THROW_WITH_ERRNO(FileIOError(message), info.call_, info.err_);
    }
}

void File::throwOnNotFound(const NotFoundInfo& info)
{
    throwOnError(info);
    const std::string message = formatString("File '%s' does not exist or is not accessible.\n%s",
                                             info.filename_.string().c_str(),
                                             info.message_);
    GMX_THROW_WITH_ERRNO(InvalidInputError(message), info.call_, info.err_);
}

// static
bool File::exists(const std::filesystem::path& filename, NotFoundHandler onNotFound)
{
    if (filename.empty())
    {
        return false;
    }
    FILE* test = std::fopen(filename.string().c_str(), "r");
    if (test == nullptr)
    {
        const bool   wasError = (errno != ENOENT && errno != ENOTDIR);
        NotFoundInfo info(filename, "The file could not be opened.", "fopen", wasError, errno);
        onNotFound(info);
        return false;
    }
    else
    {
        std::fclose(test);
        auto status = std::filesystem::status(filename);
        if (status.permissions() == std::filesystem::perms::unknown)
        {
            NotFoundInfo info(filename, "File information could not be read.", "permissions", true, 0);
            onNotFound(info);
            return false;
        }
        if (!std::filesystem::is_regular_file(filename))
        {
            NotFoundInfo info(filename, "The file is not a regular file.", "notafile", true, 0);
            onNotFound(info);
            return false;
        }
        return true;
    }
}

} // namespace gmx
