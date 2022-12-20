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
/*! \libinternal \file
 * \brief
 * Declares functions for OS-independent path handling.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_PATH_H
#define GMX_UTILITY_PATH_H

#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace gmx
{

/*! \brief
 * Split PATH environment variable into search paths
 *
 * \param[in] pathEnv String to split.
 * \returns vector of filesystem paths to search.
 */
std::vector<std::filesystem::path> splitPathEnvironment(const std::string& pathEnv);
//! Get collection of possible executable paths.
std::vector<std::filesystem::path> getSystemExecutablePaths();
//! Strip source prefix from path.
std::string stripSourcePrefix(const char* path);
//! Concatenate before extension
std::filesystem::path concatenateBeforeExtension(const std::filesystem::path& path,
                                                 const std::string&           addition);
//! Remove extension from file path.
std::filesystem::path stripExtension(const std::filesystem::path& path);
//! Check if file extension of \p path without final '.' matches \p extension.
bool extensionMatches(const std::filesystem::path& path, std::string_view extension);

class File
{
public:
    struct NotFoundInfo
    {
        NotFoundInfo(const std::filesystem::path& filename,
                     const char*                  message,
                     const char*                  call,
                     bool                         wasError,
                     int                          err) :
            filename_(filename), message_(message), call_(call), wasError_(wasError), err_(err)
        {
        }

        const std::filesystem::path& filename_;
        const char*                  message_;
        const char*                  call_;
        bool                         wasError_;
        int                          err_;
    };

    static void              returnFalseOnError(const NotFoundInfo& info);
    static void              throwOnError(const NotFoundInfo& info);
    [[noreturn]] static void throwOnNotFound(const NotFoundInfo& info);

    typedef void (*NotFoundHandler)(const NotFoundInfo& info);

    /*! \brief
     * Checks whether a file exists and is a regular file.
     *
     * \param[in] filename    Path to the file to check.
     * \param[in] onNotFound  Function to call when the file does not
     *     exists or there is an error accessing it.
     * \returns   `true` if \p filename exists and is accessible.
     *
     * Does not throw, unless onNotFound throws.
     */
    static bool exists(const std::filesystem::path& filename, NotFoundHandler onNotFound);

private:
    // Disallow instantiation.
    File();
};

} // namespace gmx

#endif
