/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Implements gmx::DirectoryEnumerator.
 *
 * \author Erik Lindahl (original C implementation)
 * \author Teemu Murtola <teemu.murtola@gmail.com> (C++ wrapper + errno handling)
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/directoryenumerator.h"

#include "config.h"

#include <cerrno>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#if HAVE_DIRENT_H
#    include <dirent.h>
#endif
#if GMX_NATIVE_WINDOWS
#    include <io.h>
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
class DirectoryEnumerator::Impl
{
public:
    static Impl* init(const std::filesystem::path& dirname, bool bThrow)
    {
        if (!std::filesystem::is_directory(dirname)
            || (std::filesystem::is_symlink(dirname)
                && !std::filesystem::is_directory(std::filesystem::read_symlink(dirname))))
        {
            if (bThrow)
            {
                const int         code    = errno;
                const std::string message = formatString("Failed to list files in directory '%s'",
                                                         dirname.string().c_str());
                GMX_THROW_WITH_ERRNO(FileIOError(message), "opendir", code);
            }
            return nullptr;
        }
        return new Impl(dirname);
    }
    explicit Impl(const std::filesystem::path& path) : iter(path) {}
    ~Impl() {}

    std::optional<std::filesystem::path> nextFile()
    {
        if (iter == std::filesystem::directory_iterator())
        {
            return std::nullopt;
        }
        if (!std::filesystem::is_regular_file(*iter) && !std::filesystem::is_directory(*iter))
        {
            ++iter;
            return nextFile();
        }
        auto filename = iter->path().filename();
        ++iter;
        return std::optional(filename);
    }

private:
    std::filesystem::directory_iterator iter;
};

/********************************************************************
 * DirectoryEnumerator
 */

// static
std::vector<std::filesystem::path>
DirectoryEnumerator::enumerateFilesWithExtension(const std::filesystem::path& dirname,
                                                 const std::string&           extension,
                                                 bool                         bThrow)
{
    std::vector<std::filesystem::path> result;
    DirectoryEnumerator                dir(dirname, bThrow);
    auto                               nextName = dir.nextFile();
    while (nextName.has_value())
    {
        if (debug)
        {
            std::fprintf(debug, "dir '%s' file '%s'\n", dirname.string().c_str(), nextName->string().c_str());
        }
        // TODO: What about case sensitivity?
        if (endsWith(nextName.value().string(), extension))
        {
            result.emplace_back(nextName.value());
        }
        nextName = dir.nextFile();
    }

    std::sort(result.begin(), result.end());
    return result;
}


DirectoryEnumerator::DirectoryEnumerator(const std::filesystem::path& dirname, bool bThrow) :
    impl_(nullptr)
{
    GMX_RELEASE_ASSERT(!dirname.empty(), "Attempted to open empty/null directory path");
    impl_.reset(Impl::init(dirname, bThrow));
}

DirectoryEnumerator::~DirectoryEnumerator() {}

std::optional<std::filesystem::path> DirectoryEnumerator::nextFile()
{
    if (impl_ == nullptr)
    {
        return std::nullopt;
    }
    return impl_->nextFile();
}

} // namespace gmx
