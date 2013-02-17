/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "path.h"

#include "gmx_header_config.h"

#include <errno.h>
#include <sys/stat.h>

#ifdef GMX_NATIVE_WINDOWS
#include <direct.h>
#endif

namespace
{

//! Directory separator to use when joining paths.
const char cDirSeparator = '/';
//! Directory separators to use when parsing paths.
const char cDirSeparators[] = "/\\";

} // namespace

namespace gmx
{

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

std::pair<std::string, std::string>
Path::splitToPathAndFilename(const std::string &path)
{
    size_t pos = path.find_last_of(cDirSeparators);
    if (pos == std::string::npos)
    {
        return std::make_pair(std::string(), path);
    }
    return std::make_pair(path.substr(0, pos), path.substr(pos+1));
}


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
