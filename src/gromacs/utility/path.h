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

#include <string>
#include <utility>
#include <vector>

namespace gmx
{

class Path
{
    public:
        static bool containsDirectory(const std::string &path);
        static bool isAbsolute(const char *path);
        static bool isAbsolute(const std::string &path);
        static bool isEquivalent(const std::string &path1,
                                 const std::string &path2);

        static std::string join(const std::string &path1,
                                const std::string &path2);
        static std::string join(const std::string &path1,
                                const std::string &path2,
                                const std::string &path3);
        static std::string normalize(const std::string &path);
        static std::string getParentPath(const std::string &path);
        static std::string getFilename(const std::string &path);
        static bool hasExtension(const std::string &path);
        static std::string stripExtension(const std::string &path);

        static bool exists(const char *path);
        static bool exists(const std::string &path);
        static std::string getWorkingDirectory();

        static void splitPathEnvironment(const std::string        &pathEnv,
                                         std::vector<std::string> *result);
        static std::vector<std::string> getExecutablePaths();

        static std::string resolveSymlinks(const std::string &path);

    private:
        // Disallow instantiation.
        Path();
};


class Directory
{
    public:
        static int create(const char *path);
        static int create(const std::string &path);
        static bool exists(const char *path);
        static bool exists(const std::string &path);

    private:
        // Disallow instantiation.
        Directory();
};

} // namespace gmx

#endif
