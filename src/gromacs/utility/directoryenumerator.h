/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares gmx::DirectoryEnumerator.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_DIRECTORYENUMERATOR_H
#define GMX_UTILITY_DIRECTORYENUMERATOR_H

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \libinternal \brief
 * Lists files in a directory.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class DirectoryEnumerator
{
    public:
        /*! \brief
         * Convenience function to list files with certain extension from a
         * directory.
         *
         * \param[in]  dirname   Path to the directory to list.
         * \param[in]  extension List files with the given extension
         *     (or suffix in file name).
         * \param[in]  bThrow    Whether failure to open the directory should throw.
         * \returns    List of files with the given extension in \p dirname.
         * \throws std::bad_alloc if out of memory.
         * \throws FileIOError if opening the directory fails and `bThrow == true`.
         * \throws FileIOError if some other I/O error occurs.
         */
        static std::vector<std::string> enumerateFilesWithExtension(
            const char *dirname, const char *extension, bool bThrow);

        /*! \brief
         * Opens a directory for listing.
         *
         * \param[in] dirname Path to the directory to list.
         * \param[in] bThrow  Whether failure to open the directory should throw.
         * \throws std::bad_alloc if out of memory.
         * \throws FileIOError if opening the directory fails and `bThrow == true`
         */
        explicit DirectoryEnumerator(const char *dirname, bool bThrow = true);
        //! \copydoc DirectoryEnumerator(const char *, bool)
        explicit DirectoryEnumerator(const std::string &dirname, bool bThrow = true);
        ~DirectoryEnumerator();

        /*! \brief
         * Gets next file in a directory.
         *
         * \param[out] filename  Name of the next file.
         * \returns `false` if there were no more files.
         * \throws  std::bad_alloc if out of memory.
         * \throws  FileIOError if listing the next file fails.
         *
         * If all files from the directory have been returned (or there are no
         * files in the directory and this is the first call), the method
         * returns `false` and \p filename is cleared.
         * Otherwise, the return value is `true` and the first/next file name
         * is returned in \p filename.
         * \p filename will not contain any path information, only the name of
         * the file.
         *
         * If `bThrow` passed to the constructor was `false` and the directory
         * was not successfully opened, the first call to this function will
         * return `false`.
         */
        bool nextFile(std::string *filename);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
