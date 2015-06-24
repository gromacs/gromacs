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
 * Declares gmx::FileOutputRedirectorInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FILEREDIRECTOR_H
#define GMX_UTILITY_FILEREDIRECTOR_H

#include <string>

#include "gromacs/utility/file.h"

namespace gmx
{

/*! \libinternal \brief
 * Allows overriding file existence checks from code that supports it.
 *
 * The calling code should take in this interface and use the methods in it
 * all file system operations that need to support this redirection.
 * By default, the code can then use defaultFileInputRedirector() in case no
 * redirection is needed.
 *
 * This allows tests to override the file existence checks without actually
 * using the file system.
 *
 * With some further refactoring of the File class, this could also support
 * redirecting input files from in-memory buffers as well, but for now the
 * current capabilities are sufficient.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class FileInputRedirectorInterface
{
    public:
        virtual ~FileInputRedirectorInterface();

        /*! \brief
         * Checks whether the provided path exists (and is a file).
         */
        virtual bool fileExists(const char *filename) const = 0;

        //! Convenience method to check file existence using an std::string path.
        bool fileExists(const std::string &filename) const
        {
            return fileExists(filename.c_str());
        }
};

/*! \libinternal \brief
 * Allows capturing `stdout` and file output from code that supports it.
 *
 * The calling code should take in this interface and use the File objects
 * it returns for all output that needs to support this redirection.
 * By default, the code can then use defaultFileOutputRedirector() in case no
 * redirection is needed.
 *
 * This allows tests to capture the file output without duplicating the
 * knowledge of which files are actually produced.  With some further
 * refactoring of the File class, this could support capturing the output into
 * in-memory buffers as well, but for now the current capabilities are
 * sufficient.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class FileOutputRedirectorInterface
{
    public:
        virtual ~FileOutputRedirectorInterface();

        /*! \brief
         * Returns a File object to use for `stdout` output.
         */
        virtual File &standardOutput() = 0;
        /*! \brief
         * Returns a File object to use for output to a given file.
         *
         * \param[in] filename  Requested file name.
         */
        virtual FileInitializer openFileForWriting(const char *filename) = 0;

        //! Convenience method to open a file using an std::string path.
        FileInitializer openFileForWriting(const std::string &filename)
        {
            return openFileForWriting(filename.c_str());
        }
};

//! \cond libapi
/*! \brief
 * Returns default implementation for FileInputRedirectorInterface.
 *
 * The returned implementation does not redirect anything, but just uses the
 * file system normally.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
FileInputRedirectorInterface &defaultFileInputRedirector();
/*! \brief
 * Returns default implementation for FileOutputRedirectorInterface.
 *
 * The returned implementation does not redirect anything, but just opens the
 * files at requested locations.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
FileOutputRedirectorInterface &defaultFileOutputRedirector();
//! \endcond

} // namespace gmx

#endif
