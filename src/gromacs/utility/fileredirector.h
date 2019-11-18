/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2019, by the GROMACS development team, led by
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
 * Declares gmx::IFileOutputRedirector.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FILEREDIRECTOR_H
#define GMX_UTILITY_FILEREDIRECTOR_H

#include <string>

#include "gromacs/utility/path.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

/*! \libinternal \brief
 * Allows overriding file existence checks from code that supports it.
 *
 * The calling code should take in this interface and use the methods in it
 * all file system operations that need to support this redirection.
 *
 * This allows tests to override the file existence checks without actually
 * using the file system.  See IFileOutputRedirector for notes on
 * a typical usage pattern.
 *
 * With some further refactoring of the File class, this could also support
 * redirecting input files from in-memory buffers as well, but for now the
 * current capabilities are sufficient.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class IFileInputRedirector
{
public:
    virtual ~IFileInputRedirector();

    /*! \brief
     * Checks whether the provided path exists (and is a file).
     *
     * The \p onNotFound can be used to influence the behavior on error
     * conditions.  Functions to pass as this parameter are provided as
     * members of gmx::File.
     */
    virtual bool fileExists(const char* filename, File::NotFoundHandler onNotFound) const = 0;

    //! Convenience method to check file existence using an std::string path.
    bool fileExists(const std::string& filename, File::NotFoundHandler onNotFound) const
    {
        return fileExists(filename.c_str(), onNotFound);
    }
};

/*! \libinternal \brief
 * Allows capturing `stdout` and file output from code that supports it.
 *
 * The calling code should take in this interface and use the stream objects
 * it returns for all output that needs to support this redirection.
 *
 * Currently, the (nearly) only purpose for this interface is for unit tests to
 * capture the file output without duplicating the knowledge of which files are
 * actually produced.  The tests can also replace actual files with in-memory
 * streams (e.g., a StringOutputStream), and test the output without actually
 * accessing the file system and managing actual files.
 *
 * As the main user for non-default implementation of this interface is tests,
 * code using this interface generally uses a pattern where the redirector is
 * initialized to defaultFileOutputRedirector(), and a separate setter is
 * provided for tests to change the default.  This allows code outside the
 * tests (and outside the code actually calling the redirector) to be written
 * as if this interface did not exist (i.e., they do not need to pass the
 * default instance).
 *
 * Also, the interface only supports text files, but can be generalized if/when
 * there is a need for binary streams (see also TextOutputStream).
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class IFileOutputRedirector
{
public:
    virtual ~IFileOutputRedirector();

    /*! \brief
     * Returns a stream to use for `stdout` output.
     */
    virtual TextOutputStream& standardOutput() = 0;
    /*! \brief
     * Returns a stream to use for output to a file at a given path.
     *
     * \param[in] filename  Requested file name.
     */
    virtual TextOutputStreamPointer openTextOutputFile(const char* filename) = 0;

    //! Convenience method to open a stream using an std::string path.
    TextOutputStreamPointer openTextOutputFile(const std::string& filename)
    {
        return openTextOutputFile(filename.c_str());
    }
};

//! \cond libapi
/*! \brief
 * Returns default implementation for IFileInputRedirector.
 *
 * The returned implementation does not redirect anything, but just uses the
 * file system normally.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
IFileInputRedirector& defaultFileInputRedirector();
/*! \brief
 * Returns default implementation for IFileOutputRedirector.
 *
 * The returned implementation does not redirect anything, but just opens the
 * files at requested locations.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
IFileOutputRedirector& defaultFileOutputRedirector();
//! \endcond

} // namespace gmx

#endif
