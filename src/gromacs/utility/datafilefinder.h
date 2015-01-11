/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::DataFileFinder and related classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_LIBFILEFINDER_H
#define GMX_UTILITY_LIBFILEFINDER_H

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class DataFileFinder;

/*! \brief
 * Search parameters for DataFileFinder.
 *
 * This class implements a named parameter idiom for DataFileFinder::findFile()
 * and DataFileFinder::openFile() to support an easily readable and
 * customizable way of searching data files.
 *
 * By default, the search first considers the current directory, followed by
 * specified and default data directories, and an exception is thrown if the
 * file could not be found.
 * To skip searching in the current directory, use includeCurrentDir().
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class DataFileOptions
{
    public:
        /*! \brief
         * Constructs default options for searching for a file with the
         * specified name.
         *
         * \param[in] filename  File name to search for.
         *
         * This constructor is not explicit to allow passing a simple string to
         * DataFileFinder methods to search for the string with the default
         * parameters.
         */
        DataFileOptions(const char *filename)
            : filename_(filename), bCurrentDir_(true), bThrow_(true)
        {
        }
        //! \copydoc DataFileOptions(const char *)
        DataFileOptions(const std::string &filename)
            : filename_(filename.c_str()), bCurrentDir_(true), bThrow_(true)
        {
        }

        //! Sets whether to search in the current (working) directory.
        DataFileOptions &includeCurrentDir(bool bInclude)
        {
            bCurrentDir_ = bInclude;
            return *this;
        }
        //! Sets whether an exception is thrown if the file could not be found.
        DataFileOptions &throwIfNotFound(bool bThrow)
        {
            bThrow_ = bThrow;
            return *this;
        }

    private:
        const char *filename_;
        bool        bCurrentDir_;
        bool        bThrow_;

        /*! \brief
         * Needed to access the members without otherwise unnecessary accessors.
         */
        friend class DataFileFinder;
};

/*! \brief
 * Information about a data file found by DataFileFinder::enumerateFiles().
 *
 * \inpublicapi
 * \ingroup module_utility
 */
struct DataFileInfo
{
    //! Initializes the structure with given values.
    DataFileInfo(const std::string &dir, const std::string &name, bool bDefault)
        : dir(dir), name(name), bFromDefaultDir(bDefault)
    {
    }

    /*! \brief
     * Directory from which the file was found.
     *
     * If the file was found from the current directory, this will be `"."`.
     * In other cases, this will be a full path (except if the user-provided
     * search path contains relative paths).
     */
    std::string dir;
    /*! \brief
     * Name of the file without any directory name.
     */
    std::string name;
    /*! \brief
     * Whether the file was found from the default directory.
     *
     * If `true`, the file was found from the default installation data
     * directory, not from the current directory or any user-provided (through
     * DataFileFinder::setSearchPathFromEnv()) location.
     * \todo
     * Consider replacing with an enum that identifies the source (current dir,
     * GMXLIB, default).
     */
    bool        bFromDefaultDir;
};

/*! \brief
 * Searches data files from a set of paths.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class DataFileFinder
{
    public:
        /*! \brief
         * Constructs a default data file finder.
         *
         * The constructed finder searches only in the directory specified by
         * the global program context (see ProgramContextInterface), and
         * optionally in the current directory.
         *
         * Does not throw.
         */
        DataFileFinder();
        ~DataFileFinder();

        /*! \brief
         * Adds search path from an environment variable.
         *
         * \param[in] envVarName  Name of the environment variable to use.
         * \throws std::bad_alloc if out of memory.
         *
         * If the specified environment variable is set, it is interpreted like
         * a `PATH` environment variable on the platform (split at appropriate
         * separators), and each path found is added to the search path this
         * finder searches.  The added paths take precedence over the default
         * directory specified by the global program context, but the current
         * directory is searched first.
         */
        void setSearchPathFromEnv(const char *envVarName);

        /*! \brief
         * Opens a data file if found.
         *
         * \param[in] options  Identifies the file to be searched for.
         * \returns The opened file handle, or `NULL` if the file could not be
         *     found and exceptions were turned off.
         * \throws  std::bad_alloc if out of memory.
         * \throws  FileIOError if
         *   - no such file can be found, and \p options specifies that an
         *     exception should be thrown, or
         *   - there is an error opening the file (note that a file is skipped
         *     during the search if the user does not have rights to open the
         *     file at all).
         *
         * See findFile() for more details.
         */
        FILE *openFile(const DataFileOptions &options) const;
        /*! \brief
         * Finds a full path to a data file if found.
         *
         * \param[in] options  Identifies the file to be searched for.
         * \returns Full path to the data file, or an empty string if the file
         *     could not be found and exceptions were turned off.
         * \throws  std::bad_alloc if out of memory.
         * \throws  FileIOError if no such file can be found, and \p options
         *     specifies that an exception should be thrown.
         *
         * Searches for a data file in the search paths configured for the
         * finder, as well as in the current directory if so required.
         * Returns the full path to the first file found.
         */
        std::string findFile(const DataFileOptions &options) const;
        /*! \brief
         * Enumerates files in the data directories.
         *
         * \param[in] options  Idenfies files to be searched for.
         * \returns Information about each found file.
         * \throws  std::bad_alloc if out of memory.
         * \throws  FileIOError if no such file can be found, and \p options
         *     specifies that an exception should be thrown.
         *
         * Enumerates all files in the data directories that have the
         * extension/suffix specified by the file name in \p options.
         * Unlike findFile() and openFile(), this only works on files that are
         * in the actual data directories, not for any entry within
         * subdirectories of those.
         * See DataFileInfo for details on what is returned for each found
         * file.
         * Files from the same directory will be returned as a continuous block
         * in the returned vector.
         */
        std::vector<DataFileInfo> enumerateFiles(const DataFileOptions &options) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
