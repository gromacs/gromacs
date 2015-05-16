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
/*! \internal \file
 * \brief
 * Implements gmx::DataFileFinder.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "datafilefinder.h"

#include <cstdlib>

#include <string>
#include <vector>

#include "buildinfo.h"
#include "gromacs/utility/directoryenumerator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * DataFileFinder::Impl
 */

class DataFileFinder::Impl
{
    public:
        static std::string getDefaultPath();

        Impl() : envName_(NULL), bEnvIsSet_(false) {}

        const char               *envName_;
        bool                      bEnvIsSet_;
        std::vector<std::string>  searchPath_;
};

std::string DataFileFinder::Impl::getDefaultPath()
{
    const InstallationPrefixInfo installPrefix
        = getProgramContext().installationPrefix();
    if (!isNullOrEmpty(installPrefix.path))
    {
        const char *const dataPath
            = installPrefix.bSourceLayout ? "share" : DATA_INSTALL_DIR;
        return Path::join(installPrefix.path, dataPath, "top");
    }
    return std::string();
}

/********************************************************************
 * DataFileFinder
 */

DataFileFinder::DataFileFinder()
    : impl_(NULL)
{
}

DataFileFinder::~DataFileFinder()
{
}

void DataFileFinder::setSearchPathFromEnv(const char *envVarName)
{
    if (!impl_.get())
    {
        impl_.reset(new Impl());
    }
    impl_->envName_ = envVarName;
    const char *const lib = getenv(envVarName);
    if (lib != NULL)
    {
        impl_->bEnvIsSet_ = true;
        Path::splitPathEnvironment(lib, &impl_->searchPath_);
    }
}

FILE *DataFileFinder::openFile(const DataFileOptions &options) const
{
    // TODO: There is a small race here, since there is some time between
    // the exists() calls and actually opening the file.  It would be better
    // to leave the file open after a successful exists() if the desire is to
    // actually open the file.
    std::string filename = findFile(options);
    if (filename.empty())
    {
        return NULL;
    }
#if 0
    if (debug)
    {
        fprintf(debug, "Opening library file %s\n", fn);
    }
#endif
    return File::openRawHandle(filename, "r");
}

std::string DataFileFinder::findFile(const DataFileOptions &options) const
{
    if (options.bCurrentDir_ && Path::exists(options.filename_))
    {
        return options.filename_;
    }
    if (impl_.get())
    {
        std::vector<std::string>::const_iterator i;
        for (i = impl_->searchPath_.begin(); i != impl_->searchPath_.end(); ++i)
        {
            // TODO: Deal with an empty search path entry more reasonably.
            std::string testPath = Path::join(*i, options.filename_);
            // TODO: Consider skipping directories.
            if (Path::exists(testPath))
            {
                return testPath;
            }
        }
    }
    const std::string &defaultPath = Impl::getDefaultPath();
    if (!defaultPath.empty())
    {
        std::string testPath = Path::join(defaultPath, options.filename_);
        if (Path::exists(testPath))
        {
            return testPath;
        }
    }
    if (options.bThrow_)
    {
        const char *const envName   = (impl_.get() ? impl_->envName_ : NULL);
        const bool        bEnvIsSet = (impl_.get() ? impl_->bEnvIsSet_ : false);
        std::string       message(
                formatString("Library file '%s' not found", options.filename_));
        if (options.bCurrentDir_)
        {
            message.append(" in current dir nor");
        }
        if (bEnvIsSet)
        {
            message.append(formatString(" in your %s path nor", envName));
        }
        message.append(" in the default directories.\nThe following paths were searched:");
        if (options.bCurrentDir_)
        {
            message.append("\n  ");
            message.append(Path::getWorkingDirectory());
            message.append(" (current dir)");
        }
        if (impl_.get())
        {
            std::vector<std::string>::const_iterator i;
            for (i = impl_->searchPath_.begin(); i != impl_->searchPath_.end(); ++i)
            {
                message.append("\n  ");
                message.append(*i);
            }
        }
        if (!defaultPath.empty())
        {
            message.append("\n  ");
            message.append(defaultPath);
            message.append(" (default)");
        }
        if (!bEnvIsSet && envName != NULL)
        {
            message.append(
                    formatString("\nYou can set additional directories to search "
                                 "with the %s path variable.", envName));
        }
        GMX_THROW(FileIOError(message));
    }
    return std::string();
}

std::vector<DataFileInfo>
DataFileFinder::enumerateFiles(const DataFileOptions &options) const
{
    // TODO: Consider if not being able to list one of the directories should
    // really be a fatal error. Or alternatively, check somewhere else that
    // paths in GMXLIB are valid.
    std::vector<DataFileInfo>                result;
    std::vector<std::string>::const_iterator i;
    if (options.bCurrentDir_)
    {
        std::vector<std::string> files
            = DirectoryEnumerator::enumerateFilesWithExtension(
                        ".", options.filename_, false);
        for (i = files.begin(); i != files.end(); ++i)
        {
            result.push_back(DataFileInfo(".", *i, false));
        }
    }
    if (impl_.get())
    {
        std::vector<std::string>::const_iterator j;
        for (j = impl_->searchPath_.begin(); j != impl_->searchPath_.end(); ++j)
        {
            std::vector<std::string> files
                = DirectoryEnumerator::enumerateFilesWithExtension(
                            j->c_str(), options.filename_, false);
            for (i = files.begin(); i != files.end(); ++i)
            {
                result.push_back(DataFileInfo(*j, *i, false));
            }
        }
    }
    const std::string &defaultPath = Impl::getDefaultPath();
    if (!defaultPath.empty())
    {
        std::vector<std::string> files
            = DirectoryEnumerator::enumerateFilesWithExtension(
                        defaultPath.c_str(), options.filename_, false);
        for (i = files.begin(); i != files.end(); ++i)
        {
            result.push_back(DataFileInfo(defaultPath, *i, true));
        }
    }
    if (result.empty() && options.bThrow_)
    {
        // TODO: Print the search path as is done in findFile().
        std::string message(
                formatString("Could not find any files ending on '%s' in the "
                             "current directory or the GROMACS library search path",
                             options.filename_));
        GMX_THROW(FileIOError(message));
    }
    return result;
}

} // namespace gmx
