/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Implements gmx::DataFileFinder.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/datafilefinder.h"

#include <cstdlib>

#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gromacs/utility/directoryenumerator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "buildinfo.h"

namespace gmx
{

/********************************************************************
 * DataFileFinder::Impl
 */

class DataFileFinder::Impl
{
public:
    static std::filesystem::path getDefaultPath();

    Impl() : envName_(nullptr), bEnvIsSet_(false) {}

    const char*                        envName_;
    bool                               bEnvIsSet_;
    std::vector<std::filesystem::path> searchPath_;
};

std::filesystem::path DataFileFinder::Impl::getDefaultPath()
{
    const InstallationPrefixInfo installPrefix = getProgramContext().installationPrefix();
    if (!installPrefix.path_.empty())
    {
        const char* const dataPath = installPrefix.sourceLayoutTreeLike_ ? "share" : GMX_INSTALL_GMXDATADIR;
        return std::filesystem::path(installPrefix.path_).append(dataPath).append("top");
    }
    return std::filesystem::path();
}

/********************************************************************
 * DataFileFinder
 */

DataFileFinder::DataFileFinder() : impl_(nullptr) {}

DataFileFinder::~DataFileFinder() {}

void DataFileFinder::setSearchPathFromEnv(const char* envVarName)
{
    if (!impl_)
    {
        impl_ = std::make_unique<Impl>();
    }
    impl_->envName_       = envVarName;
    const char* const lib = getenv(envVarName);
    if (!isNullOrEmpty(lib))
    {
        auto&                              path        = impl_->searchPath_; // convenience
        auto                               defaultPath = impl_->getDefaultPath();
        std::vector<std::filesystem::path> tmpPath     = splitPathEnvironment(lib);
        std::set<std::filesystem::path>    pathsSeen;
        pathsSeen.insert(defaultPath);
        for (auto& d : tmpPath)
        {
            if (!pathsSeen.count(d))
            {
                path.push_back(d);
                pathsSeen.insert(d);
            }
        }
        impl_->bEnvIsSet_ = true;
    }
}

FilePtr DataFileFinder::openFile(const DataFileOptions& options) const
{
    // TODO: There is a small race here, since there is some time between
    // the exists() calls and actually opening the file.  It would be better
    // to leave the file open after a successful exists() if the desire is to
    // actually open the file.
    auto filename = findFile(options);
    if (filename.empty())
    {
        return nullptr;
    }
    return TextInputFile::openRawHandle(filename);
}

std::filesystem::path DataFileFinder::findFile(const DataFileOptions& options) const
{
    if (options.bCurrentDir_ && std::filesystem::exists(options.filename_))
    {
        return options.filename_;
    }
    if (impl_ != nullptr)
    {
        for (auto path : impl_->searchPath_)
        {
            // TODO: Deal with an empty search path entry more reasonably.
            path.append(options.filename_.string());
            // TODO: Consider skipping directories.
            if (std::filesystem::exists(path))
            {
                return path;
            }
        }
    }
    const auto& defaultPath = Impl::getDefaultPath();
    if (!defaultPath.empty())
    {
        auto testPath = std::filesystem::path(defaultPath).append(options.filename_.string());
        if (std::filesystem::exists(testPath))
        {
            return testPath;
        }
    }
    if (options.bThrow_)
    {
        const char* const envName   = (impl_ != nullptr ? impl_->envName_ : nullptr);
        const bool        bEnvIsSet = (impl_ != nullptr ? impl_->bEnvIsSet_ : false);
        std::string       message(
                formatString("Library file '%s' not found", options.filename_.string().c_str()));
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
            message.append(std::filesystem::current_path().string());
            message.append(" (current dir)");
        }
        if (impl_ != nullptr)
        {
            for (const auto& path : impl_->searchPath_)
            {
                message.append("\n  ");
                message.append(path.string());
            }
        }
        if (!defaultPath.empty())
        {
            message.append("\n  ");
            message.append(defaultPath.string());
            message.append(" (default)");
        }
        if (!bEnvIsSet && envName != nullptr)
        {
            message.append(
                    formatString("\nYou can set additional directories to search "
                                 "with the %s path variable.",
                                 envName));
        }
        GMX_THROW(FileIOError(message));
    }
    return std::string();
}

std::vector<DataFileInfo> DataFileFinder::enumerateFiles(const DataFileOptions& options) const
{
    // TODO: Consider if not being able to list one of the directories should
    // really be a fatal error. Or alternatively, check somewhere else that
    // paths in GMXLIB are valid.
    std::vector<DataFileInfo> result;
    if (options.bCurrentDir_)
    {
        auto files = DirectoryEnumerator::enumerateFilesWithExtension(
                std::filesystem::current_path(), options.filename_.string(), false);
        for (const auto& file : files)
        {
            result.emplace_back(".", file, false);
        }
    }
    if (impl_ != nullptr)
    {
        for (const auto& path : impl_->searchPath_)
        {
            auto files = DirectoryEnumerator::enumerateFilesWithExtension(
                    path, options.filename_.string(), false);
            for (const auto& file : files)
            {
                result.emplace_back(path, file, false);
            }
        }
    }
    const auto& defaultPath = Impl::getDefaultPath();
    if (!defaultPath.empty())
    {
        auto files = DirectoryEnumerator::enumerateFilesWithExtension(
                defaultPath, options.filename_.string(), false);
        for (const auto& file : files)
        {
            result.emplace_back(defaultPath, file, true);
        }
    }
    if (result.empty() && options.bThrow_)
    {
        // TODO: Print the search path as is done in findFile().
        std::string message(
                formatString("Could not find any files ending on '%s' in the "
                             "current directory or the GROMACS library search path",
                             options.filename_.string().c_str()));
        GMX_THROW(FileIOError(message));
    }
    return result;
}

} // namespace gmx
