/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal
 * \file
 * \brief Implements gmx::ExtendableForceField.
 */
#include "gmxpre.h"

#include "extendableforcefield.h"

#include <utility>

#include "gromacs/utility/datafilefinder.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "fflibutil.h"

namespace gmx
{

ExtendableForceField::ExtendableForceField(const std::filesystem::path& forceFieldDir) :
    ExtendableForceField(forceFieldDir, {})
{
}

ExtendableForceField::ExtendableForceField(const std::filesystem::path&          forceFieldDir,
                                           ArrayRef<const std::filesystem::path> addonDirs) :
    directories_{ forceFieldDir }
{
    directories_.insert(directories_.end(), addonDirs.begin(), addonDirs.end());

    if (forceFieldDir.empty())
    {
        throw InvalidInputError("No force field directory was specified.");
    }
    if (getLibraryFileFinder().findFile(DataFileOptions(forceFieldDir).throwIfNotFound(false)).empty())
    {
        throw InvalidInputError(formatString("The force field directory '%s' does not exist.",
                                             forceFieldDir.string().c_str()));
    }
}

ArrayRef<const std::filesystem::path> ExtendableForceField::searchDirectories() const
{
    return directories_;
}

const std::filesystem::path& ExtendableForceField::forceFieldName() const
{
    return directories_.front();
}

std::optional<std::filesystem::path> ExtendableForceField::findSinglePath(const std::filesystem::path& fileName) const
{
    GMX_ASSERT(fileName.is_relative() && fileName.filename() == fileName,
               "fileName must be a single-component relative file name");

    for (const auto& dir : directories_)
    {
        if (!dir.empty() && fflib_fexist(dir / fileName))
        {
            return dir / fileName;
        }
    }
    return std::nullopt;
}

std::filesystem::path ExtendableForceField::findRequiredPath(const std::filesystem::path& fileName) const
{
    if (auto found = findSinglePath(fileName))
    {
        return *found;
    }
    throw InvalidInputError(formatString("File '%s' not found in force field directory '%s'",
                                         fileName.string().c_str(),
                                         directories_.front().string().c_str()));
}

std::vector<std::filesystem::path> ExtendableForceField::findAllPaths(const std::filesystem::path& fileName) const
{
    GMX_ASSERT(fileName.is_relative() && fileName.filename() == fileName,
               "fileName must be a single-component relative file name");

    std::vector<std::filesystem::path> result;
    for (const auto& dir : directories_)
    {
        if (!dir.empty() && fflib_fexist(dir / fileName))
        {
            result.push_back(dir / fileName);
        }
    }
    return result;
}

std::vector<std::filesystem::path> ExtendableForceField::findAllPathsEndingWith(std::string_view fileEnd,
                                                                                OnEmptySearch onEmpty) const
{
    std::vector<std::filesystem::path> result;
    for (const auto& dir : directories_)
    {
        if (!dir.empty())
        {
            // bFatalError is always false: the class throws exceptions instead of
            // exiting through the legacy error handler.
            auto found = fflib_search_file_end(dir, std::string(fileEnd).c_str(), false);
            result.insert(result.end(), found.begin(), found.end());
        }
    }
    if (result.empty() && onEmpty == OnEmptySearch::ThrowError)
    {
        std::string message = formatString(
                "Could not find any files ending on '%s' "
                "in the force field directory '%s'",
                std::string(fileEnd).c_str(),
                directories_.front().string().c_str());
        if (directories_.size() > 1)
        {
            message += formatString(" or the %zu addon directories: %s",
                                    directories_.size() - 1,
                                    formatAndJoin(directories_.begin() + 1,
                                                  directories_.end(),
                                                  ", ",
                                                  [](const std::filesystem::path& dir)
                                                  { return "'" + dir.string() + "'"; })
                                            .c_str());
        }
        throw InvalidInputError(message);
    }
    return result;
}

} // namespace gmx
