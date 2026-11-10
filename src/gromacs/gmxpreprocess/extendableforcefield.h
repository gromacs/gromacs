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
#ifndef GMX_GMXPREPROCESS_EXTENDABLEFORCEFIELD_H
#define GMX_GMXPREPROCESS_EXTENDABLEFORCEFIELD_H

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \brief A force field directory, optionally extended by addon directories.
 *
 * Encapsulates where force field data files (residue topologies, hydrogen
 * and termini databases, ion topologies, ...) are searched for, so that
 * code that reads such files can focus on parsing and topology handling
 * instead of repeating search mechanics.
 *
 * The search order is deterministic: the base force field directory comes
 * first, followed by the addon directories in the order they were given
 * to the constructor. A file found in the base force field therefore always
 * takes precedence over an addon file of the same name, and an addon's
 * relative order among the addons defines its priority.
 *
 * The base force field directory must be given and must exist; empty addon
 * directory entries are allowed and never match anything.
 */
class ExtendableForceField
{
public:
    /*! \brief What to do when a search finds no files at all. */
    enum class OnEmptySearch
    {
        ReturnEmptyList, /**< Return an empty list. */
        ThrowError       /**< Throw gmx::InvalidInputError. */
    };

    /*! \brief Construct for the given base force field directory.
     *
     * \p forceFieldDir may be a directory name relative to the GROMACS
     * data directories or a full path.
     *
     * \throws gmx::InvalidInputError if \p forceFieldDir is empty or no
     *         such directory can be found.
     */
    explicit ExtendableForceField(const std::filesystem::path& forceFieldDir);

    /*! \brief Construct for a base force field directory extended by addons.
     *
     * \p forceFieldDir and the \p addonDirs entries may be directory names
     * relative to the GROMACS data directories or full paths.
     *
     * \throws gmx::InvalidInputError if \p forceFieldDir is empty or no
     *         such directory can be found.
     */
    ExtendableForceField(const std::filesystem::path&          forceFieldDir,
                         ArrayRef<const std::filesystem::path> addonDirs);

    //! The directories that are searched, base force field first.
    //! References into *this, and is invalid after *this is destroyed.
    ArrayRef<const std::filesystem::path> searchDirectories() const;

    //! The base force field directory, as given to the constructor.
    const std::filesystem::path& forceFieldName() const;

    //! \p fileName in the first search directory that contains it.
    /*
     * \p fileName must be a relative, single-component file name
     * ("file.itp", not "dir/file.itp").
     */
    std::optional<std::filesystem::path> findSinglePath(const std::filesystem::path& fileName) const;

    //! Like findSinglePath(), but throws gmx::InvalidInputError if not found.
    std::filesystem::path findRequiredPath(const std::filesystem::path& fileName) const;

    //! Every search directory that contains \p fileName, base force field first.
    /*
     * \p fileName must be a relative, single-component file name
     * ("file.itp", not "dir/file.itp").
     */
    std::vector<std::filesystem::path> findAllPaths(const std::filesystem::path& fileName) const;

    //! All files ending on \p fileEnd in the search directories, base force field first.
    std::vector<std::filesystem::path>
    findAllPathsEndingWith(std::string_view fileEnd,
                           OnEmptySearch    onEmpty = OnEmptySearch::ReturnEmptyList) const;

private:
    std::vector<std::filesystem::path> directories_;
};

} // namespace gmx

#endif
