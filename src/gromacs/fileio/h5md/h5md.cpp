/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief I/o interface to H5MD HDF5 files.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#include "gmxpre.h"

#include "h5md.h"

#include "config.h"

#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/utility/exceptions.h"

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "h5md_attribute.h"
#    include "h5md_error.h"
#    include "h5md_group.h"
#    include "h5md_guard.h"
#    include "h5md_util.h"
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")
#else
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")
#endif // GMX_USE_HDF5

namespace gmx
{

H5md::H5md(const std::filesystem::path& fileName, const H5mdFileMode mode)
{
#if GMX_USE_HDF5
    /* Disable automatic HDF5 error output, e.g. when items are not found. Explicit H5EPrint2() will
     * still print error messages. */
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    switch (mode)
    {
        case H5mdFileMode::Write:
            file_ = H5Fcreate(
                    fileName.string().c_str(), H5F_ACC_TRUNC, H5Pcreate(H5P_FILE_CREATE), H5P_DEFAULT);
            break;
        case H5mdFileMode::Read:
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
        default: throw NotImplementedError("Appending to H5MD is not implemented yet.");
    }
    gmx::throwUponInvalidHid(file_, "Cannot open H5MD file.");

    filemode_ = mode;

#else
    GMX_UNUSED_VALUE(fileName);
    GMX_UNUSED_VALUE(mode);
    throw NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

H5md::~H5md()
{
#if GMX_USE_HDF5
    if (handleIsValid(file_))
    {
        // Do not throw exceptions when flushing from the destructor.
        flush(false);
        H5Fclose(file_);
    }

    /* Do not throw, if GMX_USE_HDF5 is false, in the destructor. */

#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
hid_t H5md::fileid() const
{
#if GMX_USE_HDF5
    return file_;
#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::flush(bool throwExceptionUponError)
{
#if GMX_USE_HDF5
    if (throwExceptionUponError)
    {
        GMX_ASSERT(handleIsValid(file_), "Cannot flush an invalid H5MD file.");
        gmx::throwUponH5mdError(H5Fflush(file_, H5F_SCOPE_LOCAL) < 0, "Error flushing H5MD.");
    }
    else
    {
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }

#else
    GMX_UNUSED_VALUE(throwExceptionUponError);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setAuthor(const std::string& authorName)
{
#if GMX_USE_HDF5
    const auto [authorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/author"));
    setAttribute(authorGroup, "name", authorName.c_str());
#else
    GMX_UNUSED_VALUE(authorName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::author()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/author"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/author"));
        return getAttribute(group, "name");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setCreatorProgramName(const std::string& creatorName)
{
#if GMX_USE_HDF5
    const auto [creatorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/creator"));
    setAttribute(creatorGroup, "name", creatorName.c_str());
#else
    GMX_UNUSED_VALUE(creatorName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::creatorProgramName()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/creator"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/creator"));
        return getAttribute(group, "name");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setCreatorProgramVersion(const std::string& version)
{
#if GMX_USE_HDF5
    const auto [creatorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/creator"));
    setAttribute(creatorGroup, "version", version.c_str());
#else
    GMX_UNUSED_VALUE(version);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::creatorProgramVersion()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/creator"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/creator"));
        return getAttribute(group, "version");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
