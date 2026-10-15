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
/*! \internal \file
 * \brief Implements CP2K support information and version-dependent MPI requirements.
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "qmmmversion.h"

#include <array>
#include <charconv>
#include <string>
#include <string_view>
#include <system_error>

#if GMX_CP2K
#    include <libcp2k.h>
#endif

namespace gmx
{

std::string qmmmDescription()
{
#if GMX_CP2K
    std::array<char, 100> cp2kVersion{};
    cp2k_get_version(cp2kVersion.data(), 100);
    return cp2kVersion[0] != '\0' ? std::string("enabled (") + cp2kVersion.data() + ')' : "enabled";
#else
    return "disabled";
#endif
}

bool qmmmRequiresMpiThreadMultiple(std::string_view description)
{
    // Disabled CP2K support imposes no MPI threading requirement.
    if (description == "disabled")
    {
        return false;
    }
    /* Expect the description produced by qmmmDescription().
     * If the format or the version cannot be recognized, conservatively return true.
     */
    constexpr std::string_view prefix = "enabled (CP2K version ";
    if (description.substr(0, prefix.size()) != prefix || description.back() != ')')
    {
        return true;
    }
    // Strip the description wrapper, leaving the version and any development marker.
    description.remove_prefix(prefix.size());
    description.remove_suffix(1);
    if (description.empty())
    {
        return true;
    }

    // Parse the required major.minor components; reject missing or overflowing numbers.
    unsigned int major = 0, minor = 0;
    const char*  end         = description.data() + description.size();
    const auto   majorResult = std::from_chars(description.data(), end, major);
    if (majorResult.ec != std::errc{} || majorResult.ptr == end || *majorResult.ptr != '.')
    {
        return true;
    }
    const auto minorResult = std::from_chars(majorResult.ptr + 1, end, minor);
    if (minorResult.ec != std::errc{})
    {
        return true;
    }
    // Validate an optional patch component, although it does not affect the requirement.
    std::string_view suffix(minorResult.ptr, end - minorResult.ptr);
    if (!suffix.empty() && suffix.front() == '.')
    {
        unsigned int patch       = 0;
        const auto   patchResult = std::from_chars(suffix.data() + 1, end, patch);
        if (patchResult.ec != std::errc{})
        {
            return true;
        }
        suffix = std::string_view(patchResult.ptr, end - patchResult.ptr);
    }
    // Only an empty suffix or CP2K's exact development marker is recognized.
    const bool isDevelopmentVersion = suffix == " (Development Version)";
    if (!suffix.empty() && !isDevelopmentVersion)
    {
        return true;
    }
    // CP2K master gained the requirement before advancing its version string from 2026.2.
    return major >= 2027 || (major == 2026 && minor == 2 && isDevelopmentVersion)
           || (major == 2026 && minor > 2);
}

} // namespace gmx
