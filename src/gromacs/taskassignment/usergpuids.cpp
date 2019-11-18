/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Defines routines for handling user-specified GPU IDs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "usergpuids.h"

#include <cctype>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \brief Parse a GPU ID specifier string into a container.
 *
 * \param[in]   gpuIdString  String like "013" or "0,1,3" typically
 *                           supplied by the user.
 *                           Must contain only unique decimal digits, or only decimal
 *                           digits separated by comma delimiters. A terminal
 *                           comma is accceptable (and required to specify a
 *                           single ID that is larger than 9).
 *
 * \returns  A vector of numeric IDs extracted from \c gpuIdString.
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',').
 */
static std::vector<int> parseGpuDeviceIdentifierList(const std::string& gpuIdString)
{
    std::vector<int> digits;
    auto             foundCommaDelimiters = gpuIdString.find(',') != std::string::npos;
    if (!foundCommaDelimiters)
    {
        for (const auto& c : gpuIdString)
        {
            if (std::isdigit(c) == 0)
            {
                GMX_THROW(InvalidInputError(
                        formatString("Invalid character in GPU ID string: \"%c\"\n", c)));
            }
            // Convert each character in the token to an integer
            digits.push_back(c - '0');
        }
    }
    else
    {
        if (gpuIdString[0] == ',')
        {
            GMX_THROW(InvalidInputError("Invalid use of leading comma in GPU ID string"));
        }
        std::istringstream ss(gpuIdString);
        std::string        token;
        digits.reserve(gpuIdString.length());
        token.reserve(gpuIdString.length());
        while (std::getline(ss, token, ','))
        {
            // Convert the whole token to an integer
            if (token.empty())
            {
                GMX_THROW(InvalidInputError("Invalid use of comma in GPU ID string"));
            }
            digits.push_back(std::stoi(token));
        }
    }
    return digits;
}

std::vector<int> parseUserGpuIdString(const std::string& gpuIdString)
{
    // An optional comma is used to separate GPU IDs assigned to the
    // same type of task, which will be useful for any nodes that have
    // more than ten GPUs.

    auto digits = parseGpuDeviceIdentifierList(gpuIdString);

    // Check and enforce that no duplicate IDs are allowed
    for (size_t i = 0; i != digits.size(); ++i)
    {
        for (size_t j = i + 1; j != digits.size(); ++j)
        {
            if (digits[i] == digits[j])
            {
                GMX_THROW(
                        InvalidInputError(formatString("The string of available GPU device IDs "
                                                       "'%s' may not contain duplicate device IDs",
                                                       gpuIdString.c_str())));
            }
        }
    }
    return digits;
}

std::vector<int> makeGpuIdsToUse(const gmx_gpu_info_t& gpuInfo, const std::string& gpuIdsAvailableString)
{
    auto             compatibleGpus  = getCompatibleGpus(gpuInfo);
    std::vector<int> gpuIdsAvailable = parseUserGpuIdString(gpuIdsAvailableString);

    if (gpuIdsAvailable.empty())
    {
        return compatibleGpus;
    }

    std::vector<int> gpuIdsToUse;
    gpuIdsToUse.reserve(gpuIdsAvailable.size());
    std::vector<int> availableGpuIdsThatAreIncompatible;
    for (const auto& availableGpuId : gpuIdsAvailable)
    {
        bool availableGpuIsCompatible = false;
        for (const auto& compatibleGpuId : compatibleGpus)
        {
            if (availableGpuId == compatibleGpuId)
            {
                availableGpuIsCompatible = true;
                break;
            }
        }
        if (availableGpuIsCompatible)
        {
            gpuIdsToUse.push_back(availableGpuId);
        }
        else
        {
            // Prepare data for an error message about all incompatible available GPU IDs.
            availableGpuIdsThatAreIncompatible.push_back(availableGpuId);
        }
    }
    if (!availableGpuIdsThatAreIncompatible.empty())
    {
        auto message = "You requested mdrun to use GPUs with IDs " + gpuIdsAvailableString
                       + ", but that includes the following incompatible GPUs: "
                       + formatAndJoin(availableGpuIdsThatAreIncompatible, ",", StringFormatter("%d"))
                       + ". Request only compatible GPUs.";
        GMX_THROW(InvalidInputError(message));
    }
    return gpuIdsToUse;
}

std::vector<int> parseUserTaskAssignmentString(const std::string& gpuIdString)
{
    // Implement any additional constraints here that need to be imposed

    return parseGpuDeviceIdentifierList(gpuIdString);
}

std::vector<int> makeGpuIds(ArrayRef<const int> compatibleGpus, size_t numGpuTasks)
{
    std::vector<int> gpuIdsToUse;

    gpuIdsToUse.reserve(numGpuTasks);

    auto currentGpuId = compatibleGpus.begin();
    for (size_t i = 0; i != numGpuTasks; ++i)
    {
        GMX_ASSERT(!compatibleGpus.empty(),
                   "Must have compatible GPUs from which to build a list of GPU IDs to use");
        gpuIdsToUse.push_back(*currentGpuId);
        ++currentGpuId;
        if (currentGpuId == compatibleGpus.end())
        {
            // Wrap around and assign tasks again.
            currentGpuId = compatibleGpus.begin();
        }
    }
    std::sort(gpuIdsToUse.begin(), gpuIdsToUse.end());
    return gpuIdsToUse;
}

std::string makeGpuIdString(const std::vector<int>& gpuIds, int totalNumberOfTasks)
{
    auto resultGpuIds = makeGpuIds(gpuIds, totalNumberOfTasks);
    return formatAndJoin(resultGpuIds, ",", StringFormatter("%d"));
}

void checkUserGpuIds(const gmx_gpu_info_t&   gpu_info,
                     const std::vector<int>& compatibleGpus,
                     const std::vector<int>& gpuIds)
{
    bool        foundIncompatibleGpuIds = false;
    std::string message =
            "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (const auto& gpuId : gpuIds)
    {
        if (std::find(compatibleGpus.begin(), compatibleGpus.end(), gpuId) == compatibleGpus.end())
        {
            foundIncompatibleGpuIds = true;
            message += gmx::formatString("    GPU #%d: %s\n", gpuId,
                                         getGpuCompatibilityDescription(gpu_info, gpuId));
        }
    }
    if (foundIncompatibleGpuIds)
    {
        GMX_THROW(InconsistentInputError(message));
    }
}

} // namespace gmx
