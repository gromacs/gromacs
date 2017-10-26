/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

std::vector<int>
parseUserGpuIds(const std::string &gpuIdString)
{
    // An optional comma is used to separate GPU IDs assigned to the
    // same type of task, which will be useful for any nodes that have
    // more than ten GPUs.

    std::vector<int> digits;
    auto             foundCommaDelimiters = gpuIdString.find(',') != std::string::npos;
    if (!foundCommaDelimiters)
    {
        for (const auto &c : gpuIdString)
        {
            if (std::isdigit(c) == 0)
            {
                GMX_THROW(InvalidInputError(formatString("Invalid character in GPU ID string: \"%c\"\n", c)));
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

std::vector<int>
makeGpuIds(ArrayRef<const int> compatibleGpus,
           size_t              numGpuTasks)
{
    std::vector<int> gpuIdsToUse;

    gpuIdsToUse.reserve(numGpuTasks);

    auto currentGpuId = compatibleGpus.begin();
    for (size_t i = 0; i != numGpuTasks; ++i)
    {
        GMX_ASSERT(!compatibleGpus.empty(), "Must have compatible GPUs from which to build a list of GPU IDs to use");
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

std::string
makeGpuIdString(const std::vector<int> &gpuIds,
                int                     totalNumberOfTasks)
{
    auto resultGpuIds = makeGpuIds(gpuIds, totalNumberOfTasks);
    return formatAndJoin(resultGpuIds, ",", StringFormatter("%d"));
}

void checkUserGpuIds(const gmx_gpu_info_t   &gpu_info,
                     const std::vector<int> &compatibleGpus,
                     const std::vector<int> &gpuIds)
{
    bool        foundIncompatibleGpuIds = false;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (const auto &gpuId : gpuIds)
    {
        if (std::find(compatibleGpus.begin(), compatibleGpus.end(), gpuId) == compatibleGpus.end())
        {
            foundIncompatibleGpuIds = true;
            message                += gmx::formatString("    GPU #%d: %s\n",
                                                        gpuId,
                                                        getGpuCompatibilityDescription(gpu_info, gpuId));
        }
    }
    if (foundIncompatibleGpuIds)
    {
        GMX_THROW(InconsistentInputError(message));
    }
}

} // namespace
