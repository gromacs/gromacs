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
#include <functional>
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

void checkUserTaskAssignmentGpuIds(const gmx_gpu_info_t     &gpu_info,
                                   const std::vector<int>   &compatibleGpus,
                                   const GpuTaskAssignments &gpuTaskAssignmentsOnNode)
{
    bool        foundIncompatibleGpuIds = false;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (const auto &gpuTaskAssignmentOnRank : gpuTaskAssignmentsOnNode)
    {
        for (const auto &gpuTaskAssignment : gpuTaskAssignmentOnRank)
        {
            int gpuId = gpuTaskAssignment.deviceId_;
            if (std::find(compatibleGpus.begin(), compatibleGpus.end(), gpuId) == compatibleGpus.end())
            {
                foundIncompatibleGpuIds = true;
                message                += gmx::formatString("    GPU #%d: %s\n",
                                                            gpuId,
                                                            getGpuCompatibilityDescription(gpu_info, gpuId));
            }
        }
    }
    if (foundIncompatibleGpuIds)
    {
        GMX_THROW(InconsistentInputError(message));
    }
}

namespace
{

/*! \brief Append to \c assignment all the GPU IDs for tasks of type
 * \c taskType that are in \c taskMapping.
 *
 * \throws InvalidInputError  If the user input is malformed.
 *         std::bad_alloc     If out of memory.
 */
void
addAssignmentsForTaskType(GpuTaskAssignment *assignment,
                          const std::string taskMapping,
                          GpuTask taskType)
{
    GMX_ASSERT(assignment, "Need valid task assignment to extend");

    auto gpuIds = splitDelimitedString(taskMapping, ',');
    for (const auto &gpuId : gpuIds)
    {
        size_t firstDigitIndex = gpuId.find_first_of("0123456789");
        if (gpuId.compare(0, firstDigitIndex, "G") != 0 &&
            gpuId.compare(0, firstDigitIndex, "GPU") != 0)
        {
            GMX_THROW(InvalidInputError
                      (formatString("Only device assignments of the form 'G0', 'GPU1', etc. are supported, not '%s'.", gpuId.c_str())));
        }
        // This should be the end of the string
        size_t nextNonDigitIndex = gpuId.find_first_not_of("0123456789", firstDigitIndex);
        if (nextNonDigitIndex != std::string::npos)
        {
            GMX_THROW(InvalidInputError
                      (formatString("Only device assignments of the form 'G0', 'GPU1', etc. are supported, not '%s'.", gpuId.c_str())));
        }
        std::string digits = gpuId.substr(firstDigitIndex);
        assignment->emplace_back(GpuTaskMapping{taskType, std::stoi(digits)});
    }
}

/*! \brief Return the vector of tasks on this rank described by \c rankAssignment.
 *
 * \throws InvalidInputError  If the user input is malformed.
 *         std::bad_alloc     If out of memory.
 */
GpuTaskAssignment
parseRankAssignment(const std::string &rankAssignment)
{
    GpuTaskAssignment assignmentsForRank;

    auto taskAssignments = splitDelimitedString(rankAssignment, ';');
    for (const auto &taskAssignment : taskAssignments)
    {
        auto colonIndex = taskAssignment.find(':');
        if (colonIndex == std::string::npos)
        {
            GMX_THROW(InvalidInputError
                      (formatString("Found no colon separating task type from assigned devices in '%s'.", taskAssignment.c_str())));
        }
        GpuTask taskType;
        if (std::string::npos != taskAssignment.find("NB", 0, colonIndex))
        {
            taskType = GpuTask::Nonbonded;
        }
        else if (std::string::npos != taskAssignment.find("PME", 0, colonIndex))
        {
            taskType = GpuTask::Pme;
        }
        else
        {
            GMX_THROW(InvalidInputError
                      (formatString("The only supported task types are 'NB' and 'PME', use one of those, rather than '%s'.", taskAssignment.c_str())));
        }
        std::string gpuIdString = taskAssignment.substr(colonIndex+1);
        addAssignmentsForTaskType(&assignmentsForRank, gpuIdString, taskType);
    }

    return assignmentsForRank;
}

/*! \brief Return the contents of \c input in all upper case, having
 * removed all whitespace characters.
 *
 * \throws std::bad_alloc     If out of memory.
 */
std::string cleanInput(const std::string &input)
{
    std::string output;
    output.reserve(input.size());
    for (const char &c : input)
    {
        if (!std::isspace(c))
        {
            output.push_back(std::toupper(c));
        }
    }
    return output;
}

}   // namespace

// TODO This code makes a lot of strings, vectors of strings and
// substrings. It would be nicer to be able to use a string_view
// class, but we don't have that yet. Fortunately, it only runs once
// per simulation, and only if the user makes an explicit task
// assignment.
GpuTaskAssignments
parseUserTaskAssignment(const std::string &userTaskAssignmentString)
{
    GpuTaskAssignments assignments;
    if (userTaskAssignmentString.empty())
    {
        return assignments;
    }
    try
    {
        std::string cleanedAssignment = cleanInput(userTaskAssignmentString);
        auto theStart = cleanedAssignment.begin();
        auto theEnd = cleanedAssignment.end();

        // Find matching pairs of brackets, and then the assigment
        // string between them.
        auto currentStart = theStart;
        while (currentStart != theEnd)
        {
            if (*currentStart != '[')
            {
                GMX_THROW(InvalidInputError
                          (formatString("Found text '%s' that was not part of a rank assignment.",
                                        std::string(currentStart, theEnd).c_str())));
            }
            auto leftBracketIt = currentStart;
            auto rightBracketIt = std::find(leftBracketIt+1, theEnd, ']');
            if (rightBracketIt == theEnd)
            {
                GMX_THROW(InvalidInputError
                          (formatString("Found text '%s' that does not include a rank assignment terminator ']'.",
                                        std::string(leftBracketIt+1, theEnd).c_str())));
            }
            assignments.push_back(parseRankAssignment(std::string(leftBracketIt+1, rightBracketIt)));
            
            currentStart = rightBracketIt + 1;
        }
    }
    catch (InvalidInputError &ex)
    {
        ex.prependContext(formatString("While parsing task assignment string '%s', an error was encountered:", userTaskAssignmentString.c_str()));
        throw;
    }
    return assignments;
}

} // namespace
