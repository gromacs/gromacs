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
 * \brief Definitions of utility functions to read and write index groups in H5MD files.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_indexgroup.h"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/functions.h"
#include "gromacs/utility/stringutil.h"


namespace gmx
{

//! \brief Name of dataset for the atom ranges in an index group
constexpr char c_atomRangesDataSetName[] = "ranges";
//! \brief Name of dataset for the indices of a plain index group
constexpr char c_simpleIndicesDataSetName[] = "indices";
//! \brief Number of elements in a row of the atom ranges dataset: [start, stop, step]
constexpr int c_atomRangeRowSize = 3;

namespace
{

/*! \brief Sanity check for the input index group container and index groups to write
 *
 * \note Sanity check:
 * - If any group name is empty or contains the HDF5 path separator
 * - If any index group already exists in the H5MD file
 * - If any duplicated entries in the input index groups
 * - If any index is negative
 *
 * \param[in] baseContainer The HDF5 container to write the index groups into.
 * \param[in] indexGroups   The vector of \c IndexGroup objects to write.
 *
 * \throws H5mdError If \p indexGroups is empty, if any group name is empty or contains the
 *                   HDF5 path separator, if a group of that name already exists in
 *                   \p baseContainer or is repeated within \p indexGroups, or if any index
 *                   of a group is negative.
 */
void indexGroupArraySanityCheck(const hid_t baseContainer, ArrayRef<const IndexGroup> indexGroups)
{
    GMX_H5MD_THROW_UPON_ERROR(indexGroups.empty(), "Cannot write an empty index group array.");
    for (auto indexGroupIt = indexGroups.begin(); indexGroupIt != indexGroups.end(); ++indexGroupIt)
    {
        const IndexGroup& indexGroup = *indexGroupIt;
        // A name that HDF5 reads as a path would silently create nested groups, which the
        // reader then cannot enumerate, so reject such names rather than lose the group.
        GMX_H5MD_THROW_UPON_ERROR(indexGroup.name.empty(),
                                  "Cannot write an index group with an empty name.");
        GMX_H5MD_THROW_UPON_ERROR(
                indexGroup.name.find('/') != std::string::npos,
                formatString("Cannot write the index group named '%s': the name must not contain "
                             "'/', which HDF5 interprets as a path separator.",
                             indexGroup.name.c_str()));
        // Raise if group already exists in the file
        GMX_H5MD_THROW_UPON_ERROR(
                objectExists(baseContainer, indexGroup.name.c_str()),
                formatString("Cannot write the duplicated group named: %s.", indexGroup.name.c_str()));
        // Raise if the same name appears more than once in the input
        GMX_H5MD_THROW_UPON_ERROR(std::any_of(indexGroups.begin(),
                                              indexGroupIt,
                                              [&indexGroup](const IndexGroup& previous)
                                              { return previous.name == indexGroup.name; }),
                                  formatString("Found duplicated group in the input named: %s.",
                                               indexGroup.name.c_str()));
        GMX_H5MD_THROW_UPON_ERROR(
                std::any_of(indexGroup.particleIndices.begin(),
                            indexGroup.particleIndices.end(),
                            [](const int index) { return index < 0; }),
                formatString("Cannot write the index group named '%s': particle indices must not "
                             "be negative.",
                             indexGroup.name.c_str()));
    }
}


/*! \brief Parse the input indices into plain indices and atom ranges.
 *
 * The inspiration is drawn from \ref prettyPrintListAsRange() in stringutil.cpp.
 *
 * \note Indices remain 0-based. Atom ranges are Python-style [start, stop, step] where stop is
 * exclusive: the sequence is start, start+step, ..., < stop.
 * Atom ranges are stored as a [N, 3] dataset with each row standing for [start, stop, step].
 *
 * \param[in] inputIndices The input indices to parse, in any order. Only runs that ascend in
 *                         this order become ranges, so ascending input is stored most compactly.
 * \returns The plain indices, and the atom ranges stored as a flat vector of
 *          [start1, stop1, step1, start2, stop2, step2, ...].
 */
std::pair<std::vector<int>, std::vector<int>> parseIndices(ArrayRef<const int> inputIndices)
{
    std::vector<int> outputIndices;
    std::vector<int> outputRanges;
    Index            i = 0;
    while (i < inputIndices.ssize())
    {
        if (i + 1 < inputIndices.ssize())
        {
            const int step = inputIndices[i + 1] - inputIndices[i];
            Index     end  = i + 1;

            // Extend the arithmetic sequence as far as possible
            while (end + 1 < inputIndices.ssize() && inputIndices[end + 1] - inputIndices[end] == step)
            {
                ++end;
            }

            // Sequence of at least 3 elements with positive step -> store as range [start, stop, step]
            if (end - i >= 2 && step > 0)
            {
                outputRanges.push_back(inputIndices[i]);          // 0-based start (inclusive)
                outputRanges.push_back(inputIndices[end] + step); // 0-based stop (exclusive)
                outputRanges.push_back(step);
                i = end + 1;
                continue;
            }
        }

        // Single element or short/non-positive-step sequence -> store as plain index (0-based)
        outputIndices.push_back(inputIndices[i]);
        ++i;
    }
    return { std::move(outputIndices), std::move(outputRanges) };
}

} // namespace

std::vector<std::string> getNamesOfIndexGroups(const hid_t baseContainer)
{
    std::vector<std::string> groupNames;
    auto callback = [](hid_t group, const char* name, const H5L_info_t* /*info*/, void* opData) noexcept -> herr_t
    {
        std::string atomRangesPath    = joinStrings({ name, c_atomRangesDataSetName }, "/");
        std::string simpleIndicesPath = joinStrings({ name, c_simpleIndicesDataSetName }, "/");
        if (objectExists(group, atomRangesPath.c_str()) || objectExists(group, simpleIndicesPath.c_str()))
        {
            static_cast<std::vector<std::string>*>(opData)->push_back(name);
        }
        return H5_ITER_CONT;
    };
    const herr_t status =
            H5Literate(baseContainer, H5_INDEX_NAME, H5_ITER_INC, nullptr, callback, &groupNames);
    GMX_H5MD_THROW_UPON_ERROR(status < 0, "Failed to list the index groups of the base container.");
    return groupNames;
}

void writeIndexGroups(const hid_t baseContainer, ArrayRef<const IndexGroup> indexGroups)
{
    indexGroupArraySanityCheck(baseContainer, indexGroups);
    for (const IndexGroup& indexGroup : indexGroups)
    {
        const auto [group, groupGuard] =
                makeH5mdGroupGuard(createGroup(baseContainer, indexGroup.name.c_str()));

        const auto [indices, ranges] = parseIndices(indexGroup.particleIndices);

        /* A group with no particles is still given an `indices` data set, of length zero, so
         * that it is found when the container is enumerated and reads back as a genuinely
         * empty group. Holding no data set at all is how a string-based selection is
         * recognised, so an empty index group must not be written that way. */
        if (!indices.empty() || ranges.empty())
        {
            const H5mdFixedDataSet<int> dataset = H5mdDataSetBuilder<int>(group, c_simpleIndicesDataSetName)
                                                          .withDimension({ indices.size() })
                                                          .build();
            dataset.writeData(indices);
        }

        if (!ranges.empty())
        {
            const size_t                numRanges = ranges.size() / c_atomRangeRowSize;
            const H5mdFixedDataSet<int> dataset =
                    H5mdDataSetBuilder<int>(group, c_atomRangesDataSetName)
                            .withDimension({ numRanges, c_atomRangeRowSize })
                            .build();
            dataset.writeData(ranges);
        }
    }
}

void writeIndexGroup(const hid_t baseContainer, const IndexGroup& indexGroup)
{
    ArrayRef<const IndexGroup> indexGroups = arrayRefFromArray(&indexGroup, 1);
    writeIndexGroups(baseContainer, indexGroups);
}

std::vector<IndexGroup> readIndexGroups(const hid_t baseContainer)
{
    const std::vector<std::string> names = getNamesOfIndexGroups(baseContainer);
    std::vector<IndexGroup>        indexGroups;
    indexGroups.reserve(names.size());
    for (const auto& groupName : names)
    {
        indexGroups.push_back(readIndexGroup(baseContainer, groupName.c_str()));
    }
    return indexGroups;
}

IndexGroup readIndexGroup(const hid_t baseContainer, const char* groupName)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(baseContainer, groupName));

    const bool       hasRanges  = objectExists(group, c_atomRangesDataSetName);
    const bool       hasIndices = objectExists(group, c_simpleIndicesDataSetName);
    std::vector<int> retIndices;

    // The writer gives every index group at least one data set, even an empty group
    GMX_H5MD_THROW_UPON_ERROR(!hasRanges && !hasIndices,
                              formatString("Group %s holds neither the %s nor the %s dataset, so "
                                           "it is not an index group.",
                                           groupName,
                                           c_simpleIndicesDataSetName,
                                           c_atomRangesDataSetName));

    // Read plain indices dataset, 0-based. It is read first, so it can be read straight into
    // the output rather than through an intermediate buffer.
    if (hasIndices)
    {
        const H5mdFixedDataSet<int> dataset(group, c_simpleIndicesDataSetName);
        retIndices.resize(dataset.numValues());
        dataset.readData(retIndices);
        // As for the start of a range, a plain index is an offset and cannot be negative.
        GMX_H5MD_THROW_UPON_ERROR(
                std::any_of(retIndices.begin(), retIndices.end(), [](const int i) { return i < 0; }),
                formatString("Invalid negative index in the %s dataset of group: %s.",
                             c_simpleIndicesDataSetName,
                             groupName));
    }

    // Expand [N, 3] ranges dataset: each row is [start, stop, step], 0-based with exclusive stop
    if (hasRanges)
    {
        const H5mdFixedDataSet<int> dataset(group, c_atomRangesDataSetName);
        // The layout is dictated by the file, so check it before indexing row-wise
        const DataSetDims& rangeDims = dataset.dims();
        GMX_H5MD_THROW_UPON_ERROR(
                rangeDims.size() != 2 || rangeDims[1] != c_atomRangeRowSize,
                formatString("The '%s' dataset of group '%s' must be a [N, %d] dataset.",
                             c_atomRangesDataSetName,
                             groupName,
                             c_atomRangeRowSize));
        std::vector<int> ranges(dataset.numValues());
        dataset.readData(ranges);

        /* Validate every row and total up how many indices they expand to before appending
         * any, so that the output is grown once rather than repeatedly. A single range can
         * cover the whole system, so this is the difference between one allocation and many. */
        size_t numRangeIndices = 0;
        for (Index i = 0; i < gmx::ssize(ranges); i += c_atomRangeRowSize)
        {
            const int start = ranges[i];
            const int stop  = ranges[i + 1];
            const int step  = ranges[i + 2];

            GMX_H5MD_THROW_UPON_ERROR(
                    step <= 0,
                    formatString("Invalid step %d in atom ranges dataset of group: %s.", step, groupName));
            GMX_H5MD_THROW_UPON_ERROR(
                    start < 0,
                    formatString("Invalid negative start %d in atom ranges dataset of group: %s.",
                                 start,
                                 groupName));
            // Prevent the overflow of indices before expanding the range.
            GMX_H5MD_THROW_UPON_ERROR(
                    stop > std::numeric_limits<int>::max() - step + 1,
                    formatString("Invalid stop %d with step %d in atom ranges dataset of group: "
                                 "%s: the range would reach beyond the largest index that can be "
                                 "represented.",
                                 stop,
                                 step,
                                 groupName));

            if (stop > start)
            {
                const int span = stop - start;
                numRangeIndices += gmx::divideRoundUp(span, step);
            }
        }
        retIndices.reserve(retIndices.size() + numRangeIndices);

        for (Index i = 0; i < gmx::ssize(ranges); i += c_atomRangeRowSize)
        {
            const int start = ranges[i];
            const int stop  = ranges[i + 1];
            const int step  = ranges[i + 2];
            for (int j = start; j < stop; j += step)
            {
                retIndices.push_back(j); // already 0-based
            }
        }
    }

    // The group is the set of all indices in both data sets: sort it and keep each index once,
    // whether the file lists it in overlapping ranges, in both data sets or repeatedly.
    std::sort(retIndices.begin(), retIndices.end());
    retIndices.erase(std::unique(retIndices.begin(), retIndices.end()), retIndices.end());
    return { groupName, std::move(retIndices) };
}

} // namespace gmx
