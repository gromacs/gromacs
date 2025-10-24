/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

/*! \brief Declarations of H5md utility functions.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#ifndef GMX_FILEIO_H5MD_UTIL_H
#define GMX_FILEIO_H5MD_UTIL_H

#include <hdf5.h>

#include <cstring>

#include <optional>
#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Dimensions of HDF5 data sets and spaces.
 */
using DataSetDims = std::vector<hsize_t>;

/*! \brief Return whether object \p name exists in the given HDF5 \p container.
 *
 * \param[in] container Handle to container in which to check.
 * \param[in] name      Name of object to look for.
 *
 * \returns True if the object exists, otherwise false.
 */
inline bool objectExists(const hid_t container, const char* name) noexcept
{
    // Return value: <0 = error, 0 = does not exist, >0 = exists
    return H5Oexists_by_name(container, name, H5P_DEFAULT) > 0;
}

/*! \brief Return whether HDF5 \p handle is valid.
 *
 * \param[in] handle HDF5 handle to check.
 *
 * \returns True if the handle is valid, otherwise false.
 */
inline bool handleIsValid(const hid_t handle) noexcept
{
    // Return value: <0 = error, 0 = invalid, >0 = valid
    return H5Iis_valid(handle) > 0;
}

/*! \brief Return the full name (HDF-path) of an HDF5 \p handle.
 *
 * \param[in] handle HDF5 handle to get the name of.
 * \returns The full name of the handle, or \c std::nullopt if the handle is invalid.
 */
inline std::optional<std::string> getHandlePath(const hid_t handle)
{
    int size = static_cast<int>(H5Iget_name(handle, nullptr, 0));

    if (size > 0)
    {
        std::vector<char> buffer(size + 1, '\0');
        H5Iget_name(handle, buffer.data(), size + 1);
        return std::string(buffer.data());
    }
    else
    {
        return std::nullopt;
    }
}

/*! \brief Return the base name (last component of the HDF-path) of an HDF5 \p handle.
 *
 * \param[in] handle HDF5 handle to get the base name of.
 * \returns The base name of the handle, the full name if no '/' found in its paths,
 *          or \c std::nullopt on error.
 */
inline std::optional<std::string> getHandleBaseName(const hid_t handle)
{
    std::optional<std::string> fullName = getHandlePath(handle);
    if (!fullName.has_value())
    {
        return std::nullopt;
    }
    size_t pos = fullName.value().find_last_of('/');
    if (pos != std::string::npos)
    {
        // Return empty string for root group
        return fullName.value().substr(pos + 1);
    }
    else
    {
        // If no '/' found, return the full name
        return fullName;
    }
}


/*! \brief Iterate over a list of objects and calculate the maximum string length via a callback function.
 *
 * \tparam    Iterator The type of the iterator of the source strings.
 * \tparam    Callback The type of the callback function to convert iterator to C-string.
 * \param[in] begin    Iterator to the beginning of the source vector of strings.
 * \param[in] end      Iterator to the end of the source vector of strings.
 * \param[in] cStringFromIterator The callback function to convert iterator to C-string.
 * \returns The maximum length of the strings and the number of strings in the vector.
 */
template<typename Iterator, typename Callback>
std::pair<size_t, size_t> estimateBufferSize(Iterator begin, Iterator end, Callback cStringFromIterator)
{
    size_t maxLength = 0;
    size_t count     = 0;
    for (auto it = begin; it != end; ++it)
    {
        maxLength = std::max(maxLength, std::strlen(cStringFromIterator(it)));
        ++count;
    }
    return std::make_pair(maxLength, count);
}

/*! \brief Pack a string buffer from a range of iterators and a callback function.
 *
 * Caller that knows about their problem domain can pass in a buffer
 * whose *capacity* meets that needed. Otherwise we fall back on the
 * standard vector approach of resizing logarithmically
 *
 * Registering a callback function avoids the H5md interface adding overloads and
 * explicit specializations for multiple string types and multiple collection types
 *
 * \tparam    Iterator      The type of the iterator to write.
 * \tparam    Callback      The type of the callback function to convert iterator to C-string.
 * \param[in] begin         The beginning iterator of the data to write.
 * \param[in] end           The ending iterator of the data to write.
 * \param[in] packingBuffer The buffer to write the strings into.
 * \param[in] maxStrLength  The maximum length of the strings in the attribute.
 * \param[in] cStringFromIterator The callback function to convert iterator to C-string.
 * \returns The packed buffer containing the strings.
 */
template<typename Iterator, typename Callback>
std::vector<char> packBufferViaIterator(Iterator            begin,
                                        Iterator            end,
                                        std::vector<char>&& packingBuffer,
                                        const size_t        maxStrLength,
                                        Callback            cStringFromIterator)
{
    const size_t maxStrLengthWithNull = maxStrLength + 1;
    packingBuffer.resize(0);
    for (auto it = begin; it != end; ++it)
    {
        const char*  cString               = cStringFromIterator(it);
        const size_t numCharactersInString = std::strlen(cString);
        auto         newBufferEnd =
                packingBuffer.insert(packingBuffer.end(), cString, cString + numCharactersInString)
                + numCharactersInString;

        GMX_ASSERT(numCharactersInString <= maxStrLengthWithNull, "Error in string handling");
        packingBuffer.insert(newBufferEnd, maxStrLengthWithNull - numCharactersInString, '\0');
    }
    return packingBuffer;
}

} // namespace gmx

#endif // GMX_FILEIO_H5MD_UTIL_H
