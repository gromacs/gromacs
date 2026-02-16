
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

/*! \brief Definitions of utility functions for working with H5MD attributes.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_attribute.h"

#include <hdf5.h>

#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"


// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

void setStringAttributeByBuffer(const hid_t          container,
                                const char*          attributeName,
                                const size_t         numberOfStrings,
                                const int            maxStrLength,
                                ArrayRef<const char> buffer)
{
    GMX_H5MD_THROW_UPON_ERROR(buffer.size() < numberOfStrings * (maxStrLength + 1),
                              formatString("Buffer size is too small for attribute: %s", attributeName));

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedSizeString(maxStrLength + 1));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(H5Tget_class(dataType) != H5T_STRING,
                              formatString("Data type for attribute is not a string: %s", attributeName));

    DataSetDims dims{ hsize_t(numberOfStrings) };
    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, dims.data(), nullptr));

    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(
            H5Acreate(container, attributeName, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(attribute,
                                    formatString("Failed to create attribute: %s", attributeName));

    GMX_H5MD_THROW_UPON_ERROR(
            H5Awrite(attribute, dataType, buffer.data()) < 0,
            formatString("Failed to write vector of strings attribute: %s", attributeName));
}

template<typename ValueType>
std::optional<ValueType> getAttribute(const hid_t container, const char* attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName, H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<ValueType>(dataType),
                              formatString("Type mismatch when reading attribute: %s", attributeName));

    ValueType value{};
    GMX_H5MD_THROW_UPON_ERROR(H5Aread(attribute, dataType, &value) < 0,
                              formatString("Failed to read attribute: %s", attributeName));
    return value;
}

template<>
std::optional<std::string> getAttribute<std::string>(const hid_t container, const char* attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName, H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<std::string>(dataType),
                              formatString("Type mismatch when reading attribute: %s", attributeName));

    size_t            stringSize = H5Tget_size(dataType);
    std::vector<char> strData(stringSize);
    GMX_H5MD_THROW_UPON_ERROR(H5Aread(attribute, dataType, strData.data()) < 0,
                              formatString("Failed to read string attribute: %s", attributeName));

    std::string values(strData.data());
    return values;
}


template<typename ValueType>
std::optional<std::vector<ValueType>> getAttributeVector(const hid_t container, const char* attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName, H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<ValueType>(dataType),
                              formatString("Type mismatch when reading attribute: %s", attributeName));

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Aget_space(attribute));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataSpace, formatString("Failed to get data space for attribute: %s", attributeName));

    // Setup the size of the vector
    DataSetDims dims = { 0 };
    GMX_H5MD_THROW_UPON_ERROR(H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr) < 0,
                              formatString("Failed to get dimensions for attribute: %s", attributeName));
    const size_t nelems = dims[0];
    if (nelems == 0)
    {
        return std::vector<ValueType>{};
    }

    // Read the data
    std::vector<ValueType> values(dims[0]);
    GMX_H5MD_THROW_UPON_ERROR(H5Aread(attribute, dataType, values.data()) < 0,
                              formatString("Failed to read vector attribute: %s", attributeName));
    return values;
}

template<>
std::optional<std::vector<std::string>> getAttributeVector<std::string>(const hid_t container,
                                                                        const char* attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName, H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<std::string>(dataType),
                              formatString("Type mismatch when reading attribute: %s", attributeName));

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Aget_space(attribute));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataSpace, formatString("Failed to get data space for attribute: %s", attributeName));

    // Setup the size of the vector
    DataSetDims dims = { 0 };
    GMX_H5MD_THROW_UPON_ERROR(H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr) < 0,
                              formatString("Failed to get dimensions for attribute: %s", attributeName));
    const size_t nelems = dims[0];
    if (nelems == 0)
    {
        return std::vector<std::string>{};
    }

    // Set up a buffer and read the data
    size_t                   stringSize = H5Tget_size(dataType);
    std::vector<std::string> values(nelems);
    std::vector<char>        buffer(nelems * stringSize);

    GMX_H5MD_THROW_UPON_ERROR(
            H5Aread(attribute, dataType, buffer.data()) < 0,
            formatString("Failed to read vector of strings attribute: %s", attributeName));
    for (size_t i = 0; i < nelems; i++)
    {
        values[i] = std::string(buffer.data() + i * stringSize,
                                strnlen(buffer.data() + i * stringSize, stringSize));
    }
    return values;
}

template<typename ValueType>
void setAttribute(const hid_t container, const char* attributeName, const ValueType& value)
{
    // Initialize the data space with a scalar type
    auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataSpace, formatString("Failed to create data space for attribute: %s", attributeName));

    // NOTE: Throw if the attribute already exists (!5205)
    GMX_H5MD_THROW_UPON_ERROR(H5Aexists(container, attributeName) > 0,
                              formatString("Attribute already exists: %s", attributeName));

    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate(
            container, attributeName, hdf5DataTypeFor<ValueType>(), dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(attribute,
                                    formatString("Failed to create attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(H5Awrite(attribute, hdf5DataTypeFor<ValueType>(), &value) < 0,
                              formatString("Failed to write attribute: %s", attributeName));
}

void setAttribute(const hid_t container, const char* attributeName, const char* value)
{
    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedSizeString(strlen(value) + 1));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataType, formatString("Failed to get data type for attribute: %s", attributeName));
    GMX_H5MD_THROW_UPON_ERROR(H5Tget_class(dataType) != H5T_STRING,
                              formatString("Data type for attribute is not a string: %s", attributeName));

    // Initialize the data space with a scalar type
    auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataSpace, formatString("Failed to create data space for attribute: %s", attributeName));

    // NOTE: Throw if the attribute already exists (!5205)
    GMX_H5MD_THROW_UPON_ERROR(H5Aexists(container, attributeName) > 0,
                              formatString("Attribute already exists: %s", attributeName));

    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(
            H5Acreate2(container, attributeName, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(attribute,
                                    formatString("Failed to create attribute: %s", attributeName));

    // Write the attribute
    GMX_H5MD_THROW_UPON_ERROR(H5Awrite(attribute, dataType, value) < 0,
                              formatString("Failed to write string attribute: %s", attributeName));
}

void setAttribute(const hid_t container, const char* attributeName, const std::string& value)
{
    setAttribute(container, attributeName, value.c_str());
}

template<typename ValueType>
void setAttributeVector(const hid_t container, const char* attributeName, ArrayRef<const ValueType> values)
{
    DataSetDims dims{ static_cast<hsize_t>(values.size()) };
    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, dims.data(), nullptr));

    // NOTE: Throw if the attribute already exists (!5205)
    GMX_H5MD_THROW_UPON_ERROR(H5Aexists(container, attributeName) > 0,
                              formatString("Attribute already exists: %s", attributeName));
    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate(
            container, attributeName, hdf5DataTypeFor<ValueType>(), dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    GMX_H5MD_THROW_UPON_INVALID_HID(attribute,
                                    formatString("Failed to create attribute: %s", attributeName));

    // Vector of numerical values
    GMX_H5MD_THROW_UPON_ERROR(H5Awrite(attribute, hdf5DataTypeFor<ValueType>(), values.data()) < 0,
                              formatString("Failed to write vector attribute: %s", attributeName));
}

void setAttributeVector(const hid_t container, const char* attributeName, ArrayRef<const std::string> values)
{
    std::vector<char> buffer;
    setAttributeStringVector(container, attributeName, std::move(buffer), values.begin(), values.end());
}

void setAttributeVector(const hid_t container, const char* attributeName, ArrayRef<const char* const> values)
{
    std::vector<char> buffer;
    setAttributeStringVector(container, attributeName, std::move(buffer), values.begin(), values.end());
}

template<typename Iterator>
std::vector<char> setAttributeStringVector(const hid_t         container,
                                           const char*         attributeName,
                                           std::vector<char>&& buffer,
                                           Iterator            begin,
                                           Iterator            end)
{
    size_t maxStrLength = 0;
    size_t strCount     = 0;
    if constexpr (std::is_same_v<typename std::iterator_traits<Iterator>::value_type, const char*>)
    {
        std::tie(maxStrLength, strCount) =
                estimateBufferSize(begin, end, [](const auto& str) { return *str; });
    }
    else if constexpr (std::is_same_v<typename std::iterator_traits<Iterator>::value_type, std::string>)
    {
        std::tie(maxStrLength, strCount) =
                estimateBufferSize(begin, end, [](const auto& str) { return str->c_str(); });
    }
    else
    {
        throw FileIOError(formatString("Unsupported string type for attribute: %s", attributeName));
    }
    size_t expectedSize = (maxStrLength + 1) * strCount;
    GMX_H5MD_THROW_UPON_ERROR(
            expectedSize == 0,
            formatString("Cannot write empty string vector attribute: %s", attributeName));
    if (buffer.empty() || buffer.size() < expectedSize)
    {
        // Only enlarge the buffer when empty or the current size is smaller than needed
        buffer.reserve(expectedSize);
    }

    if constexpr (std::is_same_v<typename std::iterator_traits<Iterator>::value_type, const char*>)
    {
        buffer = packBufferViaIterator(
                begin, end, std::move(buffer), maxStrLength, [](const auto& str) { return *str; });
    }
    else if constexpr (std::is_same_v<typename std::iterator_traits<Iterator>::value_type, std::string>)
    {
        buffer = packBufferViaIterator(begin,
                                       end,
                                       std::move(buffer),
                                       maxStrLength,
                                       [](const auto& str) { return str->c_str(); });
    }
    else
    {
        throw FileIOError(formatString("Unsupported string type for attribute: %s", attributeName));
    }
    setStringAttributeByBuffer(container, attributeName, strCount, maxStrLength, buffer);
    return std::move(buffer);
}

/// @cond DO_NOT_DOCUMENT
template std::optional<int32_t>  getAttribute(const hid_t, const char*);
template std::optional<int64_t>  getAttribute(const hid_t, const char*);
template std::optional<uint32_t> getAttribute(const hid_t, const char*);
template std::optional<uint64_t> getAttribute(const hid_t, const char*);
template std::optional<float>    getAttribute(const hid_t, const char*);
template std::optional<double>   getAttribute(const hid_t, const char*);

template void setAttribute(const hid_t, const char*, const int32_t&);
template void setAttribute(const hid_t, const char*, const int64_t&);
template void setAttribute(const hid_t, const char*, const uint32_t&);
template void setAttribute(const hid_t, const char*, const uint64_t&);
template void setAttribute(const hid_t, const char*, const float&);
template void setAttribute(const hid_t, const char*, const double&);

template std::optional<std::vector<int32_t>>  getAttributeVector(const hid_t, const char*);
template std::optional<std::vector<int64_t>>  getAttributeVector(const hid_t, const char*);
template std::optional<std::vector<uint32_t>> getAttributeVector(const hid_t, const char*);
template std::optional<std::vector<uint64_t>> getAttributeVector(const hid_t, const char*);
template std::optional<std::vector<float>>    getAttributeVector(const hid_t, const char*);
template std::optional<std::vector<double>>   getAttributeVector(const hid_t, const char*);

template void setAttributeVector(const hid_t, const char*, ArrayRef<const int32_t>);
template void setAttributeVector(const hid_t, const char*, ArrayRef<const int64_t>);
template void setAttributeVector(const hid_t, const char*, ArrayRef<const uint32_t>);
template void setAttributeVector(const hid_t, const char*, ArrayRef<const uint64_t>);
template void setAttributeVector(const hid_t, const char*, ArrayRef<const float>);
template void setAttributeVector(const hid_t, const char*, ArrayRef<const double>);

template std::vector<char> setAttributeStringVector(const hid_t,
                                                    const char*,
                                                    std::vector<char>&&,
                                                    ArrayRef<const std::string>::const_iterator,
                                                    ArrayRef<const std::string>::const_iterator);
template std::vector<char> setAttributeStringVector(const hid_t,
                                                    const char*,
                                                    std::vector<char>&&,
                                                    ArrayRef<const char*>::const_iterator,
                                                    ArrayRef<const char*>::const_iterator);
/// @endcond

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
