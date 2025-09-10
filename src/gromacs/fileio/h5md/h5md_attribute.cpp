
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

#include "h5md_attribute.h"

#include <hdf5.h>

#include <cstring>

#include "gromacs/utility/basedefinitions.h"

#include "h5md_dataset.h"
#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_type.h"
#include "h5md_util.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

template<typename ValueType>
std::optional<ValueType> getAttribute(const hid_t container, const std::string& attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName.c_str(), H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(!valueTypeIsDataType<ValueType>(dataType),
                       "Type mismatch when reading attribute: " + attributeName);

    ValueType value{};
    throwUponH5mdError(H5Aread(attribute, dataType, &value) < 0,
                       "Failed to read attribute: " + attributeName);
    return value;
}

template<>
std::optional<std::string> getAttribute<std::string>(const hid_t container, const std::string& attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName.c_str(), H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(!valueTypeIsDataType<std::string>(dataType),
                       "Type mismatch when reading attribute: " + attributeName);

    size_t            stringSize = H5Tget_size(dataType);
    std::vector<char> strData(stringSize);
    throwUponH5mdError(H5Aread(attribute, dataType, strData.data()) < 0,
                       "Failed to read string attribute: " + attributeName);

    std::string values(strData.data());
    return values;
}


template<typename ValueType>
std::optional<std::vector<ValueType>> getAttributeVector(const hid_t container, const std::string& attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName.c_str(), H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(!valueTypeIsDataType<ValueType>(dataType),
                       "Type mismatch when reading attribute: " + attributeName);

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Aget_space(attribute));
    throwUponInvalidHid(dataSpace, "Failed to get data space for attribute: " + attributeName);

    // Setup the size of the vector
    DataSetDims dims = { 0 };
    throwUponH5mdError(H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr) < 0,
                       "Failed to get dimensions for attribute: " + attributeName);
    std::vector<ValueType> values(dims[0]);

    // Read the data
    throwUponH5mdError(H5Aread(attribute, dataType, values.data()) < 0,
                       "Failed to read vector attribute: " + attributeName);

    return values;
}

template<>
std::optional<std::vector<std::string>> getAttributeVector<std::string>(const hid_t container,
                                                                        const std::string& attributeName)
{
    const auto [attribute, attributeGuard] =
            makeH5mdAttributeGuard(H5Aopen(container, attributeName.c_str(), H5P_DEFAULT));
    if (!handleIsValid(attribute))
    {
        return std::nullopt;
    }

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_native_type(H5Aget_type(attribute), H5T_DIR_DEFAULT));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(!valueTypeIsDataType<std::string>(dataType),
                       "Type mismatch when reading attribute: " + attributeName);

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Aget_space(attribute));
    throwUponInvalidHid(dataSpace, "Failed to get data space for attribute: " + attributeName);

    // Setup the size of the vector
    DataSetDims dims = { 0 };
    throwUponH5mdError(H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr) < 0,
                       "Failed to get dimensions for attribute: " + attributeName);
    const size_t nelems = dims[0];

    // Set up a buffer and read the data
    size_t                   stringSize = H5Tget_size(dataType);
    std::vector<std::string> values(nelems);
    std::vector<char>        buffer(nelems * stringSize);

    throwUponH5mdError(H5Aread(attribute, dataType, buffer.data()) < 0,
                       "Failed to read vector of strings attribute: " + attributeName);
    for (size_t i = 0; i < nelems; i++)
    {
        values[i] = std::string(buffer.data() + i * stringSize,
                                strnlen(buffer.data() + i * stringSize, stringSize));
    }
    return values;
}

template<typename ValueType>
void setAttribute(const hid_t container, const std::string& attributeName, const ValueType& value)
{
    // Initialize the data space with a scalar type
    auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    throwUponInvalidHid(dataSpace, "Failed to create data space for attribute: " + attributeName);

    // NOTE: Throw if the attribute already exists (!5205)
    throwUponH5mdError(H5Aexists(container, attributeName.c_str()) > 0,
                       "Attribute already exists: " + attributeName);

    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate(
            container, attributeName.c_str(), hdf5DataTypeFor<ValueType>(), dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    throwUponInvalidHid(attribute, "Failed to create attribute: " + attributeName);
    throwUponH5mdError(H5Awrite(attribute, hdf5DataTypeFor<ValueType>(), &value) < 0,
                       "Failed to write attribute: " + attributeName);
}

template<>
void setAttribute<const char*>(const hid_t container, const std::string& attributeName, const char* const& value)
{
    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedSizeString(strlen(value) + 1));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(H5Tget_class(dataType) != H5T_STRING,
                       "Data type for attribute is not a string: " + attributeName);

    // Initialize the data space with a scalar type
    auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    throwUponInvalidHid(dataSpace, "Failed to create data space for attribute: " + attributeName);

    // NOTE: Throw if the attribute already exists (!5205)
    throwUponH5mdError(H5Aexists(container, attributeName.c_str()) > 0,
                       "Attribute already exists: " + attributeName);

    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate2(
            container, attributeName.c_str(), dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    throwUponInvalidHid(attribute, "Failed to create attribute: " + attributeName);

    // Write the attribute
    throwUponH5mdError(H5Awrite(attribute, dataType, value) < 0,
                       "Failed to write string attribute: " + attributeName);
}

template<>
void setAttribute<std::string>(const hid_t container, const std::string& attributeName, const std::string& value)
{
    setAttribute<const char*>(container, attributeName, value.c_str());
}


template<typename ValueType>
void setAttributeVector(const hid_t                   container,
                        const std::string&            attributeName,
                        const std::vector<ValueType>& values)
{
    DataSetDims dims{ static_cast<hsize_t>(values.size()) };
    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, dims.data(), nullptr));

    // NOTE: Throw if the attribute already exists (!5205)
    throwUponH5mdError(H5Aexists(container, attributeName.c_str()) > 0,
                       "Attribute already exists: " + attributeName);
    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate(
            container, attributeName.c_str(), hdf5DataTypeFor<ValueType>(), dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    throwUponInvalidHid(attribute, "Failed to create attribute: " + attributeName);

    // Vector of numerical values
    throwUponH5mdError(H5Awrite(attribute, hdf5DataTypeFor<ValueType>(), values.data()) < 0,
                       "Failed to write vector attribute: " + attributeName);
}

template<>
void setAttributeVector<const char*>(const hid_t                     container,
                                     const std::string&              attributeName,
                                     const std::vector<const char*>& values)
{
    // Get the maximum string length from the vector of strings
    const int maxStringLength = values.empty()
                                        ? 1
                                        : strlen(*std::max_element(values.begin(),
                                                                   values.end(),
                                                                   [](const char* a, const char* b)
                                                                   { return strlen(a) < strlen(b); }))
                                                  + 1;

    const auto [dataType, dataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedSizeString(maxStringLength));
    throwUponInvalidHid(dataType, "Failed to get data type for attribute: " + attributeName);
    throwUponH5mdError(H5Tget_class(dataType) != H5T_STRING,
                       "Data type for attribute is not a string: " + attributeName);

    DataSetDims dims{ static_cast<hsize_t>(values.size()) };
    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, dims.data(), nullptr));

    // NOTE: Throw if the attribute already exists (!5205)
    throwUponH5mdError(H5Aexists(container, attributeName.c_str()) > 0,
                       "Attribute already exists: " + attributeName);
    const auto [attribute, attributeGuard] = makeH5mdAttributeGuard(H5Acreate(
            container, attributeName.c_str(), dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT));
    throwUponInvalidHid(attribute, "Failed to create attribute: " + attributeName);

    // Copy the strings into the buffer
    std::vector<char> buffer(values.size() * maxStringLength);
    for (size_t i = 0; i < values.size(); i++)
    {
        std::strncpy(buffer.data() + i * maxStringLength, values[i] ? values[i] : "", maxStringLength);
    }
    throwUponH5mdError(H5Awrite(attribute, dataType, buffer.data()) < 0,
                       "Failed to write vector of strings attribute: " + attributeName);
}

template<>
void setAttributeVector<std::string>(const hid_t                     container,
                                     const std::string&              attributeName,
                                     const std::vector<std::string>& values)
{
    // Convert the vector of strings to a vector of const char* for HDF5 compatibility
    std::vector<const char*> strData(values.size());
    for (size_t i = 0; i < values.size(); ++i)
    {
        strData[i] = values[i].c_str();
    }
    setAttributeVector<const char*>(container, attributeName, strData);
}

template std::optional<int32_t> getAttribute(const hid_t container, const std::string& attributeName);
template std::optional<int64_t> getAttribute(const hid_t container, const std::string& attributeName);
template std::optional<uint32_t> getAttribute(const hid_t container, const std::string& attributeName);
template std::optional<uint64_t> getAttribute(const hid_t container, const std::string& attributeName);
template std::optional<float> getAttribute(const hid_t container, const std::string& attributeName);
template std::optional<double> getAttribute(const hid_t container, const std::string& attributeName);

template std::optional<std::vector<int32_t>>  getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);
template std::optional<std::vector<int64_t>>  getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);
template std::optional<std::vector<uint32_t>> getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);
template std::optional<std::vector<uint64_t>> getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);
template std::optional<std::vector<float>>    getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);
template std::optional<std::vector<double>>   getAttributeVector(const hid_t        container,
                                                                 const std::string& attributeName);

template void setAttribute(const hid_t container, const std::string& attributeName, const int32_t& value);
template void setAttribute(const hid_t container, const std::string& attributeName, const int64_t& value);
template void setAttribute(const hid_t container, const std::string& attributeName, const uint32_t& value);
template void setAttribute(const hid_t container, const std::string& attributeName, const uint64_t& value);
template void setAttribute(const hid_t container, const std::string& attributeName, const float& value);
template void setAttribute(const hid_t container, const std::string& attributeName, const double& value);

template void setAttributeVector(const hid_t                 container,
                                 const std::string&          attributeName,
                                 const std::vector<int32_t>& values);
template void setAttributeVector(const hid_t                 container,
                                 const std::string&          attributeName,
                                 const std::vector<int64_t>& values);
template void setAttributeVector(const hid_t                  container,
                                 const std::string&           attributeName,
                                 const std::vector<uint32_t>& values);
template void setAttributeVector(const hid_t                  container,
                                 const std::string&           attributeName,
                                 const std::vector<uint64_t>& values);
template void setAttributeVector(const hid_t               container,
                                 const std::string&        attributeName,
                                 const std::vector<float>& values);
template void setAttributeVector(const hid_t                container,
                                 const std::string&         attributeName,
                                 const std::vector<double>& values);


} // namespace gmx

CLANG_DIAGNOSTIC_RESET
