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
 * \brief Defines the writing of a KeyValueTree to H5MD.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_kvtree_writer.h"

#include <cstdint>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! \brief Reserved key of a property whose value is written as the HDF5 data set.
constexpr char c_valueKey[] = "value";

//! \brief Reserved key of a property that holds the dimensions of its data set.
constexpr char c_shapeKey[] = "shape";

//! \brief Reserved keys of a property that are written as HDF5 attributes.
constexpr std::array<const char*, 5> c_attributeKeys{ "unit",
                                                      "description",
                                                      "doi",
                                                      "url",
                                                      "version" };

//! \brief Create the scalar data set \p name of \p ValueType in \p container.
template<typename ValueType>
hid_t createScalarDataSet(const hid_t container, const char* name)
{
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    const hid_t dataSet                    = H5Dcreate(
            container, name, hdf5DataTypeFor<ValueType>(), dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    GMX_H5MD_THROW_UPON_INVALID_HID(dataSet,
                                    formatString("Error creating scalar data set `%s/%s`",
                                                 getHandlePath(container).value_or("<unknown>").c_str(),
                                                 name));

    return dataSet;
}

//! \brief Create the scalar data set \p name in \p container for a string of up to \p maxSize bytes.
hid_t createScalarFixedSizeStringDataSet(const hid_t container, const char* name, const hsize_t maxSize)
{
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(hdf5DataTypeForFixedLengthString(maxSize));
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Screate(H5S_SCALAR));
    const hid_t dataSet =
            H5Dcreate(container, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    GMX_H5MD_THROW_UPON_INVALID_HID(
            dataSet,
            formatString("Error creating scalar fixed-size string data set `%s/%s`",
                         getHandlePath(container).value_or("<unknown>").c_str(),
                         name));

    return dataSet;
}

//! \brief Write \p value into the scalar \p dataSet.
template<typename ValueType>
void writeValueToScalarDataSet(const hid_t dataSet, const ValueType& value)
{
    GMX_H5MD_THROW_UPON_ERROR(
            H5Dwrite(dataSet, hdf5DataTypeFor<ValueType>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0,
            formatString("Cannot write value into scalar data set `%s`",
                         getHandlePath(dataSet).value_or("<unknown>").c_str()));
}

//! \brief Write \p value as the scalar data set \p name of \p ValueType in \p container.
template<typename ValueType>
void writeScalarDataSet(const hid_t container, const std::string& name, const ValueType& value)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(createScalarDataSet<ValueType>(container, name.c_str()));
    writeValueToScalarDataSet(dataSet, value);
}

//! \brief Write \p value into the scalar fixed-size string \p dataSet.
void writeFixedSizeStringToScalarDataSet(const hid_t dataSet, const std::string& value)
{
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
    GMX_H5MD_THROW_UPON_ERROR(H5Dwrite(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.data()) < 0,
                              formatString("Cannot write string into scalar data set `%s`",
                                           getHandlePath(dataSet).value_or("<unknown>").c_str()));
}

/*! \brief Return the elements of \p array as \p ValueType.
 *
 * \throws H5mdError If an element is not of \p ValueType.
 */
template<typename ValueType>
std::vector<ValueType> arrayValuesAs(const KeyValueTreeArray& array)
{
    std::vector<ValueType> values;
    values.reserve(array.values().size());
    for (const KeyValueTreeValue& value : array.values())
    {
        GMX_H5MD_THROW_UPON_ERROR(!value.isType<ValueType>(),
                                  "The elements of an array must all be of one type");
        values.push_back(value.cast<ValueType>());
    }
    return values;
}

//! \brief Write \p values as the data set \p name with dimensions \p dims in \p container.
template<typename ValueType>
void writeArrayDataSet(const hid_t                   container,
                       const std::string&            name,
                       const DataSetDims&            dims,
                       const std::vector<ValueType>& values)
{
    H5mdFixedDataSet<ValueType> dataSet(
            H5mdDataSetBuilder<ValueType>(container, name).withDimension(dims).build());
    dataSet.writeData(makeConstArrayRef(values));
}

/*! \brief Return the dimensions of the data set of \p property, read from its "shape".
 *
 * Without a "shape", no dimensions are returned.
 *
 * \throws H5mdError If "shape" is not an array of positive integers.
 */
DataSetDims shapeOf(const KeyValueTreeObject& property, const std::string& key)
{
    if (!property.keyExists(c_shapeKey))
    {
        return {};
    }
    const KeyValueTreeValue& shape = property[c_shapeKey];
    GMX_H5MD_THROW_UPON_ERROR(!shape.isArray(),
                              formatString("The shape of '%s' must be an array", key.c_str()));

    DataSetDims dims;
    for (const KeyValueTreeValue& dimension : shape.asArray().values())
    {
        const int64_t size = dimension.isType<int>()       ? dimension.cast<int>()
                             : dimension.isType<int64_t>() ? dimension.cast<int64_t>()
                                                           : 0;
        GMX_H5MD_THROW_UPON_ERROR(
                size <= 0, formatString("The shape of '%s' must hold positive integers", key.c_str()));
        dims.push_back(static_cast<hsize_t>(size));
    }
    return dims;
}

/*! \brief Write a KeyValueTree value to HDF5 as a dataset.
 *
 * Supported scalar types: bool, int, int64_t, float, double, std::string.
 * Supported array element types: std::string, int, int64_t, float, double.
 * bool is stored as int (0/1) since HDF5 has no native bool type.
 * An empty array creates a zero-length string dataset.
 * Unrecognized types are silently skipped.
 *
 * \param[in] container   HDF5 group handle in which to create the dataset
 * \param[in] datasetName Name of the dataset to create
 * \param[in] value       Value from KeyValueTree
 * \param[in] shape       Dimensions of an array \p value, row-major; empty for one dimension
 *
 * \throws H5mdError If \p shape is given for a value that is not an array, if the product of
 *                   \p shape differs from the number of values, if an array mixes element
 *                   types, or if the data set cannot be created or written.
 */
void writeValueAsDataset(hid_t                    container,
                         const std::string&       datasetName,
                         const KeyValueTreeValue& value,
                         ArrayRef<const hsize_t>  shape)
{
    GMX_H5MD_THROW_UPON_ERROR(
            !shape.empty() && !value.isArray(),
            formatString("Only an array can be given a shape, which '%s' is not", datasetName.c_str()));

    if (value.isType<bool>())
    {
        // HDF5 has no native bool type; store as int (0/1)
        writeScalarDataSet(container, datasetName, static_cast<int>(value.cast<bool>()));
    }
    else if (value.isType<int>())
    {
        writeScalarDataSet(container, datasetName, value.cast<int>());
    }
    else if (value.isType<int64_t>())
    {
        writeScalarDataSet(container, datasetName, value.cast<int64_t>());
    }
    else if (value.isType<double>())
    {
        writeScalarDataSet(container, datasetName, value.cast<double>());
    }
    else if (value.isType<float>())
    {
        writeScalarDataSet(container, datasetName, value.cast<float>());
    }
    else if (value.isType<std::string>())
    {
        const std::string& str = value.cast<std::string>();
        // +1 for the null terminator
        const hsize_t maxLen = str.length() + 1;

        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                createScalarFixedSizeStringDataSet(container, datasetName.c_str(), maxLen));
        writeFixedSizeStringToScalarDataSet(dataSet, str);
    }
    else if (value.isArray())
    {
        const auto&       arr  = value.asArray();
        const DataSetDims dims = shape.empty() ? DataSetDims{ arr.values().size() }
                                               : DataSetDims(shape.begin(), shape.end());
        const hsize_t     numValues =
                std::accumulate(dims.begin(), dims.end(), hsize_t{ 1 }, std::multiplies<>());
        GMX_H5MD_THROW_UPON_ERROR(
                numValues != arr.values().size(),
                formatString("The shape of '%s' holds %llu values, but the array holds %zu",
                             datasetName.c_str(),
                             static_cast<unsigned long long>(numValues),
                             arr.values().size()));

        if (arr.values().empty())
        {
            // Element type is unknown; create a zero-length string dataset as a placeholder
            H5mdDataSetBuilder<std::string>(container, datasetName)
                    .withDimension({ 0 })
                    .withMaxStringLength(1)
                    .build();
        }
        else if (arr.values()[0].isType<std::string>())
        {
            const std::vector<std::string> stringVec = arrayValuesAs<std::string>(arr);

            // +1 for the null terminator; initialise to 1 so the buffer is never zero-length
            size_t maxLen = 1;
            for (const auto& str : stringVec)
            {
                maxLen = std::max(maxLen, str.length() + 1);
            }

            H5mdFixedDataSet<std::string> dataset(H5mdDataSetBuilder<std::string>(container, datasetName)
                                                          .withDimension(dims)
                                                          .withMaxStringLength(maxLen)
                                                          .build());
            dataset.writeData(makeConstArrayRef(stringVec));
        }
        else if (arr.values()[0].isType<int>())
        {
            writeArrayDataSet(container, datasetName, dims, arrayValuesAs<int>(arr));
        }
        else if (arr.values()[0].isType<int64_t>())
        {
            writeArrayDataSet(container, datasetName, dims, arrayValuesAs<int64_t>(arr));
        }
        else if (arr.values()[0].isType<float>())
        {
            writeArrayDataSet(container, datasetName, dims, arrayValuesAs<float>(arr));
        }
        else if (arr.values()[0].isType<double>())
        {
            writeArrayDataSet(container, datasetName, dims, arrayValuesAs<double>(arr));
        }
    }
}

/*! \brief Write a KeyValueTree value to HDF5 as an attribute.
 *
 * Supported scalar types: bool, int, int64_t, float, double, std::string.
 * Supported array element types: std::string, int, int64_t, float, double.
 * bool is stored as int (0/1) since HDF5 has no native bool type.
 * An empty array and unrecognized types are silently skipped.
 *
 * \param[in] container HDF5 group or dataset handle on which to set the attribute
 * \param[in] attrName  Name of the attribute to create
 * \param[in] value     Value from KeyValueTree
 *
 * \throws H5mdError If an array mixes element types, or if the attribute cannot be created or
 *                   written.
 */
void writeValueAsAttribute(hid_t container, const std::string& attrName, const KeyValueTreeValue& value)
{
    if (value.isType<bool>())
    {
        // HDF5 has no native bool type; store as int (0/1)
        setAttribute(container, attrName.c_str(), static_cast<int>(value.cast<bool>()));
    }
    else if (value.isType<int>())
    {
        setAttribute(container, attrName.c_str(), value.cast<int>());
    }
    else if (value.isType<int64_t>())
    {
        setAttribute(container, attrName.c_str(), value.cast<int64_t>());
    }
    else if (value.isType<double>())
    {
        setAttribute(container, attrName.c_str(), value.cast<double>());
    }
    else if (value.isType<float>())
    {
        setAttribute(container, attrName.c_str(), value.cast<float>());
    }
    else if (value.isType<std::string>())
    {
        setAttribute(container, attrName.c_str(), value.cast<std::string>());
    }
    else if (value.isArray())
    {
        const auto& arr = value.asArray();
        if (arr.values().empty())
        {
            // HDF5 attributes cannot hold an empty array of unknown type, so nothing is written
            return;
        }
        else if (arr.values()[0].isType<std::string>())
        {
            setAttributeVector(
                    container, attrName.c_str(), makeConstArrayRef(arrayValuesAs<std::string>(arr)));
        }
        else if (arr.values()[0].isType<int>())
        {
            setAttributeVector(container, attrName.c_str(), makeConstArrayRef(arrayValuesAs<int>(arr)));
        }
        else if (arr.values()[0].isType<int64_t>())
        {
            setAttributeVector(
                    container, attrName.c_str(), makeConstArrayRef(arrayValuesAs<int64_t>(arr)));
        }
        else if (arr.values()[0].isType<float>())
        {
            setAttributeVector(container, attrName.c_str(), makeConstArrayRef(arrayValuesAs<float>(arr)));
        }
        else if (arr.values()[0].isType<double>())
        {
            setAttributeVector(container, attrName.c_str(), makeConstArrayRef(arrayValuesAs<double>(arr)));
        }
    }
}

/*! \brief Write a single KeyValueTree property to H5MD.
 *
 * The property is written as set out in the file documentation.
 *
 * \param[in] container HDF5 group handle
 * \param[in] key       Property key / name of the HDF5 dataset or group to create
 * \param[in] value     Property value from KeyValueTree
 *
 * \throws H5mdError If an object holds a key that is not allowed, if an object with a "value"
 *                   holds a nested object, if its "value" is of a type
 *                   that is not supported, if an array mixes element types, if its "shape" is
 *                   invalid or does not match its "value", or if writing to HDF5 fails.
 */
void writeKeyValueTreePropertyToH5md(hid_t container, const std::string& key, const KeyValueTreeValue& value)
{
    if (value.isObject())
    {
        const auto& obj           = value.asObject();
        const bool  isValueObject = obj.keyExists(c_valueKey);

        if (isValueObject)
        {
            writeValueAsDataset(container, key, obj[c_valueKey], shapeOf(obj, key));
        }
        else
        {
            makeH5mdGroupGuard(createGroup(container, key.c_str()));
        }
        // Nothing is created for a "value" of an unsupported type, so the object may not exist
        const hid_t object = H5Oopen(container, key.c_str(), H5P_DEFAULT);
        GMX_H5MD_THROW_UPON_INVALID_HID(
                object,
                formatString("Cannot open `%s/%s`; its \"value\" may be of an unsupported type",
                             getHandlePath(container).value_or("<unknown>").c_str(),
                             key.c_str()));
        const auto [objectHandle, objGuard] = makeH5mdGuard(object, H5Oclose);

        for (const auto& prop : obj.properties())
        {
            const std::string& propKey = prop.key();
            const auto&        propVal = prop.value();

            if (isValueObject && (propKey == c_valueKey || propKey == c_shapeKey))
            {
                // Already written as the data set
                continue;
            }
            else if (propVal.isObject())
            {
                GMX_H5MD_THROW_UPON_ERROR(
                        isValueObject,
                        formatString("Property `%s` holds a value, so it cannot hold the nested "
                                     "object '%s'",
                                     getHandlePath(objectHandle).value_or(key).c_str(),
                                     propKey.c_str()));
                writeKeyValueTreePropertyToH5md(objectHandle, propKey, propVal);
            }
            else if (std::find(c_attributeKeys.begin(), c_attributeKeys.end(), propKey)
                     != c_attributeKeys.end())
            {
                // Written on the data set or group that holds the property
                writeValueAsAttribute(objectHandle, propKey, propVal);
            }
            else
            {
                GMX_THROW(H5mdError(formatString(
                        "Key '%s' of `%s` is not allowed: besides nested objects, a property can "
                        "only hold '%s', '%s' and the attribute keys %s",
                        propKey.c_str(),
                        getHandlePath(objectHandle).value_or(key).c_str(),
                        c_valueKey,
                        c_shapeKey,
                        joinStrings(c_attributeKeys.begin(), c_attributeKeys.end(), ", ").c_str())));
            }
        }
    }
}

} // namespace

void writeKeyValueTreeToH5md(hid_t container, const KeyValueTreeObject& tree)
{
    for (const auto& prop : tree.properties())
    {
        writeKeyValueTreePropertyToH5md(container, prop.key(), prop.value());
    }
}

} // namespace gmx
