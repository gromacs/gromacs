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

/*! \brief Definitions of H5md data set base class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_datasetbase.h"

#include <algorithm>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"


// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

namespace
{

//! \brief Return the number of dimensions of data set with \p dataSetHandle or 0 if it is invalid.
static hsize_t getNumDims(const hid_t dataSetHandle) noexcept
{
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSetHandle));
    return handleIsValid(dataSpace) ? H5Sget_simple_extent_ndims(dataSpace) : 0;
}

//! \brief Helper function to open a data set and throw with a unique error message if it fails.
static hid_t openDataSet(const hid_t container, const char* name)
{
    const hid_t handle = H5Dopen(container, name, H5P_DEFAULT);
    GMX_H5MD_THROW_UPON_INVALID_HID(handle, formatString("Cannot open data set with name %s.", name));
    return handle;
}

/*! \brief Verify that the data set is consistent for its templated \c ValueType.
 *
 * Verifies:
 * 1. That the base primitive of ValueType matches the internal HDF5 data set \p storedDataType.
 * 2. That the primitive data set dimensions \p dims can store values of the templated \c ValueType.
 *
 * \param[in] storedDataType Data type of data set.
 * \param[in] dims           Primitive dimensions of data set.
 *
 * \throws H5mdError if there is an inconsistency.
 */
template<typename ValueType>
static void verifyDataSetConsistency(const hid_t storedDataType, const DataSetDims& dims)
{
    if constexpr (std::is_same_v<ValueType, BasicVector<float>>)
    {
        GMX_H5MD_THROW_UPON_ERROR(dims.empty() || dims.back() != DIM,
                                  "Could not open data set: inner dimension of data set must = 3 "
                                  "for BasicVector<float>");
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<float>(storedDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
    else if constexpr (std::is_same_v<ValueType, BasicVector<double>>)
    {
        GMX_H5MD_THROW_UPON_ERROR(dims.empty() || dims.back() != DIM,
                                  "Could not open data set: inner dimension of data set must = 3 "
                                  "for BasicVector<double>");
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<double>(storedDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
    else
    {
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<ValueType>(storedDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
}

//! \brief Return whether the number of data points in \p values matches the selections.
//
// Gets the number of points in the hyperslab selection of \p memoryDataSpace
// and \p fileDataSpace, then checks this against the number of \p values.
//
// \tparam ValueType Compiled type for values.
// \param[in] values Buffer to check the size of
// \param[in] memoryDataSpace Selection of space in memory to check against
// \param[in] fileDataSpace   Selection of space in HDF5 file to check against.
// \param[in] dataSet         Data set for the operation.
//
// If \p memoryDataSpace or \p fileDataSpace is `H5S_ALL` (or any other
// invalid handle) the data space is obtained from the supplied \p dataSet.
template<typename ValueType>
static bool checkBufferSize(const ArrayRef<const ValueType> values,
                            const hid_t                     memoryDataSpace,
                            const hid_t                     fileDataSpace,
                            const hid_t                     dataSet)
{
    hssize_t numValues = gmx::ssize(values);
    if constexpr (std::is_same_v<ValueType, BasicVector<float>>
                  || std::is_same_v<ValueType, BasicVector<double>>)
    {
        numValues *= DIM;
    }

    const auto numPointsIn = [&](const hid_t selectionDataSpace)
    {
        if (handleIsValid(selectionDataSpace))
        {
            return H5Sget_select_npoints(selectionDataSpace);
        }
        else
        {
            const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
            return H5Sget_select_npoints(dataSpace);
        }
    };

    return numValues == numPointsIn(memoryDataSpace) && numValues == numPointsIn(fileDataSpace);
}

template<typename ValueType>
static bool reclaimMemory(const ArrayRef<ValueType> values,
                          const hid_t               memoryDataType,
                          const hid_t               dataSpace,
                          const hid_t               dataSet)
{
    // Only try to reclaim memory if any was allocated, otherwise
    // H5Dvlen_reclaim always returns an error
    if (values.empty())
    {
        return true;
    }

    if (handleIsValid(dataSpace))
    {
        return H5Dvlen_reclaim(memoryDataType, dataSpace, H5P_DEFAULT, values.data()) >= 0;
    }
    else if (dataSpace == H5S_ALL)
    {
        // If the data space is H5S_ALL we free the whole set
        const auto [fullDataSpace, fullDataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
        return H5Dvlen_reclaim(memoryDataType, fullDataSpace, H5P_DEFAULT, values.data()) >= 0;
    }
    else
    {
        GMX_THROW(H5mdError("Invalid data space handle: cannot reclaim memory"));
    }
}

static bool readFixedLengthStringsFromDataSet(const hid_t           dataSet,
                                              const hid_t           nativeDataType,
                                              const hid_t           memoryDataSpace,
                                              const hid_t           fileDataSpace,
                                              ArrayRef<std::string> stringValues)
{
    // Get the maximum string size (including the terminating '\0')
    const size_t      maxStringLength = H5Tget_size(nativeDataType);
    std::vector<char> readBuffer(stringValues.size() * maxStringLength);

    const bool readWasSuccessful =
            H5Dread(dataSet, nativeDataType, memoryDataSpace, fileDataSpace, H5P_DEFAULT, readBuffer.data())
            >= 0;

    if (readWasSuccessful)
    {
        int i = 0;
        // HDF5 packs the data for fixed-length strings evenly, so iterate
        // over all strings using maxStringLength as the window.
        for (auto startOfString = readBuffer.cbegin(); startOfString < readBuffer.cend();
             startOfString += maxStringLength)
        {
            const auto endOfString = std::find(startOfString, startOfString + maxStringLength, '\0');
            stringValues[i].assign(startOfString, endOfString);
            ++i;
        }
    }

    return readWasSuccessful;
}

static bool readVariableLengthStringsFromDataSet(const hid_t           dataSet,
                                                 const hid_t           nativeDataType,
                                                 const hid_t           memoryDataSpace,
                                                 const hid_t           fileDataSpace,
                                                 ArrayRef<std::string> stringValues)
{
    // For variable-length strings the HDF5 read operation expects a pointer-to-char-pointers,
    // each of which it will allocate memory for and then read the string data into.
    // We use a scope guard to reclaim any memory that has been allocated after processing.
    std::vector<char*> readBufferPointers(stringValues.size(), nullptr);
    const auto         readBufferPointersGuard = sg::make_scope_guard(
            [&]()
            {
                GMX_H5MD_THROW_UPON_ERROR(
                        !reclaimMemory(makeArrayRef(readBufferPointers), nativeDataType, memoryDataSpace, dataSet),
                        "Cannot reclaim memory after reading variable-size strings");
            });

    const bool readWasSuccessful = H5Dread(dataSet,
                                           nativeDataType,
                                           memoryDataSpace,
                                           fileDataSpace,
                                           H5P_DEFAULT,
                                           readBufferPointers.data())
                                   >= 0;
    if (readWasSuccessful)
    {
        for (int i = 0; i < gmx::ssize(stringValues); ++i)
        {
            stringValues[i].assign(readBufferPointers[i]);
        }
    }
    return readWasSuccessful;
}

static bool writeFixedLengthStringsToDataSet(const hid_t                 dataSet,
                                             const hid_t                 nativeDataType,
                                             const hid_t                 memoryDataSpace,
                                             const hid_t                 fileDataSpace,
                                             ArrayRef<const std::string> stringsToWrite)
{
    // Get the maximum string size (including the terminating '\0')
    const size_t maxStringLength = H5Tget_size(nativeDataType);

    std::vector<char> writeBuffer(stringsToWrite.size() * maxStringLength);
    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        // maxStringLength includes room for the terminating '\0' character, so
        // strncpy will always have room for all normal characters and
        // perhaps also write some null characters, relying on the
        // default initialization of the vector above to provide the
        // null for strings with maxStringLength-1 normal characters.
        std::strncpy(writeBuffer.data() + (i * maxStringLength),
                     stringsToWrite[i].c_str(),
                     maxStringLength - 1);
        GMX_ASSERT(writeBuffer[((i + 1) * maxStringLength) - 1] == '\0',
                   "String must be null terminated");
    }

    return H5Dwrite(
                   dataSet, nativeDataType, memoryDataSpace, fileDataSpace, H5P_DEFAULT, writeBuffer.data())
           >= 0;
}

static bool writeVariableLengthStringsToDataSet(const hid_t                 dataSet,
                                                const hid_t                 nativeDataType,
                                                const hid_t                 memoryDataSpace,
                                                const hid_t                 fileDataSpace,
                                                ArrayRef<const std::string> stringsToWrite)
{
    std::vector<const char*> writeBufferPointers(stringsToWrite.size(), nullptr);
    for (hsize_t i = 0; i < stringsToWrite.size(); ++i)
    {
        writeBufferPointers[i] = stringsToWrite[i].c_str();
    }

    return H5Dwrite(dataSet,
                    nativeDataType,
                    memoryDataSpace,
                    fileDataSpace,
                    H5P_DEFAULT,
                    writeBufferPointers.data())
           >= 0;
}

} // namespace

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(const hid_t dataSetHandle) :
    dataSet_{ dataSetHandle },
    storedDataType_{ H5Dget_type(dataSet_) },
    nativeDataType_{ H5Tget_native_type(storedDataType_, H5T_DIR_DEFAULT) },
    numDims_{ getNumDims(dataSet_) }
{
    GMX_H5MD_THROW_UPON_INVALID_HID(dataSet_, "Invalid handle to data set.");
    verifyDataSetConsistency<ValueType>(storedDataType_, dims());
}

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(const hid_t container, const char* name) :
    H5mdDataSetBase(openDataSet(container, name))
{
    // openDataSet is used above to throw with a unique error message
    // when opening a data set fails
}

template<typename ValueType>
H5mdDataSetBase<ValueType>::~H5mdDataSetBase() noexcept = default;

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(H5mdDataSetBase<ValueType>&&) noexcept = default;

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::id() const
{
    return dataSet_;
}

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::storedDataType() const
{
    return storedDataType_;
}

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::nativeDataType() const
{
    return nativeDataType_;
}

template<typename ValueType>
DataSetDims H5mdDataSetBase<ValueType>::dims() const
{
    DataSetDims dataSetDims(numDims_, 0);
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet_));
    GMX_H5MD_THROW_UPON_ERROR(H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr) < 0,
                              "Could not read dimensions of data set");
    return dataSetDims;
}

template<typename ValueType>
bool H5mdDataSetBase<ValueType>::read(const ArrayRef<ValueType> values,
                                      const hid_t               memoryDataSpace,
                                      const hid_t               fileDataSpace) const
{
    GMX_ASSERT(checkBufferSize(makeConstArrayRef(values), memoryDataSpace, fileDataSpace, dataSet_),
               "Number of points in container of values to read into must "
               "be equal to used hyperslab selection in memory and file");

    if constexpr (std::is_same_v<ValueType, std::string>)
    {
        if (H5Tis_variable_str(this->nativeDataType()) > 0)
        {
            return readVariableLengthStringsFromDataSet(
                    this->id(), this->nativeDataType(), memoryDataSpace, fileDataSpace, values);
        }
        else
        {
            return readFixedLengthStringsFromDataSet(
                    this->id(), this->nativeDataType(), memoryDataSpace, fileDataSpace, values);
        }
    }
    else
    {
        return H5Dread(dataSet_, nativeDataType_, memoryDataSpace, fileDataSpace, H5P_DEFAULT, values.data())
               >= 0;
    }
}

template<typename ValueType>
bool H5mdDataSetBase<ValueType>::write(const ArrayRef<const ValueType> values,
                                       const hid_t                     memoryDataSpace,
                                       const hid_t                     fileDataSpace) const
{
    GMX_ASSERT(checkBufferSize(values, memoryDataSpace, fileDataSpace, dataSet_),
               "Number of points in container of values to read into must "
               "be equal to used hyperslab selection in memory and file");

    if constexpr (std::is_same_v<ValueType, std::string>)
    {
        if (H5Tis_variable_str(this->nativeDataType()) > 0)
        {
            return writeVariableLengthStringsToDataSet(
                    this->id(), this->nativeDataType(), memoryDataSpace, fileDataSpace, values);
        }
        else
        {
            return writeFixedLengthStringsToDataSet(
                    this->id(), this->nativeDataType(), memoryDataSpace, fileDataSpace, values);
        }
    }
    else
    {
        return H5Dwrite(dataSet_, nativeDataType_, memoryDataSpace, fileDataSpace, H5P_DEFAULT, values.data())
               >= 0;
    }
}

template class H5mdDataSetBase<int32_t>;
template class H5mdDataSetBase<int64_t>;
template class H5mdDataSetBase<float>;
template class H5mdDataSetBase<double>;
template class H5mdDataSetBase<BasicVector<float>>;
template class H5mdDataSetBase<BasicVector<double>>;
template class H5mdDataSetBase<std::string>;

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
