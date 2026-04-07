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

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_type.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
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
 * 1. That the base primitive of ValueType matches the internal HDF5 data set \p nativeDataType.
 * 2. That the primitive data set dimensions \p dims can store values of the templated \c ValueType.
 *
 * \param[in] nativeDataType Native data type of data set.
 * \param[in] dims           Primitive dimensions of data set.
 *
 * \throws gmx::FileIOError if there is an inconsitency.
 */
template<typename ValueType>
static void verifyDataSetConsistency(const hid_t nativeDataType, const DataSetDims& dims)
{
    if constexpr (std::is_same_v<ValueType, gmx::BasicVector<float>>)
    {
        GMX_H5MD_THROW_UPON_ERROR(dims.empty() || dims.back() != DIM,
                                  "Could not open data set: inner dimension of data set must = 3 "
                                  "for BasicVector<float>");
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<float>(nativeDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
    else if constexpr (std::is_same_v<ValueType, gmx::BasicVector<double>>)
    {
        GMX_H5MD_THROW_UPON_ERROR(dims.empty() || dims.back() != DIM,
                                  "Could not open data set: inner dimension of data set must = 3 "
                                  "for BasicVector<double>");
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<double>(nativeDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
    else
    {
        GMX_H5MD_THROW_UPON_ERROR(!valueTypeIsDataType<ValueType>(nativeDataType),
                                  "Could not open data set: compiled type parameter does not match "
                                  "the primitive type of the data set");
    }
}

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(const hid_t dataSetHandle) :
    dataSet_{ dataSetHandle },
    dataType_{ H5Dget_type(dataSet_) },
    nativeDataType_{ H5Tget_native_type(dataType_, H5T_DIR_DEFAULT) },
    numDims_{ getNumDims(dataSet_) }
{
    GMX_H5MD_THROW_UPON_INVALID_HID(dataSet_, "Invalid handle to data set.");
    verifyDataSetConsistency<ValueType>(nativeDataType_, dims());
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
hid_t H5mdDataSetBase<ValueType>::dataType() const
{
    return dataType_;
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

template class H5mdDataSetBase<int32_t>;

template class H5mdDataSetBase<int64_t>;

template class H5mdDataSetBase<float>;

template class H5mdDataSetBase<double>;

template class H5mdDataSetBase<gmx::BasicVector<float>>;

template class H5mdDataSetBase<gmx::BasicVector<double>>;

template class H5mdDataSetBase<std::string>;

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
