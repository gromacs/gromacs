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

/*! \brief Declarations of H5md data type utility routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#ifndef GMX_FILEIO_H5MD_TYPE_H
#define GMX_FILEIO_H5MD_TYPE_H

#include <hdf5.h>

#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

// A note regarding HDF5 native types:
// We here get the native types through defined constants in the library of the form H5T_NATIVE_INT32_g
// for 32-bit integers, and so on. Note the _g suffix! Most of the HDF5 documentation uses the same
// define without that _g suffix to obtain the types. However, these are in fact *macros* which
// first call H5open() to initialize the library, and then return the variable with the _g suffix.
//
// Example:
// #define H5T_NATIVE_INT32 (H5open(), H5T_NATIVE_INT32_g)
//
// We should only ever be using this module when the library has been initialized, so the call to H5open
// is not necessary: we only want the type! This is why the functions below use the _g suffix form,
// while most documentation does not.


/*! \brief Return a handle to a native HDF5 type corresponding to the templated value type.
 *
 * The data type referenced by the handle is immutable and does not need to be closed.
 *
 * \tparam ValueType Type of value.
 * \returns Handle to native HDF5 type.
 */
template<typename ValueType>
hid_t hdf5DataTypeFor() noexcept;

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<int32_t>() noexcept
{
    return H5T_NATIVE_INT32_g;
}

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<int64_t>() noexcept
{
    return H5T_NATIVE_INT64_g;
}

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<float>() noexcept
{
    return H5T_NATIVE_FLOAT_g;
}

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<double>() noexcept
{
    return H5T_NATIVE_DOUBLE_g;
}

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<uint32_t>() noexcept
{
    return H5T_NATIVE_UINT32_g;
}

//! \copydoc hdf5DataTypeFor
template<>
inline hid_t hdf5DataTypeFor<uint64_t>() noexcept
{
    return H5T_NATIVE_UINT64_g;
}

/*! \brief Return a handle to a HDF5 data type for fixed-size strings.
 *
 * The string type is created with a UTF8 character set and null termination, so stored
 * strings may be shorter than the maximum length but a fixed amount of memory is always
 * allocated. Data sets of fixed size strings are contiguous in memory and can be compressed.
 *
 * The returned handle must be closed with a call to H5Tclose to avoid leaking resources.
 *
 * \param[in] maxStringLength Fixed maximum size of string.
 *
 * \throws gmx::FileIOError if \p maxStringLength is not a positive value.
 */
inline hid_t hdf5DataTypeForFixedSizeString(const hsize_t maxStringLength)
{
    const hid_t dataType = H5Tcopy(H5T_C_S1_g);
    H5Tset_cset(dataType, H5T_CSET_UTF8);
    throwUponH5mdError(H5Tset_size(dataType, maxStringLength) < 0,
                       "Invalid fixed-size for string type");
    return dataType;
}

/*! \brief Return true if the template value type matches the HDF5 data type.
 *
 * \tparam ValueType Type of value to check.
 * \param[in] dataType Handle to HDF5 data type to check against.
 * \returns True if the types match, otherwise false.
 */
template<typename ValueType>
bool valueTypeIsDataType(const hid_t dataType) noexcept
{
    const hid_t valueType = hdf5DataTypeFor<ValueType>();
    return H5Tequal(valueType, dataType) > 0;
}

/*! \copydoc valueTypeIsDataType()
 * \brief Specialization of valueTypeIsDataType() for char* type.
 */
template<>
inline bool valueTypeIsDataType<char*>(const hid_t dataType) noexcept
{
    return H5Tget_class(dataType) == H5T_STRING;
}

/*! \copydoc valueTypeIsDataType()
 * \brief Specialization of valueTypeIsDataType() for const char* type.
 */
template<>
inline bool valueTypeIsDataType<const char*>(const hid_t dataType) noexcept
{
    return H5Tget_class(dataType) == H5T_STRING;
}

/*! \copydoc valueTypeIsDataType()
 * \brief Specialization of valueTypeIsDataType() for std::string type.
 */
template<>
inline bool valueTypeIsDataType<std::string>(const hid_t dataType) noexcept
{
    return H5Tget_class(dataType) == H5T_STRING;
}

/*! \brief Return true if the template value type matches the HDF5 data type.
 *
 * Helper function to allow for ValueType inference from a buffer.
 *
 * \tparam ValueType Type of value to check.
 * \param[in] dataType Handle to HDF5 data type to check against.
 * \param[in] valueBuffer Unused buffer which allows the compiler to infer the value type.
 * \returns True if the types match, otherwise false.
 */
template<typename ValueType>
bool valueTypeIsDataType(const hid_t dataType, const ArrayRef<const ValueType> gmx_unused valueBuffer) noexcept
{
    return valueTypeIsDataType<ValueType>(dataType);
}

} // namespace gmx

#endif // GMX_FILEIO_H5MD_TYPE_H
